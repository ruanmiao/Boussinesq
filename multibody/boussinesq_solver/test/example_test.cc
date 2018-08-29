#include <memory>
#include <vector>
#include <limits>

#include <gtest/gtest.h>
#include <gflags/gflags.h>

#include "drake/common/find_resource.h"
#include "drake/geometry/mesh_query/vtk_io.h"
#include "drake/geometry/mesh_query/mesh_query.h"
#include "drake/multibody/shapes/geometry.h"
#include "drake/multibody/boussinesq_solver/jacobian_H_matrix.h"
#include "drake/multibody/boussinesq_solver/compliance_matrix.h"
#include "drake/multibody/boussinesq_solver/objects_contact_model.h"
#include "drake/multibody/boussinesq_solver/test_helper.h"
#include "drake/solvers/moby_lcp_solver.h"

#include <iostream>
#include <fstream>

#include <fmt/format.h>
#include <fmt/ostream.h>

#define PRINT_VAR(a) std::cout << #a": " << a << std::endl;

namespace drake {
namespace multibody {
namespace boussinesq_solver {
namespace {

using Eigen::AngleAxisd;
using Eigen::Translation3d;
using Eigen::Isometry3d;
using Eigen::Vector3d;

using geometry::mesh_query::Mesh;
using geometry::mesh_query::LoadMeshFromObj;
using geometry::mesh_query::OutputMeshToVTK;
using geometry::mesh_query::AppendCellCenteredVectorFieldToVTK;
using geometry::mesh_query::AppendCellCenteredVectorFieldToVTK;
using geometry::mesh_query::AppendNodeCenteredVectorFieldToVTK;
using geometry::mesh_query::AppendNodeCenteredScalarFieldToVTK;
using geometry::PenetrationAsTrianglePair;
using geometry::mesh_query::OutputSegmentsToVTK;

int RunSpherePlaneModel(
    double penetration, double sigma,
    double young_modulus_star_sphere,
    double young_modulus_star_plane) {

  // NOTE: Only run meshes 1-4. The rest are garbage.

  // Load mesh for a sphere.
  const bool flip_normals = true;
  std::unique_ptr<Mesh<double>> sphere = LoadMeshFromObj(
      "drake/multibody/boussinesq_solver/test/Mesh_1/sphere.obj", flip_normals);
  sphere->mesh_index = 0;

  std::unique_ptr<Mesh<double>> plane = LoadMeshFromObj(
      "drake/multibody/boussinesq_solver/test/Mesh_1/plane.obj", flip_normals);
  plane->mesh_index = 1;

  const double radius = 1;
  const double z_WSo = radius - penetration;


  // Place sphere a "penetration" distance below z = 0.
  // Apply an arbirary rotation for testing.
  Isometry3d X_WSphere{Translation3d{Vector3d(0, 0, z_WSo)}};
//  (void) z_WSo;
//  Isometry3d X_WSphere{Translation3d{Vector3d(0, 0, 0)}};
  X_WSphere.linear() = MatrixX<double>::Identity(3, 3);

  // The top of the plane is placed at z = 0
  Isometry3d X_WPlane{Translation3d{Vector3d(0, 0, 0.0)}};
  X_WPlane.linear() = MatrixX<double>::Identity(3, 3);

  std::unique_ptr<BoussinesqContactModelResults<double>>
      boussinesq_results_by_force =
      CalcContactSpatialForceBetweenMeshes(
          *sphere, X_WSphere, young_modulus_star_sphere,
          *plane, X_WPlane, young_modulus_star_plane,
          sigma);

  double young_modulus_star = 1.0 /
      (1.0 / young_modulus_star_sphere + 1.0 / young_modulus_star_plane);

  double force_Hertz = 4.0 / 3.0 * young_modulus_star *
      std::sqrt(radius) * pow(penetration, 3.0 / 2.0);



  auto MakeDeformationVectorField = [](const Mesh<double>& patch, const VectorX<double> u_normal) {
    std::vector<Vector3<double>> u;
    for (int i = 0; i < u_normal.size(); ++i) {
      u.push_back(patch.node_normals_G[i] * u_normal[i]);
    }
    return u;
  };


  const int num_nodes_A = boussinesq_results_by_force->object_A_patch->points_G.size();
  const int num_nodes_B = boussinesq_results_by_force->object_B_patch->points_G.size();
  const int num_nodes = num_nodes_A + num_nodes_B;
  (void) num_nodes;


  std::ofstream patch_file("sphere_patch.vtk");
  OutputMeshToVTK(patch_file, boussinesq_results_by_force->object_A_patch->points_G, boussinesq_results_by_force->object_A_patch->triangles, X_WSphere);
  patch_file << "POINT_DATA " << num_nodes_A << std::endl;
  AppendNodeCenteredScalarFieldToVTK(patch_file, "deformation", boussinesq_results_by_force->deformation_patch_A);
  AppendNodeCenteredScalarFieldToVTK(patch_file, "normal_stress", boussinesq_results_by_force->pressure_patch_A);
  AppendNodeCenteredVectorFieldToVTK(
      patch_file, "normals", boussinesq_results_by_force->object_A_patch->node_normals_G, X_WSphere);
  auto u_A = MakeDeformationVectorField(
      *boussinesq_results_by_force->object_A_patch, boussinesq_results_by_force->deformation_patch_A);
  AppendNodeCenteredVectorFieldToVTK(patch_file, "u", u_A, X_WSphere);
  patch_file.close();

  patch_file.open("plane_patch.vtk", std::ios::out);
  OutputMeshToVTK(patch_file, boussinesq_results_by_force->object_B_patch->points_G, boussinesq_results_by_force->object_B_patch->triangles, X_WPlane);
  patch_file << "POINT_DATA " << num_nodes_B << std::endl;
  AppendNodeCenteredScalarFieldToVTK(patch_file, "deformation", boussinesq_results_by_force->deformation_patch_B);
  AppendNodeCenteredScalarFieldToVTK(patch_file, "normal_stress", boussinesq_results_by_force->pressure_patch_B);
  AppendNodeCenteredVectorFieldToVTK(
      patch_file, "normals", boussinesq_results_by_force->object_B_patch->node_normals_G, X_WPlane);
  auto u_B = MakeDeformationVectorField(
      *boussinesq_results_by_force->object_B_patch, boussinesq_results_by_force->deformation_patch_B);
  AppendNodeCenteredVectorFieldToVTK(patch_file, "u", u_B, X_WPlane);
  patch_file.close();

  // Compute some statistics and print them out to the std::out
  double max_uA_from_phi = 0;
  for (int i = 0; i < num_nodes_A; i++) {
    double u = u_A[i](2);
    if (u > max_uA_from_phi)
      max_uA_from_phi = u;
  }
  double max_uB_from_phi = 0;
  for (int i = 0; i < num_nodes_B; i++) {
    double u = u_B[i](2);
    if (-u > max_uB_from_phi)
      max_uB_from_phi = -u;
  }
  PRINT_VAR(max_uA_from_phi);
  PRINT_VAR(max_uB_from_phi);

  const auto& results = *boussinesq_results_by_force->contact_results;
  std::vector<Vector3d> pointsA(results.size());
  std::transform(results.begin(), results.end(), pointsA.begin(),
                 [](const PenetrationAsTrianglePair<double>& pair) {
                   return pair.p_WoAs_W;
                 });

  std::vector<Vector3d> pointsB(results.size());
  std::transform(results.begin(), results.end(), pointsB.begin(),
                 [](const PenetrationAsTrianglePair<double>& pair) {
                   return pair.p_WoBs_W;
                 });

  std::ofstream file("pairs.vtk");
  OutputSegmentsToVTK(file, pointsA, pointsB);


  VectorX<double> kkt_twice(2*boussinesq_results_by_force->kkt_multipliers.size());
  VectorX<double> phi0_twice(2*boussinesq_results_by_force->phi0.size());
  VectorX<double> phi_u_twice(2*boussinesq_results_by_force->phi0.size());
  VectorX<double> phi_twice(2*boussinesq_results_by_force->phi0.size());
  for (int i=0;i<boussinesq_results_by_force->kkt_multipliers.size(); ++i) {
    kkt_twice(2*i) = boussinesq_results_by_force->kkt_multipliers(i);
    kkt_twice(2*i+1) = boussinesq_results_by_force->kkt_multipliers(i);

    phi0_twice(2*i) = boussinesq_results_by_force->phi0(i);
    phi0_twice(2*i+1) = boussinesq_results_by_force->phi0(i);

    phi_u_twice(2*i) = boussinesq_results_by_force->phi_u(i);
    phi_u_twice(2*i+1) = boussinesq_results_by_force->phi_u(i);

    phi_twice(2*i) = phi0_twice(2*i) + phi_u_twice(2*i);
    phi_twice(2*i+1) = phi0_twice(2*i + 1) + phi_u_twice(2*i + 1);
  }
  file << "POINT_DATA " << kkt_twice.rows() << std::endl;
  AppendNodeCenteredScalarFieldToVTK(file, "kkt_multipliers", kkt_twice);
  AppendNodeCenteredScalarFieldToVTK(file, "phi0", phi0_twice);
  AppendNodeCenteredScalarFieldToVTK(file, "phi_u", phi_u_twice);  // phi_u = W*pi
  AppendNodeCenteredScalarFieldToVTK(file, "phi", phi_twice);      // phi_u = phi0 + W*pi = phi0 + phi_u
  file.close();

  PRINT_VAR(force_Hertz);
  return 0;
}

GTEST_TEST(ExampleTest, SpherePlane) {

//  std::vector<double> indentations{0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.1};
  std::vector<double> indentations{0.1};

//  std::vector<double> vec_young_modulus_star_sphere{1000, 1, 1};
//  std::vector<double> vec_young_modulus_star_plane{1, 1000, 1};
  std::vector<double> vec_young_modulus_star_sphere{1.0};
  std::vector<double> vec_young_modulus_star_plane{1000.0};

  int num_indentations = indentations.size();
  int num_youngs = vec_young_modulus_star_plane.size();

  for (int i = 0; i < num_indentations; i++) {

    double indentation = indentations[i];
    std::cout << "*********** Indentation: "
              << indentation
              << "*****************."
              << std::endl;

    for(int i_param = 0; i_param < num_youngs; i_param++) {

      double young_modulus_star_sphere = vec_young_modulus_star_sphere[i_param];
      double young_modulus_star_plane = vec_young_modulus_star_plane[i_param];

      std::cout << "----------E_sphere: "
                << young_modulus_star_sphere
                << "        E_plane: "
                   << young_modulus_star_plane
                << " -------------------------"
                << std::endl;

      RunSpherePlaneModel(indentation, 0.1 /* sigma */,
                          young_modulus_star_sphere, young_modulus_star_plane);
    }



    std::cout << std::endl;
  }

}

}
}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake





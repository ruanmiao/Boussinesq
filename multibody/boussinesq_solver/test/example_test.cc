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
    std::ofstream& file_force_sphere, std::ofstream& file_force_plane,
    std::ofstream& file_penetration, std::ofstream& file_force_hertz) {
  const bool flip_normals = true;

  // Load mesh for a sphere.
  std::unique_ptr<Mesh<double>> sphere = LoadMeshFromObj(
      "drake/multibody/boussinesq_solver/test/Mesh_2/sphere.obj", flip_normals);
  sphere->mesh_index = 0;

  std::unique_ptr<Mesh<double>> plane = LoadMeshFromObj(
      "drake/multibody/boussinesq_solver/test/Mesh_2/plane.obj", flip_normals);
  plane->mesh_index = 1;

  const double radius = 1;
  const double z_WSo = radius - penetration;

  const double young_modulus_star_A = 1.0;
  const double young_modulus_star_B = 1.0;
//  const double young_modulus_star_B = std::numeric_limits<double>::infinity();

  // Place sphere a "penetration" distance below z = 0.
  // Apply an arbirary rotation for testing.
  Isometry3d X_WSphere{Translation3d{Vector3d(0, 0, z_WSo)}};
  X_WSphere.linear() = MatrixX<double>::Identity(3, 3);

  // The top of the plane is placed at z = 0
  Isometry3d X_WPlane{Translation3d{Vector3d(0, 0, 0.0)}};
  X_WPlane.linear() = MatrixX<double>::Identity(3, 3);


  std::unique_ptr<BoussinesqContactModelResults<double>> boussinesq_results =
      CalcContactSpatialForceBetweenMeshes(
          *sphere, X_WSphere, young_modulus_star_A,
          *plane, X_WPlane, young_modulus_star_B,
          sigma);


  Vector3<double> FA = boussinesq_results->F_Ao_W.translational();
  Vector3<double> FB = boussinesq_results->F_Bo_W.translational();

  double force_Hertz = 4.0 / 3.0 /
      (1 / young_modulus_star_A + 1 / young_modulus_star_B) *
      std::sqrt(radius) * pow(penetration, 3.0 / 2.0);
  double ratio_FA_sphere_over_Hertz = FA(2) / force_Hertz;
  double ratio_FB_sphere_over_Hertz = FB(2) / force_Hertz;

  file_force_sphere << fmt::format("{:.8f}, {:.8f}, {:.8f}\n",
                                   FA(0), FA(1), FA(2));

  file_force_plane << fmt::format("{:.8f}, {:.8f}, {:.8f}\n",
                                  FB(0), FB(1), FB(2));

  file_force_hertz << fmt::format("{:.8f}, {:.8f}, {:.8f}\n",
                                  force_Hertz,
                                  ratio_FA_sphere_over_Hertz,
                                  ratio_FB_sphere_over_Hertz);
  file_penetration << fmt::format(
      "{:.8f}, {:d}, {:d}\n",
      penetration,
      boussinesq_results->object_A_patch->triangles.size(),
      boussinesq_results->object_B_patch->triangles.size());


//  std::ofstream patch_file("sphere_patch.vtk");
//  OutputMeshToVTK(patch_file, boussinesq_results->object_A_patch->points_G, boussinesq_results->object_A_patch->triangles, X_WSphere);
//  AppendNodeCenteredScalarFieldToVTK(patch_file, "normals_stress", boussinesq_results->pressure_patch_A);
//  patch_file.close();
//
//  patch_file.open("plane_patch.vtk", std::ios::out);
//  OutputMeshToVTK(patch_file, boussinesq_results->object_B_patch->points_G, boussinesq_results->object_B_patch->triangles, X_WPlane);
//  AppendNodeCenteredScalarFieldToVTK(patch_file, "normals_stress", boussinesq_results->pressure_patch_B);
//  patch_file.close();

  auto MakeDeformationVectorField = [](const Mesh<double>& patch, const VectorX<double> u_normal) {
    std::vector<Vector3<double>> u;
    for (int i = 0; i < u_normal.size(); ++i) {
      u.push_back(patch.node_normals_G[i] * u_normal[i]);
    }
    return u;
  };

  std::ofstream patch_file("sphere_patch.vtk");
  OutputMeshToVTK(patch_file, boussinesq_results->object_A_patch->points_G, boussinesq_results->object_A_patch->triangles, X_WSphere);
  patch_file << "POINT_DATA " << boussinesq_results->object_A_patch->points_G.size() << std::endl;
  AppendNodeCenteredScalarFieldToVTK(patch_file, "deformation", boussinesq_results->deformation_patch_A);
  AppendNodeCenteredScalarFieldToVTK(patch_file, "normal_stress", boussinesq_results->pressure_patch_A);
  AppendNodeCenteredVectorFieldToVTK(
      patch_file, "normals", boussinesq_results->object_A_patch->node_normals_G, X_WSphere);
  auto u_A = MakeDeformationVectorField(
      *boussinesq_results->object_A_patch, boussinesq_results->deformation_patch_A);
  AppendNodeCenteredVectorFieldToVTK(patch_file, "u", u_A, X_WSphere);
  patch_file.close();

  patch_file.open("plane_patch.vtk", std::ios::out);
  OutputMeshToVTK(patch_file, boussinesq_results->object_B_patch->points_G, boussinesq_results->object_B_patch->triangles, X_WPlane);
  patch_file << "POINT_DATA " << boussinesq_results->object_B_patch->points_G.size() << std::endl;
  AppendNodeCenteredScalarFieldToVTK(patch_file, "deformation", boussinesq_results->deformation_patch_B);
  AppendNodeCenteredScalarFieldToVTK(patch_file, "normal_stress", boussinesq_results->pressure_patch_B);
  AppendNodeCenteredVectorFieldToVTK(
      patch_file, "normals", boussinesq_results->object_B_patch->node_normals_G, X_WPlane);
  auto u_B = MakeDeformationVectorField(
      *boussinesq_results->object_B_patch, boussinesq_results->deformation_patch_B);
  AppendNodeCenteredVectorFieldToVTK(patch_file, "u", u_B, X_WPlane);
  patch_file.close();


  const auto& results = *boussinesq_results->contact_results;
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
  VectorX<double> kkt_twice(2*boussinesq_results->kkt_multipliers.size());
  for (int i=0;i<boussinesq_results->kkt_multipliers.size(); ++i) {
    kkt_twice(2*i) = boussinesq_results->kkt_multipliers(i);
    kkt_twice(2*i+1) = boussinesq_results->kkt_multipliers(i);
  }
  AppendNodeCenteredScalarFieldToVTK(file, "kkt_multipliers", kkt_twice);
  file.close();


  return 0;
}




GTEST_TEST(ExampleTest, SpherePlane) {

  bool control_run = false;
  if (control_run) {
    std::vector<double> penetrations{0.001, 0.002, 0.003, 0.004, 0.005, 0.006,
                                     0.008, 0.009, 0.01, 0.015, 0.016, 0.018,
                                     0.019, 0.02, 0.021, 0.022, 0.024, 0.025,
                                     0.03, 0.04, 0.05, 0.06, 0.08, 0.1, 0.12,
                                     0.15, 0.16, 0.18};

    int size_data = penetrations.size();
    std::ofstream file_force_sphere("force_sphere.txt");
    std::ofstream file_force_plane("force_plane.txt");
    std::ofstream file_force_hertz("force_hertz.txt");
    std::ofstream file_penetration("penetration.txt");
    // std::ofstream file_penetration("patch_size.txt");

    file_force_sphere << "element size approximates 0.1. Sigma > 0.2 "
                      << std::endl;
    file_force_plane << "element size approximates 0.1. Sigma > 0.2 "
                     << std::endl;
    file_penetration << "element size approximates 0.1. Sigma > 0.2 "
                     << std::endl;
    file_force_hertz << "element size approximates 0.1. Sigma > 0.2 "
                     << std::endl;

//  for (int i = 0; i < size_data; i++) {
//    RunSpherePlaneModel(penetrations[i], 0.0,
//                        file_force_sphere, file_force_plane,
//                        file_penetration, file_force_hertz);
//  }

    for (int i = 0; i < size_data; i++) {
      RunSpherePlaneModel(penetrations[i], 0.2 + penetrations[i] / 2.0,
                          file_force_sphere, file_force_plane,
                          file_penetration, file_force_hertz);
    }

    file_force_sphere.close();
    file_force_plane.close();
    file_penetration.close();
    file_force_hertz.close();
  }
  std::ofstream file_force_sphere("force_sphere.txt");
  std::ofstream file_force_plane("force_plane.txt");
  std::ofstream file_force_hertz("force_hertz.txt");
  std::ofstream file_penetration("penetration.txt");

  RunSpherePlaneModel(0.1, 0.1,
                      file_force_sphere, file_force_plane,
                      file_penetration, file_force_hertz);
  file_force_sphere.close();
  file_force_plane.close();
  file_penetration.close();
  file_force_hertz.close();
}

}
}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake





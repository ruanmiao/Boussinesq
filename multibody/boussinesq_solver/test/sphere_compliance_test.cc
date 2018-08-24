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


int SolveCaTimesF() {


  // Load mesh, strip and curved strip.
//  const bool flip_normals = false;
//  std::unique_ptr<Mesh<double>> sphere = LoadMeshFromObj(
//      "drake/multibody/boussinesq_solver/test/Mesh_4/strip_curved.obj", flip_normals);
//  sphere->mesh_index = 0;
//
//  std::unique_ptr<Mesh<double>> plane = LoadMeshFromObj(
//      "drake/multibody/boussinesq_solver/test/Mesh_4/strip.obj", flip_normals);
//  plane->mesh_index = 1;

  // Load mesh for a sphere.
  const bool flip_normals = true;
  std::unique_ptr<Mesh<double>> sphere = LoadMeshFromObj(
      "drake/multibody/boussinesq_solver/test/Mesh_4/sphere.obj", flip_normals);
  sphere->mesh_index = 0;

  std::unique_ptr<Mesh<double>> plane = LoadMeshFromObj(
      "drake/multibody/boussinesq_solver/test/Mesh_4/plane.obj", flip_normals);
  plane->mesh_index = 1;


  const double radius = 1;
  const double penetration = 0.1;
  const double z_WSo = radius - penetration;
  const double sigma = 0.1;

  const double young_modulus_star_sphere = 1.0;
  const double young_modulus_star_plane = 1.0;

  double young_modulus_star = 1.0;

//  double young_modulus_star = 1.0 /
//      (1.0 / young_modulus_star_sphere + 1.0 / young_modulus_star_plane);

  double force_Hertz = 4.0 / 3.0 * young_modulus_star *
      std::sqrt(radius) * pow(penetration, 3.0 / 2.0);

  PRINT_VAR(force_Hertz);

  // Place sphere a "penetration" distance below z = 0.
  // Apply an arbirary rotation for testing.
  Isometry3d X_WSphere{Translation3d{Vector3d(0, 0, z_WSo)}};
  X_WSphere.linear() = MatrixX<double>::Identity(3, 3);

  // The top of the plane is placed at z = 0
  Isometry3d X_WPlane{Translation3d{Vector3d(0, 0, 0.0)}};
  X_WPlane.linear() = MatrixX<double>::Identity(3, 3);


  auto owned_results =
      std::make_unique<std::vector<PenetrationAsTrianglePair<double>>>();
  std::vector<PenetrationAsTrianglePair<double>>& results = *owned_results;


  results =
      geometry::mesh_query::MeshToMeshQuery(X_WSphere, *sphere, X_WPlane, *plane, sigma);

  auto patches =
      geometry::mesh_query::MakeLocalPatchMeshes(&results, *sphere, *plane);
  std::unique_ptr<Mesh<double>> object_A_patch = std::move(patches.first);
  std::unique_ptr<Mesh<double>> object_B_patch = std::move(patches.second);

  const int num_nodes_A = object_A_patch->points_G.size();
  const int num_nodes_B = object_B_patch->points_G.size();
  const int num_nodes = num_nodes_A + num_nodes_B;

  PRINT_VAR(num_nodes_A);
  PRINT_VAR(num_nodes_B);


  // Use LCP to solve the deformation and the pressure, and force
  const MatrixX<double> C_A = CalcComplianceMatrix(
      object_A_patch->points_G,
      object_A_patch->triangles,
      1 / ( young_modulus_star_sphere * M_PI));

  const MatrixX<double> C_B = CalcComplianceMatrix(
      object_B_patch->points_G,
      object_B_patch->triangles,
      1 / (young_modulus_star_plane * M_PI));

  VectorX <double> area_matrix = VectorX<double>::Zero(num_nodes);
  VectorX <double> area_inv_matrix = VectorX<double>::Zero(num_nodes);
  for(int i = 0; i < num_nodes_A; i++) {
    area_matrix(i) = object_A_patch->node_areas(i);
    area_inv_matrix(i) = 1.0 / object_A_patch->node_areas(i);
  }
  for(int i = 0; i < num_nodes_B; i++) {
    area_matrix(num_nodes_A + i) =
        object_B_patch->node_areas(i);
    area_inv_matrix(num_nodes_A + i) =
        1.0 / object_B_patch->node_areas(i);
  }

  const VectorX<double> area_A
      = area_matrix.head(num_nodes_A);
  const VectorX<double> area_B
      = area_matrix.tail(num_nodes_B);

  const VectorX<double> area_inv_A
      = area_inv_matrix.head(num_nodes_A);
  const VectorX<double> area_inv_B
      = area_inv_matrix.tail(num_nodes_B);

  const MatrixX<double> Ca_A = C_A * area_inv_A.asDiagonal();
  const MatrixX<double> Ca_B = C_B * area_inv_B.asDiagonal();


  // Test 1. h0 defined be the distance between the 2 surfaces in the z-axis direction

  VectorX<double> h0_A_normal_to_plane(num_nodes_A);
  for (int i_gap = 0; i_gap < num_nodes_A; i_gap++) {
    const Vector3d pos = X_WSphere * object_A_patch->points_G.at(i_gap);
    h0_A_normal_to_plane(i_gap) = pos(2);
  }

  VectorX<double> h0_B_normal_to_plane(num_nodes_B);
  for (int i_gap = 0; i_gap < num_nodes_B; i_gap++) {
    const Vector3d pos = X_WPlane * object_B_patch->points_G.at(i_gap);
    double r = sqrt(pow(pos(0), 2) + pow(pos(1), 2));
    h0_B_normal_to_plane(i_gap) = z_WSo - sqrt(pow(radius, 2) - std::min(pow(r, 2), pow(radius, 2)));
  }

  // Calculate the h0 based on Hook's law/

  double ratio_A = 1.0;
  double ratio_B = 1.0;


//  double ratio_A = (1/young_modulus_star_sphere) / (1/young_modulus_star);
//  double ratio_B = (1/young_modulus_star_plane) / (1/young_modulus_star);

//  (void) ratio_A;
//  (void) ratio_B;

  VectorX<double> h0_A_ratio_normal_to_plane = h0_A_normal_to_plane * ratio_A;
  VectorX<double> h0_B_ratio_normal_to_plane = h0_B_normal_to_plane * ratio_B;

  //Result by solving h0 + Ca * f with LCP, with h0 defined by Hook's law
  solvers::MobyLCPSolver<double> moby_LCP_solver;
  VectorX<double> force_sol_A(num_nodes_A);
  VectorX<double> force_sol_B(num_nodes_B);

  bool solved = moby_LCP_solver.SolveLcpLemke(
      Ca_A, h0_A_ratio_normal_to_plane, &force_sol_A);
  DRAKE_DEMAND(solved);
  solved = moby_LCP_solver.SolveLcpLemke(
      Ca_B, h0_B_ratio_normal_to_plane, &force_sol_B);
  DRAKE_DEMAND(solved);

  double force_A_test1 = force_sol_A.sum();
  double force_B_test1 = force_sol_B.sum();

  std::cout
      << "Test1: Result by solving h0 + Ca * f with LCP, with h0 defined by the "
      << "distance between the two surfaces in the z-axis direction normal "
      << "to the plane"
      << std::endl;
  PRINT_VAR(force_A_test1);
  PRINT_VAR(force_B_test1);
  std::cout << std::endl;



  // Test 2. h0 defined be the distance between the 2 surfaces in the radius direction
  VectorX<double> h0_A_normal_to_sphere(num_nodes_A);
  for (int i_gap = 0; i_gap < num_nodes_A; i_gap++) {
    const Vector3d pos = X_WSphere * object_A_patch->points_G.at(i_gap);
    h0_A_normal_to_sphere(i_gap) = radius * pos(2) / (z_WSo - pos(2));
  }

  VectorX<double> h0_B_normal_to_sphere(num_nodes_B);
  for (int i_gap = 0; i_gap < num_nodes_B; i_gap++) {
    const Vector3d pos = X_WPlane * object_B_patch->points_G.at(i_gap);
    double r = sqrt(pow(pos(0), 2) + pow(pos(1), 2) + pow(pos(2), 2));
    h0_B_normal_to_sphere(i_gap) = r - radius;
  }

//  (void) ratio_A;
//  (void) ratio_B;

  VectorX<double> h0_A_ratio_normal_to_sphere = h0_A_normal_to_sphere * ratio_A;
  VectorX<double> h0_B_ratio_normal_to_sphere = h0_B_normal_to_sphere * ratio_B;

  //Results by solving h0 + Cp with LCP, with h0 defined by Hook's law

  solved = moby_LCP_solver.SolveLcpLemke(
      Ca_A, h0_A_ratio_normal_to_sphere, &force_sol_A);
  DRAKE_DEMAND(solved);
  solved = moby_LCP_solver.SolveLcpLemke(
      Ca_B, h0_B_ratio_normal_to_sphere, &force_sol_B);
  DRAKE_DEMAND(solved);

  double force_A_test2 = force_sol_A.sum();
  double force_B_test2 = force_sol_B.sum();

  std::cout
      << "Test2: Result by solving h0 + Ca * f with LCP, with h0 defined by the "
      << "distance between the two surfaces in the radius direction normal "
      << "to the sphere surface"
      << std::endl;
  PRINT_VAR(force_A_test2);
  PRINT_VAR(force_B_test2);
  std::cout << std::endl;




//  std::ofstream data_file("compliance_sphere.txt");
//
//  data_file << "Test1: solve h0 + Ca * f with LCP, with h0 defined by the "
//            << "distance between the two surfaces in the z-axis direction normal "
//            << "to the plane."
//            << std::endl;

  return 0;

}



//GTEST_TEST(ExampleTest, ComplianceWithAreaInv) {
//
//  SolveCaTimesF();
//
//}


}
}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake

int main(int argc, char* argv[]) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  return drake::multibody::boussinesq_solver::SolveCaTimesF();
}




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
    std::ofstream& file_force_sphere, std::ofstream& file_force_plane,
    std::ofstream& file_penetration, std::ofstream& file_force_hertz) {

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
      "drake/multibody/boussinesq_solver/test/Mesh_3/sphere.obj", flip_normals);
  sphere->mesh_index = 0;

  std::unique_ptr<Mesh<double>> plane = LoadMeshFromObj(
      "drake/multibody/boussinesq_solver/test/Mesh_3/plane.obj", flip_normals);
  plane->mesh_index = 1;




  // Flip the two meshes
// const bool flip_normals = true;
//  std::unique_ptr<Mesh<double>> plane = LoadMeshFromObj(
//      "drake/multibody/boussinesq_solver/test/Mesh_3/sphere.obj", flip_normals);
//  plane->mesh_index = 0;
//
//  std::unique_ptr<Mesh<double>> sphere = LoadMeshFromObj(
//      "drake/multibody/boussinesq_solver/test/Mesh_3/plane.obj", flip_normals);
//  sphere->mesh_index = 1;


  const double radius = 1;
  const double z_WSo = radius - penetration;

  const double young_modulus_star_sphere = 1.0;
  const double young_modulus_star_plane = 1000.0;

  // Place sphere a "penetration" distance below z = 0.
  // Apply an arbirary rotation for testing.
  Isometry3d X_WSphere{Translation3d{Vector3d(0, 0, z_WSo)}};
//  (void) z_WSo;
//  Isometry3d X_WSphere{Translation3d{Vector3d(0, 0, 0)}};
  X_WSphere.linear() = MatrixX<double>::Identity(3, 3);

  // The top of the plane is placed at z = 0
  Isometry3d X_WPlane{Translation3d{Vector3d(0, 0, 0.0)}};
  X_WPlane.linear() = MatrixX<double>::Identity(3, 3);


  // Flip the two mesh, now X_WPlane is for Object B, which is now the spere.
//  Isometry3d X_WPlane{Translation3d{Vector3d(0, 0, z_WSo)}};
//  X_WPlane.linear() = MatrixX<double>::Identity(3, 3);
//
//  // The top of the plane is placed at z = 0
//  Isometry3d X_WSphere{Translation3d{Vector3d(0, 0, 0.0)}};
//  X_WSphere.linear() = MatrixX<double>::Identity(3, 3);



  std::unique_ptr<BoussinesqContactModelResults<double>>
      boussinesq_results_by_force =
      CalcContactSpatialForceBetweenMeshes(
          *sphere, X_WSphere, young_modulus_star_sphere,
          *plane, X_WPlane, young_modulus_star_plane,
          sigma);



//  std::unique_ptr<BoussinesqContactModelResults<double>>
//      boussinesq_results_by_pressure =
//      CalcContactSpatialForceBetweenMeshesByPressure(
//          *sphere, X_WSphere, young_modulus_star_sphere,
//          *plane, X_WPlane, young_modulus_star_plane,
//          sigma);


  Vector3<double> FA_by_force = boussinesq_results_by_force->F_Ao_W.translational();
  Vector3<double> FB_by_force = boussinesq_results_by_force->F_Bo_W.translational();

//  Vector3<double> FA_by_pressure = boussinesq_results_by_pressure->F_Ao_W.translational();
//  Vector3<double> FB_by_pressure = boussinesq_results_by_pressure->F_Bo_W.translational();

  double young_modulus_star = 1.0 /
      (1.0 / young_modulus_star_sphere + 1.0 / young_modulus_star_plane);

  double force_Hertz = 4.0 / 3.0 * young_modulus_star *
      std::sqrt(radius) * pow(penetration, 3.0 / 2.0);

  double ratio_FA_sphere_over_Hertz = FA_by_force(2) / force_Hertz;
  double ratio_FB_sphere_over_Hertz = FB_by_force(2) / force_Hertz;

  file_force_sphere << "Result solved when kkt-multiplier is force: "
                       << std::endl;
  file_force_sphere << fmt::format("{:.8f}, {:.8f}, {:.8f}\n",
                                   FA_by_force(0),
                                   FA_by_force(1), FA_by_force(2));
//  file_force_sphere << "Result solved when kkt-multiplier is pressure: "
//                    << std::endl;
//  file_force_sphere << fmt::format("{:.8f}, {:.8f}, {:.8f}\n",
//                                   FA_by_pressure(0),
//                                   FA_by_pressure(1), FA_by_pressure(2));


  file_force_plane << "Result solved when kkt-multiplier is force: "
                    << std::endl;
  file_force_plane << fmt::format("{:.8f}, {:.8f}, {:.8f}\n",
                                  FB_by_force(0), FB_by_force(1), FB_by_force(2));
//  file_force_plane << "Result solved when kkt-multiplier is pressure: "
//                    << std::endl;
//  file_force_plane << fmt::format("{:.8f}, {:.8f}, {:.8f}\n",
//                                  FB_by_pressure(0),
//                                  FB_by_pressure(1), FB_by_pressure(2));

  file_force_hertz << fmt::format("{:.8f}, {:.8f}, {:.8f}\n",
                                  force_Hertz,
                                  ratio_FA_sphere_over_Hertz,
                                  ratio_FB_sphere_over_Hertz);
  file_penetration << fmt::format(
      "{:.8f}, {:d}, {:d}\n",
      penetration,
      boussinesq_results_by_force->object_A_patch->triangles.size(),
      boussinesq_results_by_force->object_B_patch->triangles.size());


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



// Deformation by kkt
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
  AppendNodeCenteredScalarFieldToVTK(file, "phi_u", phi_u_twice);
  AppendNodeCenteredScalarFieldToVTK(file, "phi", phi_twice);
  file.close();




#if 0


  // Use LCP to solve the deformation and the pressure, and force
  const MatrixX<double> C_A = CalcComplianceMatrix(
      boussinesq_results_by_force->object_A_patch->points_G,
      boussinesq_results_by_force->object_A_patch->triangles,
      1 / ( young_modulus_star_sphere * M_PI));

  const MatrixX<double> C_B = CalcComplianceMatrix(
      boussinesq_results_by_force->object_B_patch->points_G,
      boussinesq_results_by_force->object_B_patch->triangles,
      1 / (young_modulus_star_plane * M_PI));

  MatrixX<double> area_matrix = MatrixX<double>::Zero(num_nodes, num_nodes);
  for(int i = 0; i < num_nodes_A; i++) {
    area_matrix(i, i) = boussinesq_results_by_force->object_A_patch->node_areas(i);
  }
  for(int i = 0; i < num_nodes_B; i++) {
    area_matrix(num_nodes_A + i, num_nodes_A + i) =
        boussinesq_results_by_force->object_B_patch->node_areas(i);
  }

  const MatrixX<double> area_A
      = area_matrix.block(0, 0, num_nodes_A, num_nodes_A);
  const MatrixX<double> area_B
      = area_matrix.block(num_nodes_A, num_nodes_A, num_nodes_B, num_nodes_B);
  const MatrixX<double> Ca_A = C_A * area_A.inverse();
  const MatrixX<double> Ca_B = C_B * area_B.inverse();






  VectorX<double> h0_gap_A(num_nodes_A);
  for (int i_gap = 0; i_gap < num_nodes_A; i_gap++) {
    const Vector3d pos = X_WSphere * boussinesq_results_by_force->object_A_patch->points_G.at(i_gap);
    h0_gap_A(i_gap) = pos(2);
  }

  VectorX<double> h0_gap_B(num_nodes_B);
  for (int i_gap = 0; i_gap < num_nodes_B; i_gap++) {
    const Vector3d pos = X_WPlane * boussinesq_results_by_force->object_B_patch->points_G.at(i_gap);
    double r = sqrt(pow(pos(0), 2) + pow(pos(1), 2));
    h0_gap_B(i_gap) = z_WSo - sqrt(pow(radius, 2) - std::min(pow(r, 2), pow(radius, 2)));
  }


  // Calculate the h0 based on Hook's law/

  double ratio_A = (1/young_modulus_star_sphere) / (1/young_modulus_star);
  double ratio_B = (1/young_modulus_star_plane) / (1/young_modulus_star);

//  (void) ratio_A;
//  (void) ratio_B;

  VectorX<double> h0_gap_A_ratio = h0_gap_A * ratio_A;
  VectorX<double> h0_gap_B_ratio = h0_gap_B * ratio_B;



  //Result by solving h0 + Ca * f with LCP, with h0 defined by Hook's law

  solvers::MobyLCPSolver<double> moby_LCP_solver;
  VectorX<double> force_sol_A(num_nodes_A);
  VectorX<double> force_sol_B(num_nodes_B);


  bool solved = moby_LCP_solver.SolveLcpLemke(
      Ca_A, h0_gap_A_ratio, &force_sol_A);
  DRAKE_DEMAND(solved);
  solved = moby_LCP_solver.SolveLcpLemke(
      Ca_B, h0_gap_B_ratio, &force_sol_B);
  DRAKE_DEMAND(solved);


  double force_A = force_sol_A.sum();
  double force_B = force_sol_B.sum();

  std::cout
      << "Result by solving h0 + Ca * f with LCP, with h0 defined by Hook's law"
         << std::endl;
  PRINT_VAR(force_A);
  PRINT_VAR(force_B);
  std::cout << std::endl;



  //Results by solving h0 + Cp with LCP, with h0 defined by Hook's law

  VectorX<double> pressure_sol_A(num_nodes_A);
  VectorX<double> pressure_sol_B(num_nodes_B);

  solved = moby_LCP_solver.SolveLcpLemke(
      C_A, h0_gap_A_ratio, &pressure_sol_A);
  DRAKE_DEMAND(solved);
  solved = moby_LCP_solver.SolveLcpLemke(
      C_B, h0_gap_B_ratio, &pressure_sol_B);
  DRAKE_DEMAND(solved);

  VectorX<double> u_A_test(num_nodes_A);
  VectorX<double> u_B_test(num_nodes_B);


  u_A_test = C_A * pressure_sol_A;
  u_B_test = C_B * pressure_sol_B;

  double max_uA = 0;
  for (int i = 0; i < num_nodes_A; i++) {
    if (u_A_test(i) > max_uA) max_uA = u_A_test(i);
  }
  double max_uB = 0;
  for (int i = 0; i < num_nodes_B; i++) {
    if (u_B_test(i) > max_uB) max_uB = u_B_test(i);
  }


  double total_force_A = CalcForceOverMesh(
      boussinesq_results_by_force->object_A_patch->points_G,
      boussinesq_results_by_force->object_A_patch->triangles,
      pressure_sol_A);

  double total_force_B = CalcForceOverMesh(
      boussinesq_results_by_force->object_B_patch->points_G,
      boussinesq_results_by_force->object_B_patch->triangles,
      pressure_sol_B);


  std::cout
      << "Results by solving h0 + Cp with LCP"
      << std::endl;
  PRINT_VAR(max_uA);
  PRINT_VAR(max_uB)
  PRINT_VAR(total_force_A);
  PRINT_VAR(total_force_B);
  std::cout << std::endl;



#endif




#if 0
  // Test of phi = Hu

//  solved = moby_LCP_solver.SolveLcpLemke(
//      C_A, h0_gap_A, &pressure_sol_A);
//  DRAKE_DEMAND(solved);
//  solved = moby_LCP_solver.SolveLcpLemke(
//      C_B, h0_gap_B, &pressure_sol_B);
//  DRAKE_DEMAND(solved);

  u_A_test = C_A * pressure_sol_A;
  u_B_test = C_B * pressure_sol_B;

  VectorX<double> u_test(num_nodes_A + num_nodes_B);
  u_test.head(num_nodes_A) = u_A_test;
  u_test.tail(num_nodes_B) = u_B_test;

  VectorX<double> phi_deformation = boussinesq_results_by_force->H * (-u_test);
  VectorX<double> phi = boussinesq_results_by_force->phi0 + phi_deformation;

  std::ofstream phi_file("phi_distance.txt");
  phi_file << "phi, phi0, phi_deformation, ratio phi / phi0" << std::endl;
  for (int i = 0; i < phi.rows(); i++) {
    phi_file << fmt::format("{:.8f}, {:.8f}, {:.8f}, {:.8f}\n",
                            phi(i),
                            boussinesq_results_by_force->phi0(i),
                            phi_deformation(i),
    fabs(phi(i) / boussinesq_results_by_force->phi0(i)));
  }
  phi_file.close();


#endif

  PRINT_VAR(force_Hertz);


  return 0;
}


//void SolveCaTimesF(const std::vector<Vector3<double>>& object_A_points_G,
//                   const std::vector<Vector3<double>>& object_B_points_G,
//                   const std::vector<Vector3<int>>& object_A_triangles,
//                   const std::vector<Vector3<int>>& object_B_triangles,
//                   const VectorX<double>& object_A_node_areas,
//                   const VectorX<double>& object_B_node_areas,
//                   const Eigen::Isometry3d& X_WA,
//                   const Eigen::Isometry3d& X_WB,
//                   double young_modulus_star_A,
//                   double young_modulus_star_B,
//                   double radius = 1.0,
//                   double z_WSo = 0.9) {
//
//
//  // Use LCP to solve the deformation and the pressure, and force
//  const MatrixX<double> C_A = CalcComplianceMatrix(
//      object_A_points_G,
//      object_A_triangles,
//      1 / ( young_modulus_star_A * M_PI));
//
//  const MatrixX<double> C_B = CalcComplianceMatrix(
//      object_B_points_G,
//      object_B_triangles,
//      1 / (young_modulus_star_B * M_PI));
//
//  int num_nodes_A = object_A_points_G.size();
//  int num_nodes_B = object_B_points_G.size();
//  int num_nodes = num_nodes_A + num_nodes_B;
//
//  VectorX<double> area_matrix = VectorX<double>::Zero(num_nodes);
//  VectorX<double> area_matrix_inv = VectorX<double>::Zero(num_nodes);
//
//  for(int i = 0; i < num_nodes_A; i++) {
//    area_matrix(i) = object_A_node_areas(i);
//    area_matrix_inv(i) = 1.0 / object_A_node_areas(i);
//  }
//  for(int i = 0; i < num_nodes_B; i++) {
//    area_matrix(num_nodes_A + i) =
//        object_B_node_areas(i);
//    area_matrix_inv(num_nodes_A + i) =
//        object_B_node_areas(i);
//  }
//
//  const VectorX<double> area_A
//      = area_matrix.head(num_nodes_A);
//  const VectorX<double> area_B
//      = area_matrix.tail(num_nodes_B);
//
//  const VectorX<double> area_A_inv
//      = area_matrix_inv.head(num_nodes_A);
//  const VectorX<double> area_B_inv
//      = area_matrix_inv.tail(num_nodes_B);
//
//  const MatrixX<double> Ca_A = C_A * area_A_inv.asDiagonal();
//  const MatrixX<double> Ca_B = C_B * area_B_inv.asDiagonal();
//
//
//
//  VectorX<double> h0_gap_A(num_nodes_A);
//  for (int i_gap = 0; i_gap < num_nodes_A; i_gap++) {
//    const Vector3d pos = X_WA * object_A_points_G.at(i_gap);
//    h0_gap_A(i_gap) = pos(2);
//  }
//
//  VectorX<double> h0_gap_B(num_nodes_B);
//  for (int i_gap = 0; i_gap < num_nodes_B; i_gap++) {
//    const Vector3d pos = X_WB * object_B_points_G.at(i_gap);
//    double r = sqrt(pow(pos(0), 2) + pow(pos(1), 2));
//    h0_gap_B(i_gap) = z_WSo - sqrt(pow(radius, 2) - std::min(pow(r, 2), pow(radius, 2)));
//  }
//
//
//  // Calculate the h0 based on Hook's law
//  double young_modulus_star = 1.0 /
//      (1.0 / young_modulus_star_A + 1.0 / young_modulus_star_B);
//
//  double ratio_A = (1/young_modulus_star_A) / (1/young_modulus_star);
//  double ratio_B = (1/young_modulus_star_B) / (1/young_modulus_star);
//
////  (void) ratio_A;
////  (void) ratio_B;
//
//  VectorX<double> h0_gap_A_ratio = h0_gap_A * ratio_A;
//  VectorX<double> h0_gap_B_ratio = h0_gap_B * ratio_B;
//
//
//
//  //Result by solving h0 + Ca * f with LCP, with h0 defined by Hook's law
//
//  solvers::MobyLCPSolver<double> moby_LCP_solver;
//  VectorX<double> force_sol_A(num_nodes_A);
//  VectorX<double> force_sol_B(num_nodes_B);
//
//
//  bool solved = moby_LCP_solver.SolveLcpLemke(
//      Ca_A, h0_gap_A_ratio, &force_sol_A);
//  DRAKE_DEMAND(solved);
//  solved = moby_LCP_solver.SolveLcpLemke(
//      Ca_B, h0_gap_B_ratio, &force_sol_B);
//  DRAKE_DEMAND(solved);
//
//
//  double force_A = force_sol_A.sum();
//  double force_B = force_sol_B.sum();
//
//  std::cout
//      << "Result by solving h0 + Ca * f with LCP, with h0 defined by Hook's law"
//         << std::endl;
//  PRINT_VAR(force_A);
//  PRINT_VAR(force_B);
//  std::cout << std::endl;
//
//
//
//  //Results by solving h0 + Cp with LCP, with h0 defined by Hook's law
//
//  VectorX<double> pressure_sol_A(num_nodes_A);
//  VectorX<double> pressure_sol_B(num_nodes_B);
//
//  solved = moby_LCP_solver.SolveLcpLemke(
//      C_A, h0_gap_A_ratio, &pressure_sol_A);
//  DRAKE_DEMAND(solved);
//  solved = moby_LCP_solver.SolveLcpLemke(
//      C_B, h0_gap_B_ratio, &pressure_sol_B);
//  DRAKE_DEMAND(solved);
//
//  VectorX<double> u_A_test(num_nodes_A);
//  VectorX<double> u_B_test(num_nodes_B);
//
//
//  u_A_test = C_A * pressure_sol_A;
//  u_B_test = C_B * pressure_sol_B;
//
//  double max_uA = 0;
//  for (int i = 0; i < num_nodes_A; i++) {
//    if (u_A_test(i) > max_uA) max_uA = u_A_test(i);
//  }
//  double max_uB = 0;
//  for (int i = 0; i < num_nodes_B; i++) {
//    if (u_B_test(i) > max_uB) max_uB = u_B_test(i);
//  }
//
//
//  double total_force_A = CalcForceOverMesh(
//      object_A_points_G,
//      object_A_triangles,
//      pressure_sol_A);
//
//  double total_force_B = CalcForceOverMesh(
//      object_B_points_G,
//      object_B_triangles,
//      pressure_sol_B);
//
//
//  std::cout
//      << "Results by solving h0 + Cp with LCP"
//      << std::endl;
//  PRINT_VAR(max_uA);
//  PRINT_VAR(max_uB)
//  PRINT_VAR(total_force_A);
//  PRINT_VAR(total_force_B);
//  std::cout << std::endl;
//
//
//
//}



GTEST_TEST(ExampleTest, ComplianceWithAreaInv) {
//  const bool flip_normals = true;

//  // Load mesh for a sphere.
//  std::unique_ptr<Mesh<double>> sphere = LoadMeshFromObj(
//      "drake/multibody/boussinesq_solver/test/Mesh_3/sphere.obj", flip_normals);
//  sphere->mesh_index = 0;
//
//  std::unique_ptr<Mesh<double>> plane = LoadMeshFromObj(
//      "drake/multibody/boussinesq_solver/test/Mesh_3/plane.obj", flip_normals);
//  plane->mesh_index = 1;
//
//  const double young_modulus_star_A = 10000000.0;
//  const double young_modulus_star_B = 1.0;
//  double z_WSo = 0.9;
//  double sigma = 0.1;
//
//  Isometry3d X_WSphere{Translation3d{Vector3d(0, 0, z_WSo)}};
//  X_WSphere.linear() = MatrixX<double>::Identity(3, 3);
//
//  // The top of the plane is placed at z = 0
//  Isometry3d X_WPlane{Translation3d{Vector3d(0, 0, 0.0)}};
//  X_WPlane.linear() = MatrixX<double>::Identity(3, 3);

//  auto owned_results =
//      std::make_unique<std::vector<PenetrationAsTrianglePair<double>>>();
//  std::vector<PenetrationAsTrianglePair<double>>& results = *owned_results;
//
//
//  results =
//      geometry::mesh_query::MeshToMeshQuery(X_WSphere, *sphere, X_WPlane, *plane, sigma);
//
//  auto patches =
//      geometry::mesh_query::MakeLocalPatchMeshes(&results, *sphere, *plane);
//  std::unique_ptr<Mesh<double>> object_A_patch = std::move(patches.first);
//  std::unique_ptr<Mesh<double>> object_B_patch = std::move(patches.second);
//
//  SolveCaTimesF(object_A_patch->points_G, object_B_patch->points_G,
//                object_A_patch->triangles, object_B_patch->triangles,
//                object_A_patch->node_areas, object_B_patch->node_areas,
//                X_WSphere, X_WPlane,
//                young_modulus_star_A, young_modulus_star_B);

}




GTEST_TEST(ExampleTest, SpherePlane) {


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





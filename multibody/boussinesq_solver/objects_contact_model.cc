#include "drake/multibody/boussinesq_solver/objects_contact_model.h"

#include <memory>
#include <vector>
#include <iostream>
#include <string>
#include <chrono>
#include <utility>

#include "drake/multibody/boussinesq_solver/math_helper.h"
#include "drake/multibody/boussinesq_solver/jacobian_H_matrix.h"
#include "drake/multibody/boussinesq_solver/compliance_matrix.h"
#include "drake/solvers/moby_lcp_solver.h"

#define PRINT_VAR(a) std::cout << #a": " << a << std::endl;

namespace drake {
namespace multibody {
namespace boussinesq_solver {

using Eigen::AngleAxisd;
using Eigen::Translation3d;
using Eigen::Isometry3d;
using Eigen::Vector3d;

using geometry::mesh_query::Mesh;
using geometry::PenetrationAsTrianglePair;
using geometry::mesh_query::MakeLocalPatchMeshes;
using clock = std::chrono::steady_clock;

// B is the one at the bottom
std::unique_ptr<BoussinesqContactModelResults<double>>
CalcContactSpatialForceBetweenMeshes(
    const Mesh<double>& object_A,
    const Isometry3d& X_WA,
    const double young_modulus_star_A,
    const Mesh<double>& object_B,
    const Isometry3d& X_WB,
    const double young_modulus_star_B,
    double sigma) {
  auto owned_results =
      std::make_unique<std::vector<PenetrationAsTrianglePair<double>>>();
  std::vector<PenetrationAsTrianglePair<double>>& results = *owned_results;

  const clock::time_point start_query = clock::now();

  results = MeshToMeshQuery(X_WA, object_A, X_WB, object_B, sigma);

  const clock::time_point end_query = clock::now();
  double querying_time = std::chrono::duration<double>(end_query
      - start_query).count();
  PRINT_VAR(querying_time);

  // This call creates the two patches on each mesh and updates "results" so
  // that the triangle indexes in each pair are "local indexes" to the patch
  // meshes.
  auto patches = MakeLocalPatchMeshes(&results, object_A, object_B);
  std::unique_ptr<Mesh<double>> object_A_patch = std::move(patches.first);
  std::unique_ptr<Mesh<double>> object_B_patch = std::move(patches.second);

  const int num_phis = results.size();
  const int num_nodes_A = object_A_patch->points_G.size();
  const int num_nodes_B = object_B_patch->points_G.size();
  const int num_nodes = num_nodes_A + num_nodes_B;

  const MatrixX<double> H = CalcJacobianHMatrix(
      results, object_A_patch->points_G, object_B_patch->points_G,
      object_A_patch->mesh_index, object_B_patch->mesh_index,
      young_modulus_star_A, young_modulus_star_B);

  DRAKE_DEMAND(H.rows() == num_phis);
  DRAKE_DEMAND(H.cols() == num_nodes);

  // Calculate combined Young Modulus

  double young_modulus_star = 1.0 /
      (1.0 / young_modulus_star_A + 1.0 / young_modulus_star_B);

  (void) young_modulus_star;

  VectorX<double> area_matrix = VectorX<double>::Zero(num_nodes);
  VectorX<double> area_matrix_inv = VectorX<double>::Zero(num_nodes);

  for(int i = 0; i < num_nodes_A; i++) {
    area_matrix(i) = object_A_patch->node_areas(i);
    area_matrix_inv(i) = 1.0 / object_A_patch->node_areas(i);
  }
  for(int i = 0; i < num_nodes_B; i++) {
    area_matrix(num_nodes_A + i) = object_B_patch->node_areas(i);
    area_matrix_inv(num_nodes_A + i) = 1.0 / object_B_patch->node_areas(i);
  }


  const clock::time_point start_Compliance = clock::now();

  const MatrixX<double> C_A = CalcComplianceMatrix(
      object_A_patch->points_G, object_A_patch->triangles,
      1.0 / (young_modulus_star_A * M_PI));

  const MatrixX<double> C_B = CalcComplianceMatrix(
      object_B_patch->points_G, object_B_patch->triangles,
      1.0 / (young_modulus_star_B * M_PI));

  MatrixX<double> C = MatrixX<double>::Zero(num_nodes, num_nodes);
  C.block(0, 0, num_nodes_A, num_nodes_A) = C_A;
  C.block(num_nodes_A, num_nodes_A, num_nodes_B, num_nodes_B) = C_B;


  const clock::time_point end_Compliance = clock::now();
  double compliance_time
      = std::chrono::duration<double>(end_Compliance - start_Compliance).count();
  PRINT_VAR(compliance_time);


  // TODO: (mengyao) Tobe changed back
  MatrixX<double> H_trans = H.transpose();

  // Delassus operator.
  //  const MatrixX<double> W = H * C * H_trans;
  MatrixX<double> C_times_area_inv = C * area_matrix_inv.asDiagonal();
  DRAKE_ASSERT(C_times_area_inv.rows() == num_nodes);
  DRAKE_ASSERT(C_times_area_inv.cols() == num_nodes);
  MatrixX<double> W = H * C_times_area_inv * H_trans;

  VectorX<double> phi0(num_phis);
  for (int i = 0; i < num_phis; ++i) {
    phi0(i) = results[i].signed_distance;
  }
  VectorX<double> kkt_multipliers(num_phis);

  // Solve QP problem as an LCP (expensive solution).
  const clock::time_point start_LCP = clock::now();
  solvers::MobyLCPSolver<double> moby_LCP_solver;
  const bool solved = moby_LCP_solver.SolveLcpLemke(W, phi0, &kkt_multipliers);
  DRAKE_DEMAND(solved);
  const clock::time_point end_LCP = clock::now();
  double LCP_time = std::chrono::duration<double>(end_LCP - start_LCP).count();
  PRINT_VAR(LCP_time);

  VectorX<double> phi_deformation = H * C_times_area_inv * H_trans * kkt_multipliers;
  VectorX<double> pressure = area_matrix_inv.asDiagonal() * H_trans * kkt_multipliers;
  VectorX<double> u_deformation = C * pressure;


  SpatialForce<double> F_Ao_W;
  F_Ao_W.SetZero();
  SpatialForce<double> F_Bo_W;
  F_Bo_W.SetZero();

  // WEIGHTED NORMALS FORMULATION.
#if 0
  // Computation of forces directly from multipliers
  auto sign = [](double x) {
    // std::copysign(x, y):
    //   Returns a value with the magnitude of x and the sign of y.
    return std::copysign(1.0, x);
  };
  for (int k = 0; k < num_phis; ++k) {
    const auto& query = results[k];

    // Undeformed surface normals.
    const auto& nA0_W = query.normal_A_W;
    const auto& nB0_W = query.normal_B_W;

    // Approximation to the deformed surface normal on A.
    const Vector3<double> nA_W =
        (young_modulus_star_B * nA0_W -
         young_modulus_star_A * nB0_W).normalized();

    Vector3<double> n_AtoB_W = (query.p_WoBs_W - query.p_WoAs_W).normalized();
    const Vector3<double> that = -sign(query.signed_distance) * n_AtoB_W;

    F_Ao_W.translational() += that.dot(nA_W) * kkt_multipliers(k) * nA_W;
  }

  // Since nB_W = -nA_W
  F_Bo_W.translational() += -F_Ao_W.translational();
#endif

  // ORIGINAL FORMULATION using the normals on the undeformed surfaces.
  // When kkt multipliers are forces:
  VectorX<double> force_on_nodes = H_trans * kkt_multipliers;

  for (int i = 0; i < num_nodes_A; ++i) {
    F_Ao_W.translational() +=
        force_on_nodes(i) * object_A_patch->node_normals_G[i];
  }

  for (int i = 0; i < num_nodes_B; ++i) {
    F_Bo_W.translational() +=
        force_on_nodes(num_nodes_A + i) * object_B_patch->node_normals_G[i];
  }

  // I am done using the patches.
  std::unique_ptr<BoussinesqContactModelResults<double>> boussinesq_results =
      std::make_unique<BoussinesqContactModelResults<double>>();

  boussinesq_results->F_Ao_W = F_Ao_W;
  boussinesq_results->F_Bo_W = F_Bo_W;
  boussinesq_results->object_A_patch = std::move(object_A_patch);  // now object_A_patch is nullptr.
  boussinesq_results->object_B_patch = std::move(object_B_patch);

  boussinesq_results->pressure_patch_A =
      pressure.segment(0, num_nodes_A);
  boussinesq_results->pressure_patch_B =
      pressure.segment(num_nodes_A, num_nodes_B);
  boussinesq_results->deformation_patch_A =
      u_deformation.segment(0, num_nodes_A);
  boussinesq_results->deformation_patch_B =
      u_deformation.segment(num_nodes_A, num_nodes_B);

  boussinesq_results->contact_results = std::move(owned_results);
  boussinesq_results->phi0 = phi0;
  boussinesq_results->phi_u = phi_deformation;
  boussinesq_results->H = H;
  boussinesq_results->kkt_multipliers = kkt_multipliers;

  return std::move(boussinesq_results);
}









// The function by doing pressure rather than force (solve phi = phi0 + HCH' pi)
// B is the one at the bottom
std::unique_ptr<BoussinesqContactModelResults<double>>
CalcContactSpatialForceBetweenMeshesByPressure(
    const Mesh<double>& object_A,
    const Isometry3d& X_WA,
    const double young_modulus_star_A,
    const Mesh<double>& object_B,
    const Isometry3d& X_WB,
    const double young_modulus_star_B,
    double sigma,
    bool press_in) {
  std::cout << "Results by solving pressure" << std::endl;

  auto owned_results =
      std::make_unique<std::vector<PenetrationAsTrianglePair<double>>>();
  std::vector<PenetrationAsTrianglePair<double>>& results = *owned_results;

  const clock::time_point start_query = clock::now();

  results = MeshToMeshQuery(X_WA, object_A, X_WB, object_B, sigma);

  const clock::time_point end_query = clock::now();
  double querying_time = std::chrono::duration<double>(end_query
                                                           - start_query).count();
  PRINT_VAR(querying_time);

  // This call creates the two patches on each mesh and updates "results" so
  // that the triangle indexes in each pair are "local indexes" to the patch
  // meshes.
  auto patches = MakeLocalPatchMeshes(&results, object_A, object_B);
  std::unique_ptr<Mesh<double>> object_A_patch = std::move(patches.first);
  std::unique_ptr<Mesh<double>> object_B_patch = std::move(patches.second);

  const int num_phis = results.size();
  const int num_nodes_A = object_A_patch->points_G.size();
  const int num_nodes_B = object_B_patch->points_G.size();
  const int num_nodes = num_nodes_A + num_nodes_B;


  const MatrixX<double> H = CalcJacobianHMatrix(
      results, object_A_patch->points_G, object_B_patch->points_G,
      object_A_patch->mesh_index, object_B_patch->mesh_index,
      young_modulus_star_A, young_modulus_star_B);


  PRINT_VAR(object_A_patch->points_G.size());
  PRINT_VAR(object_B_patch->points_G.size());


  DRAKE_DEMAND(H.rows() == num_phis);
  DRAKE_DEMAND(H.cols() == num_nodes);

  // Calculate combined Young Modulus

  double young_modulus_star = 1 /
      (1 / young_modulus_star_A + 1 / young_modulus_star_B);

  (void) young_modulus_star;



//  const MatrixX<double> C_A = CalcComplianceMatrix(
//      object_A_patch->points_G, object_A_patch->triangles,
//      1 / (young_modulus_star * M_PI));
//
//  const MatrixX<double> C_B = CalcComplianceMatrix(
//      object_B_patch->points_G, object_B_patch->triangles,
//      1 / (young_modulus_star * M_PI));




  MatrixX<double> area_matrix = MatrixX<double>::Zero(num_nodes, num_nodes);

  for(int i = 0; i < num_nodes_A; i++) {
    area_matrix(i, i) = object_A_patch->node_areas(i);
  }
  for(int i = 0; i < num_nodes_B; i++) {
    area_matrix(num_nodes_A + i, num_nodes_A + i) = object_B_patch->node_areas(i);
  }


  const clock::time_point start_Compliance = clock::now();

  const MatrixX<double> C_A = CalcComplianceMatrix(
      object_A_patch->points_G, object_A_patch->triangles,
      1 / (young_modulus_star_A * M_PI));

  const MatrixX<double> C_B = CalcComplianceMatrix(
      object_B_patch->points_G, object_B_patch->triangles,
      1 / (young_modulus_star_B * M_PI));

  MatrixX<double> C = MatrixX<double>::Zero(num_nodes, num_nodes);
  C.block(0, 0, num_nodes_A, num_nodes_A) = C_A;
  C.block(num_nodes_A, num_nodes_A, num_nodes_B, num_nodes_B) = C_B;


  const clock::time_point end_Compliance = clock::now();
  double compliance_time
      = std::chrono::duration<double>(end_Compliance - start_Compliance).count();
  PRINT_VAR(compliance_time);




  // TODO: (mengyao) Tobe changed back
  MatrixX<double> H_trans = H.transpose();

  // Delassus operator.
  const MatrixX<double> W = H * C * H_trans;



  VectorX<double> phi0(num_phis);
  for (int i = 0; i < num_phis; ++i) {
    phi0(i) = results[i].signed_distance;
  }

  VectorX<double> kkt_multipliers(num_phis);


  const clock::time_point start_LCP = clock::now();


  solvers::MobyLCPSolver<double> moby_LCP_solver;
  const bool solved = moby_LCP_solver.SolveLcpLemke(W, phi0, &kkt_multipliers);
  DRAKE_DEMAND(solved);

  const clock::time_point end_LCP = clock::now();
  double LCP_time = std::chrono::duration<double>(end_LCP - start_LCP).count();
  PRINT_VAR(LCP_time);



  // Estimate pressure from KKT multipliers.
  // Integrate pressure to get the spatial force.
  // CODE FOR INTEGRATING FORCES HERE.


  VectorX<double> phi_deformation = H * C * H_trans * kkt_multipliers;


// TODO: (mengyao) modify it later

  VectorX<double> pressure = H_trans * kkt_multipliers;

  VectorX<double> u_deformation = C * pressure;


  VectorX<double> phi = phi0 + phi_deformation;
  std::ofstream phi_file("phi_distance_object_contact.txt");
  phi_file << "phi, phi0, phi_deformation, ratio phi / phi0" << std::endl;
  for (int i = 0; i < phi.rows(); i++) {
    phi_file << fmt::format("{:.8f}, {:.8f}, {:.8f}, {:.8f}\n",
                            phi(i),
                            phi0(i),
                            phi_deformation(i),
                            fabs(phi(i) / phi0(i)));
  }
  phi_file.close();




//  for (int i = 0; i < pressure.rows(); i++) {
//    if(pressure(i) > 0.0) {
//      pressure(i) = 0.0;
//    }
//  }



//  PRINT_VAR(phi_deformation.rows());
//  PRINT_VAR((C * phi_deformation).rows());
//  PRINT_VAR(C.rows());
//  PRINT_VAR(C.cols());


  SpatialForce<double> F_Ao_W;
  F_Ao_W.SetZero();
  SpatialForce<double> F_Bo_W;
  F_Bo_W.SetZero();


  (void) press_in;

  int num_triangles_A = object_A_patch->triangles.size();
  int num_triangles_B = object_B_patch->triangles.size();

  // TODO: If works, should refine the robustness with the "press_in"
  for (int i = 0; i < num_triangles_A; ++i) {

    Vector3<int>& triangle = object_A_patch->triangles[i];
    const Vector3d& p1 = object_A_patch->points_G[triangle(0)];
    const Vector3d& p2 = object_A_patch->points_G[triangle(1)];
    const Vector3d& p3 = object_A_patch->points_G[triangle(2)];

    double sign = 1.0;
    double pressure1 = pressure(triangle(0)) * sign;
    double pressure2 = pressure(triangle(1)) * sign;
    double pressure3 = pressure(triangle(2)) * sign;

    const double tol_limits = std::numeric_limits<double>:: epsilon();
    (void) tol_limits;
    Vector3d area_vec =
        CalcTriangleArea(p1, p2, p3);

    if ((pressure1 <= 0.0) &&
        (pressure2 <= 0.0) &&
        (pressure3 <= 0.0)) {
      F_Ao_W.translational() +=
          (pressure1 + pressure2 + pressure3) *
                  (object_A_patch->face_normals_G[i])
              * area_vec.norm() / 3.0;

    }
    else if ((pressure1 >= 0.0) &&
        (pressure2 >= 0.0) &&
        (pressure3 >= 0.0))
    {
      continue;
    }
    // If two nodes are with negtive values:
    else if ((((pressure1 >= 0.0) && (pressure2 >= 0.0))) ||
        (((pressure2 >= 0.0) && (pressure3 >= 0.0))) ||
        (((pressure1 >= 0.0) && (pressure3 >= 0.0)))) {
      if ((pressure1 >= 0.0) && (pressure2 >= 0.0)) {
        F_Ao_W.translational() +=
            (object_A_patch->node_normals_G[triangle(2)]) *
            CalcInterpolatedVolumeExcludingTwoNode(
            p1, pressure1, p2, pressure2, p3, pressure3);
      }
      else if ((pressure2 >= 0.0) && (pressure3 >= 0.0)) {
        F_Ao_W.translational() +=
            (object_A_patch->node_normals_G[triangle(0)]) *
            CalcInterpolatedVolumeExcludingTwoNode(
                p2, pressure2, p3, pressure3, p1, pressure1);
      }
      else if ((pressure1 >= 0.0) && (pressure3 >= 0.0)) {
        F_Ao_W.translational() +=
            (object_A_patch->node_normals_G[triangle(1)]) *
            CalcInterpolatedVolumeExcludingTwoNode(
                p3, pressure3, p1, pressure1, p2, pressure2);
      }
    }
    // If only one of the nodes is negtive
    else {
      if (pressure1 >= 0.0) {
        F_Ao_W.translational() +=
            (object_A_patch->node_normals_G[triangle(1)] +
                object_A_patch->node_normals_G[triangle(2)]) / 2.0 *
            CalcInterpolatedVolumeExcludingOneNode(
                p2, pressure2, p3, pressure3, p1, pressure1);
      }
      else if (pressure2 >= 0.0) {
        F_Ao_W.translational() +=
            (object_A_patch->node_normals_G[triangle(0)] +
                object_A_patch->node_normals_G[triangle(2)]) / 2.0 *
            CalcInterpolatedVolumeExcludingOneNode(
                p3, pressure3, p1, pressure1, p2, pressure2);
      }
      else {
        F_Ao_W.translational() +=
            (object_A_patch->node_normals_G[triangle(0)] +
                object_A_patch->node_normals_G[triangle(1)]) / 2.0 *
            CalcInterpolatedVolumeExcludingOneNode(
                p1, pressure1, p2, pressure2, p3, pressure3);
      }
    }
  }

  for (int i = 0; i < num_triangles_B; ++i) {

    Vector3<int>& triangle = object_B_patch->triangles[i];
    const Vector3d& p1 = object_B_patch->points_G[triangle(0)];
    const Vector3d& p2 = object_B_patch->points_G[triangle(1)];
    const Vector3d& p3 = object_B_patch->points_G[triangle(2)];

    double sign = 1.0;
    double pressure1 = pressure(num_nodes_A + triangle(0)) * sign;
    double pressure2 = pressure(num_nodes_A + triangle(1)) * sign;
    double pressure3 = pressure(num_nodes_A + triangle(2)) * sign;

    const double tol_limits = std::numeric_limits<double>:: epsilon();
    (void) tol_limits;
    Vector3d area_vec =
        CalcTriangleArea(p1, p2, p3);

    if ((pressure1 <= 0.0) &&
        (pressure2 <= 0.0) &&
        (pressure3 <= 0.0)) {
      F_Bo_W.translational() +=
          (pressure1 + pressure2 + pressure3) *
              (object_B_patch->face_normals_G[i])
              * area_vec.norm() / 3.0;
    }
    else if ((pressure1 >= 0.0) &&
        (pressure2 >= 0.0) &&
        (pressure3 >= 0.0))
    {
      continue;
    }
      // If two nodes are with negtive values:
    else if ((((pressure1 >= 0.0) && (pressure2 >= 0.0))) ||
        (((pressure2 >= 0.0) && (pressure3 >= 0.0))) ||
        (((pressure1 >= 0.0) && (pressure3 >= 0.0)))) {
      if ((pressure1 >= 0.0) && (pressure2 >= 0.0)) {
        F_Bo_W.translational() +=
            (object_B_patch->node_normals_G[triangle(2)]) *
                CalcInterpolatedVolumeExcludingTwoNode(
                    p1, pressure1, p2, pressure2, p3, pressure3);
      }
      else if ((pressure2 >= 0.0) && (pressure3 >= 0.0)) {
        F_Bo_W.translational() +=
            (object_B_patch->node_normals_G[triangle(0)]) *
                CalcInterpolatedVolumeExcludingTwoNode(
                    p2, pressure2, p3, pressure3, p1, pressure1);
      }
      else if ((pressure1 >= 0.0) && (pressure3 >= 0.0)) {
        F_Bo_W.translational() +=
            (object_B_patch->node_normals_G[triangle(1)]) *
                CalcInterpolatedVolumeExcludingTwoNode(
                    p3, pressure3, p1, pressure1, p2, pressure2);
      }
    }
      // If only one of the nodes is negtive
    else {
      if (pressure1 >= 0.0) {
        F_Bo_W.translational() +=
            (object_B_patch->node_normals_G[triangle(1)] +
                object_B_patch->node_normals_G[triangle(2)]) / 2.0 *
                CalcInterpolatedVolumeExcludingOneNode(
                    p2, pressure2, p3, pressure3, p1, pressure1);
      }
      else if (pressure2 >= 0.0) {
        F_Bo_W.translational() +=
            (object_B_patch->node_normals_G[triangle(0)] +
                object_B_patch->node_normals_G[triangle(2)]) / 2.0 *
                CalcInterpolatedVolumeExcludingOneNode(
                    p3, pressure3, p1, pressure1, p2, pressure2);
      }
      else {
        F_Bo_W.translational() +=
            (object_B_patch->node_normals_G[triangle(0)] +
                object_B_patch->node_normals_G[triangle(1)]) / 2.0 *
                CalcInterpolatedVolumeExcludingOneNode(
                    p1, pressure1, p2, pressure2, p3, pressure3);
      }
    }
  }




//  for (int i = 0; i < num_nodes_A; ++i) {
//
//    // TODO: (mengyao: to be clean-up)
//    double area = object_A_patch->node_areas[i];
//
//    F_Ao_W.translational() += pressure(i) * area
//        * object_A_patch->node_normals_G[i];
//
////    F_Ao_W.translational() += pressure(i) * area
////        * (-Vector3d::UnitZ());
//  }
//
//  for (int i = 0; i < num_nodes_B; ++i) {
//
//    // TODO: (mengyao: to be clean-up)
//    double area = object_B_patch->node_areas[i];
//    F_Bo_W.translational() += pressure(num_nodes_A + i) * area
//        * object_B_patch->node_normals_G[i];
//  }







  // I am done using the patches.
  std::unique_ptr<BoussinesqContactModelResults<double>> boussinesq_results =
      std::make_unique<BoussinesqContactModelResults<double>>();

  boussinesq_results->F_Ao_W = F_Ao_W;
  boussinesq_results->F_Bo_W = F_Bo_W;
  boussinesq_results->object_A_patch = std::move(object_A_patch);  // now object_A_patch is nullptr.
  boussinesq_results->object_B_patch = std::move(object_B_patch);

  boussinesq_results->pressure_patch_A =
      pressure.segment(0, num_nodes_A);
  boussinesq_results->pressure_patch_B =
      pressure.segment(num_nodes_A, num_nodes_B);
  boussinesq_results->deformation_patch_A =
      u_deformation.segment(0, num_nodes_A);
  boussinesq_results->deformation_patch_B =
      u_deformation.segment(num_nodes_A, num_nodes_B);

  boussinesq_results->contact_results = std::move(owned_results);
  boussinesq_results->phi0 = phi0;
  boussinesq_results->H = H;


  PRINT_VAR(F_Ao_W.translational());
  PRINT_VAR(F_Bo_W.translational());
  std::cout << std::endl;


  boussinesq_results->kkt_multipliers = kkt_multipliers;

  return std::move(boussinesq_results);
}























}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake





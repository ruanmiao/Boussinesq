#include "drake/multibody/boussinesq_solver/compliance_matrix.h"

#include <gtest/gtest.h>

#include "drake/common/eigen_types.h"
#include "drake/common/test_utilities/eigen_matrix_compare.h"
#include "drake/multibody/boussinesq_solver/test_helper.h"
#include "drake/solvers/moby_lcp_solver.h"

namespace drake {
namespace multibody {
namespace boussinesq_solver {
namespace {

using Eigen::Vector2d;
using Eigen::Vector3d;

/// The expected values for this test are the results by running Matlab. The
/// precision (15 digits) of the expected values if the same as the "long"
/// in Matlab
/// The "expected" results is solved numerically by Matlab LCP Solver with
/// its tolerance being 1e-8.
/// The tolerance here is set accordance to the Matlab function tolerance.
GTEST_TEST(ComplianceMatrixTest, ShpereContact) {
  const double r_sphere = 1.0;
  const double E_modulus = 1.0;
  const double indent_ratio = 0.2;
  const double z0_sphere = r_sphere * (1 - indent_ratio);
  const double r_contact = sqrt(r_sphere * r_sphere - z0_sphere * z0_sphere);
  const double r_mesh = r_contact * 1.1;
  const int elements_per_r = 8;

  std::pair<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3i>>
      mesh_data = MeshCircle(
          MatrixX<double>::Zero(2, 1), r_mesh, elements_per_r);

  const std::vector<Eigen::Vector3d>& points_in_mesh = mesh_data.first;
  const std::vector<Eigen::Vector3i>& triangles_in_mesh = mesh_data.second;

  int num_nodes = points_in_mesh.size();

  MatrixX<double> compliance = CalcComplianceMatrix(
      points_in_mesh, triangles_in_mesh, 1 / (E_modulus * M_PI));

  VectorX<double> h0_gap(num_nodes);
  for (int i_gap = 0; i_gap < num_nodes; i_gap++) {
    const Vector3d pos = points_in_mesh.at(i_gap);
    h0_gap(i_gap) =
        z0_sphere - sqrt(pow(r_sphere, 2) -
            std::min(pow(pos(0), 2) + pow(pos(1), 2), pow(r_sphere, 2)));
  }

  solvers::MobyLCPSolver<double> moby_LCP_solver;
  VectorX<double> pressure_sol(num_nodes);

  bool solved = moby_LCP_solver.SolveLcpLemke(
      compliance, h0_gap, &pressure_sol);
  OutputMeshToVTK(points_in_mesh, triangles_in_mesh, pressure_sol);

  double total_force = CalcForceOverMesh(
      points_in_mesh, triangles_in_mesh, pressure_sol);

  double expected_force = 0.117810975771820;
  EXPECT_NEAR(total_force, expected_force, 10 * 1e-8);
  EXPECT_GT(moby_LCP_solver.get_num_pivots(), 0);
  EXPECT_TRUE(solved);
}


}  // namespace
}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake

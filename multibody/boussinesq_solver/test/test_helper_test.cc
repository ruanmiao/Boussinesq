#include "drake/multibody/boussinesq_solver/test_helper.h"

#include <gtest/gtest.h>

#include "drake/common/eigen_types.h"
#include "drake/common/test_utilities/eigen_matrix_compare.h"

namespace drake {
namespace multibody {
namespace boussinesq_solver {
namespace {

using Eigen::Vector2d;

/// The expected values for this files are the results by running Matlab. The
/// precision (15 digits) of the expected values if the same as the "long"
/// in Matlab
GTEST_TEST(GetPressureIntegrandR, Square3x3) {
  const double length = 1;
  const Vector2d p1(0.0, 0.0);
  const Vector2d p2(length, 0.0);
  const Vector2d p3(0.0, length);
  const int nodes_per_edge = 3;

  const std::pair<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3i>>&
      mesh_data = MeshSquare(p1, p2, p3, nodes_per_edge, nodes_per_edge);

  VectorX<double> pressures = GetPressureIntegrandR(mesh_data.first);
  VectorX<double> expected(9);
  expected << 0.0, length / (nodes_per_edge - 1), length,
              length / (nodes_per_edge - 1),
              sqrt(2 * length / (nodes_per_edge - 1) * length /
                  (nodes_per_edge - 1)),
              sqrt(5 * length / (nodes_per_edge - 1) * length /
                  (nodes_per_edge - 1)),
              length,
              sqrt(5 * length / (nodes_per_edge - 1) * length /
                  (nodes_per_edge - 1)),
              sqrt(2 * length * length);
EXPECT_TRUE(
      CompareMatrices(
          pressures, expected,
          10 * std::numeric_limits<double>::epsilon(),
          MatrixCompareType::absolute));
}

GTEST_TEST(MeshCircle, NumPerR3) {
  const Vector2d origin(0.0, 0.0);
  const double radius = 1.0;
  int num_pr = 3;
  std::pair<std::vector<Eigen::Vector3d>,
            std::vector<Eigen::Vector3i>> mesh_data = MeshCircle(
                origin, radius, num_pr);
  const std::vector<Eigen::Vector3d>& points_in_mesh = mesh_data.first;
  const std::vector<Eigen::Vector3i>& triangles_in_mesh = mesh_data.second;

  const int num_nodes = points_in_mesh.size();
  VectorX<double> z_values = VectorX<double>::Zero(num_nodes);

  OutputMeshToVTK(points_in_mesh, triangles_in_mesh, z_values, false);
}

GTEST_TEST(MeshSphere, General1) {
  const Vector3<double> origin(0.0, 0.0, 0.8);
  const double radius = 1.0;
  const double half_sector = M_PI / 2;
  int num_pr = 4;
  const std::pair<std::vector<Eigen::Vector3d>,
            std::vector<Eigen::Vector3i>>& mesh_data = MeshSphere(
      origin, radius, num_pr, half_sector);

  const std::vector<Eigen::Vector3d>& points_in_mesh = mesh_data.first;
  const std::vector<Eigen::Vector3i>& triangles_in_mesh = mesh_data.second;

  const int num_nodes = points_in_mesh.size();
  VectorX<double> z_values(num_nodes);
  for (int i_node = 0; i_node < num_nodes; i_node++) {
    const VectorX<double>& pos = points_in_mesh[i_node];
    z_values[i_node] = pos(2);
  }
  OutputMeshToVTK(points_in_mesh, triangles_in_mesh, z_values);
}

}  // namespace
}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake

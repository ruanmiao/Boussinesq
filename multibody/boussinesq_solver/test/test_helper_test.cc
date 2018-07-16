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
  MatrixX<int> expected_M(24, 3);
  expected_M << 2, 3, 1, 3, 4, 1, 4, 5, 1,
                5, 6, 1, 6, 7, 1, 7, 2, 1,
                8, 9, 2, 9, 2, 3, 9, 10, 3,
                10, 11, 3, 11, 3, 4, 11, 12, 4,
                12, 13, 4, 13, 4, 5, 13, 14, 5,
                14, 15, 5, 15, 5, 6, 15, 16, 6,
                16, 17, 6, 17, 6, 7, 17, 18, 7,
                18, 19, 7, 19, 7, 2, 19, 8, 2;
  MatrixX<int> expected = expected_M - MatrixX<int>::Ones(24,3);
  std::vector<Eigen::Vector3i> tri_vectors = mesh_data.second;
  MatrixX<int> results(24, 3);
  for (int it = 0; it < 24; it++) {
    results.row(it) = tri_vectors[it];
  }
  EXPECT_TRUE(
      CompareMatrices(
          results, expected,
          10 * std::numeric_limits<double>::epsilon(),
          MatrixCompareType::absolute));
}

}  // namespace
}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake

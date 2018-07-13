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

}  // namespace
}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake

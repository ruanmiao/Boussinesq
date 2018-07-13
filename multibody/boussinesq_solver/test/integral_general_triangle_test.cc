#include "drake/multibody/boussinesq_solver/integral_general_triangle.h"

#include <gtest/gtest.h>

#include "drake/common/eigen_types.h"
#include "drake/common/test_utilities/eigen_matrix_compare.h"

namespace drake {
namespace multibody {
namespace boussinesq_solver {
namespace {

using Eigen::Vector2d;
using Eigen::Vector3d;

/// The expected values for this files are the results by running Matlab. The
/// precision (15 digits) of the expected values if the same as the "long"
/// in Matlab
GTEST_TEST(IntegralGeneralTriangleTest, Clockwise) {
  Vector2d p1, p2, p3;
  p1 << 0.0, 2.0;
  p2 << 1.0, 1.0;
  p3 << 0.0, 0.5;

  Vector3d res = CalGenralTriangleCompliance(p1, p2, p3);
  Vector3d expected_res;
  expected_res << 0.184823056895932, 0.207292259729693, 0.256305654879497;
  EXPECT_TRUE(
      CompareMatrices(
          res,
          expected_res,
          10 * std::numeric_limits<double>::epsilon(),
          MatrixCompareType::absolute));
}

GTEST_TEST(IntegralGeneralTriangleTest, EdgeAlignWithX) {
  Vector2d p1, p2, p3;
  p1 << 0.0, 2.0;
  p2 << 2.0, 0.0;
  p3 << 1.0, 0.0;

  Vector3d res = CalGenralTriangleCompliance(p1, p2, p3);
  Vector3d expected_res;
  expected_res << 0.251061663754629, 0.237592641787778, 0.282610891162499;
  EXPECT_TRUE(
      CompareMatrices(
          res,
          expected_res,
          10 * std::numeric_limits<double>::epsilon(),
          MatrixCompareType::absolute));
}

GTEST_TEST(IntegralGeneralTriangleTest, VertixAtOrigin) {
  Vector2d p1, p2, p3;
  p1 << 0.0, 0.0;
  p2 << -2.0, 0.0;
  p3 << 0.0, -1.0;

  Vector3d res = CalGenralTriangleCompliance(p1, p2, p3);
  Vector3d expected_res;
  expected_res << 0.860817881928008, 0.372163576385602, 0.488654305542406;
  EXPECT_TRUE(
      CompareMatrices(
          res,
          expected_res,
          10 * std::numeric_limits<double>::epsilon(),
          MatrixCompareType::absolute));
}

}  // namespace
}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake

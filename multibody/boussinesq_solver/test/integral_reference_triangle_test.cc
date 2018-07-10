#include "drake/multibody/boussinesq_solver/integral_reference_triangle.h"

#include <gtest/gtest.h>

#include "drake/common/eigen_types.h"
#include "drake/common/test_utilities/eigen_matrix_compare.h"

namespace drake {
namespace multibody {
namespace boussinesq_solver {
namespace {

using Eigen::Vector2d;
using Eigen::Vector3d;

GTEST_TEST(IntegralReferenceTriangleTest, CoLinear) {
  Vector2d p1, p2;
  p1 << 1.0, 0.0;
  p2 << 2.0, 0.0;

  Vector3d res = CalcIntegralReferenceTriangle(p1, p2);
  Vector3d expected_res;
  expected_res << 0.0, 0.0, 0.0;
  EXPECT_TRUE(
      CompareMatrices(
          res,
          expected_res,
          10 * std::numeric_limits<double>::epsilon(),
          MatrixCompareType::absolute));
}

GTEST_TEST(IntegralReferenceTriangleTest, Clockwise) {
  Vector2d p1, p2;
  p1 << -1.0, 1.0;
  p2 << 2.0, 0.0;

  Vector3d res = CalcIntegralReferenceTriangle(p1, p2);
  Vector3d expected_res;
  // Note: these expected values were previously computed with a prototype in
  // Matlab.
  expected_res << -0.382048107824599,  -0.560357885846892, -2.063121765365275;
  EXPECT_TRUE(
      CompareMatrices(
          res,
          expected_res,
          10 * std::numeric_limits<double>::epsilon(),
          MatrixCompareType::absolute));
}

GTEST_TEST(IntegralReferenceTriangleTest, LargeOffset) {
  Vector2d p1, p2;
  p1 << -1.0, -1.0;
  p2 << 2.0, 0.0;

  Vector3d res = CalcIntegralReferenceTriangle(p1, p2);
  Vector3d expected_res;
  // Note: these expected values were previously computed with a prototype in
  // Matlab.
  expected_res << 0.382048107824599,  -0.560357885846892, 2.063121765365274;
  EXPECT_TRUE(
      CompareMatrices(
      res,
      expected_res,
      10 * std::numeric_limits<double>::epsilon(),
      MatrixCompareType::absolute));
}

GTEST_TEST(IntegralReferenceTriangleTest, P1Origin) {
  Vector2d p1, p2;
  p1 << 0.0, 0.0;
  p2 << 2.0, 0.0;

  Vector3d res = CalcIntegralReferenceTriangle(p1, p2);
  Vector3d expected_res;
  // Note: these expected values were previously computed with a prototype in
  // Matlab.
  expected_res << 0.0, 0.0, 0.0;
  EXPECT_TRUE(
      CompareMatrices(
          res,
          expected_res,
          10 * std::numeric_limits<double>::epsilon(),
          MatrixCompareType::absolute));
}

GTEST_TEST(IntegralReferenceTriangleTest, AlignY) {
  Vector2d p1, p2;
  p1 << 0.0, 2.0;
  p2 << 0.0, 1.0;

  Vector3d res = CalcIntegralReferenceTriangle(p1, p2);
  Vector3d expected_res;
  expected_res << 0.0, 0.0, 0.0;
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

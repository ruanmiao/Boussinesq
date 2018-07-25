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

/// The expected values for this files are the results by running Matlab. The
/// precision (15 digits) of the expected values if the same as the "long"
/// in Matlab
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

GTEST_TEST(SplitedGeneralTriangleTest, ColinearOndEdgep1p2) {
  Vector2d p1, p2;
  const double r = 1.0;
  p1 << r / 2 * cos(M_PI * 2 / 3), r / 2 * sin(M_PI * 2 / 3);
  p2 << r * cos(M_PI * 2 / 3), r * sin(M_PI * 2 / 3);

  Vector3d res = CalcIntegralReferenceTriangle(p1, -p2);
  Vector3d expected_res;
  expected_res << 0.0, 0.0, 0.0;
  EXPECT_TRUE(
      CompareMatrices(
          res,
          expected_res,
          10 * std::numeric_limits<double>::epsilon(),
          MatrixCompareType::absolute));
}

/// The expected values for this files are the results by running Matlab. The
/// precision (6 digits) of the expected values if the same as the tol. of the
/// Integration function in Matlab
GTEST_TEST(IntegralReferenceTriangle3DTest, GeneralCase1) {
  Vector2d p1, p2;
  p1 << 2.0, 2.0;
  p2 << 0.5, 1.5;

  Vector3d res = CalcIntegralReferenceTriangle(p1, p2, 0.8);
  Vector3d expected_res;
  expected_res << 0.465063869455651, 0.673390142473203, 0.649887021833660;
  EXPECT_TRUE(
      CompareMatrices(
          res,
          expected_res, 10 * 1e-5,
          MatrixCompareType::absolute));
}

/// The expected values for this files are the results by running Matlab. The
/// precision (6 digits) of the expected values if the same as the tol. of the
/// Integration function in Matlab
GTEST_TEST(IntegralReferenceTriangle3DTest, General2) {
  Vector2d p1, p2;
  p1 << -2, 0;
  p2 << -1, 0;

  Vector3d res = CalcIntegralReferenceTriangle(p1, p2, 0.8);
  Vector3d expected_res;
  expected_res << 0.0, 0.0, 0.0;
  EXPECT_TRUE(
      CompareMatrices(
          res,
          expected_res, 10 * 1e-5,
          MatrixCompareType::absolute));
}


}  // namespace
}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake

#include "drake/multibody/boussinesq_solver/integral_general_triangle.h"

#include <gtest/gtest.h>

#include "drake/common/eigen_types.h"
#include "drake/common/test_utilities/eigen_matrix_compare.h"
#include "drake/multibody/boussinesq_solver/integral_reference_triangle.h"

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

  Vector3d res = CalcGeneralTriangleCompliance(p1, p2, p3);
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

  Vector3d res = CalcGeneralTriangleCompliance(p1, p2, p3);
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

  Vector3d res = CalcGeneralTriangleCompliance(p1, p2, p3);
  Vector3d expected_res;
  expected_res << 0.860817881928008, 0.372163576385602, 0.488654305542406;
  EXPECT_TRUE(
      CompareMatrices(
          res,
          expected_res,
          10 * std::numeric_limits<double>::epsilon(),
          MatrixCompareType::absolute));
}

GTEST_TEST(IntegralGeneralTriangleTest, ColinearOndEdge1) {
  Vector2d p1, p2, p3;
  const double r = 1.0;
  p1 << r / 2 * cos(M_PI * 2 / 3), r / 2 * sin(M_PI * 2 / 3);
  p2 << r * cos(M_PI * 2 / 3), r * sin(M_PI * 2 / 3);
  p3 << r * cos(M_PI * 5 / 6), r * sin(M_PI * 5 / 6);

  Vector3d res = CalcGeneralTriangleCompliance(p1, p2, p3);
  Vector3d expected_res;
  expected_res << 0.058157920941819, 0.049345918612694, 0.049622290224941;

  EXPECT_TRUE(
      CompareMatrices(
          res,
          expected_res,
          10 * std::numeric_limits<double>::epsilon(),
          MatrixCompareType::absolute));
}

GTEST_TEST(IntegralGeneralTriangleZNonZeroTest, ConstantPressure) {
  Vector2d p1, p2, p3;
  p1 << 5.0, 7.0;
  p2 << 7.0, 6.0;
  p3 << 5.5, 5.5;

  double zA = 0.8;
  Vector3d res = CalcGeneralTriangleCompliance(zA, p1, p2, p3);
  double res_interpolate = res(0) + res(1) + res(2);

  const Vector3<double>& I_12 = CalcIntegralReferenceTriangle(p1, p2, zA);
  const Vector3<double>& I_23 = CalcIntegralReferenceTriangle(p2, p3, zA);
  const Vector3<double>& I_31 = CalcIntegralReferenceTriangle(p3, p1, zA);

  double expected_res = fabs(I_12(2) + I_23(2) + I_31(2));

  EXPECT_NEAR(
      res_interpolate, expected_res,
      10 * std::numeric_limits<double>::epsilon());
}

GTEST_TEST(IntegralGeneralTriangleZNonZeroTest, LinearPressure) {
  Vector2d p1, p2, p3;
  p1 << 5.0, 7.0;
  p2 << 7.0, 6.0;
  p3 << 5.5, 5.5;

  double zA = 0.8;
  Vector3d res = CalcGeneralTriangleCompliance(zA, p1, p2, p3);
  double res_interpolate = res(0) * p1(0) + res(1) * p2(0) + res(2) * p3(0);

  const Vector3<double>& I_12 = CalcIntegralReferenceTriangle(p1, p2, zA);
  const Vector3<double>& I_23 = CalcIntegralReferenceTriangle(p2, p3, zA);
  const Vector3<double>& I_31 = CalcIntegralReferenceTriangle(p3, p1, zA);

  double expected_res = fabs(I_12(0) + I_23(0) + I_31(0));

  EXPECT_NEAR(
      res_interpolate, expected_res,
      10 * std::numeric_limits<double>::epsilon());
}

GTEST_TEST(IntegralGeneralTriangleZNonZeroTest, OnePointAtOrigin) {
  Vector2d p1, p2, p3;
  p1 << 0.0, 0.0;
  p2 << 7.0, 6.0;
  p3 << 5.5, 5.5;

  double zA = 0.8;
  Vector3d res = CalcGeneralTriangleCompliance(zA, p1, p2, p3);
  double res_interpolate = res(0) + res(1) + res(2);

  const Vector3<double>& I_12 = CalcIntegralReferenceTriangle(p1, p2, zA);
  const Vector3<double>& I_23 = CalcIntegralReferenceTriangle(p2, p3, zA);
  const Vector3<double>& I_31 = CalcIntegralReferenceTriangle(p3, p1, zA);

  double expected_res = fabs(I_12(2) + I_23(2) + I_31(2));

  EXPECT_NEAR(
      res_interpolate, expected_res,
      10 * std::numeric_limits<double>::epsilon());
}

GTEST_TEST(IntegralGeneralTriangleZNonZeroTest, OneEdgeAlignedWithX) {
  Vector2d p1, p2, p3;
  p1 << 0.0, 0.0;
  p2 << 7.0, 0.0;
  p3 << 5.5, 0.5;

  double zA = 0.8;
  Vector3d res = CalcGeneralTriangleCompliance(zA, p1, p2, p3);
  double res_interpolate = res(0) * p1(0) + res(1) * p2(0) + res(2) * p3(0);

  const Vector3<double>& I_12 = CalcIntegralReferenceTriangle(p1, p2, zA);
  const Vector3<double>& I_23 = CalcIntegralReferenceTriangle(p2, p3, zA);
  const Vector3<double>& I_31 = CalcIntegralReferenceTriangle(p3, p1, zA);

  double expected_res = fabs(I_12(0) + I_23(0) + I_31(0));

  EXPECT_NEAR(
      res_interpolate, expected_res,
      10 * std::numeric_limits<double>::epsilon());
}

GTEST_TEST(IntegralGeneralTriangleZNonZeroTest, OneEdgeAlignedWithXZ0) {
  Vector2d p1, p2, p3;
  p1 << 0.0, 0.0;
  p2 << 7.0, 0.0;
  p3 << 5.5, 0.5;

  double zA = 0.0;
  Vector3d res = CalcGeneralTriangleCompliance(zA, p1, p2, p3);
  double res_interpolate = res(0) * p1(0) + res(1) * p2(0) + res(2) * p3(0);

  const Vector3<double>& I_12 = CalcIntegralReferenceTriangle(p1, p2, zA);
  const Vector3<double>& I_23 = CalcIntegralReferenceTriangle(p2, p3, zA);
  const Vector3<double>& I_31 = CalcIntegralReferenceTriangle(p3, p1, zA);

  double expected_res = fabs(I_12(0) + I_23(0) + I_31(0));

  EXPECT_NEAR(
      res_interpolate, expected_res,
      10 * std::numeric_limits<double>::epsilon());
}

}  // namespace
}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake

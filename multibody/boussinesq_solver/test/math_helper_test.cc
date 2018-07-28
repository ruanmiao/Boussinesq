#include "drake/multibody/boussinesq_solver/math_helper.h"

#include <gtest/gtest.h>

#include "drake/common/eigen_types.h"

namespace drake {
namespace multibody {
namespace boussinesq_solver {
namespace {

using Eigen::Vector2d;
using Eigen::Vector3d;

GTEST_TEST(TriangleOrientationTest, ClockWise) {
  const Vector2d p1(0.0, 1.0);
  const Vector2d p2(1.0, 0.0);
  const Vector2d p3(0.0, 0.0);
  const int orientation = CalcTriangleOrientation(p1, p2, p3);
  EXPECT_EQ(orientation, -1);
}

GTEST_TEST(TriangleOrientationTest, CoLinear) {
  const Vector2d p1(1.0, 1.0);
  const Vector2d p2(2.0, 2.0);
  const Vector2d p3(0.0, 0.0);
  const int orientation = CalcTriangleOrientation(p1, p2, p3);
  EXPECT_EQ(orientation, 0);
}

GTEST_TEST(TriangleOrientationTest, CounterClockWise) {
  const Vector2d p1(1.0, 0.0);
  const Vector2d p2(0.0, 1.0);
  const Vector2d p3(0.0, 0.0);
  const int orientation = CalcTriangleOrientation(p1, p2, p3);
  EXPECT_EQ(orientation, 1);
}

GTEST_TEST(TriangleOrientationTest, Colinear2over3PI) {
  Vector2d p1, p2, p3;
  const double r = 1.0;
  p1 << r / 2 * cos(M_PI * 2 / 3), r / 2 * sin(M_PI * 2 / 3);
  p2 << r * cos(M_PI * 2 / 3), r * sin(M_PI * 2 / 3);
  p3(0.0, 0.0);
  const int orientation = CalcTriangleOrientation(p1, p2, p3);
  EXPECT_EQ(orientation, 0);
}

GTEST_TEST(TriangleAreaTest, AreaValue) {
  const Vector2d p1(0.0, 0.0);
  const Vector2d p2(1.0, 0.0);
  const Vector2d p3(0.0, 1.0);
  const double area = CalcTriangleArea(p1, p2, p3);
  const double expected_area = 0.5;
  EXPECT_NEAR(area, expected_area, 10 * std::numeric_limits<double>::epsilon());
}

GTEST_TEST(TriangleAreaTest, AreaSign) {
  const Vector2d p1(0.0, 0.0);
  const Vector2d p2(1.0, 0.3);
  const Vector2d p3(-0.2, 1.2);
  const double area = CalcTriangleArea(p1, p2, p3);
  EXPECT_GT(area, 0.0);
}

/// The expected values for this test are the results by running Matlab. The
/// precision (15 digits) of the expected values if the same as the "long"
/// in Matlab
GTEST_TEST(IntegralJm0nN1Test, PosToNeg) {
  const double theta_0 = M_PI / 3;
  const double theta_f = -M_PI / 3;
  const double Jmn = CalcIntegralJ0minus1(theta_0, theta_f);
  const double expected_Jmn = -2.633915793849633;
  EXPECT_NEAR(Jmn, expected_Jmn, 10 * std::numeric_limits<double>::epsilon());
}

/// The expected values for this test are the results by running Matlab. The
/// precision (15 digits) of the expected values if the same as the "long"
/// in Matlab
GTEST_TEST(IntegralJm0nN1Test, PiOverTwo) {
  const double theta_0 = M_PI / 2;
  const double theta_f = 0;
  const double Jmn = CalcIntegralJ0minus1(theta_0, theta_f);
  const double expected_Jmn = -38.025003373828866;
  EXPECT_NEAR(Jmn, expected_Jmn, 10 * std::numeric_limits<double>::epsilon());
}

GTEST_TEST(IntegralJm1nN2Test, PosToNeg) {
  const double theta_0 = M_PI / 3;
  const double theta_f = -M_PI / 3;
  const double Jmn = CalcIntegralJ1minus2(theta_0, theta_f);
  const double expected_Jmn = 0;
  EXPECT_NEAR(Jmn, expected_Jmn, 10 * std::numeric_limits<double>::epsilon());
}

/// The expected values for this test are the results by running Matlab. The
/// precision (15 digits) of the expected values if the same as the "long"
/// in Matlab
GTEST_TEST(IntegralJm1nN2Test, PiOverTwo) {
  const double theta_0 = M_PI / 2;
  const double theta_f = 0;
  const double Jmn = CalcIntegralJ1minus2(theta_0, theta_f);
  const double expected_Jmn = -1.633123935319537e+16;
  EXPECT_NEAR(
      Jmn, expected_Jmn, 10e+16 * std::numeric_limits<double>::epsilon());
}

/// The expected values for this test are the results by running Matlab. The
/// precision (15 digits) of the expected values if the same as the "long"
/// in Matlab
GTEST_TEST(IntegralIm0nN1P1Test, NonTrivialAlpha) {
  const double theta_0 = 2 * M_PI / 5;
  const double theta_f = -M_PI / 3;
  const double Jmn = CalcIntegralI0minus1P1(theta_0, theta_f, 0.111);
  const double expected_Jmn = -3.151400810889063;
  EXPECT_NEAR(
      Jmn, expected_Jmn, 10 * std::numeric_limits<double>::epsilon());
}

/// The expected values for this test are the results by running Matlab. The
/// precision (15 digits) of the expected values if the same as the "long"
/// in Matlab
GTEST_TEST(IntegralIm1nN2P1Test, NonTrivialAlpha) {
  const double theta_0 = -M_PI / 3;
  const double theta_f = M_PI / 4;
  const double Jmn = CalcIntegralI1minus2P1(theta_0, theta_f, 0.111);
  const double expected_Jmn = -0.583448860928930;
  EXPECT_NEAR(
      Jmn, expected_Jmn, 10 * std::numeric_limits<double>::epsilon());
}

/// The expected values for this test are the results by running Matlab. The
/// precision (15 digits) of the expected values if the same as the "long"
/// in Matlab
GTEST_TEST(IntegralIm2nN1PN1Test, NonTrivialAlpha) {
  const double theta_0 = M_PI * 2 / 5;
  const double theta_f = -M_PI / 4;
  const double Jmn = CalcIntegralI2minus1Pminus1(theta_0, theta_f, 0.15);
  const double expected_Jmn = -1.073471020695750;
  EXPECT_NEAR(
      Jmn, expected_Jmn, 10 * std::numeric_limits<double>::epsilon());
}

/// The expected values for this test are the results by running Matlab. The
/// precision (15 digits) of the expected values if the same as the "long"
/// in Matlab
GTEST_TEST(IntegralIm1n0PN1Test, NonTrivialAlpha) {
  const double theta_0 = M_PI * 2 / 5;
  const double theta_f = -M_PI / 4;
  const double Jmn = CalcIntegralI10Pminus1(theta_0, theta_f, 0.15);
  const double expected_Jmn = -0.401394895991949;
  EXPECT_NEAR(
      Jmn, expected_Jmn, 10 * std::numeric_limits<double>::epsilon());
}

GTEST_TEST(CalcTransformationFromTriangleFrameTest, GeneralCase) {
  const Vector3d p1(6.0, 4.0, 2.0);
  const Vector3d p2(7.0, 9.0, 1.8);
  const Vector3d p3(5.0, 8.0, 2.5);
  const Vector3d xA(4.0, 3.0, -1.0);

  const Vector3d u12 = p2 - p1;
  const Vector3d u23 = p3 - p2;
  const Vector3d u31 = p1 - p3;

  const Eigen::Isometry3d X_WT = CalcTransformationFromTriangleFrame(
      p1, p2, p3, xA);
  const Eigen::Isometry3d X_TW = X_WT.inverse();

  const Eigen::Vector3d xA_T = X_TW * xA;
  const Eigen::Vector3d p1_T = X_TW * p1;
  const Eigen::Vector3d p2_T = X_TW * p2;
  const Eigen::Vector3d p3_T = X_TW * p3;

  const Vector3d u12_T = p2_T - p1_T;
  const Vector3d u23_T = p3_T - p2_T;
  const Vector3d u31_T = p1_T - p3_T;

  EXPECT_NEAR(
      xA_T(0), 0.0, 10 * std::numeric_limits<double>::epsilon());
  EXPECT_NEAR(
      xA_T(1), 0.0, 10 * std::numeric_limits<double>::epsilon());

  EXPECT_NEAR(
      p1_T(2), 0.0, 10 * std::numeric_limits<double>::epsilon());
  EXPECT_NEAR(
      p2_T(2), 0.0, 10 * std::numeric_limits<double>::epsilon());
  EXPECT_NEAR(
      p3_T(2), 0.0, 10 * std::numeric_limits<double>::epsilon());

  EXPECT_NEAR(
      u12.norm(), u12_T.norm(), 10 * std::numeric_limits<double>::epsilon());
  EXPECT_NEAR(
      u23.norm(), u23_T.norm(), 10 * std::numeric_limits<double>::epsilon());
  EXPECT_NEAR(
      u31.norm(), u31_T.norm(), 10 * std::numeric_limits<double>::epsilon());
}


}  // namespace
}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake

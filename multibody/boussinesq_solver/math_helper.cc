#include "drake/multibody/boussinesq_solver/math_helper.h"

#include <limits>
#include "drake/common/drake_assert.h"

#include <iostream>
#define PRINT_VAR(a) std::cout << #a": " << a << std::endl;

namespace drake {
namespace multibody {
namespace boussinesq_solver {

double CalcTriangleArea(const Vector2<double>& p1,
                        const Vector2<double>& p2,
                        const Vector2<double>& p3) {
  const Vector2<double> u1 = p2 - p1;
  const Vector2<double> u2 = p3 - p1;
  return (u1(0) * u2(1) - u1(1) * u2(0)) / 2.0;
}

Vector3<double> CalcTriangleArea(const Vector3<double>& p1,
                        const Vector3<double>& p2,
                        const Vector3<double>& p3) {
  const Vector3<double> u1 = p2 - p1;
  const Vector3<double> u2 = p3 - p1;
  return u1.cross(u2)/2.0;
}

int CalcTriangleOrientation(const Vector2<double>& p1,
                            const Vector2<double>& p2,
                            const Vector2<double>& p3) {
  const Vector2<double> u1 = p2 - p1;
  const Vector2<double> u2 = p3 - p1;
  const double signed_area = (u1(0) * u2(1) - u1(1) * u2(0)) / 2.0;

  const double tol = 10 * std::numeric_limits<double>::epsilon();
  int is_clockwise = 0;
  if (signed_area > tol) is_clockwise = 1;
  if (signed_area < -tol) is_clockwise = -1;
  return is_clockwise;
}

double CalcIntegralJ0minus1(double theta_0, double theta_f) {
  double Jmn = log(fabs(1 / cos(theta_f) + tan(theta_f))) -
      log(fabs(1 / (cos(theta_0)) + tan(theta_0)));
  return Jmn;
}

double CalcIntegralJ1minus2(double theta_0, double theta_f) {
  double Jmn = 1 / (cos(theta_f)) - 1 / cos(theta_0);
  return Jmn;
}

double CalcIntegralJ00(double theta_0, double theta_f) {
  return theta_f - theta_0;
}

double CalcIntegralI0minus1P1(double theta_0, double theta_f, double alpha) {
  if (fabs(alpha - 0) < 10 * std::numeric_limits<double>::epsilon()) {
    return  CalcIntegralJ0minus1(theta_0, theta_f);
  }
  double sol0_1 = alpha * asin(alpha * sin(theta_0));
  double solf_1 = alpha * asin(alpha * sin(theta_f));

  double sol0_2 = sqrt(1 - pow(alpha, 2)) *
      atanh((sqrt(2 * (1 / pow(alpha, 2) - 1)) * sin(theta_0))
                / sqrt(cos(2 * theta_0) + 2 * (1 / pow(alpha, 2) - 1) + 1));

  double solf_2 = sqrt(1 - pow(alpha, 2)) *
      atanh((sqrt(2 * (1 / pow(alpha, 2) - 1)) * sin(theta_f))
                / sqrt(cos(2 * theta_f) + 2 * (1 / pow(alpha, 2) - 1) + 1));
  return  solf_1 + solf_2 - sol0_1 - sol0_2;
}

double CalcIntegralI1minus2P1(double theta_0, double theta_f, double alpha) {
  if (fabs(alpha - 0) < 10 * std::numeric_limits<double>::epsilon()) {
    return  CalcIntegralJ1minus2(theta_0, theta_f);
  }
  double sol0_1 = alpha * asinh(
      (alpha * sqrt(1 - pow(alpha, 2)) * cos(theta_0)) / (pow(alpha, 2) - 1));
  double solf_1 = alpha * asinh(
      (alpha * sqrt(1 - pow(alpha, 2)) * cos(theta_f)) / (pow(alpha, 2) - 1));

  double sol0_2 = sqrt(pow(alpha, 2) * pow(cos(theta_0), 2) - pow(alpha, 2) + 1)
      / cos(theta_0);
  double solf_2 = sqrt(pow(alpha, 2) * pow(cos(theta_f), 2) - pow(alpha, 2) + 1)
      / cos(theta_f);

  return solf_1 + solf_2 - sol0_1 - sol0_2;
}

double CalcIntegralI2minus1Pminus1(
    double theta_0, double theta_f, double alpha) {
  double one_minus_a2 = 1 - pow(alpha, 2);

  double sol0_1 = - 1 / alpha * asin(alpha * sin(theta_0));
  double solf_1 = - 1 / alpha * asin(alpha * sin(theta_f));

  double sol0_2 = sqrt(1/one_minus_a2) *
      atanh((sqrt(2 * (1 / pow(alpha, 2) - 1)) * sin(theta_0))
                / sqrt(cos(2 * theta_0) + 2 * (1 / pow(alpha, 2) - 1) + 1));
  double solf_2 = sqrt(1/one_minus_a2) *
      atanh((sqrt(2 * (1 / pow(alpha, 2) - 1)) * sin(theta_f))
                / sqrt(cos(2 * theta_f) + 2 * (1 / pow(alpha, 2) - 1) + 1));

  return solf_1 + solf_2 - sol0_1 - sol0_2;
}

double CalcIntegralI10Pminus1(double theta_0, double theta_f, double alpha) {
  double sol0 = asinh((alpha * sqrt(1 - pow(alpha, 2)) * cos(theta_0)) /
      (pow(alpha, 2) - 1)) / alpha;
  double solf = asinh((alpha * sqrt(1 - pow(alpha, 2)) * cos(theta_f)) /
      (pow(alpha, 2) - 1)) / alpha;
  return solf - sol0;
}

Eigen::Isometry3d CalcTransformationFromTriangleFrame(
    const Eigen::Vector3d& p1, const Eigen::Vector3d& p2,
    const Eigen::Vector3d& p3, const Eigen::Vector3d& xA) {
  const Eigen::Vector3d u1 = p2 - p1;
  const Eigen::Vector3d u2 = p3 - p1;
  const Eigen::Vector3d area = u1.cross(u2);

  if (fabs(area.norm()) < 10 * std::numeric_limits<double>::epsilon()) {
    PRINT_VAR(p1);
    PRINT_VAR(p2);
    PRINT_VAR(p3);
    PRINT_VAR(xA);
  }
  DRAKE_ASSERT(fabs(area.norm()) > 10 * std::numeric_limits<double>::epsilon());

  const Eigen::Vector3d z_T = area / area.norm();
  const Eigen::Vector3d x_T = u1 / u1.norm();
  const Eigen::Vector3d y_T = z_T.cross(x_T);

  Eigen::Matrix3d R_WT;
  R_WT.col(0) = x_T;
  R_WT.col(1) = y_T;
  R_WT.col(2) = z_T;

  const Eigen::Vector3d T0_W = xA + z_T.dot(p1 - xA) * z_T;
  Eigen::Isometry3d X_WT;
  X_WT.linear() = R_WT;
  X_WT.translation() = T0_W;

  return X_WT;
}

}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake

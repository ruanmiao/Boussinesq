#include "drake/multibody/boussinesq_solver/integral_reference_triangle.h"

#include <limits>

#include "drake/multibody/boussinesq_solver/math_helper.h"

namespace drake {
namespace multibody {
namespace boussinesq_solver {

Vector3<double> CalcIntegralReferenceTriangle(const Vector2<double>& p1,
                                              const Vector2<double>& p2,
                                              double zA) {
  (void) zA;
  Vector3<double> results;

  const double psi_offset_from_x = atan2(p1(1), p1(0));
  const double r1 = p1.norm();
  const double r2 = p2.norm();

  const double kTolerance = 10 * std::numeric_limits<double>::epsilon();
  if (r1 < kTolerance || r2 < kTolerance) {
    return Vector3<double>::Zero();
  }

  const Eigen::MatrixXd origin = Eigen::MatrixXd::Zero(2, 1);
  const double signed_theta =
      acos(p1.dot(p2) / (p1.norm() * p2.norm())) *
          CalcTriangleOrientation(p1, p2, origin);

  if (fabs(signed_theta) < kTolerance) {
    return Vector3<double>::Zero();
  }

  const double a = (r2 * cos(signed_theta) - r1) / (r2 * sin(signed_theta));
  const double beta = sqrt(r1 * r1 / (1 + a * a));
  const double phi_integral_offset = atan(a);

  const double theta_0 = phi_integral_offset;
  const double theta_f = phi_integral_offset + signed_theta;

  Matrix2<double> rotation_matrix_phi;
  rotation_matrix_phi << cos(phi_integral_offset), sin(phi_integral_offset),
      -sin(phi_integral_offset), cos(phi_integral_offset);

  Matrix2<double> rotation_matrix_psi;
  rotation_matrix_psi << cos(psi_offset_from_x), -sin(psi_offset_from_x),
      sin(psi_offset_from_x), cos(psi_offset_from_x);

  Vector2<double> integral_without_rotation;
  integral_without_rotation << CalcIntegralJm0nN1(theta_0, theta_f),
      CalcIntegralJm1nN2(theta_0, theta_f);

  results.head(2) = beta * beta /
      2 * rotation_matrix_psi * rotation_matrix_phi * integral_without_rotation;
  results(2) = beta * CalcIntegralJm0nN1(theta_0, theta_f);

  return results;
}



}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake

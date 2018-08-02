#include "drake/multibody/boussinesq_solver/integral_reference_triangle.h"

#include <limits>

#include "drake/common/drake_assert.h"
#include "drake/multibody/boussinesq_solver/math_helper.h"

namespace drake {
namespace multibody {
namespace boussinesq_solver {

// Reference: M. Carley, Potential Integrals on Triangles, arXiv:1201.4938v1,
// [math.NA], 24 Jan 2012.
// Notes: The integral of the integrants below over the triangle surrounded
// by O, p₁, p₂, where O is the origin.
// results(0) =  ∫∫ (x/R) dxdy
// results(1) =  ∫∫ (y/R) dxdy
// results(2) =  ∫∫ (1/R) dxdy
// results(0) = Rot(ψ) Rot(-φ) [∫∫ (r*cos(τ + φ)/R) r dr dτ]
// results(1)                 [∫∫ (r*sin(τ + φ)/R) r dr dτ]
// The interval of the above integration is [τ ∈ [θ₀, θf]],
// where θ₀ = φ, θf = φ + Θ, r₁ = ||p₁||, r₂ = ||p₂||,
// a = (r₂ * cos(Θ) - r₁) / (r₂ * sin(Θ))
// φ = atan(a), ψ = atan2(yₚ₂, xₚ₁)
// β = √{r₁²/(1+a²)}
// Θ is the signed angle from Op₁ to Op₂.
//
// The results of ∫∫ (r*cos(τ + φ)/R) r dr dτ, ∫∫ (r*sin(τ + φ)/R) r dr dτ,
// and ∫∫ (1/R) r dr dτ can be looked up from Table 1 and Table 4 in the
// reference paper.
// Special note: the value of Jmn when m = 0, n = -1 is wrong in Table 4. It
// should be -1/a ln|secθ + tanθ| instead.
Vector3<double> CalcIntegralReferenceTriangle(const Vector2<double>& p1,
                                              const Vector2<double>& p2) {
  Vector3<double> results;

  const double psi_offset_from_x = atan2(p1(1), p1(0));
  const double r1 = p1.norm();
  const double r2 = p2.norm();

  const double kTolerance = 10 * std::numeric_limits<double>::epsilon();
  if (r1 < kTolerance || r2 < kTolerance) {
    return Vector3<double>::Zero();
  }

  const Eigen::MatrixXd origin = Eigen::MatrixXd::Zero(2, 1);
  const double dot_product =
      (p1(0) * p2(0) + p1(1) * p2(1)) / (p1.norm() * p2.norm());

  double a_signed_theta = acos(dot_product);
  if (dot_product > 1.0 || dot_product < -1.0) a_signed_theta = 0;

  const double orientation = CalcTriangleOrientation(p1, p2, origin);
  const double signed_theta = a_signed_theta * orientation;

  if (fabs(signed_theta) < kTolerance ||
      fabs(M_PI - signed_theta) < kTolerance) {
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
  integral_without_rotation <<
                            CalcIntegralJ0minus1(theta_0, theta_f),
      CalcIntegralJ1minus2(theta_0, theta_f);

  results.head(2) = beta * beta /
      2 * rotation_matrix_psi * rotation_matrix_phi * integral_without_rotation;
  results(2) = beta * CalcIntegralJ0minus1(theta_0, theta_f);

  return results;
}

Vector3<double> CalcIntegralReferenceTriangle(const Vector2<double>& p1,
                                              const Vector2<double>& p2,
                                              double zA) {
  if (pow(zA, 2) < 10 * std::numeric_limits<double>::epsilon()) {
    return  CalcIntegralReferenceTriangle(p1, p2);
  }

  Vector3<double> results;

  const double psi_offset_from_x = atan2(p1(1), p1(0));
  const double r1 = p1.norm();
  const double r2 = p2.norm();

  const double kTolerance = 10 * std::numeric_limits<double>::epsilon();
  if (r1 < kTolerance || r2 < kTolerance) {
    return Vector3<double>::Zero();
  }

  const Eigen::MatrixXd origin = Eigen::MatrixXd::Zero(2, 1);
  const double dot_product =
      (p1(0) * p2(0) + p1(1) * p2(1)) / (p1.norm() * p2.norm());

  double a_signed_theta = acos(dot_product);
  if (dot_product > 1.0 || dot_product < -1.0) a_signed_theta = 0;

  const double orientation = CalcTriangleOrientation(p1, p2, origin);
  const double signed_theta = a_signed_theta * orientation;

  if (fabs(signed_theta) < kTolerance ||
      fabs(M_PI - signed_theta) < kTolerance) {
    return Vector3<double>::Zero();
  }

  const double a = (r2 * cos(signed_theta) - r1) / (r2 * sin(signed_theta));
  const double beta = sqrt((pow(r1, 2) + pow(zA, 2) * (1 + pow(a, 2))) /
      (1 + pow(a, 2)));
  const double phi_integral_offset = atan(a);

  const double alpha = sqrt(pow(zA, 2) / pow(beta, 2));
  const double alpha_prime = sqrt(1 - pow(alpha, 2));

  const double theta_0 = phi_integral_offset;
  const double theta_f = phi_integral_offset + signed_theta;

  const double delta_0 = sqrt(1 - pow(alpha * sin(theta_0), 2));
  const double delta_f = sqrt(1 - pow(alpha * sin(theta_f), 2));

  Matrix2<double> rotation_matrix_phi;
  rotation_matrix_phi << cos(phi_integral_offset), sin(phi_integral_offset),
      -sin(phi_integral_offset), cos(phi_integral_offset);

  Matrix2<double> rotation_matrix_psi;
  rotation_matrix_psi << cos(psi_offset_from_x), -sin(psi_offset_from_x),
      sin(psi_offset_from_x), cos(psi_offset_from_x);

  Vector2<double> integral_without_rotation;

  integral_without_rotation(0) = pow(beta, 2) * alpha_prime / 2 *
      CalcIntegralI0minus1P1(theta_0, theta_f, alpha) +
      pow(zA, 2) * alpha_prime / 2 *
          CalcIntegralI2minus1Pminus1(theta_0, theta_f, alpha) -
      pow(zA, 2) / 4 * ((sin(theta_f) *
          log((delta_f + alpha_prime) / (delta_f - alpha_prime))) -
          (sin(theta_0) *
              log((delta_0 + alpha_prime) / (delta_0 - alpha_prime))));

  integral_without_rotation(1) = pow(beta, 2) * alpha_prime / 2 *
      CalcIntegralI1minus2P1(theta_0, theta_f, alpha) -
      pow(zA, 2) * alpha_prime / 2 *
          CalcIntegralI10Pminus1(theta_0, theta_f, alpha) +
      pow(zA, 2) / 4 * ((cos(theta_f) *
          log((delta_f + alpha_prime) / (delta_f - alpha_prime))) -
          (cos(theta_0) *
              log((delta_0 + alpha_prime) / (delta_0 - alpha_prime))));

  DRAKE_ASSERT(fabs(delta_0 - alpha_prime) >
      std::numeric_limits<double>::epsilon());
  DRAKE_ASSERT(fabs(delta_f - alpha_prime) >
      std::numeric_limits<double>::epsilon());

  DRAKE_ASSERT(fabs(integral_without_rotation(0)) <
      std::numeric_limits<double>::infinity());
  DRAKE_ASSERT(fabs(integral_without_rotation(1)) <
      std::numeric_limits<double>::infinity());

  results.head(2) =
      rotation_matrix_psi * rotation_matrix_phi * integral_without_rotation;
  results(2) = beta * CalcIntegralI0minus1P1(theta_0, theta_f, alpha) -
      fabs(zA) * CalcIntegralJ00(theta_0, theta_f);

  return results;
}



}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake

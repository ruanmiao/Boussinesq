#include "drake/multibody/boussinesq_solver/math_helper.h"

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

int CalcTriangleOrientation(const Vector2<double>& p1,
                            const Vector2<double>& p2,
                            const Vector2<double>& p3) {
  const Vector2<double> u1 = p2 - p1;
  const Vector2<double> u2 = p3 - p1;
  const double signed_area = (u1(0) * u2(1) - u1(1) * u2(0)) / 2.0;

  int is_clockwise = 0;
  if (signed_area > 0) is_clockwise = 1;
  if (signed_area < 0) is_clockwise = -1;
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

}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake

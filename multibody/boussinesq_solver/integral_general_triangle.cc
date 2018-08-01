#include "drake/multibody/boussinesq_solver/integral_general_triangle.h"

#include <limits>

#include "drake/multibody/boussinesq_solver/integral_reference_triangle.h"
#include "drake/multibody/boussinesq_solver/math_helper.h"

namespace drake {
namespace multibody {
namespace boussinesq_solver {

Vector3<double> CalcGeneralTriangleCompliance(const Vector2<double>& p1,
                                              const Vector2<double>& p2,
                                              const Vector2<double>& p3,
                                              double k_const) {
  return CalcGeneralTriangleCompliance(0.0, p1, p2, p3, k_const);
}

Vector3<double> CalcGeneralTriangleCompliance(
    double zA,
    const Vector2<double>& p1,
    const Vector2<double>& p2,
    const Vector2<double>& p3,
    double k_const) {

  if (fabs(zA) < std::numeric_limits<double>::epsilon()) {
    zA = 0.0;
  }
  const int is_counter_cw = CalcTriangleOrientation(p1, p2, p3);

  const Vector3<double>& I_12 = CalcIntegralReferenceTriangle(p1, p2, zA);
  const Vector3<double>& I_23 = CalcIntegralReferenceTriangle(p2, p3, zA);
  const Vector3<double>& I_31 = CalcIntegralReferenceTriangle(p3, p1, zA);

  Vector3<double> I = I_12 + I_23 + I_31;
  if (is_counter_cw < 0) { I = -I; }

  Matrix2<double> M;
  M.col(0) = p1 - p3;
  M.col(1) = p2 - p3;
  Matrix2<double> Mi = M.inverse();
  Vector2<double> m = -Mi * p3;

  Matrix3<double> K_abc;
  Vector3<double> ka, kb, kc;
  ka << Mi(0, 0), Mi(1, 0), (-Mi(0, 0) - Mi(1, 0));
  kb << Mi(0, 1), Mi(1, 1), (-Mi(0, 1) - Mi(1, 1));
  kc << m(0), m(1), (1 - m(0) - m(1));
  K_abc.col(0) = ka;
  K_abc.col(1) = kb;
  K_abc.col(2) = kc;

  const Vector3<double> compliance = k_const * K_abc * I;
  return compliance;
}


}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake

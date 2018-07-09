#include "drake/multibody/boussinesq_solver/integral_reference_triangle.h"

namespace drake {
namespace multibody {
namespace boussinesq_solver {

int TriangleOrientation(const Vector2<double>& p1,
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

}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake

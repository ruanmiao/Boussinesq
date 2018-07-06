#include "drake/multibody/boussinesq_solver/triangle_area.h"

namespace drake {
namespace multibody {
namespace boussinesq_solver {

double TriangleArea(const Vector2<double>& p1,
                    const Vector2<double>& p2,
                    const Vector2<double>& p3) {
  const Vector2<double> u1 = p2 - p1;
  const Vector2<double> u2 = p3 - p1;
  return (u1(0) * u2(1) - u1(1) * u2(0)) / 2.0;
}

}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake

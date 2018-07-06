#pragma once

#include "drake/common/eigen_types.h"

namespace drake {
namespace multibody {
namespace boussinesq_solver {

/// This method computes the area of the triangle defined by points p1, p2 and
/// p3. The return value is positive if the points are oriented counterclockwise
/// and negative otherwise.
/// @param p1 First vertex.
/// @param p2 Second vertex.
/// @param p3 Third vertex.
/// @returns The are of the triangle defined by p1, p2 and p3.
double TriangleArea(const Vector2<double>& p1,
                    const Vector2<double>& p2,
                    const Vector2<double>& p3);

}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake

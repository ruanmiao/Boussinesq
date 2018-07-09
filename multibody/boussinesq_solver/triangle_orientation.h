#pragma once

#include "drake/common/eigen_types.h"

namespace drake {
namespace multibody {
namespace boussinesq_solver {

/// This method check the orientation of the three points p1, p2, p3. Return 1
/// if the three points are ordered counter-clockwise. -1 if clockwise. 0 if
/// the three are co-linear.
/// @param p1 First vertex.
/// @param p2 Second vertex.
/// @param p3 Third vertex.
/// @returns The orientation of the triangle defined by p1, p2 and p3.
int TriangleOrientation(const Vector2<double>& p1,
                    const Vector2<double>& p2,
                    const Vector2<double>& p3);
}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake

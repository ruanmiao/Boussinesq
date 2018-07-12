#pragma once

#include "drake/common/eigen_types.h"

namespace drake {
namespace multibody {
namespace boussinesq_solver {

/// This method returns a 3x1 vector containing the compliance of a triangle
/// with respect to the origin. The triangle is defined by the three input
/// 2D vertices p1, p2, and p3. The last input k_const is the
/// compliance-related constant.
/// Furthermore, regardless the orientation of the triangle, the result will
/// be the same.
/// @param p1 First vertex.
/// @param p2 Second vertex.
/// @param p3 Third vertex.
/// @param k_const The compliance-related constant.
/// @returns 3D vector, with V(0), V(1), and V(2) each storing the compliance
/// of the three vertices p1, p2, and p3 in the same order with respect to the
/// origin.
Vector3<double> CalGenralTriangleCompliance(
    const Vector2<double>& p1,
    const Vector2<double>& p2,
    const Vector2<double>& p3,
    double k_const = 1);


}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake

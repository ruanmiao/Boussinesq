#pragma once

#include "drake/common/eigen_types.h"

namespace drake {
namespace multibody {
namespace boussinesq_solver {

/// This method returns a vector containing the integrals of the integrands
/// x/R, y/R, and 1/R over a triangle defined by vertices p1, p2 and the origin
/// O. The triangle formed by vertices p1, p2 and O, in that order, might form
/// a clockwise oriented as well as a counter-clockwise oriented triangle.
/// @param p1 First vertex.
/// @param p2 Second vertex.
/// @param zA The signed distance from xA to the plane of the triangle .
/// @returns V(0), V(1), and V(2) each storing one integral result for : x/R,
/// y/R, and 1/R, in that order. If p1, p2 and O are in counter-clockwise order,
/// the last integral on 1/R, V(2), is positive. Otherwise V(2) is negative.
Vector3<double> CalcIntegralReferenceTriangle(
    const Vector2<double> &p1, const Vector2<double> &p2, const double zA = 0);

}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake

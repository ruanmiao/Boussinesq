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
/// @returns 3D vector C, with C(0), C(1), and C(2) each storing the compliance
/// of the three vertices p1, p2, and p3 in the same order with respect to the
/// field point.
/// For example, let p be the vector, with p(i) being the pressure at each node
/// pi. Then the deformation of the field point induced by the triangular
/// element consist of p1, p2, and p3 can be interpolated as u = Cp.
Vector3<double> CalcGeneralTriangleCompliance(
    const Vector2<double>& p1,
    const Vector2<double>& p2,
    const Vector2<double>& p3,
    double k_const = 1);

/// This method returns a 3x1 vector containing the compliance of a triangle
/// with respect to the field point, with its projection at the origin in the
/// plane of the triangle, zA is the signed distance from the field point to
/// the plane. The triangle is defined by the three input
/// 2D vertices p1, p2, and p3. The last input k_const is the
/// compliance-related constant.
/// Furthermore, regardless the orientation of the triangle, the result will
/// be the same.
/// @param zA Singed distance from the field point to the plane
/// @param p1 First vertex.
/// @param p2 Second vertex.
/// @param p3 Third vertex.
/// @param k_const The compliance-related constant.
/// @returns 3D vector C, with C(0), C(1), and C(2) each storing the compliance
/// of the three vertices p1, p2, and p3 in the same order with respect to the
/// field point.
/// For example, let p be the vector, with p(i) being the pressure at each node
/// pi. Then the deformation of the field point induced by the triangular
/// element consist of p1, p2, and p3 can be interpolated as u = Cp.
Vector3<double> CalcGeneralTriangleCompliance(
    double zA,
    const Vector2<double>& p1,
    const Vector2<double>& p2,
    const Vector2<double>& p3,
    double k_const = 1);

}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake

#pragma once

#include <vector>

#include "drake/common/eigen_types.h"

namespace drake {
namespace multibody {
namespace boussinesq_solver {

/// This method returns a 1xn row vector containing the compliance of all the
/// vertices with respect to a vertex (node_A) in the mesh, where n is the
/// number of vertices in the mesh.
/// @param points_in_mesh Each entry stores the position of a node.
/// @param triangles_in_mesh Each entry stores the three indexes of the nodes
/// forming the triangle
/// @param node_A Reference point.
/// @param k_const The compliance-related constant.
/// @returns 1xn vector/row matrix, with v(i) being the compliance from point
/// i to point A
MatrixX<double> CalcRowComplianceMatrix(
    const std::vector<Eigen::Vector3d>& points_in_mesh,
    const std::vector<Eigen::Vector3i>& triangles_in_mesh,
    const Vector3<double>& node_A,
    double k_const = 1);


}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake

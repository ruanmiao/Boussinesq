#pragma once

#include <vector>

#include "drake/common/eigen_types.h"

namespace drake {
namespace multibody {
namespace boussinesq_solver {

/// This method returns an nxn matrix containing the compliance of all the
/// vertices with respect to all vertices in the mesh, where n is the
/// number of vertices in the mesh.
/// @param points_in_mesh Each entry stores the position of a node.
/// @param triangles_in_mesh Each entry stores the three indexes of the nodes
/// forming the triangle
/// @param k_const The compliance-related constant.
/// @returns nxn matrix, with M(i,j) being the compliance from point
/// j to point i.
MatrixX<double> CalcComplianceMatrix(
    const std::vector<Eigen::Vector3d>& points_in_mesh,
    const std::vector<Eigen::Vector3i>& triangles_in_mesh,
    double k_const = 1.0);



}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake

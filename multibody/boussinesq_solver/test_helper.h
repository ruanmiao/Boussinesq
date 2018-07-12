#pragma once

#include <utility>
#include <vector>

#include "drake/common/eigen_types.h"

namespace drake {
namespace multibody {
namespace boussinesq_solver {

/// Given the 2D Parallelogram area defined by the 4 vertices, with the number
/// of elements for each edge specified.
//  TODO: The layout of the mesh.
/// @param p1 First vertex.
/// @param p2 Second vertex.
/// @param p3 Third vertex.
/// @param p4 Third vertex.
/// @returns The first vector in the return is the location of all nodes in
/// the mesh. The second vector in the return is the node indexes of all the
/// triangles in the mesh.
std::pair<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3i>>
MeshSquare(
    const Vector2<double>& p1, const Vector2<double>& p2,
    const Vector2<double>& p3, const Vector2<double>& p4,
    int num_px = 2, int num_py = 2);

}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake

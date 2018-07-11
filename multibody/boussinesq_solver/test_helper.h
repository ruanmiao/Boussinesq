#pragma once

#include <utility>
#include <vector>

#include "drake/common/eigen_types.h"

namespace drake {
namespace multibody {
namespace boussinesq_solver {

/// Given the 2D Parallelogram area defined by the 4 vertices, with the number
/// of elements for each edge specified. The layout of the mesh looks like
/// the following:
/// @param p1 First vertex.
/// @param p2 Second vertex.
/// @param p3 Third vertex.
/// @returns The are of the triangle defined by p1, p2 and p3.
std::pair<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3i>>
MeshSquare(
    const Vector2<double>& p1, const Vector2<double>& p2,
    const Vector2<double>& p3, const Vector2<double>& p4,
    int num_px = 2, int num_py = 2);

}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake

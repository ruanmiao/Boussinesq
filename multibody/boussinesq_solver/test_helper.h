#pragma once

#include <utility>
#include <vector>

#include "drake/common/eigen_types.h"

namespace drake {
namespace multibody {
namespace boussinesq_solver {

/// Given the 2D Parallelogram area defined by the 4 vertices, with the number
/// of elements for each edge specified.
/// Layout of the 4 vertices:       p3 --- (p4)
///                                 p1 --- p2
// TODO(Mengyao): The layout of the mesh.
/// @param p1 First vertex.
/// @param p2 Second vertex.
/// @param p3 Third vertex.
/// @param num_px Number of vertices per edge in x direction.
/// @param num_py Number of vertices per edge in y direction.
/// @returns The first vector in the return is the location of all nodes in
/// the mesh. The second vector in the return is the node indexes of all the
/// triangles in the mesh.
std::pair<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3i>>
MeshSquare(
    const Vector2<double>& p1, const Vector2<double>& p2,
  const Vector2<double>& p3, int num_px = 2, int num_py = 2);

/// Return the pressure on each node given the points locations, where the
/// pressure field function is p(pᵢ) = rᵢ.
/// @param points_in_mesh: Locations of all nodes in the mesh
/// @returns The pressures for all the nodes in the mesh
Eigen::VectorXd GetPressureIntegrandR(
    const std::vector<Eigen::Vector3d>& points_in_mesh);

/// Return the pressure on each node given the points locations, where the
/// pressure field function is p(pᵢ) = xᵢ.
/// @param points_in_mesh: Locations of all nodes in the mesh
/// @returns The pressures for all the nodes in the mesh
Eigen::VectorXd GetPressureIntegrandX(
    const std::vector<Eigen::Vector3d>& points_in_mesh);


}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake

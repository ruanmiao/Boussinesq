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

/// Given the 2D Circle area defined by an origin and a radius, with the number
/// of nodes per radius (including the origin)
/// @param o Origin.
/// @param radius Radius.
/// @param num_pr Number of vertices per radius.
/// @returns The first vector in the return is the location of all nodes in
/// the mesh. The second vector in the return is the node indexes of all the
/// triangles in the mesh.
std::pair<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3i>>
MeshCircle(
    const Vector2<double>& o, double radius, int num_pr);

/// Generate the mesh on a subset of the surface of a sphere given the origin
/// and the radius of the sphere, with the sector specifying the size of the
/// meshing area. The number of Vertices per half sector (origin included) is
/// given as well.
/// @param o Origin.
/// @param radius Radius.
/// @param num_pr Number of vertices per radius.
/// @param half_sector Value of half the sector
/// @returns The first vector in the return is the location of all nodes in
/// the mesh. The second vector in the return is the node indexes of all the
/// triangles in the mesh.
std::pair<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3i>>
MeshSphere(
    const Vector3<double> &o, double radius, int num_pr, double half_sector);

/// Help to determine the index of a node in the whole mesh. Only used in the
/// MeshSquare above. The three inputs specified the three iterator indexes in
/// the three nested loops in MeshSquare.
/// @param ir Looping index related to the radius.
/// @param itheta Looping index specifies which part of the circle.
/// @param ipos Looping index specified the exact position.
/// @returns The global node index of the node in the mesh.
int FindIndexInMeshCircle(int ir, int itheta, int ipos);

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

/// Helper function outputs the mesh data to .vtk files for visualization
/// support.
/// @param points_in_mesh Poses of nodes in the mesh
/// @param triangles_in_mesh Indexes of points of triangles in the mesh
/// @param filename (with .vtk extension)
/// @param enable Output to file or disregard this function
/// @returns none
bool OutputMeshToVTK(
    const std::vector<Vector3<double>>& points_in_mesh,
                     const std::vector<Vector3<int>>& triangles_in_mesh,
                     const VectorX<double>& values,
                     bool enabled = true);

/// Helper function computes the total force applied knowing the pressure over
/// the mesh, where the mesh is in a 2D plane
/// @param points_in_mesh Poses of nodes in the mesh
/// @param triangles_in_mesh Indexes of points of triangles in the mesh
/// @param pressure Pressure at each nodes in the mesh
/// @returns Total force
double CalcForceOverMesh(const std::vector<Vector3<double>>& points_in_mesh,
                         const std::vector<Vector3<int>>& triangles_in_mesh,
                         const VectorX<double>& pressure);

/// Helper function computes the total force applied knowing the pressure over
/// the mesh of the surface of a sphere, where the mesh is in the 3D space.
/// @param points_in_mesh Poses of nodes in the mesh
/// @param triangles_in_mesh Indexes of points of triangles in the mesh
/// @param pressure Pressure at each nodes in the mesh
/// @param center Center of the sphere
/// @param z_n normal direction of the half plane
/// @returns Total force
double CalcForceOverMeshOfSphere(
    const std::vector<Vector3<double>>& points_in_mesh,
    const std::vector<Vector3<int>>& triangles_in_mesh,
    const VectorX<double>& pressure,
    const Vector3<double>& center,
    const Vector3<double>& pressure_dir = Vector3<double>::UnitZ());

}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake

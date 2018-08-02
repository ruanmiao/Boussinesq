#pragma once

#include <memory>
#include <vector>

#include "drake/common/drake_assert.h"
#include "drake/common/drake_copyable.h"
#include "drake/common/drake_throw.h"
#include "drake/common/eigen_types.h"
#include "drake/geometry/query_results/penetration_as_triangle_pair.h"

namespace drake {
namespace geometry {
namespace mesh_query {

template <typename T>
struct Mesh {
  // A unique identifier within a model containing seeral meshes.
  int mesh_index{-1};  // Invalid initialization to a negative value.

  // Mesh nodes measured and expressed in a the mesh geometry frame G.
  std::vector<Vector3<double>> points_G;
  // Mesh triangles.
  std::vector<Vector3<int>> triangles;
  // Mesh face normals expressed in the mesh geometry frame G, the same as the
  // nodes.
  std::vector<Vector3<double>> face_normals_G;

  std::vector<Vector3<double>> node_normals_G;

  // A conveniently pre-computed array such that for a node_index then
  // node_element[node_index].first is an element of which node_index is a
  // vertex.
  // node_element[node_index].second tells us the local index of node_index in
  // that element.
  // The element assigned to node_index is arbitrary.
  std::vector<std::pair<int, int>> node_element;
};

std::vector<Vector3<double>> CalcMeshFaceNormals(
    const std::vector<Vector3<double>>& points_A,
    const std::vector<Vector3<int>>& triangles);

std::vector<Vector3<double>> CalcAreaWeightedNormals(const Mesh<double>& mesh);

void FlipNormals(Mesh<double>* mesh);

template <typename T>
struct PointMeshDistance {
  DRAKE_DEFAULT_COPY_AND_MOVE_AND_ASSIGN(PointMeshDistance)

  PointMeshDistance() = default;

  int triangle_index;

  // Distance from query point Q to point on the mesh P.
  double distance;

  // Point on the mesh, measured and expressed in a frame F.
  Vector3<T> p_FP;

  // Normal at point P expressed in frame F.
  Vector3<T> normal_F;

  // Indexes on the mesh to the triangle nodes.
  Vector3<int> triangle;

  // The barycentric coordinates for point P on the mesh.
  Vector3<T> barycentric_P;

};

std::vector<PenetrationAsTrianglePair<double>> MeshToMeshQuery(
    const Isometry3<double>& X_WA, const Mesh<double>& meshA,
    const Isometry3<double>& X_WB, const Mesh<double>& meshB);

/// Computes the signed distance from the query point Q to mesh A.
bool CalcPointToMeshNegativeDistance(
    const Isometry3<double>& X_FA,
    const std::vector<Vector3<double>>& points_A,
    const std::vector<Vector3<int>>& triangles,
    const std::vector<Vector3<double>>& triangle_normals_A,
    const Vector3<double>& p_FQ,
    PointMeshDistance<double>* point_mesh_distance_ptr = nullptr);

}  // namespace mesh_query
}  // namespace geometry
}  // namespace drae

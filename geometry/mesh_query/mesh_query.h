#pragma once

#include <memory>
#include <vector>

#include "drake/common/drake_assert.h"
#include "drake/common/drake_copyable.h"
#include "drake/common/drake_throw.h"
#include "drake/common/eigen_types.h"
#include "drake/geometry/query_results/penetration_as_triangle_pair.h"
#include "drake/geometry/mesh_query/mesh.h"

namespace drake {
namespace geometry {
namespace mesh_query {

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
    const Isometry3<double>& X_WB, const Mesh<double>& meshB, double sigma = 0);

/// Computes the signed distance from the query point Q to mesh A.
/// Returns true if the point is inside the mesh.
bool CalcPointToMeshNegativeDistance(
    const Isometry3<double>& X_FA,
    const std::vector<Vector3<double>>& points_A,
    const std::vector<Vector3<int>>& triangles,
    const std::vector<Vector3<double>>& face_normals_A,
    const std::vector<Vector3<double>>& node_normals_A,
    const Vector3<double>& p_FQ,
    PointMeshDistance<double>* point_mesh_distance_ptr = nullptr);

/// Computes the signed distance from the query point Q to mesh A.
/// IMPORTANT!!!: It assumes the point is outside the mesh.
void CalcPointToMeshPositiveDistance(
    const Isometry3<double>& X_FA,
    const std::vector<Vector3<double>>& points_A,
    const std::vector<Vector3<int>>& triangles,
    const std::vector<Vector3<double>>& face_normals_A,
    const std::vector<Vector3<double>>& node_normals_A,
    const Vector3<double>& p_FQ,
    PointMeshDistance<double>* point_mesh_distance_ptr = nullptr,
        bool debug_mode = false);

std::pair<std::unique_ptr<Mesh<double>>, std::unique_ptr<Mesh<double>>>
MakeLocalPatchMeshes(
    std::vector<PenetrationAsTrianglePair<double>>* pairs,
    const Mesh<double>& meshA, const Mesh<double>& meshB);

}  // namespace mesh_query
}  // namespace geometry
}  // namespace drae

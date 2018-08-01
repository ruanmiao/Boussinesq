#include "drake/geometry/mesh_query/mesh_query.h"

#include <algorithm>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

namespace drake {
namespace geometry {
namespace mesh_query {

std::vector<PenetrationAsTrianglePair<double>> MeshToMeshQuery(
    const Isometry3<double>& X_FA,
    const std::vector<Vector3<double>>& meshA_points_A,
    const std::vector<Vector3<int>>& meshA_triangles,
    const Isometry3<double>& X_FB,
    const std::vector<Vector3<double>>& meshB_points_B,
    const std::vector<Vector3<int>>& meshB_triangles) {
  std::vector<PenetrationAsTrianglePair<double>> pairs;
  return pairs;
}

std::vector<Vector3<double>> CalcMeshFaceNormals(
    const std::vector<Vector3<double>>& points_A,
    const std::vector<Vector3<int>>& triangles) {
  const int num_triangles = triangles.size();

  std::vector<Vector3<double>> normals;
  normals.reserve(num_triangles);

  for (int triangle_index = 0;
       triangle_index < num_triangles; ++triangle_index) {
    const Vector3<int>& triangle = triangles[triangle_index];

    const Vector3<double>& p_AP1 = points_A[triangle[0]];
    const Vector3<double>& p_AP2 = points_A[triangle[1]];
    const Vector3<double>& p_AP3 = points_A[triangle[2]];

    // The orientation of the normal follows the right-hand rule.
    // It is expressed in the same frame in which the mesh points are expressed.
    const Vector3<double> u1_A = p_AP2 - p_AP1;
    const Vector3<double> u2_A = p_AP3 - p_AP1;
    const Vector3<double> normal_A = u1_A.cross(u2_A);

    normals.push_back(normal_A.normalized());
  }

  return normals;
}

bool CalcPointToMeshNegativeDistance(
    const Isometry3<double>& X_FA,
    const std::vector<Vector3<double>>& points_A,
    const std::vector<Vector3<int>>& triangles,
    const std::vector<Vector3<double>>& triangle_normals_A,
    const Vector3<double>& p_FQ,
    PointMeshDistance<double>* point_mesh_distance_ptr) {
  using std::abs;
  using std::min;

  // Point Q measured and expressed in the mesh frame A.
  const Vector3<double> p_AQ = X_FA.inverse() * p_FQ;

  // P1, P2 and P3 are counter-clockwise around plane_normal_A, which defines
  // the "positive" orientation for a triangle's area.
  auto CalcTriangleArea = [](
      const Vector3<double>& p_AP1,
      const Vector3<double>& p_AP2,
      const Vector3<double>& p_AP3,
      const Vector3<double>& plane_normal_A) {
    // The orientation of the normal follows the right-hand rule.
    // It is expressed in the same frame in which the mesh points are expressed.
    const Vector3<double> u1_A = p_AP2 - p_AP1;
    const Vector3<double> u2_A = p_AP3 - p_AP1;
    const Vector3<double> area_vector_A = u1_A.cross(u2_A);
    // By dotting with the plane's normal we get the appropriate sign for the
    // are.
    return plane_normal_A.dot(area_vector_A) / 2.0;
  };

  // Results from the first scan over triangles:
  double max_neg_dist = -std::numeric_limits<double>::infinity();
  int max_neg_dist_triangle = -1;
  Vector3<double> p_AP;
  double A1, A2, A3;

  for (size_t triangle_index = 0;
       triangle_index < triangles.size(); ++triangle_index) {
    const Vector3<int>& triangle = triangles[triangle_index];
    const Vector3<double>& normal_A = triangle_normals_A[triangle_index];

    const Vector3<double>& p_AP1 = points_A[triangle[0]];
    const Vector3<double>& p_AP2 = points_A[triangle[1]];
    const Vector3<double>& p_AP3 = points_A[triangle[2]];

    // Distance from point Q to the plane on which the triangle lies.
    const double plane_distance = normal_A.dot(p_AQ - p_AP1);

    // point is outside convex mesh. Thus we are done.
    if (plane_distance > 0) return false;

    // Save the triangle with the minimum distance.
    if (plane_distance > max_neg_dist) {
      // point on the triangle's plane.
      const Vector3<double> p_AP_local = p_AQ - plane_distance * normal_A;
      const double A1_local = CalcTriangleArea(p_AP2, p_AP3, p_AP_local, normal_A);
      if (A1_local < 0) continue; // Outside triangle.
      const double A2_local = CalcTriangleArea(p_AP3, p_AP1, p_AP_local, normal_A);
      if (A2_local < 0) continue; // Outside triangle.
      const double A3_local = CalcTriangleArea(p_AP1, p_AP2, p_AP_local, normal_A);
      if (A3_local < 0) continue; // Outside triangle.

      // If here, the projection lies inside the triangle and therefore we save
      // the data that we already computed for later.

      max_neg_dist = plane_distance;
      max_neg_dist_triangle = triangle_index;
      A1 = A1_local;
      A2 = A2_local;
      A3 = A3_local;
      p_AP = p_AP_local;
    }
  }

  // If we got here is because all distances are negative and we are inside.
  // Sanity check that.
  DRAKE_DEMAND(max_neg_dist_triangle >= 0);
  DRAKE_DEMAND(max_neg_dist <= 0);

  if (point_mesh_distance_ptr == nullptr) return true;

  //////////////////////////////////////////////////////////////////
  // We know we are inside. Compute additional information and exit:
  //////////////////////////////////////////////////////////////////

  // The signed distance closest to the boundary.
  const double distance = max_neg_dist;

  // The closest plane does correspond to the triangle when inside.
  const int triangle_index = max_neg_dist_triangle;

  // normal on the mesh
  const Vector3<double>& normal_A = triangle_normals_A[triangle_index];

  // point on the mesh.
  //const Vector3<double> p_AP = p_AQ - distance * normal_A;

  // Triangle indexes.
  const Vector3<int>& triangle = triangles[triangle_index];

  // The three areas must be positive when inside a convex mesh.
  DRAKE_DEMAND(A1 >= 0);
  DRAKE_DEMAND(A2 >= 0);
  DRAKE_DEMAND(A3 >= 0);

  // Total triangle area
  const double area = A1 + A2 + A3;

  const Vector3<double> barycentric_P(A1 / area, A2 / area, A3 / area);

  PointMeshDistance<double>& point_mesh_distance = *point_mesh_distance_ptr;

  point_mesh_distance.triangle_index = triangle_index;
  point_mesh_distance.distance = distance;
  point_mesh_distance.p_FP = X_FA * p_AP;
  point_mesh_distance.normal_F = X_FA * normal_A;
  point_mesh_distance.triangle = triangle;
  point_mesh_distance.barycentric_P = barycentric_P;

  // Confirm point is inside the mesh.
  return true;
}

}  // namespace mesh_query
}  // namespace geometry
}  // namespace drake
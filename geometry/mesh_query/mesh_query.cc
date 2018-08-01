#include "drake/geometry/mesh_query/mesh_query.h"

#include <algorithm>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

namespace drake {
namespace geometry {
namespace mesh_query {

// Tolerance used to test when the distance is close to zero. This number is
// made sligthly larger than zero to avoid divission by zero when computing a
// normalized vector between a query point Q and its projection P on the surface
// of a mesh.
const double kNearSurfaceTolerance = 10 *
    std::numeric_limits<double>::epsilon();

std::vector<PenetrationAsTrianglePair<double>> MeshToMeshQuery(
    const Isometry3<double>& X_WA, const Mesh<double>& meshA,
    const Isometry3<double>& X_WB, const Mesh<double>& meshB) {
  std::vector<PenetrationAsTrianglePair<double>> pairs;

  pairs.clear();

  auto Mesh1NodesVsMesh2Surface = [&pairs](
      const Isometry3<double>& X_WM1, const Mesh<double>& mesh1,
      const Isometry3<double>& X_WM2, const Mesh<double>& mesh2) {
    for (size_t node_index = 0; node_index < mesh1.points_G.size(); ++node_index) {
      const Vector3<double>& p_AQ = mesh1.points_G[node_index];
      const Vector3<double> p_WQ = X_WM1 * p_AQ;

      PointMeshDistance<double> point_mesh_result;

      const bool is_inside = CalcPointToMeshNegativeDistance(
          X_WM2, mesh2.points_G, mesh2.triangles, mesh2.face_normals_G, p_WQ,
          &point_mesh_result);

      if (is_inside) {
        PenetrationAsTrianglePair<double> result;

        result.signed_distance = point_mesh_result.distance;

        //////////////////////////////////////////////////////////////////////////
        // MESH A INFO
        //////////////////////////////////////////////////////////////////////////

        result.meshA_index = mesh1.mesh_index;
        result.triangle_A = mesh1.node_element[node_index].first;
        result.p_WoAs_W = p_WQ;

        const Vector3<double> p_WP = point_mesh_result.p_FP;
        const double distance = point_mesh_result.distance;

        // For now we are assuming the distance is zero. Therefore verify this.
        DRAKE_DEMAND(distance <= -kNearSurfaceTolerance);

        result.normal_A_W = (p_WQ - p_WP) / (-distance);

        // Since we assume that node_element points to local node "zero" in the
        // mesh A triangle, the barycentric coordinates are 1, 0, 0.
        result.barycentric_A =
            Vector3<double>::Unit(mesh1.node_element[node_index].second);

        //////////////////////////////////////////////////////////////////////////
        // MESH B INFO
        //////////////////////////////////////////////////////////////////////////
        result.meshB_index = mesh2.mesh_index;
        result.triangle_B = point_mesh_result.triangle_index;
        result.p_WoBs_W = point_mesh_result.p_FP;
        result.barycentric_B = point_mesh_result.barycentric_P;
        result.normal_B_W = point_mesh_result.normal_F;

        pairs.push_back(result);
      }
    }

  };

  // Scan each node on Mesh A and perform a point-mesh distance query with
  // mesh B.
  Mesh1NodesVsMesh2Surface(X_WA, meshA, X_WB, meshB);

  // Reverse roles of mesh A and B. Now scan each node on Mesh B and perform a
  // point-mesh distance query with mesh A.
  Mesh1NodesVsMesh2Surface(X_WB, meshB, X_WA, meshA);

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
    // The check is made against a small tolerance to avoid zero distances.
    if (plane_distance > -kNearSurfaceTolerance) return false;

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
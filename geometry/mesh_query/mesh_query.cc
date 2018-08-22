#include "drake/geometry/mesh_query/mesh_query.h"
#include "drake/multibody/shapes/geometry.h"

#include <algorithm>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

#include <iostream>
#define PRINT_VAR(a) std::cout << #a": " << a << std::endl;

namespace drake {
namespace geometry {
namespace mesh_query {

Vector3<double> CalcClosestPointOnTriangle(
    const Vector3<double>& p,
    const Vector3<double>& a,
    const Vector3<double>& b,
    const Vector3<double>& c);

Vector3<double> CalcClosestPointOnTriangleBarycentricCoordinates(
    const Vector3<double>& p,
    const Vector3<double>& a,
    const Vector3<double>& b,
    const Vector3<double>& c);

// Tolerance used to test when the distance is close to zero. This number is
// made sligthly larger than zero to avoid divission by zero when computing a
// normalized vector between a query point Q and its projection P on the surface
// of a mesh.
const double kNearSurfaceTolerance = 10 *
    std::numeric_limits<double>::epsilon();

std::vector<PenetrationAsTrianglePair<double>> MeshToMeshQuery(
    const Isometry3<double>& X_WA, const Mesh<double>& meshA,
    const Isometry3<double>& X_WB, const Mesh<double>& meshB,
    double sigma) {
  std::vector<PenetrationAsTrianglePair<double>> pairs;

  pairs.clear();

  auto Mesh1NodesVsMesh2Surface = [&pairs, sigma](
      const Isometry3<double>& X_WM1, const Mesh<double>& mesh1,
      const Isometry3<double>& X_WM2, const Mesh<double>& mesh2) {
    for (size_t node_index = 0; node_index < mesh1.points_G.size(); ++node_index) {
      const Vector3<double>& p_AQ = mesh1.points_G[node_index];
      const Vector3<double> p_WQ = X_WM1 * p_AQ;

      PointMeshDistance<double> point_mesh_result;

      const bool is_inside = CalcPointToMeshNegativeDistance(
          X_WM2, mesh2.points_G, mesh2.triangles,
          mesh2.face_normals_G, mesh2.node_normals_G,
          p_WQ,
          &point_mesh_result);

      if (!is_inside) {
        CalcPointToMeshPositiveDistance(
            X_WM2, mesh2.points_G, mesh2.triangles,
            mesh2.face_normals_G, mesh2.node_normals_G,
            p_WQ,
            &point_mesh_result);
      }

      bool is_inside_sigma = point_mesh_result.distance < sigma;

      if (is_inside_sigma) {
        PenetrationAsTrianglePair<double> result;

        result.signed_distance = point_mesh_result.distance;

        //////////////////////////////////////////////////////////////////////////
        // MESH A INFO
        //////////////////////////////////////////////////////////////////////////

        result.meshA_index = mesh1.mesh_index;
        result.triangleA_index = mesh1.node_element[node_index].first;
        result.triangleA = mesh1.triangles[result.triangleA_index];
        result.p_WoAs_W = p_WQ;

        //const Vector3<double> p_WP = point_mesh_result.p_FP;
        //const double distance = point_mesh_result.distance;

        // For now we are assuming the distance is zero. Therefore verify this.
        //DRAKE_DEMAND(distance <= -kNearSurfaceTolerance);

        // Frame G is the frame of the mesh1. Therefore we must transform to
        // the world frame.
        result.normal_A_W = X_WM1.linear() * mesh1.node_normals_G[node_index];

        // Since we assume that node_element points to local node "zero" in the
        // mesh A triangle, the barycentric coordinates are 1, 0, 0.
        result.barycentric_A =
            Vector3<double>::Unit(mesh1.node_element[node_index].second);

        //////////////////////////////////////////////////////////////////////////
        // MESH B INFO
        //////////////////////////////////////////////////////////////////////////
        result.meshB_index = mesh2.mesh_index;
        result.triangleB_index = point_mesh_result.triangle_index;
        result.triangleB = mesh2.triangles[result.triangleB_index];
        result.p_WoBs_W = point_mesh_result.p_FP;
        result.barycentric_B = point_mesh_result.barycentric_P;
        // Frame F IS the world frame W on output from
        // CalcPointToMeshNegativeDistance().
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
  // Mesh1NodesVsMesh2Surface(X_WB, meshB, X_WA, meshA);

  return pairs;
}

std::pair<std::unique_ptr<Mesh<double>>, std::unique_ptr<Mesh<double>>>
MakeLocalPatchMeshes(
    std::vector<PenetrationAsTrianglePair<double>>* pairs,
    const Mesh<double>& meshA, const Mesh<double>& meshB) {
  auto meshA_patch = std::make_unique<Mesh<double>>();
  auto meshB_patch = std::make_unique<Mesh<double>>();

  meshA_patch->mesh_index = meshA.mesh_index;
  meshB_patch->mesh_index = meshB.mesh_index;

  std::set<int> patchA_triangles;
  std::set<int> patchB_triangles;

  std::set<int> patchA_nodes;
  std::set<int> patchB_nodes;

  auto InsertTriangle = [](
      int triangle_index, const Vector3<int>& triangle,
      std::set<int>* patch_triangles,
      std::set<int>* patch_nodes) {
    patch_triangles->insert(triangle_index);
    for (int i = 0; i < 3; ++i) {
      patch_nodes->insert(triangle[i]);
    }
  };

  auto InsertNodeAndAdjacentTriangles = [InsertTriangle](
      int node_index, const std::vector<int>& node_triangles,
      const std::vector<Vector3<int>>& mesh_triangles,
      std::set<int>* patch_triangles,
      std::set<int>* patch_nodes) {
    patch_nodes->insert(node_index);
    for (int triangle_index : node_triangles) {
      InsertTriangle(triangle_index, mesh_triangles[triangle_index],
                     patch_triangles, patch_nodes);
    }
  };

  auto InsertTriangleAndAdjacentTriangles = [InsertNodeAndAdjacentTriangles](
      int triangle_index, const Vector3<int>& triangle,
      const std::vector<Vector3<int>>& mesh_triangles,
      const std::vector<std::vector<int>>& nodes_triangles,
      std::set<int>* patch_triangles,
      std::set<int>* patch_nodes) {
    patch_triangles->insert(triangle_index);

    for (int i = 0; i < 3; ++i) {
      const int node_index = triangle[i];
      const auto& node_triangles = nodes_triangles[node_index];
      InsertNodeAndAdjacentTriangles(
          node_index, node_triangles, mesh_triangles,
          patch_triangles, patch_nodes);
    }
  };

  // Crete the set of triangles in the patch for each mesh.
  // 1) First add the triangles directly referenced by th query pairs.
  for (const auto& pair : *pairs) {
    int triangle_index = pair.triangleA_index;
    if (pair.meshA_index == meshA.mesh_index ) {
      const auto& triangle = meshA.triangles[triangle_index];
      InsertTriangleAndAdjacentTriangles(
          triangle_index, triangle,
          meshA.triangles, meshA.node_triangles,
          &patchA_triangles, &patchA_nodes);
    } else {
      const auto& triangle = meshB.triangles[triangle_index];
      InsertTriangleAndAdjacentTriangles(
          triangle_index, triangle,
          meshB.triangles, meshB.node_triangles,
          &patchB_triangles, &patchB_nodes);
    }

    triangle_index = pair.triangleB_index;
    if (pair.meshB_index == meshA.mesh_index) {
      const auto& triangle = meshA.triangles[triangle_index];
      InsertTriangleAndAdjacentTriangles(
          triangle_index, triangle,
          meshA.triangles, meshA.node_triangles,
          &patchA_triangles, &patchA_nodes);
    } else {
      const auto& triangle = meshB.triangles[triangle_index];
      InsertTriangleAndAdjacentTriangles(
          triangle_index, triangle,
          meshB.triangles, meshB.node_triangles,
          &patchB_triangles, &patchB_nodes);
    }
  }

  auto ConvertPatchSetsToMesh = [](
      const std::set<int>& patch_nodes,
      const std::set<int>& patch_triangles,
      const Mesh<double>& mesh,
      Mesh<double>* patch_mesh) {
    std::vector<int> triangles_map(mesh.triangles.size(), -1);

    std::vector<int> nodes_map(mesh.points_G.size(), -1);
    patch_mesh->node_areas.resize(mesh.points_G.size());
    for (int node_index : patch_nodes) {
      // Define the local patch index to node_index
      const int local_node_index = patch_mesh->points_G.size();
      nodes_map[node_index] = local_node_index;
      patch_mesh->points_G.push_back(mesh.points_G[node_index]);
      patch_mesh->node_areas[local_node_index] = mesh.node_areas[node_index];

      patch_mesh->node_normals_G.push_back(mesh.node_normals_G[node_index]);
      patch_mesh->node_triangles.push_back(mesh.node_triangles[node_index]);
    }

    for (int triangle_index : patch_triangles) {
      const auto& triangle = mesh.triangles[triangle_index];

      DRAKE_DEMAND(nodes_map[triangle[0]] >= 0);
      DRAKE_DEMAND(nodes_map[triangle[1]] >= 0);
      DRAKE_DEMAND(nodes_map[triangle[2]] >= 0);

      const Vector3<int> pach_triangle(
          nodes_map.at(triangle[0]),
          nodes_map.at(triangle[1]),
          nodes_map.at(triangle[2]));

      triangles_map[triangle_index] = patch_mesh->triangles.size();
      patch_mesh->triangles.push_back(pach_triangle);
      patch_mesh->face_normals_G.push_back(mesh.face_normals_G[triangle_index]);
    }

    return triangles_map;
  };

  auto patchA_triangles_map =
      ConvertPatchSetsToMesh(patchA_nodes, patchA_triangles, meshA,
                         meshA_patch.get());

  auto patchB_triangles_map =
      ConvertPatchSetsToMesh(patchB_nodes, patchB_triangles, meshB,
                         meshB_patch.get());

  // Convert indexes from the global indexing in the original meshes to local
  // patch indexing.
  for (auto& pair : *pairs) {
    pair.triangleA_index = pair.meshA_index == meshA.mesh_index ?
                      patchA_triangles_map.at(pair.triangleA_index) :
                      patchB_triangles_map.at(pair.triangleA_index);
    DRAKE_DEMAND(pair.triangleA_index >= 0);

    if (pair.meshA_index == meshA.mesh_index) {
      pair.triangleA = meshA_patch->triangles[pair.triangleA_index];
    } else {
      pair.triangleA = meshB_patch->triangles[pair.triangleA_index];
    }

    pair.triangleB_index = pair.meshB_index == meshA.mesh_index ?
                      patchA_triangles_map.at(pair.triangleB_index) :
                      patchB_triangles_map.at(pair.triangleB_index);

    if (pair.meshB_index == meshA.mesh_index) {
      pair.triangleB = meshA_patch->triangles[pair.triangleB_index];
    } else {
      pair.triangleB = meshB_patch->triangles[pair.triangleB_index];
    }

    DRAKE_DEMAND(pair.triangleB_index >= 0);
  }

  return std::make_pair(std::move(meshA_patch), std::move(meshB_patch));
};

void CalcPointToMeshPositiveDistance(
    const Isometry3<double>& X_FA,
    const std::vector<Vector3<double>>& points_A,
    const std::vector<Vector3<int>>& triangles,
    const std::vector<Vector3<double>>& face_normals_A,
    const std::vector<Vector3<double>>& node_normals_A,
    const Vector3<double>& p_FQ,
    PointMeshDistance<double>* point_mesh_distance_ptr) {
  using std::abs;
  using std::min;

  const int num_nodes = points_A.size();
  const int num_elements = triangles.size();
  DRAKE_DEMAND(static_cast<int>(face_normals_A.size()) == num_elements);
  DRAKE_DEMAND(static_cast<int>(node_normals_A.size()) == num_nodes);

  // Point Q measured and expressed in the mesh frame A.
  const Vector3<double> p_AQ = X_FA.inverse() * p_FQ;

  double min_distance2 = std::numeric_limits<double>::infinity();
  int min_distance_triangle = -1;
  Vector3<double> min_dist_barycentric;
  Vector3<double> min_dist_p_AP;

  // "Brute force" scans all triangles and computes the distance to each.
  for (size_t triangle_index = 0;
       triangle_index < triangles.size(); ++triangle_index) {
    const Vector3<int>& triangle = triangles[triangle_index];

    const Vector3<double>& p_AP1 = points_A[triangle[0]];
    const Vector3<double>& p_AP2 = points_A[triangle[1]];
    const Vector3<double>& p_AP3 = points_A[triangle[2]];

    const Vector3<double> barycentric_P =
        CalcClosestPointOnTriangleBarycentricCoordinates(p_AQ, p_AP1, p_AP2, p_AP3);

    const Vector3<double> p_AP =
        barycentric_P[0] * p_AP1 +
        barycentric_P[1] * p_AP2 +
        barycentric_P[2] * p_AP3;

    DRAKE_DEMAND(barycentric_P[0] >=0);
    DRAKE_DEMAND(barycentric_P[1] >=0);
    DRAKE_DEMAND(barycentric_P[2] >=0);
    //PRINT_VAR(barycentric_P[0]+barycentric_P[1]+barycentric_P[2]);
    DRAKE_DEMAND(std::abs(barycentric_P[0]+barycentric_P[1]+barycentric_P[2] - 1.0) < 5e-15);



    const double distance2 = (p_AP-p_AQ).squaredNorm();

    if (distance2 < min_distance2) {
      min_distance2 = distance2;
      min_distance_triangle = triangle_index;
      min_dist_barycentric = barycentric_P;
      min_dist_p_AP = p_AP;
    }
  }

  // Triangle indexes.
  const Vector3<int>& triangle = triangles[min_distance_triangle];

  // Interpolate normal to point on P on the mesh.
  const Vector3<double>& normal_P1_A = node_normals_A[triangle[0]];
  const Vector3<double>& normal_P2_A = node_normals_A[triangle[1]];
  const Vector3<double>& normal_P3_A = node_normals_A[triangle[2]];

  const Vector3<double> normal_A =
      min_dist_barycentric[0] * normal_P1_A +
      min_dist_barycentric[1] * normal_P2_A +
      min_dist_barycentric[2] * normal_P3_A;

  PointMeshDistance<double>& point_mesh_distance = *point_mesh_distance_ptr;

  point_mesh_distance.triangle_index = min_distance_triangle;
  point_mesh_distance.distance = sqrt(min_distance2);
  point_mesh_distance.p_FP = X_FA * min_dist_p_AP;
  point_mesh_distance.normal_F = X_FA.linear() * normal_A;
  point_mesh_distance.triangle = triangle;
  point_mesh_distance.barycentric_P = min_dist_barycentric;
}

// Returns true if the point is inside the convex mesh.
bool CalcPointToMeshNegativeDistance(
    const Isometry3<double>& X_FA,
    const std::vector<Vector3<double>>& points_A,
    const std::vector<Vector3<int>>& triangles,
    const std::vector<Vector3<double>>& face_normals_A,
    const std::vector<Vector3<double>>& node_normals_A,
    const Vector3<double>& p_FQ,
    PointMeshDistance<double>* point_mesh_distance_ptr) {
  using std::abs;
  using std::min;

  const int num_nodes = points_A.size();
  const int num_elements = triangles.size();
  DRAKE_DEMAND(static_cast<int>(face_normals_A.size()) == num_elements);
  DRAKE_DEMAND(static_cast<int>(node_normals_A.size()) == num_nodes);

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
  double A1(0), A2(0), A3(0);

  for (size_t triangle_index = 0;
       triangle_index < triangles.size(); ++triangle_index) {
    const Vector3<int>& triangle = triangles[triangle_index];
    const Vector3<double>& face_normal_A = face_normals_A[triangle_index];

    const Vector3<double>& p_AP1 = points_A[triangle[0]];
    const Vector3<double>& p_AP2 = points_A[triangle[1]];
    const Vector3<double>& p_AP3 = points_A[triangle[2]];

    // Distance from point Q to the plane on which the triangle lies.
    const double plane_distance = face_normal_A.dot(p_AQ - p_AP1);

    // point is outside convex mesh. Thus we are done.
    // The check is made against a small tolerance to avoid zero distances.
    if (plane_distance > -kNearSurfaceTolerance) return false;

    // Save the triangle with the minimum distance.
    if (plane_distance > max_neg_dist) {
      // point on the triangle's plane.
      const Vector3<double> p_AP_local = p_AQ - plane_distance * face_normal_A;
      const double A1_local = CalcTriangleArea(p_AP2, p_AP3, p_AP_local, face_normal_A);
      if (A1_local < 0) continue; // Outside triangle.
      const double A2_local = CalcTriangleArea(p_AP3, p_AP1, p_AP_local, face_normal_A);
      if (A2_local < 0) continue; // Outside triangle.
      const double A3_local = CalcTriangleArea(p_AP1, p_AP2, p_AP_local, face_normal_A);
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
  //const Vector3<double>& normal_A = face_normals_A[triangle_index];

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

  // Interpolate normal to point on P on the mesh.
  const Vector3<double>& normal_P1_A = node_normals_A[triangle[0]];
  const Vector3<double>& normal_P2_A = node_normals_A[triangle[1]];
  const Vector3<double>& normal_P3_A = node_normals_A[triangle[2]];

  const Vector3<double> normal_A =
      barycentric_P[0] * normal_P1_A +
      barycentric_P[1] * normal_P2_A +
      barycentric_P[2] * normal_P3_A;

  PointMeshDistance<double>& point_mesh_distance = *point_mesh_distance_ptr;

  point_mesh_distance.triangle_index = triangle_index;
  point_mesh_distance.distance = distance;
  point_mesh_distance.p_FP = X_FA * p_AP;
  point_mesh_distance.normal_F = X_FA.linear() * normal_A;
  point_mesh_distance.triangle = triangle;
  point_mesh_distance.barycentric_P = barycentric_P;

  // Confirm point is inside the mesh.
  return true;
}

void FlipNormals(Mesh<double>* mesh) {
  for (auto& triangle : mesh->triangles) {
    std::swap(triangle[0], triangle[1]);
  }
  for (auto& normal : mesh->face_normals_G) {
    normal = -normal;
  }
  for (auto& normal : mesh->node_normals_G) {
    normal = -normal;
  }
}

/// This method computes the closest point P' on the triangle abc to point p.
/// Algorithm in Section 5.1.5 of Ericson, C., 2004. Real-time collision
/// detection. CRC Press.
Vector3<double> CalcClosestPointOnTriangle(
    const Vector3<double>& p,
    const Vector3<double>& a,
    const Vector3<double>& b,
    const Vector3<double>& c) {
  const Vector3<double> ab = b - a;
  const Vector3<double> ac = c - a;
  const Vector3<double> bc = c - b;

  // Compute parametric position s for projetion P' of P on AB,
  // P' = A + s*AB, s = snom/(snom+sdenom).
  const double snom = ab.dot(p-a);
  const double sdenom = (p-b).dot(a-b);

  // Compute parametric position t for projection P' of P on AC,
  // P' = A + t*AC, s = tnom/(tnom+tdenom).
  const double tnom = ac.dot(p-a);
  const double tdenom = (p-c).dot(a-c);

  if (snom <= 0.0 && tnom <= 0.0) return a;  // Vertex region early out.

  // Compute parametric position u for projection P' of P on BC,
  // P' = B + u*BC, u = unom/(unom+udenom).
  const double unom = bc.dot(p-b);
  const double udenom = (p-c).dot(b-c);

  if (sdenom <= 0.0 && unom <= 0.0) return b;    // Vertex region early out.
  if (tdenom <= 0.0 && udenom <= 0.0) return c;  // Vertex region early out.

  // P is outside (or on) AB if the triple product [N PA PB] <= 0.
  const Vector3<double> n = (b-a).cross(c-a);
  const double vc = n.dot((a-p).cross(b-p));

  // If P is outside AB and within feature region of AB, return the projection
  // of P onto AB.
  if (vc <= 0.0 && snom >= 0.0 && sdenom >= 0.0)
    return a + snom/(snom+sdenom) * ab;

  // P is outside (or on) BC if the triple product [N PB PC] <= 0.
  const double va = n.dot((b-p).cross(c-p));

  // If P is outside BC and within feature region of BC,
  // return projection of P onto BC.
  if (va <= 0.0 && unom >= 0.0 && udenom >= 0.0)
    return b + unom/(unom+udenom)*bc;

  // P is outside (or on) CA if the triple product [N PC PA] <= 0.
  const double vb = n.dot((c-p).cross(a-p));

  // If P is outside CA and within feature region of CA, return the projection
  // of P onto CA.
  if (vb <= 0 && tnom >= 0.0 && tdenom >= 0.0)
    return a + tnom/(tnom+tdenom)*ac;

  // P must project inside the face region. Compute Q using barycentric
  // coordinates.
  const double u = va/(va+vb+vc);
  const double v = vb/(va+vb+vc);
  const double w = 1.0 - u - v;  // = vc//(va+vb+vc).
  return u * a + v * b + w * c;
}

Vector3<double> CalcClosestPointOnTriangleBarycentricCoordinates(
    const Vector3<double>& p,
    const Vector3<double>& a,
    const Vector3<double>& b,
    const Vector3<double>& c) {
  const Vector3<double> ab = b - a;
  const Vector3<double> ac = c - a;
  const Vector3<double> bc = c - b;

  // Compute parametric position s for projetion P' of P on AB,
  // P' = A + s*AB, s = snom/(snom+sdenom).
  const double snom = ab.dot(p-a);
  const double sdenom = (p-b).dot(a-b);

  // Compute parametric position t for projection P' of P on AC,
  // P' = A + t*AC, s = tnom/(tnom+tdenom).
  const double tnom = ac.dot(p-a);
  const double tdenom = (p-c).dot(a-c);

  if (snom <= 0.0 && tnom <= 0.0) return Vector3<double>(1.0, 0.0, 0.0);  // Vertex region early out.

  // Compute parametric position u for projection P' of P on BC,
  // P' = B + u*BC, u = unom/(unom+udenom).
  const double unom = bc.dot(p-b);
  const double udenom = (p-c).dot(b-c);

  if (sdenom <= 0.0 && unom <= 0.0) return Vector3<double>(0.0, 1.0, 0.0);    // Vertex region early out.
  if (tdenom <= 0.0 && udenom <= 0.0) return Vector3<double>(0.0, 0.0, 1.0);  // Vertex region early out.

  // P is outside (or on) AB if the triple product [N PA PB] <= 0.
  const Vector3<double> n = (b-a).cross(c-a);
  const double vc = n.dot((a-p).cross(b-p));

  // If P is outside AB and within feature region of AB, return the projection
  // of P onto AB.
  if (vc <= 0.0 && snom >= 0.0 && sdenom >= 0.0)
    return Vector3<double>(sdenom/(snom+sdenom), snom/(snom+sdenom), 0.0);

  // P is outside (or on) BC if the triple product [N PB PC] <= 0.
  const double va = n.dot((b-p).cross(c-p));

  // If P is outside BC and within feature region of BC,
  // return projection of P onto BC.
  if (va <= 0.0 && unom >= 0.0 && udenom >= 0.0)
    return Vector3<double>(0.0, udenom/(unom+udenom), unom/(unom+udenom));

  // P is outside (or on) CA if the triple product [N PC PA] <= 0.
  const double vb = n.dot((c-p).cross(a-p));

  // If P is outside CA and within feature region of CA, return the projection
  // of P onto CA.
  if (vb <= 0 && tnom >= 0.0 && tdenom >= 0.0)
    return Vector3<double>(tnom/(tnom+tdenom), 0.0, tdenom/(tnom+tdenom));

  // P must project inside the face region. Compute Q using barycentric
  // coordinates.
  const double u = va/(va+vb+vc);
  const double v = vb/(va+vb+vc);
  const double w = 1.0 - u - v;  // = vc//(va+vb+vc).
  return Vector3<double>(u, v, w);
}

}  // namespace mesh_query
}  // namespace geometry
}  // namespace drake
#include "drake/geometry/mesh_query/mesh.h"

#include <fstream>
#include <memory>
#include <vector>

#include <fmt/format.h>
#include <fmt/ostream.h>

#include "drake/common/drake_assert.h"
#include "drake/common/drake_copyable.h"
#include "drake/common/drake_throw.h"
#include "drake/common/eigen_types.h"
#include "drake/common/find_resource.h"
#include "drake/multibody/shapes/geometry.h"

namespace drake {
namespace geometry {
namespace mesh_query {

double CalcTriangleArea(
    const Vector3<double>& p_AP1,
    const Vector3<double>& p_AP2,
    const Vector3<double>& p_AP3) {
  // The orientation of the normal follows the right-hand rule.
  // It is expressed in the same frame in which the mesh points are expressed.
  const Vector3<double> u1_A = p_AP2 - p_AP1;
  const Vector3<double> u2_A = p_AP3 - p_AP1;
  const Vector3<double> area_vector_A = u1_A.cross(u2_A);
  // By dotting with the plane's normal we get the appropriate sign for the
  // are.
  return area_vector_A.norm() / 2.0;
};

std::unique_ptr<Mesh<double>> LoadMeshFromObj(
    const std::string& file_name, bool flip_normals) {
  const auto resource_name = FindResourceOrThrow(file_name);
  DrakeShapes::Mesh mesh_loader(resource_name, resource_name);

  auto mesh = std::make_unique<Mesh<double>>();
  mesh_loader.LoadObjFile(&mesh->points_G, &mesh->triangles);

  // Compute normals.
  mesh->face_normals_G = CalcMeshFaceNormals(mesh->points_G, mesh->triangles);
  mesh->node_normals_G = CalcAreaWeightedNormals(*mesh);

  // Flipping normals changes the indexing. There we MUST flip normals BEFORE
  // computing Mesh::node_element.
  if (flip_normals) FlipNormals(mesh.get());

  const int num_points = mesh->points_G.size();
  // Allocate and initialize to invalid index values.
  mesh->node_element.resize(num_points, std::make_pair(-1, -1));

  // Arbitrarily fill in mesh->node_element.
  const int num_elements = mesh->triangles.size();
  for (int element_index = 0; element_index < num_elements; ++element_index) {
    const auto& triangle = mesh->triangles[element_index];
    for (int i = 0; i < 3; ++i) {
      const int node_index = triangle[i];
      if (mesh->node_element[node_index].first < 0) {  // not yet initialized.
        mesh->node_element[node_index].first = element_index;
        mesh->node_element[node_index].second = i;
      }
    }
  }

  // Fill in the inverse of the connectivities.
  mesh->node_triangles.resize(num_points);
  for (int element_index = 0; element_index < num_elements; ++element_index) {
    const auto &triangle = mesh->triangles[element_index];
    for (int i = 0; i < 3; ++i) {
      const int node_index = triangle[i];
      mesh->node_triangles[node_index].push_back(element_index);
    }
  }

  return mesh;
}

std::vector<Vector3<double>> CalcAreaWeightedNormals(
    const Mesh<double>& mesh) {
  const auto& points_F = mesh.points_G;
  const auto& triangles = mesh.triangles;
  const auto& face_normals_F = mesh.face_normals_G;

  const int num_nodes = points_F.size();
  const int num_triangles = triangles.size();

  // Ensure normals were already computed.
  DRAKE_DEMAND(static_cast<int>(face_normals_F.size()) == num_triangles);

  // Compute the area for each triangle in the mesh.
  VectorX<double> triangle_area(num_triangles);
  for (int ie = 0; ie < num_triangles; ++ie) {
    const auto& triangle = triangles[ie];
    const Vector3<double>& p_FP1 = points_F[triangle[0]];
    const Vector3<double>& p_FP2 = points_F[triangle[1]];
    const Vector3<double>& p_FP3 = points_F[triangle[2]];
    triangle_area[ie] = CalcTriangleArea(p_FP1, p_FP2, p_FP3);
  }

  // Compute a version of "lumped" mass matrix by simply adding the area of
  // all triangles adjacent to a node.
  VectorX<double> node_area = VectorX<double>::Zero(num_nodes);
  std::vector<Vector3<double>> rhs(num_nodes, Vector3<double>::Zero());
  for (int triangle_index = 0;
       triangle_index < num_triangles; ++triangle_index) {
    const auto& triangle = triangles[triangle_index];
    const double area = triangle_area[triangle_index];
    for (int local_index = 0; local_index < 3; ++local_index) {
      int node_index = triangle[local_index];
      node_area[node_index] += area;
      rhs[node_index] += face_normals_F[triangle_index] * area;
    }
  }

  std::vector<Vector3<double>> node_normals_F(num_nodes);
  for (int node_index = 0;  node_index < num_nodes; ++node_index) {
    node_normals_F[node_index] = rhs[node_index] / node_area[node_index];
    node_normals_F[node_index].normalize();
  }
  return node_normals_F;
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

void OutputMeshToOBJ(
    const std::string& file_name,
    const std::vector<Vector3<double>>& points,
    const std::vector<Vector3<int>>& triangles) {
  const int num_nodes = points.size();
  const int num_tris = triangles.size();

  std::ofstream file(file_name);

  for (int i_node = 0; i_node < num_nodes; i_node++) {
    const Vector3<double> pos = points[i_node];
    file << "v " << pos[0] << " " << pos[1]
         << " " << pos[2] << std::endl;
  }

  file << std::endl;
  for (int i_tri = 0; i_tri < num_tris; i_tri++) {
    const Vector3<int> tri = triangles[i_tri];
    file << "f " << tri[0] + 1 << " " << tri[1] + 1
         << " " << tri[2] + 1 << std::endl;
  }
  file.close();
}

}  // namespace mesh_query
}  // namespace geometry
}  // namespace drae

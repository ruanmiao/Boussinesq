#pragma once

#include <fstream>
#include <memory>
#include <set>
#include <vector>

#include <fmt/format.h>
#include <fmt/ostream.h>

#include "drake/common/drake_assert.h"
#include "drake/common/drake_copyable.h"
#include "drake/common/drake_throw.h"
#include "drake/common/eigen_types.h"

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

  VectorX<double> node_areas;

  // A conveniently pre-computed array such that for a node_index then
  // node_element[node_index].first is an element of which node_index is a
  // vertex.
  // node_element[node_index].second tells us the local index of node_index in
  // that element.
  // The element assigned to node_index is arbitrary.
  std::vector<std::pair<int, int>> node_element;

  // A vector of size num_nodes = points_G.size().
  // Entry at node_index contains the set of all triangles adjacent to
  // node_index.
  std::vector<std::vector<int>> node_triangles;
};

double CalcTriangleArea(
    const Vector3<double>& p_AP1,
    const Vector3<double>& p_AP2,
    const Vector3<double>& p_AP3);

std::vector<Vector3<double>> CalcMeshFaceNormals(
    const std::vector<Vector3<double>>& points_A,
    const std::vector<Vector3<int>>& triangles);

std::vector<Vector3<double>> CalcAreaWeightedNormals(const Mesh<double>& mesh);

void FlipNormals(Mesh<double>* mesh);

void OutputMeshToOBJ(
    const std::string& file_name,
    const std::vector<Vector3<double>>& points,
    const std::vector<Vector3<int>>& triangles);

std::unique_ptr<Mesh<double>> LoadMeshFromObj(
    const std::string& file_name, bool flip_normals);

}  // namespace mesh_query
}  // namespace geometry
}  // namespace drae

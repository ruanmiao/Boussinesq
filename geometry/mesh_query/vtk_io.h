#pragma once

#include <fstream>
#include <memory>
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


void OutputMeshToVTK(
    std::ofstream& file,
    const std::vector<Vector3<double>>& points_G,
    const std::vector<Vector3<int>>& triangles,
    const Isometry3<double>& X_WG = Isometry3<double>::Identity()) {
  const int num_nodes = points_G.size();
  const int num_tris = triangles.size();

  // Header for the VTK file.
  file << "# vtk DataFile Version 3.0" << std::endl;
  file << "Visualize mesh data" << std::endl;
  file << "ASCII" << std::endl;
  file << std::endl;
  file << "DATASET UNSTRUCTURED_GRID" << std::endl;

  // Header fot he points_G data.
  file << "POINTS " << num_nodes << " double" << std::endl;

  // Write the points_G data.
  for (const auto& point_G : points_G) {
    const auto point_W = X_WG * point_G;
    file << fmt::format("{:.8f} {:.8f} {:.8f}\n",
                        point_W[0], point_W[1], point_W[2]);
  }
  file << std::endl;

  // Header for the elements (triangles in this case).
  file << "CELLS " << num_tris << " " << num_tris * 4 << std::endl;


  for (const auto& triangle : triangles) {
    file << "3 " << triangle[0] << " " << triangle[1] << " " << triangle[2] <<
         std::endl;
  }
  file << std::endl;

  // VTK needs us to spell out the element type for each element.
  // Cell type = 5 is for triangles.
  file << "CELL_TYPES " << num_tris << std::endl;
  for (int i_tri = 0; i_tri < num_tris; i_tri++) {
    file << "5" << std::endl;
  }
  file << std::endl;
}

void OutputMeshToVTK(
    const std::string& file_name,
    const std::vector<Vector3<double>>& points_G,
    const std::vector<Vector3<int>>& triangles,
    const Isometry3<double>& X_WG = Isometry3<double>::Identity()) {
  std::ofstream file(file_name);
  OutputMeshToVTK(file, points_G, triangles, X_WG);
  file.close();
}

void OutputScatteredPointsToVTK(
    std::ofstream& file,
    const std::vector<Vector3<double>>& points) {
  const int num_nodes = points.size();

  // Header for the VTK file.
  file << "# vtk DataFile Version 3.0" << std::endl;
  file << "Visualize mesh data" << std::endl;
  file << "ASCII" << std::endl;
  file << std::endl;
  file << "DATASET UNSTRUCTURED_GRID" << std::endl;

  // Header fot he points data.
  file << "POINTS " << num_nodes << " double" << std::endl;

  // Write the points data.
  for (const auto& point : points) {
    file << fmt::format("{:.8f} {:.8f} {:.8f}\n", point[0], point[1], point[2]);
  }
  file << std::endl;

  // VTK needs us to spell out the element type for each element.
  // Cell type = 1 is for a single point.
  file << "CELL_TYPES " << num_nodes << std::endl;
  for (int i_tri = 0; i_tri < num_nodes; i_tri++) {
    file << "1" << std::endl;
  }
  file << std::endl;
}

void OutputScatteredPointsToVTK(
    const std::string& file_name,
    const std::vector<Vector3<double>>& points) {
  std::ofstream file(file_name);
  OutputScatteredPointsToVTK(file, points);
  file.close();
}

void AppendCellCenteredVectorFieldToVTK(
    std::ofstream& file,
    const std::string& field_name,
    const std::vector<Vector3<double>>& vector_field_G,
    const Isometry3<double>& X_WG = Isometry3<double>::Identity()) {
  const int num_nodes = vector_field_G.size();

  file << std::endl;
  file << "CELL_DATA " << num_nodes << std::endl;
  file << "VECTORS " + field_name + " double" << std::endl;

  for (const auto& vector_G : vector_field_G) {
    const auto vector_W = X_WG.linear() * vector_G;
    file << fmt::format("{:.8f} {:.8f} {:.8f}\n",
                        vector_W[0], vector_W[1], vector_W[2]);
  }
  file << std::endl;
}

void AppendNodeCenteredVectorFieldToVTK(
    std::ofstream& file,
    const std::string& field_name,
    const std::vector<Vector3<double>>& vector_field_G,
    const Isometry3<double>& X_WG = Isometry3<double>::Identity()) {
  const int num_nodes = vector_field_G.size();

  file << std::endl;
  file << "POINT_DATA " << num_nodes << std::endl;
  file << "VECTORS " + field_name + " double" << std::endl;

  for (const auto& vector_G : vector_field_G) {
    const auto vector_W = X_WG.linear() * vector_G;
    file << fmt::format("{:.8f} {:.8f} {:.8f}\n",
                        vector_W[0], vector_W[1], vector_W[2]);
  }
  file << std::endl;
}

void OutputSegmentsToVTK(
    std::ofstream& file,
    const std::vector<Vector3<double>>& start_G,
    const std::vector<Vector3<double>>& end_G) {
  DRAKE_DEMAND(start_G.size() == end_G.size());
  const int num_segments = start_G.size();

  // Header for the VTK file.
  file << "# vtk DataFile Version 3.0" << std::endl;
  file << "Visualize mesh data" << std::endl;
  file << "ASCII" << std::endl;
  file << std::endl;
  file << "DATASET UNSTRUCTURED_GRID" << std::endl;

  // Header fot he points_G data.
  file << "POINTS " << 2 * num_segments << " double" << std::endl;

  // Write the points_G data.
  for (int i = 0; i < num_segments; ++i) {
    const auto& start = start_G[i];
    const auto& end = end_G[i];
    file << fmt::format("{:.8f} {:.8f} {:.8f}\n", start[0], start[1], start[2]);
    file << fmt::format("{:.8f} {:.8f} {:.8f}\n", end[0], end[1], end[2]);
  }
  file << std::endl;

  // Header for the elements (triangles in this case).
  file << "CELLS " << num_segments << " " << num_segments * 3 << std::endl;

  for (int i = 0; i < num_segments; ++i) {
    file << "2 " << 2*i << " " << 2*i+1 << std::endl;
  }
  file << std::endl;

  // VTK needs us to spell out the element type for each element.
  // Cell type = 3 is for segments.
  file << "CELL_TYPES " << num_segments << std::endl;
  for (int i_tri = 0; i_tri < num_segments; i_tri++) {
    file << "3" << std::endl;
  }
  file << std::endl;
}

}  // namespace mesh_query
}  // namespace geometry
}  // namespace drae

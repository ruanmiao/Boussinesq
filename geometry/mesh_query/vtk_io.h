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
    const std::string& file_name,
    const std::vector<Vector3<double>>& points,
    const std::vector<Vector3<int>>& triangles) {
  const int num_nodes = points.size();
  const int num_tris = triangles.size();

  std::ofstream file(file_name);

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

  file.close();
}

void OutputScatteredPointsToVTK(
    const std::string& file_name,
    const std::vector<Vector3<double>>& points) {
  const int num_nodes = points.size();

  std::ofstream file(file_name);

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

  file.close();
}

}  // namespace mesh_query
}  // namespace geometry
}  // namespace drae

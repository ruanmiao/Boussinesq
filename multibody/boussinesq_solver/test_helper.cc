#include "drake/multibody/boussinesq_solver/test_helper.h"

#include <fstream>
#include <iostream>
#include <string>

#include "drake/multibody/boussinesq_solver/math_helper.h"
#include "drake/common/drake_assert.h"

namespace drake {
namespace multibody {
namespace boussinesq_solver {

std::pair<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3i>>
MeshSquare(
    const Vector2<double>& p1, const Vector2<double>& p2,
    const Vector2<double>& p3, int num_px, int num_py) {

  const int num_nodes = num_px * num_py;
  const int num_tris = (num_px - 1) * (num_py - 1) * 2;
  std::vector<Vector3<double>> points(num_nodes);
  std::vector<Vector3<int>> triangles(num_tris);

  const VectorX<double>& x_positions =
      VectorX<double>::LinSpaced(num_px, p1(0), p2(0));
  const VectorX<double>& y_positions =
      VectorX<double>::LinSpaced(num_py, p1(1), p3(1));
  int i_tri = 0;

  for (int iy = 0; iy < num_py; iy++) {
    for (int ix = 0; ix < num_px; ix++) {
      const int cur_i = iy * num_px + ix;
      const Vector3<double> pos(x_positions(ix), y_positions(iy), 0);
      points[cur_i] = pos;

      if ((ix == num_px - 1) || (iy == num_py - 1)) continue;

      const Vector3<int> tri1(cur_i, cur_i + 1, cur_i + num_px);
      const Vector3<int> tri2(cur_i + 1 + num_px, cur_i + 1, cur_i + num_px);

      assert(i_tri < num_tris - 1);
      triangles[i_tri] = tri1;
      i_tri++;
      triangles[i_tri] = tri2;
      i_tri++;
    }
  }
  return std::make_pair(points, triangles);
}

std::pair<std::vector<Eigen::Vector3d>,
          std::vector<Eigen::Vector3i>> MeshCircle(
              const Vector2<double>& o, double radius, int num_pr) {
  const int num_nodes = 3 * num_pr * (num_pr - 1) + 1;
  const int num_tris = 6 * (num_pr - 1) * (num_pr - 1);

  std::vector<Vector3<double>> points(num_nodes);
  std::vector<Vector3<int>> triangles(num_tris);

  const VectorX<double>& r_positions =
      VectorX<double>::LinSpaced(num_pr, 0, radius);
  const VectorX<double>& theta_positions =
      VectorX<double>::LinSpaced(7, 0, 2 * M_PI);

  int cur_p = 0;
  int cur_t = 0;

  points[cur_p] << o(0), o(1), 0;
  cur_p++;

  for (int ir = 1; ir < num_pr; ir++) {
    for (int itheta = 0; itheta < 6; itheta++) {
      int ipos = 0;
      const double r = r_positions(ir);

      const VectorX<double>& p_positions = VectorX<double>::LinSpaced(
              ir + 1, theta_positions(itheta), theta_positions(itheta + 1));

      DRAKE_ASSERT(cur_t < num_tris);
      const Vector3<int> tri(
          FindIndexInMeshCircle(ir, itheta, ipos),
          FindIndexInMeshCircle(ir-1, itheta, ipos),
          FindIndexInMeshCircle(ir, itheta, ipos+1));
      triangles[cur_t] = tri;
      cur_t++;

      const Vector3<double> pos0(
          o(0) + r * cos(p_positions(0)),
          o(1) + r * sin(p_positions(0)), 0);
      DRAKE_ASSERT(cur_p < num_nodes);
      points[cur_p] = pos0;
      cur_p++;

      for (ipos = 1; ipos < ir; ipos++) {
        DRAKE_ASSERT(cur_t < num_tris);
        const Vector3<int> tri1(
            FindIndexInMeshCircle(ir, itheta, ipos),
            FindIndexInMeshCircle(ir - 1, itheta, ipos - 1),
            FindIndexInMeshCircle(ir - 1, itheta, ipos));
        triangles[cur_t] = tri1;
        cur_t++;

        DRAKE_ASSERT(cur_t < num_tris);
        const Vector3<int> tri2(
            FindIndexInMeshCircle(ir, itheta, ipos),
            FindIndexInMeshCircle(ir - 1, itheta, ipos),
            FindIndexInMeshCircle(ir, itheta, ipos + 1));

            triangles[cur_t] = tri2;
        cur_t++;

        const Vector3<double> pos(
            o(0) + r * cos(p_positions(ipos)),
            o(1) + r * sin(p_positions(ipos)), 0);

        DRAKE_ASSERT(cur_p < num_nodes);
        points[cur_p] = pos;
        cur_p++;
      }
    }
  }
  return std::make_pair(points, triangles);
}

std::pair<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3i>>
MeshSphere(
    const Vector3<double> &o, double radius, int num_pr, double half_sector) {
  DRAKE_ASSERT(fabs(half_sector- M_PI) > std::numeric_limits<double>::epsilon());
  const int num_nodes = 3 * num_pr * (num_pr - 1) + 1;
  const int num_tris = 6 * (num_pr - 1) * (num_pr - 1);

  std::vector<Vector3<double>> points(num_nodes);
  std::vector<Vector3<int>> triangles(num_tris);

  const VectorX<double>& sector_positions =
      VectorX<double>::LinSpaced(num_pr, 0, half_sector);

  const VectorX<double> r_positions =
      radius * sector_positions.array().sin();
  const VectorX<double> z_positions = -radius * sector_positions.array().cos();
  const VectorX<double>& theta_positions =
      VectorX<double>::LinSpaced(7, 0, 2 * M_PI);

  int cur_p = 0;
  int cur_t = 0;

  points[cur_p] << o(0), o(1), (o(2) - radius);
  cur_p++;

  for (int ir = 1; ir < num_pr; ir++) {
    for (int itheta = 0; itheta < 6; itheta++) {
      int ipos = 0;
      const double r = r_positions(ir);
      const double z = z_positions(ir);

      const VectorX<double>& p_positions = VectorX<double>::LinSpaced(
          ir + 1, theta_positions(itheta), theta_positions(itheta + 1));

      assert(cur_t < num_tris);
      const Vector3<int> tri(
          FindIndexInMeshCircle(ir, itheta, ipos),
          FindIndexInMeshCircle(ir-1, itheta, ipos),
          FindIndexInMeshCircle(ir, itheta, ipos+1));
      triangles[cur_t] = tri;
      cur_t++;

      const Vector3<double> pos0(
          o(0) + r * cos(p_positions(0)),
          o(1) + r * sin(p_positions(0)),
          o(2) + z);
      assert(cur_p < num_nodes);
      points[cur_p] = pos0;
      cur_p++;

      for (ipos = 1; ipos < ir; ipos++) {
        assert(cur_t < num_tris);
        const Vector3<int> tri1(
            FindIndexInMeshCircle(ir, itheta, ipos),
            FindIndexInMeshCircle(ir - 1, itheta, ipos - 1),
            FindIndexInMeshCircle(ir - 1, itheta, ipos));
        triangles[cur_t] = tri1;
        cur_t++;

        assert(cur_t < num_tris);
        const Vector3<int> tri2(
            FindIndexInMeshCircle(ir, itheta, ipos),
            FindIndexInMeshCircle(ir - 1, itheta, ipos),
            FindIndexInMeshCircle(ir, itheta, ipos + 1));

        triangles[cur_t] = tri2;
        cur_t++;

        const Vector3<double> pos(
            o(0) + r * cos(p_positions(ipos)),
            o(1) + r * sin(p_positions(ipos)),
            o(2) + z);

        assert(cur_p < num_nodes);
        points[cur_p] = pos;
        cur_p++;
      }
    }
  }
  return std::make_pair(points, triangles);
}

int FindIndexInMeshCircle(int ir, int itheta, int ipos) {
  const int layer_max = 3 * ir * (ir + 1) + 1;
  if (ir == 0) return 0;
  int index = 1 + 3 * ir * (ir - 1) + itheta * ir + ipos;
  if (index < layer_max) return  index;
  return index - ir * 6;
}

Eigen::VectorXd GetPressureIntegrandR(
    const std::vector<Eigen::Vector3d>& points_in_mesh) {
  const int num_nodes = points_in_mesh.size();
  Eigen::VectorXd pressures(num_nodes);

  for (int i_node = 0; i_node < num_nodes; i_node++) {
    const Vector3<double> pos = points_in_mesh[i_node];
    pressures(i_node, 0) = sqrt(pos(0) * pos(0) + pos(1) * pos(1));
  }
  return pressures;
}

Eigen::VectorXd GetPressureIntegrandX(
    const std::vector<Eigen::Vector3d>& points_in_mesh) {
  const int num_nodes = points_in_mesh.size();
  Eigen::VectorXd pressures(num_nodes);

  for (int i_node = 0; i_node < num_nodes; i_node++) {
    const Vector3<double> pos = points_in_mesh[i_node];
    pressures(i_node, 0) = pos(0);
  }
  std::cout << "cout try !!!" << std::endl;
  return pressures;
}

bool OutputMeshToVTK(
    const std::vector<Vector3<double>>& points_in_mesh,
                     const std::vector<Vector3<int>>& triangles_in_mesh,
                     const VectorX<double>& values,
                     bool enabled) {
  if (!enabled) {
    return enabled;
  }
  const int num_nodes = points_in_mesh.size();
  const int num_tris = triangles_in_mesh.size();

  std::ofstream mesh_obj;
  mesh_obj.open("mesh_data.obj");

  // output to .obj file
  for (int i_node = 0; i_node < num_nodes; i_node++) {
    const Vector3<double> pos = points_in_mesh[i_node];
    mesh_obj << "v " << pos[0] << " " << pos[1]
              << " " << pos[2] << std::endl;
  }

  mesh_obj << std::endl;
  for (int i_tri = 0; i_tri < num_tris; i_tri++) {
    const Vector3<int> tri = triangles_in_mesh[i_tri];
    mesh_obj << "f " << tri[0] + 1 << " " << tri[1] + 1
              << " " << tri[2] + 1 << std::endl;
  }
  mesh_obj.close();

  std::ofstream mesh_vtk_2d;
  mesh_vtk_2d.open("mesh_data_2D.vtk");

  std::ofstream mesh_vtk_3d;
  mesh_vtk_3d.open("mesh_data_3D.vtk");

  //  output to .vtk file.
  mesh_vtk_2d << "# vtk DataFile Version 3.0" << std::endl;
  mesh_vtk_2d << "Visualize mesh data" << std::endl;
  mesh_vtk_2d << "ASCII" << std::endl;
  mesh_vtk_2d << std::endl;
  mesh_vtk_2d << "DATASET UNSTRUCTURED_GRID" << std::endl;
  mesh_vtk_2d << "POINTS " << num_nodes << " double" << std::endl;

  mesh_vtk_3d << "# vtk DataFile Version 3.0" << std::endl;
  mesh_vtk_3d << "Visualize mesh data" << std::endl;
  mesh_vtk_3d << "ASCII" << std::endl;
  mesh_vtk_3d << std::endl;
  mesh_vtk_3d << "DATASET UNSTRUCTURED_GRID" << std::endl;
  mesh_vtk_3d << "POINTS " << num_nodes << " double" << std::endl;

  for (int i_node = 0; i_node < num_nodes; i_node++) {
    const Vector3<double> pos = points_in_mesh[i_node];
    mesh_vtk_2d << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
    mesh_vtk_3d << pos[0] << " " << pos[1]
                << " " << values[i_node] << std::endl;
  }

  mesh_vtk_2d << std::endl;
  mesh_vtk_2d << "CELLS " << num_tris << " " << num_tris * 4 << std::endl;

  mesh_vtk_3d << std::endl;
  mesh_vtk_3d << "CELLS " << num_tris << " " << num_tris * 4 << std::endl;

  for (int i_tri = 0; i_tri < num_tris; i_tri++) {
    const Vector3<int> tri = triangles_in_mesh[i_tri];
    mesh_vtk_2d << "3 " << tri[0] << " " << tri[1]
              << " " << tri[2] << std::endl;
    mesh_vtk_3d << "3 " << tri[0] << " " << tri[1]
                << " " << tri[2] << std::endl;
  }

  mesh_vtk_2d << std::endl;
  mesh_vtk_2d << "CELL_TYPES " << num_tris << std::endl;

  mesh_vtk_3d << std::endl;
  mesh_vtk_3d << "CELL_TYPES " << num_tris << std::endl;

  for (int i_tri = 0; i_tri < num_tris; i_tri++) {
    mesh_vtk_2d << "5" << std::endl;
    mesh_vtk_3d << "5" << std::endl;
  }

  mesh_vtk_2d << std::endl;
  mesh_vtk_2d << "POINT_DATA " << num_nodes << std::endl;
  mesh_vtk_2d << "SCALARS pressure double 1" << std::endl;
  mesh_vtk_2d << "LOOKUP_TABLE default" << std::endl;

  for (int i_node = 0; i_node < num_nodes; i_node++) {
    mesh_vtk_2d << values(i_node) << std::endl;
  }
  mesh_vtk_2d << std::endl;
  mesh_vtk_2d.close();

  mesh_vtk_3d << std::endl;
  mesh_vtk_3d.close();

  return true;
}

double CalcForceOverMesh(const std::vector<Vector3<double>>& points_in_mesh,
                         const std::vector<Vector3<int>>& triangles_in_mesh,
                         const VectorX<double>& pressure) {
  const int num_tris = triangles_in_mesh.size();
  double total_force = 0;

  for (int i_tri = 0; i_tri < num_tris; i_tri++) {
    const Vector3<int>& nodes = triangles_in_mesh[i_tri];
    const Vector3<double>& p1 = points_in_mesh[nodes[0]];
    const Vector3<double>& p2 = points_in_mesh[nodes[1]];
    const Vector3<double>& p3 = points_in_mesh[nodes[2]];
    Vector3<double> vector_tri_area = CalcTriangleArea(p1, p2, p3);
    double tri_area = vector_tri_area.norm();

    for (int i_node = 0; i_node < 3; i_node++) {
      total_force = total_force + tri_area / 3 * pressure[nodes[i_node]];
    }
  }
  return total_force;
}

double CalcForceOverMeshOfSphere(
    const std::vector<Vector3<double>>& points_in_mesh,
    const std::vector<Vector3<int>>& triangles_in_mesh,
    const VectorX<double>& pressure,
    const Vector3<double>& center,
    const Vector3<double>& pressure_dir) {
  const int num_tris = triangles_in_mesh.size();
  double total_force = 0;

  for (int i_tri = 0; i_tri < num_tris; i_tri++) {
    const Vector3<int>& nodes = triangles_in_mesh[i_tri];
    const Vector3<double>& p1 = points_in_mesh[nodes[0]];
    const Vector3<double>& p2 = points_in_mesh[nodes[1]];
    const Vector3<double>& p3 = points_in_mesh[nodes[2]];
    Vector3<double> vector_tri_area = CalcTriangleArea(p1, p2, p3);
    double tri_area = vector_tri_area.norm();

    for (int i_node = 0; i_node < 3; i_node++) {
      Vector3<double> dir = center - points_in_mesh[nodes[i_node]];
      const double theta_cos = dir.dot(pressure_dir)
          / (dir.norm() * pressure_dir.norm());

      total_force = total_force +
          tri_area / 3 * pressure[nodes[i_node]] * theta_cos;
    }
  }
  return total_force;
}


//void LoadObjFile(std::vector<Vector3<double>>& vertices,
//                 std::vector<Vector3<double>>& triangles,
//                 std::string obj_file_name) {
//  std::ifstream file(obj_file_name);
//  if (!file) {
//    throw std::runtime_error("Error opening file \"" + obj_file_name + "\".");
//  }
//
//  std::string path;
//  const size_t idx = obj_file_name.rfind('/');
//  if (idx != std::string::npos) {
//    path = obj_file_name.substr(0, idx + 1);
//  }
//
//  std::tinyobj::attrib_t attrib;
//  std::vector<tinyobj::shape_t> shapes;
//  std::vector<tinyobj::material_t> materials;
//  std::string err;
//
//  // Would be nice to hand off triangulation of non-triangular faces
//  // to tinyobjloader, however tinyobjloader does no checking that
//  // the triangulation is in fact correct.
//  bool do_tinyobj_triangulation = false;
//  bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &err,
//                              obj_file_name.c_str(), path.c_str(), do_tinyobj_triangulation);
//
//  // Use the boolean return value and the error string to determine
//  // if we should proceeed.
//  if (!ret || !err.empty()) {
//    throw std::runtime_error("Error parsing file \""
//                                 + obj_file_name + "\" : " + err);
//  }
//
//  // Store the vertices.
//  for (int index = 0;
//       index < static_cast<int>(attrib.vertices.size());
//       index += 3) {
//    vertices->push_back(Vector3d(attrib.vertices[index] * scale_[0],
//                                 attrib.vertices[index + 1] * scale_[1],
//                                 attrib.vertices[index + 2] * scale_[2]));
//  }
//
//  // Counter for triangles added to compensate for non-triangular faces.
//  int added_triangles = 0;
//  // Counter keeping track of indices, used to check against the number
//  // of vertices later to ensure all is correct.
//  int maximum_index = 0;
//
//  // Iterate over the shapes.
//  for (const auto& one_shape : shapes) {
//    const tinyobj::mesh_t& mesh = one_shape.mesh;
//    int index_offset = 0;
//
//    // For each face in the shape.
//    for (int face = 0;
//         face < static_cast<int>(mesh.num_face_vertices.size());
//         ++face) {
//      const int vert_count = mesh.num_face_vertices[face];
//
//      std::vector<int> indices;
//      for (int vert = 0; vert < vert_count; ++vert) {
//        // Store the vertex index.
//        int vertex_index = mesh.indices[index_offset + vert].vertex_index;
//        maximum_index = std::max(maximum_index, vertex_index);
//        indices.push_back(vertex_index);
//      }
//
//      // Make sure the face has three vertices. This can only occur if
//      // do_tiny_obj_triangulation is false and faces are non triangular.
//      if (indices.size() != 3) {
//        if (triangulate == TriangulatePolicy::kTry && indices.size() > 3) {
//          // This is a very naive triangulation.  It simply creates a fan
//          // around the 0th vertex through the loop of vertices in the polygon.
//          // The fan won't necessarily produce triangles with the best aspect
//          // ratio.
//          //
//          // The requirement that two adjacent, generated triangles have normals
//          // that deviate by no more than 30° will implicitly catch the
//          // degenerate case where applying this algorithm to a concave face
//          // would cause changes in triangle winding and, therefore, normal
//          // directions.
//          //
//          // For a polygon with n vertices indexed in the range [0, n - 1] we
//          // create triangles built on those indices in the pattern:
//          //    (0, 1, 2), (0, 2, 3), (0, 3, 4), ..., (0, n-2, n-1).
//          // The triangle consisting of (0, 1, 2) is handled below, so this
//          // starts with (0, 2, 3).
//          Vector3d lastNormal;
//          bool valid = GetNormal(*vertices, indices[0], indices[1],
//                                 indices[2], &lastNormal, kMinArea);
//          const int index_size = static_cast<int>(indices.size());
//          for (int i = 2; valid && i < index_size - 1; ++i) {
//            // OBJ file indices are 1-based; subtracting 1 makes them 0-based.
//            Vector3d testNormal;
//            valid = GetNormal(*vertices, indices[0], indices[i],
//                              indices[i + 1], &testNormal, kMinArea);
//            if (valid) {
//              double dot_product = lastNormal.dot(testNormal);
//              if (dot_product < 0) {
//                throw std::runtime_error("Trying to triangulate face number " +
//                    std::to_string(face) + " in '" + obj_file_name +
//                    "' led to bad triangles. The triangle based " +
//                    "on vertices " + std::to_string(indices[0]) +
//                    ", " + std::to_string(indices[i]) + ", and " +
//                    std::to_string(indices[i + 1]) +
//                    " (0-indexed) is wound in the opposite " +
//                    "direction from the previous triangle. " +
//                    "Consider triangulating by hand.");
//              } else if (dot_product < kCosThreshold) {
//                throw std::runtime_error("Trying to triangulate face number " +
//                    std::to_string(face) + " in '" + obj_file_name +
//                    "'.  The face is not sufficiently planar. " +
//                    "Consider triangulating by hand.");
//              }
//            }
//            lastNormal = testNormal;
//            triangles->push_back(
//                Vector3i(indices[0], indices[i], indices[i + 1]));
//          }
//          if (!valid) {
//            // Problems in computing the normal are logged in GetNormal.
//            throw std::runtime_error("Unable to triangulate face number " +
//                std::to_string(face) + " in '" + obj_file_name +
//                "'. See log for details.");
//          }
//          // A polygon of n vertices produces a fan of n - 2 triangles.
//          // One of those is a given, so we're "adding" n - 3 triangles.
//          added_triangles += index_size - 3;
//        } else {
//          throw std::runtime_error(
//              "In file \"" + obj_file_name + "\" (face " +
//                  std::to_string(face) +
//                  "). Only triangular faces are supported. However "
//                  + std::to_string(vert_count) + " indices are provided.");
//        }
//      }
//
//      triangles->push_back(Vector3i(indices[0], indices[1], indices[2]));
//
//      index_offset += vert_count;
//    }
//  }
//
//  // Verifies that the maximum index referenced when defining faces does not
//  // exceed the number of vertices (note max index starts from zero).
//  if (maximum_index >= static_cast<int>(vertices->size())) {
//    throw std::runtime_error(
//        "In file \"" + obj_file_name + "\". "
//                                       "The maximum index referenced in defining faces (" +
//            std::to_string(maximum_index) + ") "
//                                            "exceeds the number of vertices (" +
//            std::to_string(vertices->size()) + ". ");
//  }
//
//  if (added_triangles > 0) {
//    drake::log()->info("Encountered non-triangular faces in '" + obj_file_name +
//        "'. Triangulation was applied, adding " +
//        std::to_string(added_triangles) +
//        " new triangles to the mesh.");
//  }
//}

}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake

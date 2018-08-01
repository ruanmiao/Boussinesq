#include "drake/multibody/boussinesq_solver/test_helper.h"

#include <fstream>
#include <iostream>
#include <string>

#include "drake/multibody/boussinesq_solver/math_helper.h"

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

      assert(cur_t < num_tris);
      const Vector3<int> tri(
          FindIndexInMeshCircle(ir, itheta, ipos),
          FindIndexInMeshCircle(ir-1, itheta, ipos),
          FindIndexInMeshCircle(ir, itheta, ipos+1));
      triangles[cur_t] = tri;
      cur_t++;

      const Vector3<double> pos0(
          o(0) + r * cos(p_positions(0)),
          o(1) + r * sin(p_positions(0)), 0);
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
            o(1) + r * sin(p_positions(ipos)), 0);

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
                     const VectorX<double>& values) {
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




}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake

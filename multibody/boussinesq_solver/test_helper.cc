#include "drake/multibody/boussinesq_solver/test_helper.h"

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
  cur_p ++;

  for (int ir = 1; ir < num_pr; ir++) {
    for (int itheta = 0; itheta < 6; itheta++) {

      int ipos = 0;
      const double r = r_positions(ir);

      const VectorX<double>& p_positions = VectorX<double>::LinSpaced(
              ir + 1, theta_positions(itheta), theta_positions(itheta + 1));

      assert(cur_t < num_tris);
      const Vector3<int> tri(
          FindIndexInMeshCircle(ir, itheta, ipos),
          FindIndexInMeshCircle(ir, itheta, ipos+1),
          FindIndexInMeshCircle(ir-1, itheta, ipos));
      triangles[cur_t] = tri;
      cur_t ++;

      const Vector3<double> pos0(
          o(0) + r * cos(p_positions(0)),
          o(1) + r * sin(p_positions(0)), 0);
      assert(cur_p < num_nodes);
      points[cur_p] = pos0;
      cur_p ++;

      for (ipos = 1; ipos < ir; ipos ++) {

        assert(cur_t < num_tris);
        const Vector3<int> tri1(
            FindIndexInMeshCircle(ir, itheta, ipos),
            FindIndexInMeshCircle(ir - 1, itheta, ipos - 1),
            FindIndexInMeshCircle(ir - 1, itheta, ipos));
        triangles[cur_t] = tri1;
        cur_t ++;

        assert(cur_t < num_tris);
        const Vector3<int> tri2(
            FindIndexInMeshCircle(ir, itheta, ipos),
            FindIndexInMeshCircle(ir, itheta, ipos + 1),
            FindIndexInMeshCircle(ir - 1, itheta, ipos));
        triangles[cur_t] = tri2;
        cur_t ++;

        const Vector3<double> pos(
            o(0) + r * cos(p_positions(ipos)),
            o(1) + r * sin(p_positions(ipos)), 0);

        assert(cur_p < num_nodes);
        points[cur_p] = pos;
        cur_p ++;
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
  return pressures;
}

}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake

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


}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake

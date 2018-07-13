#include "drake/multibody/boussinesq_solver/compliance_matrix.h"

#include "drake/multibody/boussinesq_solver/element_compliance_matrix.h"

namespace drake {
namespace multibody {
namespace boussinesq_solver {

MatrixX<double> CalcComplianceMatrix(
    const std::vector<Eigen::Vector3d>& points_in_mesh,
    const std::vector<Eigen::Vector3i>& triangles_in_mesh,
    double k_const) {
  const int num_nodes = points_in_mesh.size();
  MatrixX<double> compliance(num_nodes, num_nodes);

  for (int i_row = 0; i_row < num_nodes; i_row++) {
    compliance.row(i_row) = CalcElementComplianceRowMatrix(points_in_mesh,
        triangles_in_mesh, points_in_mesh.at(i_row), k_const);
  }

  return compliance;
}

}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake

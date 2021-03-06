#include "drake/multibody/boussinesq_solver/row_compliance_matrix.h"

#include <limits>

#include "drake/common/drake_assert.h"
#include "drake/multibody/boussinesq_solver/integral_general_triangle.h"
#include "drake/multibody/boussinesq_solver/math_helper.h"

namespace drake {
namespace multibody {
namespace boussinesq_solver {

MatrixX<double> CalcRowComplianceMatrix(
    const std::vector<Eigen::Vector3d>& points_in_mesh,
    const std::vector<Eigen::Vector3i>& triangles_in_mesh,
    const Vector3<double>& node_A,
    double k_const) {
  const size_t num_nodes = points_in_mesh.size();
  const size_t num_tris = triangles_in_mesh.size();
  MatrixX<double> compliance = MatrixX<double>::Zero(1, num_nodes);

  for (size_t i_tri = 0; i_tri < num_tris; ++i_tri) {
    const Eigen::Vector3i& indexes = triangles_in_mesh[i_tri];
    const Eigen::Vector3d& p1 = points_in_mesh[indexes(0)];
    const Eigen::Vector3d& p2 = points_in_mesh[indexes(1)];
    const Eigen::Vector3d& p3 = points_in_mesh[indexes(2)];


    const Eigen::Isometry3d X_WT = CalcTransformationFromTriangleFrame(
        p1, p2, p3, node_A);
    const Eigen::Isometry3d X_TW = X_WT.inverse();

    const Eigen::Vector3d xA_T = X_TW * node_A;
    const Eigen::Vector3d p1_T = X_TW * p1;
    const Eigen::Vector3d p2_T = X_TW * p2;
    const Eigen::Vector3d p3_T = X_TW * p3;

    DRAKE_ASSERT(fabs(xA_T(2)) < std::numeric_limits<double>::infinity());

    Vector3<double> element_compliance = CalcGeneralTriangleCompliance(xA_T(2),
        p1_T.head(2), p2_T.head(2), p3_T.head(2), k_const);

    DRAKE_ASSERT(fabs(element_compliance.norm()) <
        std::numeric_limits<double>::infinity());

    const double c1 = compliance(indexes(0));
    compliance(indexes(0)) = c1 + element_compliance(0);
    const double c2 = compliance(indexes(1));
    compliance(indexes(1)) = c2 + element_compliance(1);
    const double c3 = compliance(indexes(2));
    compliance(indexes(2)) = c3 + element_compliance(2);
  }

  return compliance;
}
}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake

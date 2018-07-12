#include "drake/multibody/boussinesq_solver/element_compliance_matrix.h"

#include <gtest/gtest.h>

#include "drake/common/eigen_types.h"
#include "drake/common/test_utilities/eigen_matrix_compare.h"
#include "drake/multibody/boussinesq_solver/test_helper.h"

namespace drake {
namespace multibody {
namespace boussinesq_solver {
namespace {

using Eigen::Vector2d;
using Eigen::Vector3d;

GTEST_TEST(ElementComplianceMatrixTest, IntegrantR) {
  const Vector2d p1(0.0, 0.0);
  const Vector2d p2(1.0, 0.0);
  const Vector2d p3(0.0, 1.0);
  const int num_px = 3;
  const int num_py = 3;

  const std::pair<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3i>>&
      mesh_data = MeshSquare(p1, p2, p3, num_px, num_py);

  MatrixX<double> compliance_p1 = CalcElementComplianceRowMatrix(
      mesh_data.first, mesh_data.second, Vector3d(0.0, 0.0, 0.0), 1);
  VectorX<double> pressures = GetPressureIntegrandR(mesh_data.first);

  MatrixX<double> deformation_p1 = compliance_p1 * pressures;
  double results = deformation_p1(0,0);
  const double expected = 1.0;
  EXPECT_NEAR(results,
              expected,
              10 * std::numeric_limits<double>::epsilon());
}

}  // namespace
}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake

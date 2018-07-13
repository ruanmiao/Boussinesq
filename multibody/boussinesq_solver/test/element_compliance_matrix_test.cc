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

/// The expected values for this test are the results by running Matlab
/// Prototype. The precision (15 digits) of the expected values if the same as
/// the "long" in Matlab
/// The analytical solution to the result should be 1, since ∫∫ r/r drdθ=Area
/// (different from the "expected") bellow. The difference/error to the
/// analytical solution arises from the non-linear pressure field in the
/// Cardesian system and the triangular linear interpolation.
GTEST_TEST(ElementComplianceMatrixTest, IntegrantR) {
  const double length = 1;
  const Vector2d p1(0.0, 0.0);
  const Vector2d p2(length, 0.0);
  const Vector2d p3(0.0, length);
  const int num_px = 3;
  const int num_py = 3;

  const std::pair<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3i>>&
      mesh_data = MeshSquare(p1, p2, p3, num_px, num_py);

  MatrixX<double> compliance_p1 = CalcElementComplianceRowMatrix(
      mesh_data.first, mesh_data.second, Vector3d(0.0, 0.0, 0.0), 1);
  VectorX<double> pressures = GetPressureIntegrandR(mesh_data.first);

  MatrixX<double> deformation_p1 = compliance_p1 * pressures;
  double results = deformation_p1(0, 0);
  const double expected = 1.080517576210803;
  EXPECT_NEAR(results,
              expected,
              10 * std::numeric_limits<double>::epsilon());
}

/// The expected values for this test are the results by running Matlab. The
/// precision (15 digits) of the expected values if the same as the "long"
/// in Matlab
/// The "expected" results is solved numerically by Matlab functions. Since
/// the pressure filed is linear in the Cartesian coordinate, the results
/// should be the same with the "expected" solution (numerical result).
/// The tolerance here is set to the Matlab function tolerance.
GTEST_TEST(ElementComplianceMatrixTest, IntegrantX) {
  const double length = 1;
  const Vector2d p1(0.0, 0.0);
  const Vector2d p2(length, 0.0);
  const Vector2d p3(0.0, length);
  const int num_px = 3;
  const int num_py = 3;

  const std::pair<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3i>>&
      mesh_data = MeshSquare(p1, p2, p3, num_px, num_py);

  MatrixX<double> compliance_p1 = CalcElementComplianceRowMatrix(
      mesh_data.first, mesh_data.second, Vector3d(0.0, 0.0, 0.0), 1);
  VectorX<double> pressures = GetPressureIntegrandX(mesh_data.first);

  MatrixX<double> deformation_p1 = compliance_p1 * pressures;
  double results = deformation_p1(0, 0);
  const double expected = 0.647794300685803;
  EXPECT_NEAR(results,
              expected,
              10 * 7.640083832722066e-07);
}


}  // namespace
}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake

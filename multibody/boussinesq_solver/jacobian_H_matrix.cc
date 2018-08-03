#include "drake/multibody/boussinesq_solver/jacobian_H_matrix.h"

#include "drake/multibody/boussinesq_solver/compliance_matrix.h"
#include "drake/multibody/shapes/geometry.h"

namespace drake {
namespace multibody {
namespace boussinesq_solver {

using drake::geometry::PenetrationAsTrianglePair;

Eigen::MatrixXd CalcJacobianHMatrix(
    std::vector<PenetrationAsTrianglePair<double>>) {

  Eigen::MatrixXd jacobian_H_matrix;


  return jacobian_H_matrix;
}


}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake

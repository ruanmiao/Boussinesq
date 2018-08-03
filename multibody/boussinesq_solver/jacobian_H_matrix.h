#pragma once

#include <vector>

#include "drake/common/eigen_types.h"
#include "drake/geometry/mesh_query/mesh_query.h"

namespace drake {
namespace multibody {
namespace boussinesq_solver {

Eigen::MatrixXd CalcJacobianHMatrix(
    std::vector<drake::geometry::PenetrationAsTrianglePair<double>>);


}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake

#pragma once

#include <vector>

#include "drake/common/eigen_types.h"
#include "drake/geometry/mesh_query/mesh_query.h"

namespace drake {
namespace multibody {
namespace boussinesq_solver {

Eigen::MatrixXd CalcJacobianHMatrix(
    const std::vector<drake::geometry::PenetrationAsTrianglePair<double>>&
    queries,
    const std::vector<Vector3<double>>& patch_A,
    const std::vector<Vector3<double>>& patch_B,
    int patch_A_index,
    int patch_B_index);


}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake

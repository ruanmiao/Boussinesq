#include "drake/multibody/boussinesq_solver/jacobian_H_matrix.h"

#include <limits>

#include "drake/multibody/shapes/geometry.h"
#include "drake/common/drake_assert.h"

#include <iostream>
#define PRINT_VAR(a) std::cout << #a": " << a << std::endl;

namespace drake {
namespace multibody {
namespace boussinesq_solver {

using drake::geometry::PenetrationAsTrianglePair;

Eigen::MatrixXd CalcJacobianHMatrix(
    const std::vector<PenetrationAsTrianglePair<double>>& queries,
    const std::vector<Vector3<double>>& patch_A,
    const std::vector<Vector3<double>>& patch_B,
    int patch_A_index,
    int patch_B_index,
    double young_modulus_star_A,
    double young_modulus_star_B) {
  const int num_queries = queries.size();
  const int patch_A_size = patch_A.size();
  const int patch_B_size = patch_B.size();
  const int num_nodes = patch_A_size + patch_B_size;
  Eigen::MatrixXd jacobian_H_matrix =
      MatrixX<double>::Zero(num_queries, num_nodes);

  for (int i_query = 0; i_query < num_queries; ++i_query) {
    const PenetrationAsTrianglePair<double>& query = queries[i_query];

    Vector3<double> p_AtoB_W = query.p_WoBs_W - query.p_WoAs_W;

    DRAKE_ASSERT(((query.meshA_index == patch_A_index) &&
        (query.meshB_index == patch_B_index)) ||
                 ((query.meshB_index == patch_A_index) &&
                 (query.meshA_index == patch_B_index)));

    //TODO: (mengyao) Think about how to fix the case when distance is 0
    if(fabs(p_AtoB_W.norm()) < std::numeric_limits<double>::epsilon()) {
      continue;
    }
    Vector3<double> n_AtoB_W = p_AtoB_W / p_AtoB_W.norm();

    auto sign = [](double x) {
      // std::copysign(x, y):
      //   Returns a value with the magnitude of x and the sign of y.
      return std::copysign(1.0, x);
    };

    const Vector3<double> that = -sign(query.signed_distance) * n_AtoB_W;

    RowVectorX<double> S_A = RowVectorX<double>::Zero(num_nodes);
    RowVectorX<double> S_B = RowVectorX<double>::Zero(num_nodes);

    const Vector2<int> mesh_start_index(0, patch_A_size);
    for (int local_node = 0; local_node < 3; ++local_node) {
      const int nodeA_index =
          mesh_start_index[query.meshA_index] + query.triangleA[local_node];
      const int nodeB_index =
          mesh_start_index[query.meshB_index] + query.triangleB[local_node];
      DRAKE_DEMAND(nodeA_index < num_nodes);
      DRAKE_DEMAND(nodeB_index < num_nodes);
      S_A(nodeA_index) = query.barycentric_A[local_node];
      S_B(nodeB_index) = query.barycentric_B[local_node];
    }

    // The original normals on the underarmed surfaces.
    const auto& normalA0_W = query.normal_A_W;
    const auto& normalB0_W = query.normal_B_W;

    // WEIGHTED NORMALS FORMULATION.
    (void) young_modulus_star_A;
    (void) young_modulus_star_B;
#if 0
    // Approximation of the normal on the deformed contact patch on body A side.
    // The normal on B points in the opposite direction.
    const Vector3<double> normalA_W =
        (young_modulus_star_B * normalA0_W -
         young_modulus_star_A * normalB0_W).normalized();
        // Modified formulation using an approximate normal such that nB = -nA.
    jacobian_H_matrix.row(i_query) = that.dot(normalA_W) * (S_A + S_B);
#endif

    // ORIGINAL FORMULATION using the undeformed surfaces's normals
    jacobian_H_matrix.row(i_query) =
        that.dot(normalA0_W) * S_A - that.dot(normalB0_W) * S_B;
  }

  return jacobian_H_matrix;
}



}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake

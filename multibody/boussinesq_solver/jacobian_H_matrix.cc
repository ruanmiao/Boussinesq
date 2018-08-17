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

  double young_modulus_star = 1 /
      (1 / young_modulus_star_A + 1 / young_modulus_star_B);

  double ratio_A = (1/young_modulus_star_A) / (1/young_modulus_star);
  double ratio_B = (1/young_modulus_star_B) / (1/young_modulus_star);


  for (int i_query = 0; i_query < num_queries; ++i_query) {
    const PenetrationAsTrianglePair<double>& query = queries[i_query];

    Vector3<double> p_AtoB_W = query.p_WoBs_W - query.p_WoAs_W;

//    PRINT_VAR(query.meshA_index);
//    PRINT_VAR(patch_A_index);
//    PRINT_VAR(query.meshB_index);
//    PRINT_VAR(patch_B_index);

    DRAKE_ASSERT(((query.meshA_index == patch_A_index) &&
        (query.meshB_index == patch_B_index)) ||
                 ((query.meshB_index == patch_A_index) &&
                 (query.meshA_index == patch_B_index)));

    //TODO: (mengyao) Think about how to fix the case when distance is 0
    if(fabs(p_AtoB_W.norm()) < std::numeric_limits<double>::epsilon()) {
      continue;
    }
    Vector3<double> n_AtoB_W = p_AtoB_W / p_AtoB_W.norm();

    MatrixX<double> S_A = MatrixX<double>::Zero(1, num_nodes);
    MatrixX<double> S_B = MatrixX<double>::Zero(1, num_nodes);

    // Check whether this query is A to B or B to A;
    if (query.meshA_index == patch_A_index) {
      for (int i_node = 0; i_node < query.triangleA.size(); i_node++) {
        DRAKE_ASSERT(query.triangleA[i_node] < patch_A_size);
        DRAKE_ASSERT(patch_A_size + query.triangleB[i_node] < num_nodes);

        S_A(0, query.triangleA[i_node]) = query.barycentric_A[i_node];
        S_B(0, patch_A_size + query.triangleB[i_node]) = query.barycentric_B[i_node];
      }
    } else {
      for (int i_node = 0; i_node < query.triangleA.size(); i_node++) {
        DRAKE_ASSERT(query.triangleB[i_node] < patch_A_size);
        DRAKE_ASSERT(patch_A_size + query.triangleA[i_node] < num_nodes);

        S_A(0, patch_A_size + query.triangleA[i_node]) = query.barycentric_A[i_node];
        S_B(0, query.triangleB[i_node]) = query.barycentric_B[i_node];
      }
    }

    (void) ratio_A;
    (void) ratio_B;

    jacobian_H_matrix.row(i_query) = -((n_AtoB_W.dot(query.normal_B_W) * S_B -
        n_AtoB_W.dot(query.normal_A_W) * S_A));


//    jacobian_H_matrix.row(i_query) =
//        -((n_AtoB_W.dot(query.normal_B_W) * S_B * ratio_B -
//        n_AtoB_W.dot(query.normal_A_W) * S_A * ratio_A));

  }
  return jacobian_H_matrix;
}



}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake

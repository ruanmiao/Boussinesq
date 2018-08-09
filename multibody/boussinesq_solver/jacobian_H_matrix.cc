#include "drake/multibody/boussinesq_solver/jacobian_H_matrix.h"

#include <limits>

#include "drake/multibody/boussinesq_solver/compliance_matrix.h"
#include "drake/multibody/shapes/geometry.h"
#include "drake/common/drake_assert.h"

namespace drake {
namespace multibody {
namespace boussinesq_solver {

using drake::geometry::PenetrationAsTrianglePair;

Eigen::MatrixXd CalcJacobianHMatrix(
    const std::vector<PenetrationAsTrianglePair<double>>& queries,
    const std::vector<Vector3<double>>& patch_A,
    const std::vector<Vector3<double>>& patch_B,
    int patch_A_index,
    int patch_B_index) {
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

    jacobian_H_matrix.row(i_query) = (n_AtoB_W.dot(query.normal_B_W) * S_B -
        n_AtoB_W.dot(query.normal_A_W) * S_A);
  }


//  std::ofstream debug_file("examining_values.txt");
//
//  debug_file << "query - 9th: " << std::endl;
//
//  PenetrationAsTrianglePair<double> query0 = queries[8];
//  debug_file << "normal_A: " << query0.normal_A_W(0)
//                   << ", " << query0.normal_A_W(1)
//                   << ", " << query0.normal_A_W(2)
//                   << std::endl;
//  debug_file << "pos on A in world frame: " << query0.p_WoAs_W(0)
//             << ", " << query0.p_WoAs_W(1)
//             << ", " << query0.p_WoAs_W(2)
//             << std::endl;
//
//  debug_file << "triangle on A: " << query0.triangleA(0)
//             << ", " << query0.triangleA(1)
//             << ", " << query0.triangleA(2)
//             << std::endl;
//
//
//  debug_file << "normal_B: " << query0.normal_B_W(0)
//                   << ", " << query0.normal_B_W(1)
//                   << ", " << query0.normal_B_W(2)
//                   << std::endl;
//
//  debug_file << "pos on B in world frame: " << query0.p_WoBs_W(0)
//             << ", " << query0.p_WoBs_W(1)
//             << ", " << query0.p_WoBs_W(2)
//             << std::endl;
//
//
//  debug_file.close();




  return jacobian_H_matrix;
}


}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake

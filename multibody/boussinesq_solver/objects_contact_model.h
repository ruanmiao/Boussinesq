#pragma once

#include "drake/geometry/mesh_query/mesh_query.h"
#include "drake/multibody/multibody_tree/math/spatial_force.h"

namespace drake {
namespace multibody {
namespace boussinesq_solver {

// TODO (mengyao): remove the variables for debugging use
template <typename T>
struct BoussinesqContactModelResults {
  SpatialForce<T> F_Ao_W;
  SpatialForce<T> F_Bo_W;
  std::unique_ptr<geometry::mesh_query::Mesh<T>> object_A_patch;
  std::unique_ptr<geometry::mesh_query::Mesh<T>> object_B_patch;
  VectorX<double> pressure_patch_A;
  VectorX<double> pressure_patch_B;
  VectorX<double> deformation_patch_A;
  VectorX<double> deformation_patch_B;
  VectorX<double> phi0;
  MatrixX<double> H;

  std::unique_ptr<std::vector<geometry::PenetrationAsTrianglePair<double>>> contact_results;
  VectorX<double> kkt_multipliers;
};

// TODO: Complete the comment here
std::unique_ptr<BoussinesqContactModelResults<double>>
CalcContactSpatialForceBetweenMeshes(
    const geometry::mesh_query::Mesh<double>& object_A,
    const Eigen::Isometry3d& X_WA,
    const double young_modulus_star_A,
    const geometry::mesh_query::Mesh<double>& object_B,
    const Eigen::Isometry3d& X_WB,
    const double young_modulus_star_B,
    double sigma,
    bool press_in = true);



std::unique_ptr<BoussinesqContactModelResults<double>>
CalcContactSpatialForceBetweenMeshesByPressure(
    const geometry::mesh_query::Mesh<double>& object_A,
    const Eigen::Isometry3d& X_WA,
    const double young_modulus_star_A,
    const geometry::mesh_query::Mesh<double>& object_B,
    const Eigen::Isometry3d& X_WB,
    const double young_modulus_star_B,
    double sigma,
    bool press_in = true);



}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake





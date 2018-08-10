#pragma once

#include "drake/geometry/mesh_query/mesh_query.h"

namespace drake {
namespace multibody {
namespace boussinesq_solver {
namespace {

// TODO: Complete the comment here
bool CalcSpherePlaneModel();

// TODO: Complete the comment here
SpatialForce<double> CalcContactSpatialForceBetweenMeshes(
    const geometry::mesh_query::Mesh<double>& object_A,
    const Isometry3<double>& X_WA,
    const double young_modulus_star_A,
    std::string type_A,
    const geometry::mesh_query::Mesh<double>& object_B,
    const Isometry3<double>& X_WB,
    const double young_modulus_star_B, double sigma) const;

}
}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake





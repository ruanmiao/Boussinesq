#include <memory>
#include <vector>
#include <limits>

#include <gtest/gtest.h>
#include <gflags/gflags.h>

#include "drake/common/find_resource.h"
#include "drake/geometry/mesh_query/vtk_io.h"
#include "drake/geometry/mesh_query/mesh_query.h"
#include "drake/multibody/shapes/geometry.h"
#include "drake/multibody/boussinesq_solver/jacobian_H_matrix.h"
#include "drake/multibody/boussinesq_solver/compliance_matrix.h"
#include "drake/multibody/boussinesq_solver/objects_contact_model.h"
#include "drake/solvers/moby_lcp_solver.h"

#include <iostream>
#define PRINT_VAR(a) std::cout << #a": " << a << std::endl;

namespace drake {
namespace multibody {
namespace boussinesq_solver {
namespace {

using Eigen::AngleAxisd;
using Eigen::Translation3d;
using Eigen::Isometry3d;
using Eigen::Vector3d;

using geometry::mesh_query::Mesh;
using geometry::mesh_query::LoadMeshFromObj;
using geometry::mesh_query::OutputMeshToVTK;
using geometry::mesh_query::AppendCellCenteredVectorFieldToVTK;
using geometry::mesh_query::AppendCellCenteredVectorFieldToVTK;
using geometry::PenetrationAsTrianglePair;
using geometry::mesh_query::OutputSegmentsToVTK;

GTEST_TEST(ExampleTest, SpherePlane) {

  const bool flip_normals = true;

  // Load mesh for a sphere.
  std::unique_ptr<Mesh<double>> sphere = LoadMeshFromObj(
      "drake/multibody/boussinesq_solver/test/sphere.obj", flip_normals);
  sphere->mesh_index = 0;

  std::unique_ptr<Mesh<double>> plane = LoadMeshFromObj(
      "drake/multibody/boussinesq_solver/test/plane.obj", flip_normals);
  plane->mesh_index = 1;

  const double radius = 1;
  const double penetration = 0.1;
  const double z_WSo = radius - penetration;
  double sigma = 0.1;

  const double young_modulus_star_A = 1.0;
  const double young_modulus_star_B = std::numeric_limits<double>::infinity();

  // Place sphere a "penetration" distance below z = 0.
  // Apply an arbirary rotation for testing.
  Isometry3d X_WSphere{Translation3d{Vector3d(0, 0, z_WSo)}};
  X_WSphere.linear() = MatrixX<double>::Identity(3, 3);

  // The top of the plane is placed at z = 0
  Isometry3d X_WPlane{Translation3d{Vector3d(0, 0, 0.0)}};
  X_WPlane.linear() = MatrixX<double>::Identity(3, 3);


  std::unique_ptr<BoussinesqContactModelResults<double>> boussinesq_results =
      CalcContactSpatialForceBetweenMeshes(
      *sphere, X_WSphere, young_modulus_star_A,
      *plane, X_WPlane, young_modulus_star_B,
      sigma);

  PRINT_VAR(boussinesq_results->F_Ao_W);
  PRINT_VAR(boussinesq_results->F_Bo_W);

}

}
}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake





#include <memory>
#include <vector>

#include <gtest/gtest.h>
#include <gflags/gflags.h>

#include "drake/common/find_resource.h"
#include "drake/geometry/mesh_query/vtk_io.h"
#include "drake/geometry/mesh_query/mesh_query.h"
#include "drake/multibody/shapes/geometry.h"
#include "drake/multibody/boussinesq_solver/jacobian_H_matrix.h"
#include "drake/multibody/boussinesq_solver/compliance_matrix.h"
#include "drake/solvers/moby_lcp_solver.h"

#include <iostream>
#define PRINT_VAR(a) std::cout << #a": " << a << std::endl;

namespace drake {
namespace multibody {
namespace boussinesq_solver {

using Eigen::AngleAxisd;
using Eigen::Translation3d;
using Eigen::Isometry3d;
using Eigen::Vector3d;

using geometry::mesh_query::Mesh;
using geometry::mesh_query::LoadMeshFromObj;
using geometry::mesh_query::OutputMeshToVTK;
using geometry::mesh_query::AppendCellCenteredVectorFieldToVTK;
using geometry::mesh_query::AppendNodeCenteredVectorFieldToVTK;
using geometry::PenetrationAsTrianglePair;
using geometry::mesh_query::OutputSegmentsToVTK;

// A great portion of this function is copied from the
// geometry::mesh_query::Domain() implemented by Alejandro
//int DoSphereEllipsoidModel() {
//  // Load mesh for a sphere.
//  std::unique_ptr<Mesh<double>> sphere = LoadMeshFromObj(
//      "drake/geometry/mesh_query/examples/sphere.obj");
//  sphere->mesh_index = 0;
//
//  // Load mesh for an ellipsoid.
//  std::unique_ptr<Mesh<double>> ellipsoid = LoadMeshFromObj(
//      "drake/geometry/mesh_query/examples/ellipsoid.obj");
//  ellipsoid->mesh_index = 1;
//
//  // For this example normals in the obj point inward and therefore we flip
//  // them.
//  FlipNormals(sphere.get());
//  FlipNormals(ellipsoid.get());
//
//  const double radius = 1.0;
//  const double penetration = 0.05;
//  const double z_WSo = radius - penetration;
//
//  // Place sphere a "penetration" distance below z = 0.
//  // Apply an arbirary rotation for testing.
//  Isometry3d X_WSphere{Translation3d{Vector3d(0, 0, z_WSo)}};
//  X_WSphere.linear() =
//      (AngleAxisd(M_PI / 3., Vector3d::UnitX()) *
//       AngleAxisd(M_PI / 5., Vector3d::UnitY()) *
//       AngleAxisd(3. * M_PI / 8., Vector3d::UnitZ())).matrix();
//
//  // The top of the ellipsoid is placed at z = 0
//  Isometry3d X_WEllipsoid{Translation3d{Vector3d(0, 0, -1.0)}};
//  X_WEllipsoid.linear() = AngleAxisd(M_PI_2, Vector3d::UnitX()).matrix();
//
//  // Write mesh and normals to a file.
//  std::ofstream bottom_sphere_file("ellipsoid.vtk");
//  OutputMeshToVTK(bottom_sphere_file, ellipsoid->points_G, ellipsoid->triangles, X_WEllipsoid);
//  AppendCellCenteredVectorFieldToVTK(
//      bottom_sphere_file, "FaceNormals", ellipsoid->face_normals_G, X_WEllipsoid);
//  AppendNodeCenteredVectorFieldToVTK(
//      bottom_sphere_file, "NodeNormals", ellipsoid->node_normals_G, X_WEllipsoid);
//  bottom_sphere_file.close();
//
//  std::ofstream top_sphere_file("sphere.vtk");
//  OutputMeshToVTK(top_sphere_file, sphere->points_G, sphere->triangles, X_WSphere);
//  AppendCellCenteredVectorFieldToVTK(
//      top_sphere_file, "FaceNormals", sphere->face_normals_G, X_WSphere);
//  AppendNodeCenteredVectorFieldToVTK(
//      top_sphere_file, "NodeNormals", sphere->node_normals_G, X_WSphere);
//  top_sphere_file.close();
//
//  // Perform the mesh-mesh query.
//  // The triangles referenced in the pairs are "global indexes" in the original
//  // full meshes.
//  std::vector<PenetrationAsTrianglePair<double>> results = MeshToMeshQuery(
//      X_WEllipsoid, *ellipsoid,
//      X_WSphere, *sphere);
//
//  // This call creates the two patches on each mesh and updates "results" so
//  // that the triangle indexes in each pair are "local indexes" to the patch
//  // meshes.
//  auto patches = MakeLocalPatchMeshes(&results, *ellipsoid, *sphere);
//  std::unique_ptr<Mesh<double>> ellipsoid_patch = std::move(patches.first);
//  std::unique_ptr<Mesh<double>> sphere_patch = std::move(patches.second);
//
//  PRINT_VAR(sphere_patch->points_G.size());
//  PRINT_VAR(ellipsoid_patch->points_G.size());
//
//  OutputMeshToVTK("ellipsoid_patch.vtk",
//                  ellipsoid_patch->points_G, ellipsoid_patch->triangles,
//                  X_WEllipsoid);
//
//  OutputMeshToVTK("sphere_patch.vtk",
//                  sphere_patch->points_G, sphere_patch->triangles,
//                  X_WSphere);
//
//  std::vector<Vector3d> pointsA(results.size());
//  std::transform(results.begin(), results.end(), pointsA.begin(),
//                 [](const PenetrationAsTrianglePair<double>& pair) {
//                   return pair.p_WoAs_W;
//                 });
//
//  std::vector<Vector3d> pointsB(results.size());
//  std::transform(results.begin(), results.end(), pointsB.begin(),
//  [](const PenetrationAsTrianglePair<double>& pair) {
//    return pair.p_WoBs_W;
//  });
//
//  std::vector<Vector3d> normals;
//  for (const auto& result : results) {
//    normals.push_back(result.normal_A_W);
//    normals.push_back(result.normal_B_W);
//  }
//
//  {
//    std::ofstream file("pairs.vtk");
//    OutputSegmentsToVTK(file, pointsA, pointsB);
//    AppendNodeCenteredVectorFieldToVTK(
//        file, "normals", normals);
//    file.close();
//  }
//
//  // Print out patches elements just for verification.
//  std::vector<int> sphere_patch_triangles;
//  std::vector<int> ellipsoid_patch_triangles;
//  for (const auto& pair : results) {
//
//    if (pair.meshA_index == sphere->mesh_index) {
//      sphere_patch_triangles.push_back(pair.triangleA_index);
//    } else {
//      ellipsoid_patch_triangles.push_back(pair.triangleA_index);
//    }
//
//    if (pair.meshB_index == sphere->mesh_index) {
//      sphere_patch_triangles.push_back(pair.triangleB_index);
//    } else {
//      ellipsoid_patch_triangles.push_back(pair.triangleB_index);
//    }
//  }
//
//  for (auto sphere_triangle : sphere_patch_triangles) {
//    PRINT_VAR(sphere_triangle);
//  }
//
//  for (auto ellipsoid_triangle : ellipsoid_patch_triangles) {
//    PRINT_VAR(ellipsoid_triangle);
//  }
//
//
//  // Test jacobian_H_matrix:
//  const int size_sphere = sphere_patch->points_G.size();
//  const int size_ellipsoid = ellipsoid_patch->points_G.size();
//  const int E_modulus = 1;
//
//  //TODO: (mengyao) figure out a standard way to gety h0_gap
//  VectorX<double> h0_gap_sphere(size_sphere);
//  VectorX<double> h0_gap_ellipsoid(size_ellipsoid);
//
//  MatrixX<double> compliance_sphere = CalcComplianceMatrix(
//      sphere_patch->points_G, sphere_patch->triangles,
//      1 / (E_modulus * M_PI));
//  MatrixX<double> compliance_ellipsoid = CalcComplianceMatrix(
//      ellipsoid_patch->points_G, ellipsoid_patch->triangles,
//      1 / (E_modulus * M_PI));
//
//  solvers::MobyLCPSolver<double> moby_LCP_solver;
//  VectorX<double> pressure_sphere;
//  VectorX<double> pressure_ellipsoid;
//  moby_LCP_solver.SolveLcpLemke(compliance_sphere,
//                                h0_gap_sphere, &pressure_sphere);
//  moby_LCP_solver.SolveLcpLemke(compliance_ellipsoid,
//                                h0_gap_ellipsoid, &pressure_ellipsoid);
//
//  VectorX<double> u_deformation(size_sphere + size_ellipsoid);
//  u_deformation.head(size_sphere) = compliance_sphere * pressure_sphere;
//  u_deformation.tail(size_ellipsoid) =
//      compliance_ellipsoid * pressure_ellipsoid;
//
//  MatrixX<double> jacobian_H = CalcJacobianHMatrix(
//      results, sphere_patch->points_G, ellipsoid_patch->points_G,
//      sphere_patch->mesh_index, ellipsoid_patch->mesh_index);
//
//
//
//  return 0;
//}

int DoShperePlaneModel() {

  const bool flip_normals = true;

  // Load mesh for a sphere.
  std::unique_ptr<Mesh<double>> sphere = LoadMeshFromObj(
      "drake/multibody/boussinesq_solver/test/sphere.obj", flip_normals);
  sphere->mesh_index = 0;

  // TODO: (mengyao)Load mesh for an plane.
  std::unique_ptr<Mesh<double>> plane = LoadMeshFromObj(
      "drake/multibody/boussinesq_solver/test/plane.obj", flip_normals);
  plane->mesh_index = 1;

  const double radius = 1;
  const double penetration = 0.1;
  const double z_WSo = radius - penetration;

  // Place sphere a "penetration" distance below z = 0.
  // Apply an arbirary rotation for testing.
  Isometry3d X_WSphere{Translation3d{Vector3d(0, 0, z_WSo)}};
  X_WSphere.linear() = MatrixX<double>::Identity(3, 3);

  // The top of the plane is placed at z = 0
  Isometry3d X_WPlane{Translation3d{Vector3d(0, 0, 0.0)}};
  X_WPlane.linear() = MatrixX<double>::Identity(3, 3);

  // Write mesh and normals to a file.
  std::ofstream bottom_plane_file("plane.vtk");
  OutputMeshToVTK(bottom_plane_file, plane->points_G, plane->triangles, X_WPlane);
  AppendCellCenteredVectorFieldToVTK(
      bottom_plane_file, "FaceNormals", plane->face_normals_G, X_WPlane);
  AppendNodeCenteredVectorFieldToVTK(
      bottom_plane_file, "NodeNormals", plane->node_normals_G, X_WPlane);
  bottom_plane_file.close();

  std::ofstream top_sphere_file("sphere.vtk");
  OutputMeshToVTK(top_sphere_file, sphere->points_G, sphere->triangles, X_WSphere);
  AppendCellCenteredVectorFieldToVTK(
      top_sphere_file, "FaceNormals", sphere->face_normals_G, X_WSphere);
  AppendNodeCenteredVectorFieldToVTK(
      top_sphere_file, "NodeNormals", sphere->node_normals_G, X_WSphere);
  top_sphere_file.close();

  // Perform the mesh-mesh query.
  // The triangles referenced in the pairs are "global indexes" in the original
  // full meshes.
  std::vector<PenetrationAsTrianglePair<double>> results = MeshToMeshQuery(
      X_WPlane, *plane,
      X_WSphere, *sphere);

  // This call creates the two patches on each mesh and updates "results" so
  // that the triangle indexes in each pair are "local indexes" to the patch
  // meshes.
  auto patches = MakeLocalPatchMeshes(&results, *plane, *sphere);
  std::unique_ptr<Mesh<double>> plane_patch = std::move(patches.first);
  std::unique_ptr<Mesh<double>> sphere_patch = std::move(patches.second);

  PRINT_VAR(sphere_patch->points_G.size());
  PRINT_VAR(plane_patch->points_G.size());

  OutputMeshToVTK("plane_patch.vtk",
                  plane_patch->points_G, plane_patch->triangles,
                  X_WPlane);

  OutputMeshToVTK("sphere_patch.vtk",
                  sphere_patch->points_G, sphere_patch->triangles,
                  X_WSphere);

  std::vector<Vector3d> pointsA(results.size());
  std::transform(results.begin(), results.end(), pointsA.begin(),
                 [](const PenetrationAsTrianglePair<double>& pair) {
                   return pair.p_WoAs_W;
                 });

  std::vector<Vector3d> pointsB(results.size());
  std::transform(results.begin(), results.end(), pointsB.begin(),
                 [](const PenetrationAsTrianglePair<double>& pair) {
                   return pair.p_WoBs_W;
                 });

  std::vector<Vector3d> normals;
  for (const auto& result : results) {
    normals.push_back(result.normal_A_W);
    normals.push_back(result.normal_B_W);
  }

  {
    std::ofstream file("pairs.vtk");
    OutputSegmentsToVTK(file, pointsA, pointsB);
    AppendNodeCenteredVectorFieldToVTK(
        file, "normals", normals);
    file.close();
  }

  // Print out patches elements just for verification.
  std::vector<int> sphere_patch_triangles;
  std::vector<int> plane_patch_triangles;
  for (const auto& pair : results) {

    if (pair.meshA_index == sphere->mesh_index) {
      sphere_patch_triangles.push_back(pair.triangleA_index);
    } else {
      plane_patch_triangles.push_back(pair.triangleA_index);
    }

    if (pair.meshB_index == sphere->mesh_index) {
      sphere_patch_triangles.push_back(pair.triangleB_index);
    } else {
      plane_patch_triangles.push_back(pair.triangleB_index);
    }
  }

  for (auto sphere_triangle : sphere_patch_triangles) {
    PRINT_VAR(sphere_triangle);
  }

  for (auto plane_triangle : plane_patch_triangles) {
    PRINT_VAR(plane_triangle);
  }




  // Test jacobian_H_matrix:
  const int num_queries = results.size();
  const int size_sphere = sphere_patch->points_G.size();
  const int size_plane = plane_patch->points_G.size();

  VectorX<double> u_deformation =
      VectorX<double>::Zero(size_sphere + size_plane);

  std::ofstream deformation_file("deformed_distance.txt");

  deformation_file << "deformation u: "
                   << std::endl;

  for (int i_pos = 0; i_pos < size_sphere; ++i_pos) {
    Vector3d pos = X_WSphere * sphere_patch->points_G[i_pos];

    u_deformation(i_pos) = std::max(0.0,
        radius * (0.0 - pos(2)) / (z_WSo - pos(2)));

    deformation_file << fmt::format("{:.8f}\n",
                                    u_deformation(i_pos));
  }
  deformation_file << std::endl;



  MatrixX<double> jacobian_H = CalcJacobianHMatrix(
      results, sphere_patch->points_G, plane_patch->points_G,
      sphere_patch->mesh_index, plane_patch->mesh_index);

  deformation_file << "query: " << std::endl;

  PenetrationAsTrianglePair<double> query0 = results[31];

  PRINT_VAR("MeshA");
  PRINT_VAR(query0.meshA_index);
  PRINT_VAR(query0.triangleA_index);
  PRINT_VAR(query0.triangleA.transpose());
  PRINT_VAR(query0.p_WoAs_W.transpose());
  PRINT_VAR(query0.normal_A_W.transpose());
  PRINT_VAR(query0.barycentric_A.transpose());

  PRINT_VAR("MeshB");
  PRINT_VAR(query0.meshB_index);
  PRINT_VAR(query0.triangleB_index);
  PRINT_VAR(query0.triangleB.transpose());
  PRINT_VAR(query0.p_WoBs_W.transpose());
  PRINT_VAR(query0.normal_B_W.transpose());
  PRINT_VAR(query0.barycentric_B.transpose());

  deformation_file << "pos on A in world frame: " << query0.p_WoAs_W(0)
                   << ", " << query0.p_WoAs_W(1)
                   << ", " << query0.p_WoAs_W(2)
                   << std::endl;
  deformation_file << std::endl;

  deformation_file << "triangle on A: " << query0.triangleA(0)
                   << ", " << query0.triangleA(1)
                   << ", " << query0.triangleA(2)
                   << std::endl;
  deformation_file << std::endl;

  deformation_file << "barycentric_A: " << query0.barycentric_A(0)
                   << ", " << query0.barycentric_A(1)
                   << ", " << query0.barycentric_A(2)
                   << std::endl;
  deformation_file << std::endl;


  Vector3<int> triangleA = query0.triangleA;

  deformation_file << "With 0-base index convention: " << std::endl;
  for (int i_pos = 0; i_pos < 3; ++i_pos) {
    int node_index = triangleA(i_pos);
    Vector3d pos = X_WSphere * sphere_patch->points_G[node_index];
    deformation_file << "deformation information for the " << node_index + 1
                     <<  "th node on A, with index " << node_index
                     << " in the mesh (triangle A) : " << std::endl;
    deformation_file << "pos on A: " << pos(0)
                     << ", " << pos(1)
                     << ", " << pos(2)
                     << std::endl;
    deformation_file << "ration: " << (0.0 - pos(2)) / (z_WSo - pos(2))
                     << std::endl;

    deformation_file << "******** end of this node in the local coordinate ****"
                        << std::endl;
    deformation_file << std::endl;

  }


//  deformation_file << "normal_B: " << query0.normal_B_W(0)
//                   << ", " << query0.normal_B_W(1)
//                   << ", " << query0.normal_B_W(2)
//                   << std::endl;

  VectorX<double> initial_distance(num_queries);
  for(int i_query = 0; i_query < num_queries; i_query++) {
    initial_distance(i_query) = results[i_query].signed_distance;
  }

  VectorX<double> deformation_queries = jacobian_H * u_deformation;
  VectorX<double> actual_distance = initial_distance + deformation_queries;



  deformation_file << "initial_distance, deformation term, actual distance, "
                      << "actual distance / max_error: "
                      << std::endl;

  double error_max = (1 - std::sqrt(1 - std::pow(0.1 / std::sqrt(3), 2)));
  deformation_file << "max_error: "
                   << error_max
                   << std::endl;

  for (int i_query = 0; i_query < num_queries; i_query++) {
    deformation_file << fmt::format("{:.8f}, {:.8f}, {:.8f}, {:.8f} \n",
                                    initial_distance(i_query),
                                    deformation_queries(i_query),
                                    actual_distance(i_query),
                                    actual_distance(i_query) / error_max);
  }

  deformation_file.close();

  return 0;
}




GTEST_TEST(JacobianHMatrixTest, SpherePlane) {
  drake::multibody::boussinesq_solver::DoShperePlaneModel();
}




}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake





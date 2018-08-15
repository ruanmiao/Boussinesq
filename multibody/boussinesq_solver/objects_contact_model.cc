#include "drake/multibody/boussinesq_solver/objects_contact_model.h"

#include <memory>
#include <vector>
#include <iostream>
#include <string>

#include "drake/multibody/boussinesq_solver/math_helper.h"
#include "drake/multibody/boussinesq_solver/jacobian_H_matrix.h"
#include "drake/multibody/boussinesq_solver/compliance_matrix.h"
#include "drake/solvers/moby_lcp_solver.h"

#define PRINT_VAR(a) std::cout << #a": " << a << std::endl;

namespace drake {
namespace multibody {
namespace boussinesq_solver {

using Eigen::AngleAxisd;
using Eigen::Translation3d;
using Eigen::Isometry3d;
using Eigen::Vector3d;

using geometry::mesh_query::Mesh;
using geometry::PenetrationAsTrianglePair;
using geometry::mesh_query::MakeLocalPatchMeshes;

// B is the one at the bottom
std::unique_ptr<BoussinesqContactModelResults<double>>
CalcContactSpatialForceBetweenMeshes(
    const Mesh<double>& object_A,
    const Isometry3d& X_WA,
    const double young_modulus_star_A,
    const Mesh<double>& object_B,
    const Isometry3d& X_WB,
    const double young_modulus_star_B, double sigma) {

  auto owned_results =
      std::make_unique<std::vector<PenetrationAsTrianglePair<double>>>();
  std::vector<PenetrationAsTrianglePair<double>>& results = *owned_results;

  results = MeshToMeshQuery(X_WA, object_A, X_WB, object_B, sigma);

  // This call creates the two patches on each mesh and updates "results" so
  // that the triangle indexes in each pair are "local indexes" to the patch
  // meshes.
  auto patches = MakeLocalPatchMeshes(&results, object_A, object_B);
  std::unique_ptr<Mesh<double>> object_A_patch = std::move(patches.first);
  std::unique_ptr<Mesh<double>> object_B_patch = std::move(patches.second);

  const int num_phis = results.size();
  const int num_nodes_A = object_A_patch->points_G.size();
  const int num_nodes_B = object_B_patch->points_G.size();
  const int num_nodes = num_nodes_A + num_nodes_B;


  const MatrixX<double> H = CalcJacobianHMatrix(
      results, object_A_patch->points_G, object_B_patch->points_G,
      object_A_patch->mesh_index, object_B_patch->mesh_index);

  DRAKE_DEMAND(H.rows() == num_phis);
  DRAKE_DEMAND(H.cols() == num_nodes);

  const MatrixX<double> C_A = CalcComplianceMatrix(
      object_A_patch->points_G, object_A_patch->triangles,
      1 / (young_modulus_star_A * M_PI));

  const MatrixX<double> C_B = CalcComplianceMatrix(
      object_B_patch->points_G, object_B_patch->triangles,
      1 / (young_modulus_star_B * M_PI));

  MatrixX<double> C = MatrixX<double>::Zero(num_nodes, num_nodes);
  C.block(0, 0, num_nodes_A, num_nodes_A) = C_A;
  C.block(num_nodes_A, num_nodes_A, num_nodes_B, num_nodes_B) = C_B;

  // Delassus operator.
  const MatrixX<double> W = H * C * H.transpose();

  VectorX<double> phi0(num_phis);
  for (int i = 0; i < num_phis; ++i) {
    phi0(i) = results[i].signed_distance;
  }

  VectorX<double> kkt_multipliers(num_phis);
  solvers::MobyLCPSolver<double> moby_LCP_solver;
  const bool solved = moby_LCP_solver.SolveLcpLemke(W, phi0, &kkt_multipliers);
  DRAKE_DEMAND(solved);

  // Estimate pressure from KKT multipliers.
  // Integrate pressure to get the spatial force.
  // CODE FOR INTEGRATING FORCES HERE.


  VectorX<double> phi_deformation = H * C * H.transpose() * kkt_multipliers;


// TODO: (mengyao) modify it later

  VectorX<double> pressure = H.transpose() * kkt_multipliers;

  VectorX<double> u_deformation = C * pressure;
//  for (int i = 0; i < pressure.rows(); i++) {
//    if(pressure(i) < 0.0) {
//      pressure(i) = 0.0;
//    }
//  }

//  PRINT_VAR(phi_deformation.rows());
//  PRINT_VAR((C * phi_deformation).rows());
//  PRINT_VAR(C.rows());
//  PRINT_VAR(C.cols());


  SpatialForce<double> F_Ao_W;
  F_Ao_W.SetZero();
  SpatialForce<double> F_Bo_W;
  F_Bo_W.SetZero();

  for (int i = 0; i < num_nodes_A; ++i) {

    // TODO: (mengyao: to be clean-up)
    double area = object_A_patch->node_areas[i];
    F_Ao_W.translational() += pressure(i) * area
        * object_A_patch->node_normals_G[i];
  }

  for (int i = 0; i < num_nodes_B; ++i) {

    // TODO: (mengyao: to be clean-up)
    double area = object_B_patch->node_areas[i];
    F_Bo_W.translational() += pressure(num_nodes_A + i) * area
        * object_B_patch->node_normals_G[i];
  }



//  int num_triangles_A = object_A_patch->triangles.size();
//  int num_triangles_B = object_B_patch->triangles.size();
//
//  for (int i = 0; i < num_triangles_A; ++i) {
//
//    // TODO: (mengyao: to be clean-up)
//    Vector3<int>& triangle = object_A_patch->triangles[i];
//    Vector3d area_vec =
//        CalcTriangleArea(object_A_patch->points_G[triangle(0)],
//                         object_A_patch->points_G[triangle(1)],
//                         object_A_patch->points_G[triangle(2)]);
//
//    F_Ao_W.translational() +=
//        (pressure(triangle(0)) * object_A_patch->node_normals_G[triangle(0)] +
//            pressure(triangle(1)) * object_A_patch->node_normals_G[triangle(1)]
//            + pressure(triangle(2))*
//                object_A_patch->node_normals_G[triangle(2)])
//            * area_vec.norm() / 3;
//
////    F_Ao_W.translational() += (pressure(triangle(0)) + pressure(triangle(1)) +
////        pressure(triangle(2))) * area_vec.norm() / 3
////        * object_A_patch->face_normals_G[i];
//  }
//
//  for (int i = 0; i < num_triangles_B; ++i) {
//
//    // TODO: (mengyao: to be clean-up)
//    Vector3<int>& triangle = object_B_patch->triangles[i];
//    Vector3d area_vec =
//        CalcTriangleArea(object_B_patch->points_G[triangle(0)],
//                         object_B_patch->points_G[triangle(1)],
//                         object_B_patch->points_G[triangle(2)]);
//
//    F_Bo_W.translational() += (pressure(num_nodes_A + triangle(0)) *
//        object_B_patch->node_normals_G[triangle(0)]
//        + pressure(num_nodes_A + triangle(1))
//            * object_B_patch->node_normals_G[triangle(1)]
//        + pressure(num_nodes_A + triangle(2)) *
//            object_B_patch->node_normals_G[triangle(2)]) * area_vec.norm() / 3;
//
////    F_Bo_W.translational() += (pressure(num_nodes_A + triangle(0))
////        + pressure(num_nodes_A + triangle(1))
////        + pressure(num_nodes_A + triangle(2))) * area_vec.norm() / 3
////        * object_B_patch->face_normals_G[i];
//  }






  // I am done using the patches.
  std::unique_ptr<BoussinesqContactModelResults<double>> boussinesq_results =
      std::make_unique<BoussinesqContactModelResults<double>>();

  boussinesq_results->F_Ao_W = F_Ao_W;
  boussinesq_results->F_Bo_W = F_Bo_W;
  boussinesq_results->object_A_patch = std::move(object_A_patch);  // now object_A_patch is nullptr.
  boussinesq_results->object_B_patch = std::move(object_B_patch);

  boussinesq_results->pressure_patch_A =
      pressure.segment(0, num_nodes_A);
  boussinesq_results->pressure_patch_B =
      pressure.segment(num_nodes_A, num_nodes_B);
  boussinesq_results->deformation_patch_A =
      u_deformation.segment(0, num_nodes_A);
  boussinesq_results->deformation_patch_B =
      u_deformation.segment(num_nodes_A, num_nodes_B);

  boussinesq_results->contact_results = std::move(owned_results);





//  boussinesq_results->kkt_multipliers = kkt_multipliers;
  boussinesq_results->kkt_multipliers = kkt_multipliers;

  return std::move(boussinesq_results);
}



//
//#if 0
//bool CalcSpherePlaneModel() {
//  const bool flip_normals = true;
//
//  // Load mesh for a sphere.
//  std::unique_ptr<Mesh<double>> sphere = LoadMeshFromObj(
//      "drake/multibody/boussinesq_solver/test/sphere.obj", flip_normals);
//  sphere->mesh_index = 0;
//
//  std::unique_ptr<Mesh<double>> plane = LoadMeshFromObj(
//      "drake/multibody/boussinesq_solver/test/plane.obj", flip_normals);
//  plane->mesh_index = 1;
//
//  const double radius = 1;
//  const double penetration = 0.1;
//  const double z_WSo = radius - penetration;
//
//  // Place sphere a "penetration" distance below z = 0.
//  // Apply an arbirary rotation for testing.
//  Isometry3d X_WSphere{Translation3d{Vector3d(0, 0, z_WSo)}};
//  X_WSphere.linear() = MatrixX<double>::Identity(3, 3);
//
//  // The top of the plane is placed at z = 0
//  Isometry3d X_WPlane{Translation3d{Vector3d(0, 0, 0.0)}};
//  X_WPlane.linear() = MatrixX<double>::Identity(3, 3);
//
//  CalcObjectsContactModel(*sphere, X_WSphere, 1.0, "sphere",
//                          *plane, X_WPlane, 1.0, "plane");
//
//  return true;
//}
//#endif
//
//
//#if 0
//  const double radius = 1;
//  const double penetration = 0.1;
//  // TODO: the z_WSo should be deleted later
//  const double z_WSo = radius - penetration;
//
//  std::string mesh_A_name;
//  mesh_A_name.append(type_A);
//  mesh_A_name.append(".vtk");
//
//  std::string mesh_B_name;
//  mesh_B_name.append(type_B);
//  mesh_B_name.append(".vtk");
//
//  // Write mesh and normals to a file.
//  std::ofstream bottom_plane_file(mesh_A_name);
//  OutputMeshToVTK(bottom_plane_file,
//                  object_B->points_G,
//                  object_B->triangles,
//                  X_WB);
//  AppendCellCenteredVectorFieldToVTK(
//      bottom_plane_file, "FaceNormals", object_B->face_normals_G, X_WB);
//  AppendNodeCenteredVectorFieldToVTK(
//      bottom_plane_file, "NodeNormals", object_B->node_normals_G, X_WB);
//  bottom_plane_file.close();
//
//  std::ofstream top_sphere_file(mesh_B_name);
//  OutputMeshToVTK(top_sphere_file, object_A.points_G,
//                  object_A.triangles, X_WA);
//  AppendCellCenteredVectorFieldToVTK(
//      top_sphere_file, "FaceNormals", object_A->face_normals_G, X_WA);
//  AppendNodeCenteredVectorFieldToVTK(
//      top_sphere_file, "NodeNormals", object_A->node_normals_G, X_WA);
//  top_sphere_file.close();
//
//  // Perform the mesh-mesh query.
//  // The triangles referenced in the pairs are "global indexes" in the original
//  // full meshes.
//
//  //TODO:(mengyao) Try to switch the A and B here,
//  // might need to unflip the plane
//  std::vector<PenetrationAsTrianglePair<double>> results = MeshToMeshQuery(
//      X_WB, *object_B,
//      X_WA, *object_A);
//
//
//
//
//
//
//
//
//
//  // This call creates the two patches on each mesh and updates "results" so
//  // that the triangle indexes in each pair are "local indexes" to the patch
//  // meshes.
//  //TODO:(mengyao) Try to switch the A and B here,
//  // might need to unflip the plane
//  auto patches = MakeLocalPatchMeshes(&results, *object_B, *object_A);
//  std::unique_ptr<Mesh<double>> object_B_patch = std::move(patches.first);
//  std::unique_ptr<Mesh<double>> object_A_patch = std::move(patches.second);
//
//  PRINT_VAR(object_A_patch->points_G.size());
//  PRINT_VAR(object_B_patch->points_G.size());
//
//  std::string patch_A_name;
//  patch_A_name.append(type_A);
//  patch_A_name.append("_patch.vtk");
//
//  std::string patch_B_name;
//  patch_B_name.append(type_B);
//  patch_B_name.append("_patch.vtk");
//
//  OutputMeshToVTK(patch_A_name,
//                  object_B_patch->points_G, object_B_patch->triangles,
//                  X_WB);
//
//  OutputMeshToVTK(patch_B_name,
//                  object_A_patch->points_G, object_A_patch->triangles,
//                  X_WA);
//
//  std::vector<Vector3d> pointsA(results.size());
//  std::transform(results.begin(), results.end(), pointsA.begin(),
//                 [](const PenetrationAsTrianglePair<double> &pair) {
//                   return pair.p_WoAs_W;
//                 });
//
//  std::vector<Vector3d> pointsB(results.size());
//  std::transform(results.begin(), results.end(), pointsB.begin(),
//                 [](const PenetrationAsTrianglePair<double> &pair) {
//                   return pair.p_WoBs_W;
//                 });
//
//  std::vector<Vector3d> normals;
//  for (const auto &result : results) {
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
//  (void) material_A;
//  (void) material_B;
//  // Test jacobian_H_matrix:
//  const int num_queries = results.size();
//  const int size_sphere = object_A_patch->points_G.size();
//  const int size_plane = object_B_patch->points_G.size();
//
//  VectorX<double> u_deformation =
//      VectorX<double>::Zero(size_sphere + size_plane);
//
//  for (int i_pos = 0; i_pos < size_sphere; ++i_pos) {
//    Vector3d pos = X_WA * object_A_patch->points_G[i_pos];
//
//    u_deformation(i_pos) = std::max(
//        0.0, radius * (0.0 - pos(2)) / (z_WSo - pos(2)));
//  }
//
//  MatrixX<double> jacobian_H = CalcJacobianHMatrix(
//      results, object_A_patch->points_G, object_B_patch->points_G,
//      object_A_patch->mesh_index, object_B_patch->mesh_index);
//
//  VectorX<double> initial_distance(num_queries);
//  for (int i_query = 0; i_query < num_queries; i_query++) {
//    initial_distance(i_query) = results[i_query].signed_distance;
//  }
//
//  VectorX<double> deformation_queries = jacobian_H * u_deformation;
//  VectorX<double> actual_distance = initial_distance + deformation_queries;
//
//
//
//
//
//  // ********* Debug data ****************
//  // Print out patches elements just for verification.
//  std::vector<int> sphere_patch_triangles;
//  std::vector<int> plane_patch_triangles;
//  for (const auto &pair : results) {
//
//    if (pair.meshA_index == object_A->mesh_index) {
//      sphere_patch_triangles.push_back(pair.triangleA_index);
//    } else {
//      plane_patch_triangles.push_back(pair.triangleA_index);
//    }
//
//    if (pair.meshB_index == object_A->mesh_index) {
//      sphere_patch_triangles.push_back(pair.triangleB_index);
//    } else {
//      plane_patch_triangles.push_back(pair.triangleB_index);
//    }
//  }
//
//  for (auto sphere_triangle : sphere_patch_triangles) {
//    PRINT_VAR(sphere_triangle);
//  }
//
//  for (auto plane_triangle : plane_patch_triangles) {
//    PRINT_VAR(plane_triangle);
//  }
//
//  std::ofstream deformation_file("deformed_distance.txt");
//
//  deformation_file << "deformation u: "
//                   << std::endl;
//
//  for (int i_pos = 0; i_pos < size_sphere; ++i_pos) {
//    deformation_file << fmt::format("{:.8f}\n", u_deformation(i_pos));
//  }
//  deformation_file << std::endl;
//
//  deformation_file << "query: " << std::endl;
//  PenetrationAsTrianglePair<double> query0 = results[31];
//
//  PRINT_VAR("MeshA");
//  PRINT_VAR(query0.meshA_index);
//  PRINT_VAR(query0.triangleA_index);
//  PRINT_VAR(query0.triangleA.transpose());
//  PRINT_VAR(query0.p_WoAs_W.transpose());
//  PRINT_VAR(query0.normal_A_W.transpose());
//  PRINT_VAR(query0.barycentric_A.transpose());
//
//  PRINT_VAR("MeshB");
//  PRINT_VAR(query0.meshB_index);
//  PRINT_VAR(query0.triangleB_index);
//  PRINT_VAR(query0.triangleB.transpose());
//  PRINT_VAR(query0.p_WoBs_W.transpose());
//  PRINT_VAR(query0.normal_B_W.transpose());
//  PRINT_VAR(query0.barycentric_B.transpose());
//
//  deformation_file << "pos on A in world frame: " << query0.p_WoAs_W(0)
//                   << ", " << query0.p_WoAs_W(1)
//                   << ", " << query0.p_WoAs_W(2)
//                   << std::endl;
//  deformation_file << std::endl;
//
//  deformation_file << "triangle on A: " << query0.triangleA(0)
//                   << ", " << query0.triangleA(1)
//                   << ", " << query0.triangleA(2)
//                   << std::endl;
//  deformation_file << std::endl;
//
//  deformation_file << "barycentric_A: " << query0.barycentric_A(0)
//                   << ", " << query0.barycentric_A(1)
//                   << ", " << query0.barycentric_A(2)
//                   << std::endl;
//  deformation_file << std::endl;
//
//  Vector3<int> triangleA = query0.triangleA;
//
//  deformation_file << "With 0-base index convention: " << std::endl;
//  for (int i_pos = 0; i_pos < 3; ++i_pos) {
//    int node_index = triangleA(i_pos);
//    Vector3d pos = X_WA * object_A_patch->points_G[node_index];
//    deformation_file << "deformation information for the " << node_index + 1
//                     << "th node on A, with index " << node_index
//                     << " in the mesh (triangle A) : " << std::endl;
//    deformation_file << "pos on A: " << pos(0)
//                     << ", " << pos(1)
//                     << ", " << pos(2)
//                     << std::endl;
//    deformation_file << "ration: " << (0.0 - pos(2)) / (z_WSo - pos(2))
//                     << std::endl;
//
//    deformation_file << "******** end of this node in the local coordinate ****"
//                     << std::endl;
//    deformation_file << std::endl;
//
//  }
//
//  deformation_file << "initial_distance, deformation term, actual distance, "
//                   << "actual distance / max_error: "
//                   << std::endl;
//
//  double error_max = (1 - std::sqrt(1 - std::pow(0.1 / std::sqrt(3), 2)));
//  deformation_file << "max_error: "
//                   << error_max
//                   << std::endl;
//
//  for (int i_query = 0; i_query < num_queries; i_query++) {
//    deformation_file << fmt::format("{:.8f}, {:.8f}, {:.8f}, {:.8f} \n",
//                                    initial_distance(i_query),
//                                    deformation_queries(i_query),
//                                    actual_distance(i_query),
//                                    actual_distance(i_query) / error_max);
//  }
//
//  deformation_file.close();
//  // ******* Debug End *********
//
//
//  return 0;
//}
//
//#endif

}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake





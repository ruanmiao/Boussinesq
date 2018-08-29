#include <memory>
#include <vector>

#include <gflags/gflags.h>

#include "drake/common/find_resource.h"
#include "drake/geometry/mesh_query/vtk_io.h"
#include "drake/geometry/mesh_query/mesh_query.h"
#include "drake/multibody/shapes/geometry.h"

#include <iostream>
#define PRINT_VAR(a) std::cout << #a": " << a << std::endl;

namespace drake {
namespace geometry {
namespace mesh_query {
namespace {

DEFINE_double(penetration, 0.1,
              "Distance between sphere and ellipsoid. It can be negative (they are separated.)");

DEFINE_double(sigma, 0.1,
              "Queries with distance below sigma are included, including negative distance results.");

using Eigen::AngleAxisd;
using Eigen::Translation3d;
using Eigen::Isometry3d;
using Eigen::Vector3d;

int DoMain() {

  const bool flip_normals = true;

  // Load mesh for a sphere.
  std::unique_ptr<Mesh<double>> sphere = LoadMeshFromObj(
      "drake/geometry/mesh_query/examples/sphere.obj", flip_normals);
  sphere->mesh_index = 0;

  // Load mesh for an ellipsoid.
  std::unique_ptr<Mesh<double>> ellipsoid = LoadMeshFromObj(
      "drake/geometry/mesh_query/examples/ellipsoid.obj", flip_normals);
  ellipsoid->mesh_index = 1;

  const double radius = 1.0;
  const double penetration = FLAGS_penetration;
  const double z_WSo = radius - penetration;
  const double sigma = FLAGS_sigma;

  // Place sphere a "penetration" distance below z = 0.
  // Apply an arbirary rotation for testing.
  Isometry3d X_WSphere{Translation3d{Vector3d(0, 0, z_WSo)}};
  X_WSphere.linear() =
      (AngleAxisd(M_PI / 3., Vector3d::UnitX()) *
       AngleAxisd(M_PI / 5., Vector3d::UnitY()) *
       AngleAxisd(3. * M_PI / 8., Vector3d::UnitZ())).matrix();

  // The top of the ellipsoid is placed at z = 0
  Isometry3d X_WEllipsoid{Translation3d{Vector3d(0, 0, -1.0)}};
  X_WEllipsoid.linear() = AngleAxisd(M_PI_2, Vector3d::UnitX()).matrix();

  // Write mesh and normals to a file.
  std::ofstream bottom_sphere_file("ellipsoid.vtk");
  OutputMeshToVTK(bottom_sphere_file, ellipsoid->points_G, ellipsoid->triangles, X_WEllipsoid);
  AppendCellCenteredVectorFieldToVTK(
      bottom_sphere_file, "FaceNormals", ellipsoid->face_normals_G, X_WEllipsoid);
  bottom_sphere_file << "POINT_DATA " << ellipsoid->points_G.size() << std::endl;
  AppendNodeCenteredVectorFieldToVTK(
      bottom_sphere_file, "NodeNormals", ellipsoid->node_normals_G, X_WEllipsoid);
  bottom_sphere_file.close();

  std::ofstream top_sphere_file("sphere.vtk");
  OutputMeshToVTK(top_sphere_file, sphere->points_G, sphere->triangles, X_WSphere);
  AppendCellCenteredVectorFieldToVTK(
      top_sphere_file, "FaceNormals", sphere->face_normals_G, X_WSphere);
  top_sphere_file << "POINT_DATA " << sphere->points_G.size() << std::endl;
  AppendNodeCenteredVectorFieldToVTK(
      top_sphere_file, "NodeNormals", sphere->node_normals_G, X_WSphere);
  top_sphere_file.close();

  // Perform the mesh-mesh query.
  // The triangles referenced in the pairs are "global indexes" in the original
  // full meshes.
  std::vector<PenetrationAsTrianglePair<double>> results = MeshToMeshQuery(
      X_WEllipsoid, *ellipsoid,
      X_WSphere, *sphere, sigma);

  // This call creates the two patches on each mesh and updates "results" so
  // that the triangle indexes in each pair are "local indexes" to the patch
  // meshes.
  auto patches = MakeLocalPatchMeshes(&results, *ellipsoid, *sphere);
  std::unique_ptr<Mesh<double>> ellipsoid_patch = std::move(patches.first);
  std::unique_ptr<Mesh<double>> sphere_patch = std::move(patches.second);

  PRINT_VAR(sphere_patch->points_G.size());
  PRINT_VAR(ellipsoid_patch->points_G.size());

  OutputMeshToVTK("ellipsoid_patch.vtk",
                  ellipsoid_patch->points_G, ellipsoid_patch->triangles,
                  X_WEllipsoid);

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
  std::vector<int> ellipsoid_patch_triangles;
  for (const auto& pair : results) {

    if (pair.meshA_index == sphere->mesh_index) {
      sphere_patch_triangles.push_back(pair.triangleA_index);
    } else {
      ellipsoid_patch_triangles.push_back(pair.triangleA_index);
    }

    if (pair.meshB_index == sphere->mesh_index) {
      sphere_patch_triangles.push_back(pair.triangleB_index);
    } else {
      ellipsoid_patch_triangles.push_back(pair.triangleB_index);
    }
  }

  for (auto sphere_triangle : sphere_patch_triangles) {
    PRINT_VAR(sphere_triangle);
  }

  for (auto ellipsoid_triangle : ellipsoid_patch_triangles) {
    PRINT_VAR(ellipsoid_triangle);
  }

  return 0;
}

}  // namespace
}  // namespace mesh_query
}  // namespace geometry
}  // namespace drake

int main(int argc, char* argv[]) {
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  return drake::geometry::mesh_query::DoMain();
}

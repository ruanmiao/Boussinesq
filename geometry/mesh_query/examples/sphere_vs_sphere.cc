#include <memory>
#include <vector>

#include <gflags/gflags.h>

#include "drake/common/find_resource.h"
#include "drake/geometry/mesh_query/vtk_io.h"
#include "drake/geometry/mesh_query/mesh_query.h"
#include "drake/multibody/shapes/geometry.h"

namespace drake {
namespace geometry {
namespace mesh_query {
namespace {

using Eigen::Translation3d;
using Eigen::Isometry3d;
using Eigen::Vector3d;

std::unique_ptr<Mesh<double>> LoadMeshFromObj(
    const std::string& file_name) {
  const auto resource_name = FindResourceOrThrow(file_name);
  DrakeShapes::Mesh mesh_loader(resource_name, resource_name);

  auto mesh = std::make_unique<Mesh<double>>();
  mesh_loader.LoadObjFile(&mesh->points_G, &mesh->triangles);

  // Compute normals.
  mesh->face_normals_G = CalcMeshFaceNormals(mesh->points_G, mesh->triangles);
  mesh->node_normals_G = CalcAreaWeightedNormals(*mesh);

  const int num_points = mesh->points_G.size();
  // Allocate and initialize to invalid index values.
  mesh->node_element.resize(num_points, std::make_pair(-1, -1));

  // Arbitrarily fill in mesh->node_element.
  const int num_elements = mesh->triangles.size();
  for (int element_index = 0; element_index < num_elements; ++element_index) {
    const auto& triangle = mesh->triangles[element_index];
    for (int i = 0; i < 3; ++i) {
      const int node_index = triangle[i];
      if (mesh->node_element[node_index].first < 0) {  // not yet initialized.
        mesh->node_element[node_index].first = element_index;
        mesh->node_element[node_index].second = i;
      }
    }
  }

  return mesh;
}

int DoMain() {

  std::unique_ptr<Mesh<double>> top_sphere = LoadMeshFromObj(
      "drake/geometry/mesh_query/examples/sphere.obj");
  top_sphere->mesh_index = 0;

  std::unique_ptr<Mesh<double>> bottom_sphere = LoadMeshFromObj(
      "drake/geometry/mesh_query/examples/sphere.obj");
  bottom_sphere->mesh_index = 1;

  // For this example normals in the obj point inward and therefore we flip
  // them.
  FlipNormals(top_sphere.get());
  FlipNormals(bottom_sphere.get());

  // Re-compute normals just in case for verification.
  top_sphere->face_normals_G = CalcMeshFaceNormals(
      top_sphere->points_G, top_sphere->triangles);

  const double radius = 1.0;
  const double penetration = 0.1;
  const double z_WSo = radius - penetration;

  Isometry3d X_WStop{Translation3d{Vector3d(0, 0, z_WSo)}};

  Isometry3d X_WSbottom{Translation3d{Vector3d(0, 0, -1.0)}};

  // Write mesh and normals to a file.
  std::ofstream bottom_sphere_file("top_sphere.vtk");
  OutputMeshToVTK(bottom_sphere_file, bottom_sphere->points_G, bottom_sphere->triangles, X_WStop);
  AppendCellCenteredVectorFieldToVTK(
      bottom_sphere_file, "FaceNormals", bottom_sphere->face_normals_G);
  AppendNodeCenteredVectorFieldToVTK(
      bottom_sphere_file, "NodeNormals", bottom_sphere->node_normals_G);
  bottom_sphere_file.close();

  std::ofstream top_sphere_file("bottom_sphere.vtk");
  OutputMeshToVTK(top_sphere_file, top_sphere->points_G, top_sphere->triangles, X_WSbottom);
  AppendCellCenteredVectorFieldToVTK(
      top_sphere_file, "FaceNormals", top_sphere->face_normals_G);
  AppendNodeCenteredVectorFieldToVTK(
      top_sphere_file, "NodeNormals", top_sphere->node_normals_G);
  top_sphere_file.close();

  // Perform the mesh-mesh query.
  std::vector<PenetrationAsTrianglePair<double>> results = MeshToMeshQuery(
      X_WSbottom, *bottom_sphere,
      X_WStop, *top_sphere);

  std::vector<Vector3d> pointsA(results.size());
  std::transform(results.begin(), results.end(), pointsA.begin(),
                 [](const PenetrationAsTrianglePair<double>& results) {
                   return results.p_WoAs_W;
                 });

  std::vector<Vector3d> pointsB(results.size());
  std::transform(results.begin(), results.end(), pointsB.begin(),
  [](const PenetrationAsTrianglePair<double>& results) {
    return results.p_WoBs_W;
  });

#if 0
  std::vector<Vector3d> normals(2 * results.size());
  std::transform(results.begin(), results.end(), normals.begin(),
                 [](const PenetrationAsTrianglePair<double>& results) {
                   return results.normal_A_W;
                 });
  std::transform(results.begin(), results.end(), normals.begin() + results.size(),
                 [](const PenetrationAsTrianglePair<double>& results) {
                   return results.normal_B_W;
                 });
#endif

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

#if 0
  {
    std::ofstream file("pointsA.vtk");
    OutputScatteredPointsToVTK(file, pointsA);
    file.close();
  }

  {
    std::ofstream file("pointsB.vtk");
    OutputScatteredPointsToVTK(file, pointsB);
    file.close();
  }
#endif

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

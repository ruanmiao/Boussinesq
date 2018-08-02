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

  mesh->face_normals_G = CalcMeshFaceNormals(mesh->points_G, mesh->triangles);

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

  std::unique_ptr<Mesh<double>> sphere = LoadMeshFromObj(
      "drake/geometry/mesh_query/examples/sphere.obj");
  sphere->mesh_index = 0;

  std::unique_ptr<Mesh<double>> plane = LoadMeshFromObj(
      "drake/geometry/mesh_query/examples/big_triangle.obj");
  plane->mesh_index = 1;

  // For this example normals in the obj point inward and therefore we flip
  // them.
  FlipNormals(sphere.get());

  // Re-compute normals just in case for verification.
  sphere->face_normals_G = CalcMeshFaceNormals(
      sphere->points_G, sphere->triangles);

  const double radius = 1.0;
  const double penetration = 0.1;
  const double z_WSo = radius - penetration;

  Isometry3d X_WS{Translation3d{Vector3d(0, 0, z_WSo)}};

  // Write mesh and normals to a file.
  std::ofstream spere_file("sphere.vtk");
  OutputMeshToVTK(spere_file, sphere->points_G, sphere->triangles, X_WS);
  AppendCellCenteredVectorFieldToVTK(
      spere_file, "normals", sphere->face_normals_G);
  spere_file.close();

  // Perform the mesh-mesh query.
  std::vector<PenetrationAsTrianglePair<double>> results = MeshToMeshQuery(
      Isometry3d::Identity(), *plane,
      X_WS, *sphere);

  std::vector<Vector3d> pointsB(results.size());
  std::transform(results.begin(), results.end(), pointsB.begin(),
  [](const PenetrationAsTrianglePair<double>& results) {
    return results.p_WoBs_W;
  });

  std::vector<Vector3d> pointsA(results.size());
  std::transform(results.begin(), results.end(), pointsA.begin(),
                 [](const PenetrationAsTrianglePair<double>& results) {
                   return results.p_WoAs_W;
                 });
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

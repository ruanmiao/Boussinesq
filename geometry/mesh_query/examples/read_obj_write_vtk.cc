#include <memory>
#include <vector>

#include <gflags/gflags.h>

#include "drake/common/find_resource.h"
#include "drake/geometry/mesh_query/vtk_io.h"
#include "drake/multibody/shapes/geometry.h"

namespace drake {
namespace geometry {
namespace mesh_query {
namespace {

using ::DrakeShapes::Mesh;

int DoMain() {

  std::vector<Vector3<double>> points;
  std::vector<Vector3<int>> triangles;

  const auto file_name =
      FindResourceOrThrow("drake/geometry/mesh_query/examples/sphere.obj");
  Mesh mesh_loader(file_name, file_name);
  mesh_loader.LoadObjFile(&points, &triangles);

  OutputMeshToVTK("sphere.vtk", points, triangles);
  OutputScatteredPointsToVTK("sphere_points.vtk", points);
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

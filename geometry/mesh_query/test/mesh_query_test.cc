#include "drake/geometry/mesh_query/mesh_query.h"

#include "drake/common/test_utilities/eigen_matrix_compare.h"

#include <gtest/gtest.h>

namespace drake {
namespace geometry {
namespace mesh_query {
namespace {

using Eigen::Vector3d;
using Eigen::Vector3i;

std::pair<std::vector<Vector3d>, std::vector<Vector3i>> MakeCube() {
  std::vector<Vector3d> points;
  std::vector<Vector3i> triangles;

  // Bottom
  points.push_back({0.0, 0.0, 0.0});
  points.push_back({1.0, 0.0, 0.0});
  points.push_back({1.0, 1.0, 0.0});
  points.push_back({0.0, 1.0, 0.0});

  // Top
  points.push_back({0.0, 0.0, 1.0});
  points.push_back({1.0, 0.0, 1.0});
  points.push_back({1.0, 1.0, 1.0});
  points.push_back({0.0, 1.0, 1.0});

  // Bottom
  triangles.push_back({0, 2, 1});
  triangles.push_back({0, 3, 2});

  // Top
  triangles.push_back({4, 5, 6});
  triangles.push_back({6, 7, 4});

  // Front
  triangles.push_back({0, 1, 5});
  triangles.push_back({5, 4, 0});

  // Back
  triangles.push_back({3, 7, 6});
  triangles.push_back({6, 2, 3});

  // Left
  triangles.push_back({0, 4, 7});
  triangles.push_back({7, 3, 0});

  // Right
  triangles.push_back({6, 5, 1});
  triangles.push_back({1, 2, 6});

  return std::make_pair(points, triangles);
}

GTEST_TEST(VisualMaterial, CalcMeshFaceNormals) {
  const double kTolerance = 10 * std::numeric_limits<double>::epsilon();

  // Make the mesh for a simple cube.
  std::vector<Vector3d> points_W;
  std::vector<Vector3i> triangles;
  std::tie(points_W, triangles) = MakeCube();
  ASSERT_EQ(points_W.size(), 8);
  ASSERT_EQ(triangles.size(), 12);

  std::vector<Vector3<double>> normals_W =
      CalcMeshFaceNormals(points_W, triangles);
  ASSERT_EQ(normals_W.size(), 12);

  // Bottom
  EXPECT_TRUE(CompareMatrices(
      normals_W[0], -Vector3d::UnitZ(),
      kTolerance, MatrixCompareType::relative));
  EXPECT_TRUE(CompareMatrices(
      normals_W[1], -Vector3d::UnitZ(),
      kTolerance, MatrixCompareType::relative));

  // Top
  EXPECT_TRUE(CompareMatrices(
      normals_W[2], Vector3d::UnitZ(),
      kTolerance, MatrixCompareType::relative));
  EXPECT_TRUE(CompareMatrices(
      normals_W[3], Vector3d::UnitZ(),
      kTolerance, MatrixCompareType::relative));

  // Front
  EXPECT_TRUE(CompareMatrices(
      normals_W[4], -Vector3d::UnitY(),
      kTolerance, MatrixCompareType::relative));
  EXPECT_TRUE(CompareMatrices(
      normals_W[5], -Vector3d::UnitY(),
      kTolerance, MatrixCompareType::relative));

  // Back
  EXPECT_TRUE(CompareMatrices(
      normals_W[6], Vector3d::UnitY(),
      kTolerance, MatrixCompareType::relative));
  EXPECT_TRUE(CompareMatrices(
      normals_W[7], Vector3d::UnitY(),
      kTolerance, MatrixCompareType::relative));

  // Left
  EXPECT_TRUE(CompareMatrices(
      normals_W[8], -Vector3d::UnitX(),
      kTolerance, MatrixCompareType::relative));
  EXPECT_TRUE(CompareMatrices(
      normals_W[9], -Vector3d::UnitX(),
      kTolerance, MatrixCompareType::relative));

  // Right
  EXPECT_TRUE(CompareMatrices(
      normals_W[10], Vector3d::UnitX(),
      kTolerance, MatrixCompareType::relative));
  EXPECT_TRUE(CompareMatrices(
      normals_W[11], Vector3d::UnitX(),
      kTolerance, MatrixCompareType::relative));
}

}  // namespace
}  // namespace mesh_query
}  // namespace geometry
}  // namespace drake

#include "drake/geometry/mesh_query/mesh_query.h"

#include "drake/common/test_utilities/eigen_matrix_compare.h"

#include <gtest/gtest.h>

namespace drake {
namespace geometry {
namespace mesh_query {
namespace {

using Eigen::Isometry3d;
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

GTEST_TEST(MeshQueries, CalcMeshFaceNormals) {
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

struct PointToMeshQueryData {
  // Mesh nodes expressed in the mesh geometry frame G.
  std::vector<Vector3d> points_G;

  // Mesh triangles.
  std::vector<Vector3i> triangles;

  // Mesh normals expressed in the G frame.
  std::vector<Vector3d> normals_G;

  // Pose of the mesh geometry G in the world frame W.
  Isometry3d X_WG{Isometry3d::Identity()};

  // Query point in the world frame W.
  Vector3d p_WQ;
};

struct PointToMeshQueryResults {
  // `true` if the query point Q is inside the mesh.
  bool is_inside;

  // Signed distance from point Q to the mesh.
  double distance;

  // Negative implies not checked.
  int triangle_index{-1};

  Vector3d normal_W;

  Vector3d p_WP;
};

void VerifyPointToMeshQuery(
    const PointToMeshQueryData& data, PointToMeshQueryResults expected_results) {
  PointMeshDistance<double> results;
  const bool is_inside = CalcPointToMeshNegativeDistance(
      data.X_WG, data.points_G, data.triangles, data.normals_G, data.p_WQ,
      &results);
  EXPECT_EQ(is_inside, expected_results.is_inside);

  const double kTolerance = 10 * std::numeric_limits<double>::epsilon();

  EXPECT_NEAR(results.distance, expected_results.distance, kTolerance);

  // Do not verify the returned data (garage) if the point is outside.
  if (!expected_results.is_inside) return;

  if (expected_results.triangle_index >= 0) {
    EXPECT_EQ(results.triangle_index, expected_results.triangle_index);
  }

  EXPECT_TRUE(CompareMatrices(
      results.normal_F, expected_results.normal_W,
      kTolerance, MatrixCompareType::relative));

  EXPECT_TRUE(CompareMatrices(
      results.p_FP, expected_results.p_WP,
      kTolerance, MatrixCompareType::relative));
}

GTEST_TEST(MeshQueries, PointInsideCubeMesh) {
  const double kEpsilon = std::numeric_limits<double>::epsilon();

  // Make a mesh with points measured and expressed in a "geometry" frame G.
  std::vector<Vector3d> points_G;
  std::vector<Vector3i> triangles;
  std::tie(points_G, triangles) = MakeCube();
  ASSERT_EQ(points_G.size(), 8);
  ASSERT_EQ(triangles.size(), 12);
  std::vector<Vector3<double>> normals_G =
      CalcMeshFaceNormals(points_G, triangles);

  PointToMeshQueryData data;
  data.X_WG = Isometry3d::Identity();
  data.points_G = points_G;
  data.triangles = triangles;
  data.normals_G = normals_G;

  PointToMeshQueryResults expected_results;

  {
    data.p_WQ << 0.35, 0.4, 0.5;
    expected_results.is_inside = true;
    expected_results.distance = -0.35;
    expected_results.triangle_index = 8;
    expected_results.normal_W = -Vector3d::UnitX();
    expected_results.p_WP << 0.0, 0.4, 0.5;
    VerifyPointToMeshQuery(data, expected_results);
  }

  {
    data.p_WQ << 0.1, 0.45, 0.95;
    expected_results.is_inside = true;
    expected_results.distance = -0.05;
    expected_results.triangle_index = 3;
    expected_results.normal_W = Vector3d::UnitZ();
    expected_results.p_WP << 0.1, 0.45, 1.0;
    VerifyPointToMeshQuery(data, expected_results);
  }

  expected_results.triangle_index = -1;
  {
    data.p_WQ << kEpsilon, 0.5, 0.5;
    expected_results.is_inside = true;
    expected_results.distance = 0.0;
    expected_results.normal_W = -Vector3d::UnitX();
    expected_results.p_WP << 0.0, 0.5, 0.5;
    VerifyPointToMeshQuery(data, expected_results);
  }

  {
    data.p_WQ << -kEpsilon, 0.5, 0.5;
    expected_results.is_inside = false;
    VerifyPointToMeshQuery(data, expected_results);
  }

  {
    data.p_WQ << 0.1, 1.0-kEpsilon, 0.99;
    expected_results.is_inside = true;
    expected_results.distance = 0.0;
    expected_results.normal_W = Vector3d::UnitY();
    expected_results.p_WP << 0.1, 1.0, 0.99;
    VerifyPointToMeshQuery(data, expected_results);
  }

  {
    data.p_WQ << 0.1, 1.0+kEpsilon, 0.99;
    expected_results.is_inside = false;
    VerifyPointToMeshQuery(data, expected_results);
  }

  {
    data.p_WQ << 1.2, 3.0, 0.5;
    expected_results.is_inside = false;
    VerifyPointToMeshQuery(data, expected_results);
  }

  {
    data.p_WQ << -1.2, 0.5, -0.2;
    expected_results.is_inside = false;
    VerifyPointToMeshQuery(data, expected_results);
  }
}

}  // namespace
}  // namespace mesh_query
}  // namespace geometry
}  // namespace drake

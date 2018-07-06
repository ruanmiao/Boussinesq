#include "drake/multibody/boussinesq_solver/triangle_area.h"

#include <gtest/gtest.h>

#include "drake/common/eigen_types.h"

namespace drake {
namespace multibody {
namespace boussinesq_solver {
namespace {

using Eigen::Vector2d;

GTEST_TEST(TriangleAreaTest, AreaValue) {
  const Vector2d p1(0.0, 0.0);
  const Vector2d p2(1.0, 0.0);
  const Vector2d p3(0.0, 1.0);
  const double area = TriangleArea(p1, p2, p3);
  const double expected_area = 0.5;
  EXPECT_NEAR(area, expected_area, 10 * std::numeric_limits<double>::epsilon());
}

GTEST_TEST(TriangleAreaTest, AreaSign) {
  const Vector2d p1(0.0, 0.0);
  const Vector2d p2(1.0, 0.3);
  const Vector2d p3(-0.2, 1.2);
  const double area = TriangleArea(p1, p2, p3);
  EXPECT_GT(area, 0.0);
}

}  // namespace
}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake

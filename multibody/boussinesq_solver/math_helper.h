#pragma once

#include "drake/common/eigen_types.h"

namespace drake {
namespace multibody {
namespace boussinesq_solver {

/// This method computes the area of the triangle defined by points p1, p2 and
/// p3. The return value is positive if the points are oriented counterclockwise
/// and negative otherwise.
/// @param p1 First vertex.
/// @param p2 Second vertex.
/// @param p3 Third vertex.
/// @returns The are of the triangle defined by p1, p2 and p3.
double CalcTriangleArea(const Vector2<double>& p1,
                        const Vector2<double>& p2,
                        const Vector2<double>& p3);

/// This method checks the orientation of the three points p1, p2, p3.
/// It returns 1 if the three points are ordered counter-clockwise. -1 if
/// clockwise. 0 if the three are co-linear.
/// @param p1 First vertex.
/// @param p2 Second vertex.
/// @param p3 Third vertex.
/// @returns The orientation of the triangle defined by p1, p2 and p3.
int CalcTriangleOrientation(const Vector2<double>& p1,
                            const Vector2<double>& p2,
                            const Vector2<double>& p3);

/// This method returns the integral of Jₘₙ =∫ sinᵐ x cosⁿ x dx
/// in the interval x ∈ [θ₀, θf], where m = 0, n = -1;
/// @param theta_0.
/// @param theta_f.
/// @returns the integral Jmn, where m = 0, n = -1.
double CalcIntegralJ0minus1(double theta_0, double theta_f);

/// This method returns the integral of Jₘₙ =∫ sinᵐ x cosⁿ x dx
/// in the interval x ∈ [θ₀, θf], where m = 1, n = -2;
/// @param theta_0.
/// @param theta_f.
/// @returns the integral Jmn. where m = 1, n = -2.
double CalcIntegralJ1minus2(double theta_0, double theta_f);

/// This method returns the integral of Jₘₙ =∫ sinᵐ x cosⁿ x dx
/// in the interval x ∈ [θ₀, θf], where m = 0, n = 0;
/// @param theta_0.
/// @param theta_f.
/// @returns the integral Jmn. where m = 0, n = 0.
double CalcIntegralJ00(double theta_0, double theta_f);

/// This method returns the integral of Iₘₙ⁽ᵖ⁾ = ∫ sinᵐx cosⁿx Δᵖ dx
/// in the interval x ∈ [θ₀, θf], where p = 1, m = 0, n = -1;
/// Δ² = 1 - α²sin²θ
/// @param theta_0.
/// @param theta_f.
/// @returns the integral Imn^(p). where m = 0, n = -1, p = 1.
double CalcIntegralI0minus1P1(double theta_0, double theta_f, double alpha);

/// This method returns the integral of Iₘₙ⁽ᵖ⁾ = ∫ sinᵐx cosⁿx Δᵖ dx
/// in the interval x ∈ [θ₀, θf], where p = 1, m = 1, n = -2;
/// Δ² = 1 - α²sin²θ
/// @param theta_0.
/// @param theta_f.
/// @returns the integral Imn^(p). where m = 1, n = -2, p = 1.
double CalcIntegralI1minus2P1(double theta_0, double theta_f, double alpha);

/// This method returns the integral of Iₘₙ⁽ᵖ⁾ = ∫ sinᵐx cosⁿx Δᵖ dx
/// in the interval x ∈ [θ₀, θf], where p = -1, m = 2, n = -1;
/// Δ² = 1 - α²sin²θ
/// @param theta_0.
/// @param theta_f.
/// @returns the integral Imn^(p). where m = 2, n = -1, p = -1.
double CalcIntegralI2minus1Pminus1(
    double theta_0, double theta_f, double alpha);

/// This method returns the integral of Iₘₙ⁽ᵖ⁾ = ∫ sinᵐx cosⁿx Δᵖ dx
/// in the interval x ∈ [θ₀, θf], where p = -1, m = 1, n = 0;
/// Δ² = 1 - α²sin²θ
/// @param theta_0.
/// @param theta_f.
/// @returns the integral Imn^(p). where m = 1, n = 0, p = -1.
double CalcIntegralI10Pminus1(double theta_0, double theta_f, double alpha);

/// This method returns the Isometric transform from the Triangle frame to the
/// World frame. The origin of the triangle from is the projection of the field
/// point on the the plane of the triangle. The z-axis is defined by the
/// surface normal of the triangle (p1, p2, p3). THe x-axis is defined by
/// the vector from p1 to p2.
/// @param p1 The 3x1 Vector denotes the position of point 1 in the world frame
/// @param p2 The 3x1 Vector denotes the position of point 2 in the world frame
/// @param p3 The 3x1 Vector denotes the position of point 3 in the world frame
/// @param xA The 3x1 Vector denotes the position of field point
/// in the world frame
/// @returns The Isometric transform from the Triangle frame to the Wold frame.
Eigen::Isometry3d CalcTransformationFromTriangleFrame(
    const Eigen::Vector3d& p1, const Eigen::Vector3d& p2,
    const Eigen::Vector3d& p3, const Eigen::Vector3d& xA);

}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake

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
/// @returns the integral Jmn.
double CalcIntegralJm0nN1(double theta_0, double theta_f);

/// This method returns the integral of Jₘₙ =∫ sinᵐ x cosⁿ x dx
/// in the interval x ∈ [θ₀, θf], where m = 1, n = -2;
/// @param theta_0.
/// @param theta_f.
/// @returns the integral Jmn.
double CalcIntegralJm1nN2(double theta_0, double theta_f);

/// This method returns the inverse of a 2x2 matrix
/// @param matrix to be inversed
/// @returns the integral Jmn.
Matrix2<double> InverseMatrix2(Matrix2<double> matrix);


}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake

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
double TriangleArea(const Vector2<double>& p1,
                    const Vector2<double>& p2,
                    const Vector2<double>& p3);

/// This method check the orientation of the three points p1, p2, p3. Return 1
/// if the three points are ordered counter-clockwise. -1 if clockwise. 0 if
/// the three are co-linear.
/// @param p1 First vertex.
/// @param p2 Second vertex.
/// @param p3 Third vertex.
/// @returns The orientation of the triangle defined by p1, p2 and p3.
int TriangleOrientation(const Vector2<double>& p1,
                        const Vector2<double>& p2,
                        const Vector2<double>& p3 = Eigen::MatrixXd::Zero(
                            2, 1));

/// This method returns the integral of Jmn = int(sin^m * cos^n) in the interval
/// [theta_0, theta_f], where m = 0, n = -1;
/// @param theta_0 Start of the interval.
/// @param theta_f End of the interval.
/// @returns the integral Jmn.
double IntegralJm0nN1(const double& theta_0, const double& theta_f);

/// This method returns the integral of Jmn = int(sin^m * cos^n) in the interval
/// [theta_0, theta_f], where m = 1, n = -2;
/// @param theta_0 Start of the interval.
/// @param theta_f End of the interval.
/// @returns the integral Jmn.
double IntegralJm1nN2(const double& theta_0, const double& theta_f);

/// This method returns the inverse of a 2x2 matrix
/// @param matrix to be inversed
/// @returns the integral Jmn.
Matrix2<double> InverseMatrix2(Matrix2<double> matrix);


}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake

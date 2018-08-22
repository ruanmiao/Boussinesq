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


/// This method computes the vector of the area of the triangle defined by
/// points p1, p2 and p3. The norm of the vector is the area of the triangles,
/// and the direction points towards toward the surface normal of triangle
/// p1, p2, p3.
/// @param p1 First vertex.
/// @param p2 Second vertex.
/// @param p3 Third vertex.
/// @returns The are of the triangle defined by p1, p2 and p3.
Vector3<double> CalcTriangleArea(const Vector3<double>& p1,
                                 const Vector3<double>& p2,
                                 const Vector3<double>& p3);

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

/// This method return the volume of a pyramid with a trapesium base, given the
/// top and bottom length of the trapesium with its height, and the height of
/// the pyramid.
/// @param trapesium_length_1 The length of one of the parallel sides
/// @param trapesium_length_2 The length of the other parallel side
/// @param trapesium_height The height of the trapesium
/// @param pyramid_height The height of the pyramid
/// @returns The volume of the pyramid.
double CalcVolumeOfPyramidWithTrapesiumBase(
    double trapesium_length_1, double trapesium_length_2,
    double trapesium_height, double pyramid_height);

/// This method return the area of a tetrahedral given the 3 nodes at its base,
/// and its height.
/// @param p1 First node of the triangle base
/// @param p2 Second node of the triangle base
/// @param p3 Third node of the triangle base
/// @param height The height of the tetrahedral
/// @returns The volume of the tetrahedral.
double CalcVolumeOfTetrahedral(
    const Eigen::Vector3d &p1, const Eigen::Vector3d &p2,
    const Eigen::Vector3d &p3, double height);

/// This method will return the volume of the 3D geometry defined in the
/// following way: Given a base triangle p1, p2, p3, with the value on p3 has
/// its sign different from those of p1 and p2. Do an interpolation of the
/// value over the whole triangle, with z13 and z23 being the nodes on the
/// edge connecting p1, p3 and p2, p3. The 3D geometry is the one with the base
/// being the polygon p1, p2, z23, z13. The heights are the corresponding
/// values interpolated over the region.
/// @param p1_in The values of p1 and p2 are of the same sign
/// @param value_1 The value at p1.
/// @param p2_in The values of p1 and p2 are of the same sign
/// @param value_2 The value at p2.
/// @param p3_out The value of p3 is different from those of p1 and p2
/// @param value_3 The value at p3.
/// @returns The volume defined by the portion including p1 and p2.
/// The sign of the volume is the same as those of p1 and p2
double CalcInterpolatedVolumeExcludingOneNode(
    const Eigen::Vector3d& p1_in, double value_1,
    const Eigen::Vector3d& p2_in, double value_2,
    const Eigen::Vector3d& p3_out, double value_3);

/// This method will return the volume of the 3D geometry defined in the
/// following way: Given a base triangle p1, p2, p3, with the value on p3 has
/// its sign different from those of p1 and p2. Do an interpolation of the
/// value over the whole triangle, with z13 and z23 being the nodes on the
/// edge connecting p1, p3 and p2, p3. The 3D geometry is the one with the base
/// being the polygon p3, z23, z13. The heights are the corresponding
/// values interpolated over the region.
/// @param p1_out The values of p1 and p2 are of the same sign
/// @param value_1 The value at p1.
/// @param p2_out The values of p1 and p2 are of the same sign
/// @param value_2 The value at p2.
/// @param p3_in The value of p3 is different from those of p1 and p2
/// @param value_3 The value at p3.
/// @returns The volume defined by the portion including p3.
/// The sign of the volume is the same as that of p3
double CalcInterpolatedVolumeExcludingTwoNode(
    const Eigen::Vector3d& p1_out, double value_1,
    const Eigen::Vector3d& p2_out, double value_2,
    const Eigen::Vector3d& p3_in, double value_3);



}  // namespace boussinesq_solver
}  // namespace multibody
}  // namespace drake

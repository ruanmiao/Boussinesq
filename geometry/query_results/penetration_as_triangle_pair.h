#pragma once

#include "drake/common/drake_copyable.h"
#include "drake/common/eigen_types.h"
#include "drake/geometry/geometry_ids.h"

namespace drake {
namespace geometry {

/** @name     Triangle mesh penetration reported as set of point-pairs

 A unique characterization of the intersection of two penetrating meshes.
 Conceptually, the overlapping volume is bounded by two patches:
 one patch coming from each mesh which meet up at a common boundary (the
 "waterline"). Rather than reporting the volume directly, this structure
 represents the volume through point-pair samples. A point pair consists of
 one point on each patch and a collection of related data.

 For two meshes A and B, we'll define the set of points
 `{(a₁, b₁), (a₂, b₂), ..., (aₙ, bₙ)}`. For each pair `(a, b)`, where `a`
 lies on the surface of A and b lies on the surface of B, we define the
 following quantities:

   - φ, the *signed* distance from a to b: `φ = S(a, b)|a - b|`. Where
     `S(a, b) = 1 → a ⋂ B = ∅ ∧ b ⋂ A = ∅` (non-penetrating) and `
     S(a, b) = -1 → a ⋂ B = {a} ∧ b ⋂ A = {b}` (penetrating or touching).
   - A triangle from A on which point a lies, called Tᵃ.
   - The barycentric coordinates for a with respect to triangle Tᵃ: αᵃ, βᵃ,
     γᴬ. If a is a vertex of mesh A, then one of αᵃ, βᵃ, γᵃ will be one
     and the other values will be zero. In any case, they should sum to one.
   - A triangle on which point b lies, called Tᵇ.
   - The barycentric coordinates for b with respect to triangle Tᵇ: αᵇ, βᵇ, γᵇ.
   - The outward-pointing unit normal n̂ᵃ of mesh A at point a. Generally,
     `(b - a) · nᵃ ≠ φ`.
   - The outward-pointing unit normal n̂ᵇ of mesh B at point b.

 Note: this representation does not *explicitly* encode the connectivity of the
 patch. However, the connectivity can be inferred by exploring the triangles.

 <!--
    General assumptions:
      1. This struct can *only* be populated as a result of colliding two
         meshes. Attempting to do so with *any* other kind of geometric
         primitive should throw an exception.
      2. The mesh provides per-vertex normals. These normals are interpolated
         over the triangle to provide continuous normal definitions over the
         entire mesh.
 -->
 */

//@{

/** The data for a single point pair that samples the overlapping region. This
 struct doesn't specifically track which geometries these quantities are
 derived from (i.e., what is mesh A or mesh B). It relies on the greater context
 to provide that information.

 This defines the point pair (Aₛ, Bₛ) and its corresponding data.

 @tparam T The underlying scalar type. Must be a valid Eigen scalar.  */
template <typename T>
struct PenetrationAsTrianglePair {
  DRAKE_DEFAULT_COPY_AND_MOVE_AND_ASSIGN(PenetrationAsTrianglePair)

  const int meshA_index;

  /** A reference to the triangle on A on which As lies. Its lifespan is that
  of the underlying mesh. */
  const Vector3<int> triangle_A;

  /** Position vector for a point Aₛ on the surface of A measured and expressed
   in the world frame. */
  const Vector3<T> p_WoAs_W;

  /** The outward facing normal of mesh A at Aₛ expressed in the world frame. */
  const Vector3<T> normal_A_W;

  /** The barycentric coordinates of Aₛ on triangle_A.  */
  const Vector3<T> barycentric_A;

  const int meshB_index;

  /** A reference to the triangle of mesh B on which Bs lies.  */
  const Vector3<int> triangle_B;

  /** Position vector for a point Bₛ on the surface of B measured and expressed
  in the world frame. */
  const Vector3<T> p_WoBs_W;

  /** The barycentric coordinates of Bₛ on triangle_B.  */
  const Vector3<T> barycentric_B;

  /** The outward facing normal of mesh B at Bₛ expressed in the world frame. */
  const Vector3<T> normal_B_W;

  /** The signed distance from As to Bs: `φ = S(Aₛ, Bₛ)|Aₛ - Bₛ|`.  */
  const T signed_distance;

};

// @}

}  // namespace geometry
}  // namespace drake

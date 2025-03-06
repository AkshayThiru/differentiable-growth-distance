// Copyright 2025 Akshay Thirugnanam
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

/**
 * @file cylinder.h
 * @author Akshay Thirugnanam (akshay_t@berkeley.edu)
 * @date 2025-03-06
 * @brief 3D cylinder class.
 */

#ifndef DGD_GEOMETRY_3D_CYLINDER_H_
#define DGD_GEOMETRY_3D_CYLINDER_H_

#include <cassert>
#include <cmath>

#include "dgd/data_types.h"
#include "dgd/geometry/convex_set.h"

namespace dgd {

/**
 * @brief Axis-aligned 3D cylinder class.
 *
 * @note The cylinder is oriented along the x-axis.
 */
class Cylinder : public ConvexSet<3> {
 public:
  /**
   * @brief Constructs a Cylinder object.
   *
   * @param hlx    Half axis length.
   * @param radius Radius.
   * @param margin Safety margin.
   */
  Cylinder(Real hlx, Real radius, Real margin);

  ~Cylinder() {};

  Real SupportFunction(const Vec3f& n, Vec3f& sp) const final;

 private:
  const Real hlx_;    /**< Half axis length. */
  const Real radius_; /**< Radius. */
  const Real margin_; /**< Safety margin. */
};

inline Cylinder::Cylinder(Real hlx, Real radius, Real margin)
    : ConvexSet<3>(), hlx_(hlx), radius_(radius), margin_(margin) {
  assert((hlx > Real(0.0)) && (radius > Real(0.0)));
  assert(margin >= Real(0.0));
  SetInradius(std::min(hlx, radius) + margin);
}

inline Real Cylinder::SupportFunction(const Vec3f& n, Vec3f& sp) const {
  sp = margin_ * n;
  const Real k{std::sqrt(n(1) * n(1) + n(2) * n(2))};
  if (k > kEps) sp.tail<2>() += radius_ * n.tail<2>() / k;
  sp(0) += std::copysign(hlx_, n(0));
  return sp.dot(n);
}

}  // namespace dgd

#endif  // DGD_GEOMETRY_3D_CYLINDER_H_

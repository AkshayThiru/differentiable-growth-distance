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
 * @file convex_set.h
 * @author Akshay Thirugnanam (akshay_t@berkeley.edu)
 * @date 2025-02-18
 * @brief Convex set abstract class implementation.
 */

#ifndef DGD_GEOMETRY_CONVEX_SET_H_
#define DGD_GEOMETRY_CONVEX_SET_H_

#include <cassert>

#include "dgd/data_types.h"

namespace dgd {

/**
 * @brief Convex set abstract class implementing the support function.
 *
 * @attention The convex set must be a compact and solid set (i.e., closed,
 * bounded, and with a nonempty interior). Also, the origin must be in
 * the interior of the set.
 *
 * @tparam dim Dimension of the convex set (2 or 3).
 */
template <int dim>
class ConvexSet {
 protected:
  /**
   * @brief Initializes convex set properties.
   *
   * @see ConvexSet(Real inradius)
   */
  ConvexSet();

  /**
   * @brief Initializes convex set inradius.
   *
   * @param inradius Inradius for the convex set.
   * @see inradius_
   */
  ConvexSet(Real inradius);

  /**
   * @brief Inradius for the convex set.
   *
   * Radius of a ball that is centered at the origin and contained in the set.
   * Any number greater than 0 and less than the Chebyshev radius will work.
   * However, larger values can help prevent singularities in simplex
   * computations.
   */
  Real inradius_;

 public:
  /**
   * @brief Destroys the Convex Set object.
   */
  virtual ~ConvexSet() {}

  /**
   * @brief Implements the support function.
   *
   * Implements the support function for the convex set \f$C\f$, which is given
   * by: \f{align*}{
   * sv & = \max_{x \in C} \langle n, x\rangle, \\
   * sp & \in \arg \max_{x \in C} \langle n, x\rangle,
   * \f}
   * where \f$sv\f$ is the return value of the function.
   *
   * @note Safety margins (in the 2-norm) can be directly included in the
   * support function computation as: \f{align*}{
   * sv & \leftarrow sv + m, \\
   * sp & \leftarrow sp + m \cdot n,
   * \f}
   * where \f$m\f$ is the safety margin. Note that the safety margin also
   * increases the inradius.
   *
   * @param[in]  n  Normal vector with unit 2-norm.
   * @param[out] sp Support point. A point at which the maximum for the
   * support function is attained.
   * @return     Value of the support function at the normal vector.
   */
  virtual Real SupportFunction(const Vecf<dim>& n, Vecf<dim>& sp) const = 0;

  /**
   * @brief Gets the dimension of the convex set.
   *
   * @return Dimension, given by the template parameter dim.
   */
  static constexpr int Dimension();

  /**
   * @brief Gets the inradius.
   *
   * @return Inradius.
   * @see inradius_
   */
  Real GetInradius() const;

  /**
   * @brief Sets the inradius.
   *
   * @param[in] inradius Inradius (\f$> 0\f$).
   * @see inradius_
   */
  void SetInradius(Real inradius);
};

template <int dim>
inline ConvexSet<dim>::ConvexSet() : ConvexSet(kEps) {}

template <int dim>
inline ConvexSet<dim>::ConvexSet(Real inradius) : inradius_(inradius) {
  assert(inradius > Real(0.0));
}

template <int dim>
constexpr int ConvexSet<dim>::Dimension() {
  return dim;
}

template <int dim>
inline Real ConvexSet<dim>::GetInradius() const {
  return inradius_;
}

template <int dim>
inline void ConvexSet<dim>::SetInradius(Real inradius) {
  assert(inradius > Real(0.0));
  inradius_ = inradius;
}

}  // namespace dgd

#endif  // DGD_GEOMETRY_CONVEX_SET_H_

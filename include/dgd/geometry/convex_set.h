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

#include <stdexcept>

#include "dgd/data_types.h"

namespace dgd {

/**
 * @brief Support function hint struct; used internally.
 *
 * @tparam dim Dimension of the convex sets.
 */
template <int dim>
struct SupportFunctionHint {
  Vecr<dim> n_prev = Vecr<dim>::Zero();
  int idx_ws = -1;
};

/**
 * @brief Convex set abstract class implementing the support function.
 *
 * @attention The convex set must be a compact and solid set (i.e., closed,
 * bounded, and with a nonempty interior). Also, the origin must be in
 * the interior of the set. The origin is the center point of the convex set
 * in its local frame.
 *
 * @tparam dim Dimension of the convex set (2 or 3).
 */
template <int dim>
class ConvexSet {
 public:
  virtual ~ConvexSet() {}

  /**
   * @brief Implements the support function in the local frame.
   *
   * Implements the support function for the convex set \f$C\f$, which is given
   * by: \f{align*}{
   * sv & = \max_{x \in C} \langle n, x\rangle, \\
   * sp & \in \arg \max_{x \in C} \langle n, x\rangle,
   * \f}
   * where \f$sv\f$ is the return value of the function.
   *
   * @note If the normal vector n is required to have unit 2-norm, the function
   * RequireUnitNormal should return true.
   *
   * @note Safety margins (in the 2-norm) can be directly included in the
   * support function computation (when n has unit 2-norm) as: \f{align*}{
   * sv & \leftarrow sv + m, \\
   * sp & \leftarrow sp + m \cdot n,
   * \f}
   * where \f$m\f$ is the safety margin. Note that the safety margin also
   * increases the inradius.
   *
   * @attention A nonzero safety margin can result in slower convergence.
   *
   * @see RequireUnitNormal
   *
   * @param[in]     n    Normal vector.
   * @param[out]    sp   Support point. Any point at which the maximum for the
   *                     support function is attained.
   * @param[in,out] hint Additional hints.
   * @return        Support function value at the normal vector.
   */
  virtual Real SupportFunction(
      const Vecr<dim>& n, Vecr<dim>& sp,
      SupportFunctionHint<dim>* /*hint*/ = nullptr) const = 0;

  /**
   * @brief Returns the normalization requirement for the normal vector passed
   * to SupportFunction.
   *
   * If the return value is true, the argument n to SupportFunction will have
   * unit 2-norm.
   *
   * @see SupportFunction
   */
  virtual bool RequireUnitNormal() const = 0;

  /**
   * @brief Returns true if the convex set is polytopic.
   */
  virtual bool IsPolytopic() const = 0;

  static constexpr int dimension();

  Real inradius() const;

  void set_inradius(Real inradius);

 protected:
  ConvexSet();

  ConvexSet(Real inradius);

  /**
   * @brief Convex set inradius at the origin.
   *
   * Radius of a ball that is centered at the origin and contained in the set.
   * Any number greater than 0 and less than the inradius will work. However,
   * larger values can help prevent singularities in simplex computations and
   * enable faster convergence.
   */
  Real inradius_;
};

template <int dim>
inline ConvexSet<dim>::ConvexSet() : ConvexSet(kEps) {}

template <int dim>
inline ConvexSet<dim>::ConvexSet(Real inradius) : inradius_(inradius) {
  if (inradius <= 0.0) {
    throw std::domain_error("Inradius is not positive");
  }
}

template <int dim>
constexpr int ConvexSet<dim>::dimension() {
  return dim;
}

template <int dim>
inline Real ConvexSet<dim>::inradius() const {
  return inradius_;
}

template <int dim>
inline void ConvexSet<dim>::set_inradius(Real inradius) {
  if (inradius <= 0.0) {
    throw std::domain_error("Inradius is not positive");
  }
  inradius_ = inradius;
}

}  // namespace dgd

#endif  // DGD_GEOMETRY_CONVEX_SET_H_

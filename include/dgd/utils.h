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
 * @file utils.h
 * @author Akshay Thirugnanam (akshay_t@berkeley.edu)
 * @date 2025-02-28
 * @brief Utility functions.
 */

#ifndef DGD_UTILS_H_
#define DGD_UTILS_H_

#include <cmath>
#include <random>
#include <stdexcept>

#include "dgd/data_types.h"

namespace dgd {

/**
 * Math utility functions.
 */

/**
 * @brief Sets the rotation matrix using the body ZYX Euler angles.
 *
 * @param[in]  euler ZYX Euler angles, in the form (roll, pitch, yaw).
 * @param[out] rot   Rotation matrix.
 */
inline void EulerToRotation(const Vec3r& euler, Rotation3r& rot) {
  const Eigen::AngleAxis<Real> R(euler(0), Vec3r::UnitX());
  const Eigen::AngleAxis<Real> P(euler(1), Vec3r::UnitY());
  const Eigen::AngleAxis<Real> Y(euler(2), Vec3r::UnitZ());
  rot.noalias() = (Y * P * R).matrix();
}

/**
 * RNG utility functions.
 */

namespace {  // Anonymous namespace to prevent ODR violations.

std::random_device rd;
std::mt19937 generator{rd()};

inline void SetDefaultSeed() { generator.seed(5489u); }

inline void SetRandomSeed() { generator.seed(rd()); }

/**
 * @brief Returns a uniform random real number.
 */
inline Real Random(Real range_from, Real range_to) {
  if (range_from >= range_to) {
    throw std::range_error("Invalid range");
  }
  std::uniform_real_distribution<Real> dis(range_from, range_to);
  return dis(generator);
}

/**
 * @brief Returns a uniform random real number with zero mean.
 *
 * @param  range Maximum absolute value of random number.
 */
inline Real Random(Real range = 1.0) { return Random(-range, range); }

/**
 * @brief Sets a random (nonuniformly distributed) unit vector.
 */
template <int dim>
inline void RandomUnitVector(Vecr<dim>& n) {
  static_assert((dim == 2) || (dim == 3), "Invalid dimension");
  if constexpr (dim == 2) {
    n << Random(), Random();
  } else {
    n << Random(), Random(), Random();
  }
  const Real norm = n.norm();
  if (norm < kEps) {
    n = Vecr<dim>::UnitX();
  } else {
    n /= norm;
  }
}

/**
 * @name Random rotation functions
 * @brief Sets a random rotation matrix.
 */
///@{
inline void RandomRotation(Rotation2r& rot) {
  const Real angle = Random(kPi);
  rot(0, 0) = std::cos(angle);
  rot(1, 0) = std::sin(angle);
  rot(0, 1) = -rot(1, 0);
  rot(1, 1) = rot(0, 0);
}

inline void RandomRotation(Rotation3r& rot) {
  EulerToRotation(kPi * Vec3r(Random(), Random(0.5), Random()), rot);
}
///@}

/**
 * @brief Sets a random rigid body transformation matrix.
 *
 * @param[in]  range_from Lower bound of position.
 * @param[in]  range_to   Upper bound of position.
 * @param[out] tf         Transformation matrix.
 */
template <int dim>
inline void RandomRigidBodyTransform(const Vecr<dim>& range_from,
                                     const Vecr<dim>& range_to,
                                     Transformr<dim>& tf) {
  if ((range_from.array() >= range_to.array()).any()) {
    throw std::range_error("Invalid range");
  }
  Rotationr<dim> rot;
  RandomRotation(rot);
  tf.template block<dim, dim>(0, 0) = rot;
  for (int i = 0; i < dim; ++i) tf(i, dim) = Random(range_from(i), range_to(i));
  tf.template block<1, dim>(dim, 0) = Vecr<dim>::Zero().transpose();
  tf(dim, dim) = 1.0;
}

/**
 * @brief Sets a random rigid body transformation matrix with position in a cube
 * range.
 *
 * @param[in]  range_from Lower bound of position.
 * @param[in]  range_to   Upper bound of position.
 * @param[out] tf         Transformation matrix.
 */
template <int dim>
inline void RandomRigidBodyTransform(Real range_from, Real range_to,
                                     Transformr<dim>& tf) {
  RandomRigidBodyTransform<dim>(Vecr<dim>::Constant(range_from),
                                Vecr<dim>::Constant(range_to), tf);
}

}  // namespace

}  // namespace dgd

#endif  // DGD_UTILS_H_

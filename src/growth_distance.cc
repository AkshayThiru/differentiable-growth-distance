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
 * @file growth_distance.cc
 * @author Akshay Thirugnanam (akshay_t@berkeley.edu)
 * @date 2025-08-05
 * @brief Growth distance algorithm implementation.
 */

#include "dgd/growth_distance.h"

#include "dgd/data_types.h"
#include "dgd/geometry/convex_set.h"
#include "dgd/output.h"
#include "dgd/settings.h"
#include "dgd/solvers/bundle_scheme_2d.h"
#include "dgd/solvers/bundle_scheme_3d.h"
#include "dgd/solvers/solver_options.h"

namespace dgd {

namespace {

using detail::BundleScheme;
using detail::SolverType;

// Growth distance algorithm template function.
template <int dim, class C1, class C2>
inline Real GrowthDistanceImpl(const C1* set1, const Transformr<dim>& tf1,
                               const C2* set2, const Transformr<dim>& tf2,
                               const Settings& settings, Output<dim>& out,
                               bool warm_start) {
  static_assert((dim == 2) || (dim == 3), "dim must be 2 or 3");
  static_assert((C1::dimension() == dim) && (C2::dimension() == dim),
                "Convex sets are not two or three-dimensional");

  if (set1->IsPolytopic() && set2->IsPolytopic()) {
    return BundleScheme<C1, C2, SolverType::CuttingPlane, false>(
        set1, tf1, set2, tf2, settings, out, warm_start);
  } else {
    return BundleScheme<C1, C2, SolverType::TrustRegionNewton, false>(
        set1, tf1, set2, tf2, settings, out, warm_start);
  }
}

// Collision detection template function.
template <int dim, class C1, class C2>
inline bool DetectCollisionImpl(const C1* set1, const Transformr<dim>& tf1,
                                const C2* set2, const Transformr<dim>& tf2,
                                const Settings& settings, Output<dim>& out,
                                bool warm_start) {
  static_assert((dim == 2) || (dim == 3), "dim must be 2 or 3");
  static_assert((C1::dimension() == dim) && (C2::dimension() == dim),
                "Convex sets are not two or three-dimensional");

  Real gd;
  if (set1->IsPolytopic() && set2->IsPolytopic()) {
    gd = BundleScheme<C1, C2, SolverType::CuttingPlane, true>(
        set1, tf1, set2, tf2, settings, out, warm_start);
  } else {
    gd = BundleScheme<C1, C2, SolverType::TrustRegionNewton, true>(
        set1, tf1, set2, tf2, settings, out, warm_start);
  }

  return ((out.status == SolutionStatus::CoincidentCenters) ||
          ((out.status == SolutionStatus::Optimal) && (gd > Real(0.0))));
}

}  // namespace

Real GrowthDistance(const ConvexSet<2>* set1, const Transform2r& tf1,
                    const ConvexSet<2>* set2, const Transform2r& tf2,
                    const Settings& settings, Output2d& out, bool warm_start) {
  return GrowthDistanceImpl(set1, tf1, set2, tf2, settings, out, warm_start);
}

Real GrowthDistance(const ConvexSet<3>* set1, const Transform3r& tf1,
                    const ConvexSet<3>* set2, const Transform3r& tf2,
                    const Settings& settings, Output3d& out, bool warm_start) {
  return GrowthDistanceImpl(set1, tf1, set2, tf2, settings, out, warm_start);
}

bool DetectCollision(const ConvexSet<2>* set1, const Transform2r& tf1,
                     const ConvexSet<2>* set2, const Transform2r& tf2,
                     const Settings& settings, Output2d& out, bool warm_start) {
  return DetectCollisionImpl(set1, tf1, set2, tf2, settings, out, warm_start);
}

bool DetectCollision(const ConvexSet<3>* set1, const Transform3r& tf1,
                     const ConvexSet<3>* set2, const Transform3r& tf2,
                     const Settings& settings, Output3d& out, bool warm_start) {
  return DetectCollisionImpl(set1, tf1, set2, tf2, settings, out, warm_start);
}

}  // namespace dgd

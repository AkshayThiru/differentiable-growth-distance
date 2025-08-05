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
 * @file growth_distance.h
 * @author Akshay Thirugnanam (akshay_t@berkeley.edu)
 * @date 2025-02-18
 * @brief Growth distance algorithm for two and three-dimensional convex sets.
 */

#ifndef DGD_GROWTH_DISTANCE_H_
#define DGD_GROWTH_DISTANCE_H_

#include "dgd/data_types.h"
#include "dgd/geometry/convex_set.h"
#include "dgd/output.h"
#include "dgd/settings.h"

namespace dgd {

/**
 * @name Growth distance algorithm for convex sets
 * @brief Growth distance algorithm for two and three-dimensional convex sets.
 *
 * @attention When using warm-start, the following properties must be ensured:
 * The same output must be reused from the previous function call;
 * The output must not be used for other pairs of sets in between function
 * calls;
 * The order of the sets must be the same.
 *
 * @note Output from a previous collision detection call can be used to warm
 * start the growth distance algorithm.
 *
 * @param[in]     set1,set2  Convex sets.
 * @param[in]     tf1,tf2    Rigid body transformations for the sets.
 * @param[in]     settings   Settings.
 * @param[in,out] out        Output.
 * @param         warm_start Whether to use previous output for warm start.
 * @return        Growth distance lower bound.
 */
///@{
Real GrowthDistance(const ConvexSet<2>* set1, const Transform2r& tf1,
                    const ConvexSet<2>* set2, const Transform2r& tf2,
                    const Settings& settings, Output2d& out,
                    bool warm_start = false);

Real GrowthDistance(const ConvexSet<3>* set1, const Transform3r& tf1,
                    const ConvexSet<3>* set2, const Transform3r& tf2,
                    const Settings& settings, Output3d& out,
                    bool warm_start = false);
///@}

/**
 * @name Boolean collision detection algorithm for convex sets
 * @brief Collision detection algorithm for two and three-dimensional convex
 * sets.
 *
 * Returns true if the centers coincide or if the sets intersect;
 * false if a separating plane has been found or if the maximum number of
 * iterations have been reached.
 *
 * @note Output from a previous growth distance call can be used to warm
 * start the collision detection function.
 *
 * @param[in]     set1,set2  Convex sets.
 * @param[in]     tf1,tf2    Rigid body transformations for the sets.
 * @param[in]     settings   Settings.
 * @param[in,out] out        Output.
 * @param         warm_start Whether to use previous output for warm start.
 * @return        true, if the sets are colliding; false, otherwise.
 */
///@{
bool DetectCollision(const ConvexSet<2>* set1, const Transform2r& tf1,
                     const ConvexSet<2>* set2, const Transform2r& tf2,
                     const Settings& settings, Output2d& out,
                     bool warm_start = false);

bool DetectCollision(const ConvexSet<3>* set1, const Transform3r& tf1,
                     const ConvexSet<3>* set2, const Transform3r& tf2,
                     const Settings& settings, Output3d& out,
                     bool warm_start = false);
///@}

}  // namespace dgd

#endif  // DGD_GROWTH_DISTANCE_H_

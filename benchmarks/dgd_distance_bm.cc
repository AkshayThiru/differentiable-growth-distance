#include <cmath>
#include <functional>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "dgd/data_types.h"
#include "dgd/error_metrics.h"
#include "dgd/geometry/3d/ellipsoid.h"
#include "dgd/geometry/3d/frustum.h"
#include "dgd/geometry/3d/mesh.h"
#include "dgd/geometry/convex_set.h"
#include "dgd/growth_distance.h"
#include "dgd/output.h"
#include "dgd/settings.h"
#include "dgd/utils/timer.h"
#include "internal_helpers/debug.h"
#include "internal_helpers/filesystem_utils.h"
#include "internal_helpers/random_utils.h"
#include "internal_helpers/set_generator.h"

// Benchmarks to run.
// const bool dynamic_dispatch = false;
const bool cold_start = true;
const bool warm_start = true;
const bool two_dim_sets = true;
const bool three_dim_sets = true;

// Printing.
const bool print_suboptimal_run = true;
const bool exit_on_suboptimal_run = true;

// Constants.
const double position_lim = 5.0;
const double dx_max = 0.1, ang_max = dgd::kPi / 18.0;

//  Number of pairs of sets to benchmark.
const int npair = 10000;
//  Number of poses per set pair for cold-start.
const int npose_c = 1000;
//  Number of poses per set pair for warm-start.
const int npose_w = 1000;
//  Number of cold-start function calls for a given pair of sets and poses.
const int ncold = 100;
//  Number of warm-start function calls for a given pair of sets.
const int nwarm = 100;

template <class C>
using SetPtr = std::shared_ptr<C>;

template <int dim>
struct OptimalSolution {
  dgd::Vecr<dim> z1, z2, normal;
  dgd::Real gd;
  dgd::SolutionStatus status;

  void SetFromOutput(const dgd::Output<dim> out) {
    z1 = out.z1;
    z2 = out.z2;
    normal = out.normal;
    gd = out.growth_dist_ub;
    status = out.status;
  }

  void SetOutput(dgd::Output<dim>& out) const {
    out.z1 = z1;
    out.z2 = z2;
    out.normal = normal;
    out.growth_dist_ub = gd;
    out.status = status;
  }
};

template <int dim, class C1, class C2>
void ColdStart(std::function<const SetPtr<C1>()> generator1,
               std::function<const SetPtr<C2>()> generator2, dgd::Rng& rng,
               int npair, int npose) {
  double avg_solve_time = 0.0;
  double max_prim_dual_gap = 0.0;
  double max_prim_infeas_err = 0.0;
  double max_dual_infeas_err = 0.0;
  double avg_iter = 0.0;
  int nsubopt = 0;

  dgd::Timer timer;
  timer.Stop();
  dgd::Transformr<dim> tf1, tf2;
  dgd::Settings settings;
  dgd::Output<dim> out;
  for (int i = 0; i < npair; ++i) {
    const SetPtr<C1> set1 = generator1();
    const SetPtr<C2> set2 = generator2();
    for (int j = 0; j < npose; ++j) {
      dgd::internal::SetRandomRigidBodyTransforms(rng, tf1, tf2, -position_lim,
                                                  position_lim);
      // Initial call to reduce cache misses.
      dgd::GrowthDistance(set1.get(), tf1, set2.get(), tf2, settings, out);
      timer.Start();
      for (int k = 0; k < ncold; ++k) {
        dgd::GrowthDistance(set1.get(), tf1, set2.get(), tf2, settings, out);
      }
      timer.Stop();
      const auto err =
          dgd::ComputeSolutionError(set1.get(), tf1, set2.get(), tf2, out);

      if (print_suboptimal_run) {
        if (out.status == dgd::SolutionStatus::MaxIterReached) {
          dgd::internal::PrintSetup(set1.get(), tf1, set2.get(), tf2, settings,
                                    out, err);
          if (exit_on_suboptimal_run) return;
        }
      }

      avg_solve_time += timer.Elapsed() / double(ncold);
      max_prim_dual_gap = std::max(max_prim_dual_gap, err.prim_dual_gap);
      max_prim_infeas_err = std::max(max_prim_infeas_err, err.prim_infeas_err);
      max_dual_infeas_err = std::max(max_dual_infeas_err, err.dual_infeas_err);
      avg_iter += out.iter;
      nsubopt += (out.status != dgd::SolutionStatus::Optimal) &&
                 (out.status != dgd::SolutionStatus::CoincidentCenters);
    }
  }
  avg_solve_time = avg_solve_time / (npair * npose);
  avg_iter = avg_iter / (npair * npose);

  dgd::internal::PrintStatistics(avg_solve_time, max_prim_dual_gap,
                                 max_prim_infeas_err, max_dual_infeas_err,
                                 avg_iter, nsubopt);
}

template <int dim, class C1, class C2>
void WarmStart(std::function<const SetPtr<C1>()> generator1,
               std::function<const SetPtr<C2>()> generator2, dgd::Rng& rng,
               int npair, int npose) {
  double avg_solve_time = 0.0;
  double max_prim_dual_gap = 0.0;
  double max_prim_infeas_err = 0.0;
  double max_dual_infeas_err = 0.0;
  double avg_iter = 0.0;
  int nsubopt = 0;

  dgd::Timer timer;
  timer.Stop();
  dgd::Transformr<dim> tf1, tf1_c, tf2;
  dgd::Vecr<dim> dx;
  dgd::Rotationr<dim> drot;
  dgd::Settings settings;
  dgd::Output<dim> out;
  std::vector<OptimalSolution<dim>> opt_sols(nwarm);
  for (int i = 0; i < npair; ++i) {
    const SetPtr<C1> set1 = generator1();
    const SetPtr<C2> set2 = generator2();
    for (int j = 0; j < npose; ++j) {
      dgd::internal::SetRandomRigidBodyTransforms(rng, tf1, tf2, -position_lim,
                                                  position_lim);
      tf1_c = tf1;
      dgd::internal::SetRandomDisplacement(rng, dx, drot, dx_max, ang_max);
      // Initial cold-start call.
      dgd::GrowthDistance(set1.get(), tf1, set2.get(), tf2, settings, out);
      int total_iter = 0;
      timer.Start();
      for (int k = 0; k < nwarm; ++k) {
        dgd::internal::UpdateRigidBodyTransform(tf1, dx, drot);
        dgd::GrowthDistance(set1.get(), tf1, set2.get(), tf2, settings, out,
                            true);
        opt_sols[k].SetFromOutput(out);
        total_iter += out.iter;
      }
      timer.Stop();
      avg_solve_time += timer.Elapsed() / double(nwarm);
      avg_iter += total_iter / double(nwarm);

      tf1 = tf1_c;
      for (int k = 0; k < nwarm; ++k) {
        dgd::internal::UpdateRigidBodyTransform(tf1, dx, drot);
        opt_sols[k].SetOutput(out);

        const auto err =
            dgd::ComputeSolutionError(set1.get(), tf1, set2.get(), tf2, out);

        if (print_suboptimal_run) {
          if (out.status == dgd::SolutionStatus::MaxIterReached) {
            if (k == 0) {
              dgd::internal::PrintSetup<dim>(set1.get(), tf1, set2.get(), tf2,
                                             settings, out, err);
            } else {
              dgd::Output<dim> out_prev;
              opt_sols[k - 1].SetOutput(out_prev);
              dgd::internal::PrintSetup(set1.get(), tf1, set2.get(), tf2,
                                        settings, out, err, &out_prev, true);
            }
            if (exit_on_suboptimal_run) return;
          }
        }

        max_prim_dual_gap = std::max(max_prim_dual_gap, err.prim_dual_gap);
        max_prim_infeas_err =
            std::max(max_prim_infeas_err, err.prim_infeas_err);
        max_dual_infeas_err =
            std::max(max_dual_infeas_err, err.dual_infeas_err);
        nsubopt += (out.status != dgd::SolutionStatus::Optimal) &&
                   (out.status != dgd::SolutionStatus::CoincidentCenters);
      }
    }
  }
  avg_solve_time = avg_solve_time / (npair * npose);
  avg_iter = avg_iter / (npair * npose);

  dgd::internal::PrintStatistics(avg_solve_time, max_prim_dual_gap,
                                 max_prim_infeas_err, max_dual_infeas_err,
                                 avg_iter, nsubopt);
}

int main(int argc, char** argv) {
  std::string path = "./";
  if (argc < 2) {
    std::cout << "No directory specified; using current directory" << std::endl;
  } else if (argc == 2) {
    path = argv[1];
    if (!dgd::internal::IsValidDirectory(path)) return EXIT_FAILURE;
  } else {
    std::cerr << "Usage: " << argv[0] << " <asset_folder_path>" << std::endl;
    return EXIT_FAILURE;
  }

  std::vector<std::string> filenames;
  dgd::internal::GetObjFileNames(path, filenames);

  dgd::internal::ConvexSetGenerator gen(dgd::internal::ConvexSetFeatureRange{});
  gen.LoadMeshesFromObjFiles(filenames);
  if (gen.nmeshes() == 0) {
    std::cout << "No .obj files found; mesh benchmarks will not be run"
              << std::endl;
  }
  dgd::Rng rng;
  auto set_default_seed = [&gen, &rng]() -> void {
    gen.SetDefaultRngSeed();
    rng.SetDefaultSeed();
  };
  auto set_random_seed = [&gen, &rng]() -> void {
    gen.SetRandomRngSeed();
    rng.SetRandomSeed();
  };

  // Benchmarks.
  /*
  //  1. Dynamic dispatch: (frustum + ellipsoid, cold-start).
  if (dynamic_dispatch) {
    set_default_seed();
    auto generator1 = [&gen]() -> const SetPtr<dgd::ConvexSet<3>> {
      return gen.GetPrimitiveSet(dgd::internal::CurvedPrimitive3D::Frustum);
    };
    auto generator2 = [&gen]() -> const SetPtr<dgd::ConvexSet<3>> {
      return gen.GetPrimitiveSet(dgd::internal::CurvedPrimitive3D::Ellipsoid);
    };
    std::cout << "Dynamic dispatch: (frustum + ellipsoid, cold-start)"
              << std::endl;
    ColdStart<3, dgd::ConvexSet<3>, dgd::ConvexSet<3>>(generator1, generator2,
                                                       rng, npair, npose_c);
  }
  //  2. Inline call: (frustum + ellipsoid, cold-start).
  if (dynamic_dispatch) {
    set_default_seed();
    auto generator1 = [&gen]() -> const SetPtr<dgd::Frustum> {
      const auto set1 =
          gen.GetPrimitiveSet(dgd::internal::CurvedPrimitive3D::Frustum);
      return std::static_pointer_cast<dgd::Frustum>(set1);
    };
    auto generator2 = [&gen]() -> const SetPtr<dgd::Ellipsoid> {
      const auto set2 =
          gen.GetPrimitiveSet(dgd::internal::CurvedPrimitive3D::Ellipsoid);
      return std::static_pointer_cast<dgd::Ellipsoid>(set2);
    };
    std::cout << "Inline call     : (frustum + ellipsoid, cold-start)"
              << std::endl;
    ColdStart<3, dgd::Frustum, dgd::Ellipsoid>(generator1, generator2, rng,
                                               npair, npose_c);
  }

  //  3. Dynamic dispatch: (mesh + mesh, cold-start).
  if (dynamic_dispatch && (gen.nmeshes() > 0)) {
    set_default_seed();
    auto generator = [&gen]() -> const SetPtr<dgd::ConvexSet<3>> {
      return gen.GetRandomMeshSet();
    };
    std::cout << "Dynamic dispatch: (mesh + mesh, cold-start)" << std::endl;
    ColdStart<3, dgd::ConvexSet<3>, dgd::ConvexSet<3>>(generator, generator,
                                                       rng, npair, npose_c);
  }
  //  4. Inline call: (mesh + mesh, cold-start).
  if (dynamic_dispatch && (gen.nmeshes() > 0)) {
    set_default_seed();
    auto generator = [&gen]() -> const SetPtr<dgd::Mesh> {
      return std::static_pointer_cast<dgd::Mesh>(gen.GetRandomMeshSet());
    };
    std::cout << "Inline call     : (mesh + mesh, cold-start)" << std::endl;
    ColdStart<3, dgd::Mesh, dgd::Mesh>(generator, generator, rng, npair,
                                       npose_c);
  }
  */

  std::cout << std::string(50, '-') << std::endl;
  set_random_seed();
  //  5. Cold-start: (2D primitive + 2D primitive).
  if (cold_start && two_dim_sets) {
    auto generator = [&gen]() -> const SetPtr<dgd::ConvexSet<2>> {
      return gen.GetRandom2DSet();
    };
    std::cout << "Cold-start: (2D primitive + 2D primitive)" << std::endl;
    ColdStart<2, dgd::ConvexSet<2>, dgd::ConvexSet<2>>(generator, generator,
                                                       rng, npair, npose_c);
  }
  //  6. Cold-start: (3D curved primitive + 3D curved primitive).
  if (cold_start && three_dim_sets) {
    auto generator = [&gen]() -> const SetPtr<dgd::ConvexSet<3>> {
      return gen.GetRandomCurvedPrimitive3DSet();
    };
    std::cout << "Cold-start: (3D curved primitive + 3D curved primitive)"
              << std::endl;
    ColdStart<3, dgd::ConvexSet<3>, dgd::ConvexSet<3>>(generator, generator,
                                                       rng, npair, npose_w);
  }
  //  7. Cold-start: (3D primitive + 3D primitive).
  if (cold_start && three_dim_sets) {
    auto generator = [&gen]() -> const SetPtr<dgd::ConvexSet<3>> {
      return gen.GetRandomPrimitive3DSet();
    };
    std::cout << "Cold-start: (3D primitive + 3D primitive)" << std::endl;
    ColdStart<3, dgd::ConvexSet<3>, dgd::ConvexSet<3>>(generator, generator,
                                                       rng, npair, npose_w);
  }
  //  8. Cold-start: (mesh + mesh).
  if (cold_start && three_dim_sets && (gen.nmeshes() > 0)) {
    auto generator = [&gen]() -> const SetPtr<dgd::ConvexSet<3>> {
      return gen.GetRandomMeshSet();
    };
    std::cout << "Cold-start: (mesh + mesh)" << std::endl;
    ColdStart<3, dgd::ConvexSet<3>, dgd::ConvexSet<3>>(generator, generator,
                                                       rng, npair, npose_c);
  }

  std::cout << std::string(50, '-') << std::endl;
  //  9. Warm-start: (2D primitive + 2D primitive).
  if (warm_start && two_dim_sets) {
    auto generator = [&gen]() -> const SetPtr<dgd::ConvexSet<2>> {
      return gen.GetRandom2DSet();
    };
    std::cout << "Warm-start: (2D primitive + 2D primitive)" << std::endl;
    WarmStart<2, dgd::ConvexSet<2>, dgd::ConvexSet<2>>(generator, generator,
                                                       rng, npair, npose_w);
  }
  //  10. Warm-start: (3D curved primitive + 3D curved primitive).
  if (warm_start && three_dim_sets) {
    auto generator = [&gen]() -> const SetPtr<dgd::ConvexSet<3>> {
      return gen.GetRandomCurvedPrimitive3DSet();
    };
    std::cout << "Warm-start: (3D curved primitive + 3D curved primitive)"
              << std::endl;
    WarmStart<3, dgd::ConvexSet<3>, dgd::ConvexSet<3>>(generator, generator,
                                                       rng, npair, npose_w);
  }
  //  11. Warm-start: (mesh + mesh).
  if (warm_start && three_dim_sets && (gen.nmeshes() > 0)) {
    auto generator = [&gen]() -> const SetPtr<dgd::ConvexSet<3>> {
      return gen.GetRandomMeshSet();
    };
    std::cout << "Warm-start: (mesh + mesh)" << std::endl;
    WarmStart<3, dgd::ConvexSet<3>, dgd::ConvexSet<3>>(generator, generator,
                                                       rng, npair, npose_w);
  }
}

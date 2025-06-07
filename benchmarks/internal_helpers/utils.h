#ifndef DGD_BENCHMARKS_INTERNAL_HELPERS_UTILS_H_
#define DGD_BENCHMARKS_INTERNAL_HELPERS_UTILS_H_

#include <filesystem>
#include <iostream>
#include <vector>

#include "dgd/data_types.h"
#include "dgd/utils.h"

namespace dgd {

namespace internal {

// Checks if the given folder path is a valid directory.
inline bool IsValidDirectory(const std::string& path) {
  if (!std::filesystem::exists(path) || !std::filesystem::is_directory(path)) {
    std::cerr << "Error: The specified folder path is not a valid directory"
              << std::endl;
    return false;
  }
  return true;
}

// Gets the .obj file names from the given folder path.
inline void GetObjFileNames(const std::string& path,
                            std::vector<std::string>& filenames) {
  filenames.clear();
  for (const auto& entry : std::filesystem::directory_iterator(path)) {
    if (entry.is_regular_file() && entry.path().extension() == ".obj") {
      filenames.push_back(path + entry.path().filename().string());
    }
  }
}

// Sets random rigid body transforms for a pair of rigid bodies.
template <int dim>
inline void SetRandomTransform(Transformr<dim>& tf1, Transformr<dim>& tf2,
                               Real range_from, Real range_to) {
  RandomRigidBodyTransform<dim>(range_from, range_to, tf1);
  RandomRigidBodyTransform<dim>(range_from, range_to, tf2);
}

// Sets random small displacement.
inline void SetRandomScrew(Vec2r& dx, Rotation2r& drot, Real dx_max,
                           Real ang_max) {
  RandomUnitVector(dx);
  dx *= Random(dx_max);
  const Real ang = Random(ang_max);
  drot(0, 0) = std::cos(ang);
  drot(1, 0) = std::sin(ang);
  drot(0, 1) = -drot(1, 0);
  drot(1, 1) = drot(0, 0);
}

inline void SetRandomScrew(Vec3r& dx, Rotation3r& drot, Real dx_max,
                           Real ang_max) {
  RandomUnitVector(dx);
  dx *= Random(dx_max);
  EulerToRotation(ang_max * Vec3r(Random(), Random(0.5), Random()), drot);
}

// Updates the rigid body transform using the screw.
template <int dim>
inline void UpdateTransform(Transformr<dim>& tf, const Vecr<dim>& dx,
                            const Rotationr<dim>& drot) {
  tf.template block<dim, 1>(0, dim) += dx;
  tf.template block<dim, dim>(0, 0) *= drot;
}

}  // namespace internal

}  // namespace dgd

#endif  // DGD_BENCHMARKS_INTERNAL_HELPERS_UTILS_H_

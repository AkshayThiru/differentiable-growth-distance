#include "internal_helpers/set_generator.h"

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "dgd/geometry/2d/ellipse.h"
#include "dgd/geometry/2d/polygon.h"
#include "dgd/geometry/2d/rectangle.h"
#include "dgd/geometry/3d/cone.h"
#include "dgd/geometry/3d/cuboid.h"
#include "dgd/geometry/3d/cylinder.h"
#include "dgd/geometry/3d/ellipsoid.h"
#include "dgd/geometry/3d/frustum.h"
#include "dgd/geometry/3d/mesh.h"
#include "dgd/geometry/3d/polytope.h"
#include "dgd/geometry/xd/capsule.h"
#include "dgd/geometry/xd/sphere.h"
#include "dgd/graham_scan.h"
#include "dgd/utils.h"
#include "internal_helpers/mesh_loader.h"

namespace dgd {

namespace internal {

namespace {

inline Real Random(Range<Real> range) {
  return dgd::Random(range.low, range.high);
}

inline int Random(Range<int> range) {
  std::uniform_int_distribution<int> dis(range.low, range.high);
  return dis(generator);
}

template <int dim>
inline void CenterVertices(std::vector<Vecr<dim>>& vert) {
  Vecr<dim> center;
  for (const auto& v : vert) {
    center += v;
  }
  center /= static_cast<Real>(vert.size());
  for (auto& v : vert) {
    v -= center;
  }
}

}  // namespace

ConvexSetGenerator::ConvexSetGenerator(const ConvexSetFeatureRange& fr)
    : fr_(fr),
      count2_(static_cast<int>(Primitive2D::Count_)),
      ccount3_(static_cast<int>(CurvedPrimitive3D::Count_)),
      fcount3_(static_cast<int>(FlatPrimitive3D::Count_)) {
  SetDefaultSeed();
  if ((fr.polytope.nvert.low < 4) || (fr.polygon.nvert.low < 3)) {
    throw std::invalid_argument(
        "Polytopes (polygons) must have at least 4 (3) vertices.");
  }
  /*
  polygon_vert_.reserve(fr.polygon.nvert.high);
  polytope_vert_.reserve(fr.polytope.nvert.high);
  */
  meshes_.clear();
  nmeshes_ = 0;
}

void ConvexSetGenerator::SetDefaultRngSeed() const { SetDefaultSeed(); }

void ConvexSetGenerator::LoadMeshesFromObjFiles(
    const std::vector<std::string>& filenames) {
  meshes_.clear();

  MeshProperties mp{};
  for (const auto& filename : filenames) {
    try {
      mp.SetVertexMeshFromObjFile(filename);
    } catch (const std::runtime_error& e) {
      continue;
    }
    mp.SetZeroVertexCenter();

    meshes_.push_back(
        std::make_shared<Mesh>(mp.vert, mp.vgraph, 0.0, mp.inradius));
  }
  nmeshes_ = static_cast<int>(meshes_.size());
}

ConvexSetPtr<2> ConvexSetGenerator::GetPrimitiveSet(Primitive2D type) {
  switch (type) {
    case Primitive2D::Ellipse:
      return std::make_shared<Ellipse>(Random(fr_.ellipse.hlx),
                                       Random(fr_.ellipse.hly), 0.0);

      /*
      case Primitive2D::Polygon: {
        polygon_vert_.clear();

        const int nvert = Random(fr_.polygon.nvert);
        for (int i = 0; i < nvert; ++i) {
          Vec2r v{dgd::Random(), dgd::Random()};
          v = v * fr_.polygon.size / (v.lpNorm<4>() + kEps);
          polygon_vert_.push_back(v);
        }
        CenterVertices(polygon_vert_);

        return std::make_shared<Polygon>(polygon_vert_, 0.0, kEps);
      }
      */

    case Primitive2D::Rectangle:
      return std::make_shared<Rectangle>(Random(fr_.rectangle.hlx),
                                         Random(fr_.rectangle.hly), 0.0);

    case Primitive2D::Stadium:
      return std::make_shared<Stadium>(Random(fr_.capsule.hlx),
                                       Random(fr_.capsule.radius), 0.0);

    case Primitive2D::Circle:
      return std::make_shared<Circle>(Random(fr_.sphere.radius));

    default:
      throw std::invalid_argument("Invalid 2D primitive type");
  }
}

ConvexSetPtr<3> ConvexSetGenerator::GetPrimitiveSet(
    CurvedPrimitive3D type) const {
  switch (type) {
    case CurvedPrimitive3D::Capsule:
      return std::make_shared<Capsule>(Random(fr_.capsule.hlx),
                                       Random(fr_.capsule.radius), 0.0);

    case CurvedPrimitive3D::Cone:
      return std::make_shared<Cone>(Random(fr_.cone.radius),
                                    Random(fr_.cone.height), 0.0);

    case CurvedPrimitive3D::Cylinder:
      return std::make_shared<Cylinder>(Random(fr_.cylinder.hlx),
                                        Random(fr_.cylinder.radius), 0.0);

    case CurvedPrimitive3D::Ellipsoid:
      return std::make_shared<Ellipsoid>(Random(fr_.ellipsoid.hlx),
                                         Random(fr_.ellipsoid.hly),
                                         Random(fr_.ellipsoid.hlz), 0.0);

    case CurvedPrimitive3D::Frustum:
      return std::make_shared<Frustum>(Random(fr_.frustum.base_radius),
                                       Random(fr_.frustum.top_radius),
                                       Random(fr_.frustum.height), 0.0);

    case CurvedPrimitive3D::Sphere:
      return std::make_shared<Sphere>(Random(fr_.sphere.radius));

    default:
      throw std::invalid_argument("Invalid curved 3D primitive type");
  }
}

ConvexSetPtr<3> ConvexSetGenerator::GetPrimitiveSet(FlatPrimitive3D type) {
  switch (type) {
    case FlatPrimitive3D::Cuboid:
      return std::make_shared<Cuboid>(Random(fr_.cuboid.hlx),
                                      Random(fr_.cuboid.hly),
                                      Random(fr_.cuboid.hlz), 0.0);

      /*
      case FlatPrimitive3D::Polytope: {
        polytope_vert_.clear();

        if (meshes_.empty()) {
          const int nvert = Random(fr_.polytope.nvert);
          for (int i = 0; i < nvert; ++i) {
            Vec3r v{dgd::Random(), dgd::Random(), dgd::Random()};
            v = v * fr_.polytope.size / (v.lpNorm<4>() + kEps);
            polytope_vert_.push_back(v);
          }
        } else {
          const int mesh_idx = Random(Range<int>{0, nmeshes_ - 1});
          const int nvert_m = meshes_[mesh_idx]->nvertices();
          for (int i = 0; i < 4; ++i) {
            polytope_vert_.push_back(
                meshes_[mesh_idx]->vertices()[(i * nvert_m) / 3]);
          }
          for (int i = 4; i < fr_.polytope.nvert.high; ++i) {
            const int idx = Random(Range<int>{0, nvert_m - 1});
            polytope_vert_.push_back(meshes_[mesh_idx]->vertices()[idx]);
          }
        }
        CenterVertices(polytope_vert_);

        return std::make_shared<Polytope>(polytope_vert_, 0.0, kEps);
      }
      */

    default:
      throw std::invalid_argument("Invalid flat 3D primitive type");
  }
}

ConvexSetPtr<3> ConvexSetGenerator::GetRandomCurvedPrimitive3DSet() const {
  const int idx = Random(Range<int>{0, ccount3_ - 1});
  return GetPrimitiveSet(static_cast<CurvedPrimitive3D>(idx));
}

ConvexSetPtr<3> ConvexSetGenerator::GetRandomPrimitive3DSet() {
  const int idx = Random(Range<int>{0, ccount3_ + fcount3_ - 1});
  if (idx < ccount3_) {
    return GetPrimitiveSet(static_cast<CurvedPrimitive3D>(idx));
  } else {
    return GetPrimitiveSet(static_cast<FlatPrimitive3D>(idx - ccount3_));
  }
}

ConvexSetPtr<3> ConvexSetGenerator::GetRandomMeshSet(int* idx) const {
  if (meshes_.empty()) {
    throw std::runtime_error(
        "No meshes available to retrieve a random mesh set");
  }

  const int mesh_idx = Random(Range<int>{0, nmeshes_ - 1});
  if (idx) *idx = mesh_idx;
  return meshes_[mesh_idx];
}

ConvexSetPtr<2> ConvexSetGenerator::GetRandom2DSet() {
  const int idx = Random(Range<int>{0, count2_ - 1});
  return GetPrimitiveSet(static_cast<Primitive2D>(idx));
}

ConvexSetPtr<3> ConvexSetGenerator::GetRandom3DSet() {
  const int idx = Random(Range<int>{0, ccount3_ + fcount3_ + nmeshes_ - 1});
  if (idx < ccount3_ + fcount3_) {
    return GetRandomPrimitive3DSet();
  } else {
    return GetRandomMeshSet();
  }
}

}  // namespace internal

}  // namespace dgd

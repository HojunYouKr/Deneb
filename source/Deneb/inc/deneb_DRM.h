#pragma once

#include <memory>
#include <vector>

#include "deneb_element.h"


namespace deneb {
class Basis;
class BasisSurface;

class DRM {
 public:
  DRM(){};
  virtual ~DRM(){};
  virtual void GetDRMPoints(const int num_points,
                            std::vector<double>& points) const = 0;
};

class DRMVolume : public DRM {
 public:
  static std::shared_ptr<DRMVolume> GetDRM(const ElemType elemtype);

 protected:
  std::shared_ptr<Basis> approx_basis_;

 public:
  DRMVolume() : approx_basis_(nullptr){};
  virtual ~DRMVolume(){};
  std::shared_ptr<Basis> GetApproximateBasis() { return approx_basis_; };
};

// DRMVolumeTris: Poly2D, alpha-optimized points
class DRMVolumeTris : public DRMVolume {
 public:
  DRMVolumeTris();
  virtual ~DRMVolumeTris(){};
  virtual void GetDRMPoints(const int num_points,
                            std::vector<double>& points) const;
};
// DRMVolumeQuad: PolyQuadShape, Gauss-Legendre points
class DRMVolumeQuad : public DRMVolume {
 public:
  DRMVolumeQuad();
  virtual ~DRMVolumeQuad(){};
  virtual void GetDRMPoints(const int num_points,
                            std::vector<double>& points) const;
};
// DRMVolumeTets: Poly3D, alpha-optimized points
class DRMVolumeTets : public DRMVolume {
 public:
  DRMVolumeTets();
  virtual ~DRMVolumeTets(){};
  virtual void GetDRMPoints(const int num_points,
                            std::vector<double>& points) const;
};
// DRMVolumeHexa: PolyHexaShape, Gauss-Legendre points
class DRMVolumeHexa : public DRMVolume {
 public:
  DRMVolumeHexa();
  virtual ~DRMVolumeHexa(){};
  virtual void GetDRMPoints(const int num_points,
                            std::vector<double>& points) const;
};
// DRMVolumePris: PolyPrisShape, Gauss-Legendre points * alpha-optimized points
class DRMVolumePris : public DRMVolume {
 public:
  DRMVolumePris();
  virtual ~DRMVolumePris(){};
  virtual void GetDRMPoints(const int num_points,
                            std::vector<double>& points) const;
};
// DRMVolumePyra: PolyPyra, pyramid DRM-SFP points
class DRMVolumePyra : public DRMVolume {
 public:
  DRMVolumePyra();
  virtual ~DRMVolumePyra(){};
  virtual void GetDRMPoints(const int num_points,
                            std::vector<double>& points) const;
};

class DRMSurface : public DRM {
 public:
  static std::shared_ptr<DRMSurface> GetDRM(const ElemType elemtype);

 protected:
  std::shared_ptr<BasisSurface> approx_basis_;

 public:
  DRMSurface() : approx_basis_(nullptr){};
  virtual ~DRMSurface(){};
  std::shared_ptr<BasisSurface> GetApproximateBasis() { return approx_basis_; };
};

// DRMSurfaceLine: Poly2D, Gauss-Legendre points
class DRMSurfaceLine : public DRMSurface {
 public:
  DRMSurfaceLine();
  virtual ~DRMSurfaceLine(){};
  virtual void GetDRMPoints(const int num_points,
                            std::vector<double>& points) const;
};
// DRMSurfaceTris: Poly3D, Witherden-Vincent points
class DRMSurfaceTris : public DRMSurface {
 public:
  DRMSurfaceTris();
  virtual ~DRMSurfaceTris(){};
  virtual void GetDRMPoints(const int num_points,
                            std::vector<double>& points) const;
};
// DRMSurfaceQuad: Poly3D, Witherden-Vincent points
class DRMSurfaceQuad : public DRMSurface {
 public:
  DRMSurfaceQuad();
  virtual ~DRMSurfaceQuad(){};
  virtual void GetDRMPoints(const int num_points,
                            std::vector<double>& points) const;
};
}  // namespace deneb

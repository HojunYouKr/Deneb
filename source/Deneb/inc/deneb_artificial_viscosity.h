#pragma once

#include <memory>
#include <string>
#include <vector>

#define DENEB_ARTIFICIAL_VISCOSITY_NAME artificial_viscosity_global_ptr
#define DENEB_ARTIFICIAL_VISCOSITY deneb::DENEB_ARTIFICIAL_VISCOSITY_NAME
#define DENEB_ARTIFICIAL_VISCOSITY_INITIALIZE(name) \
  DENEB_ARTIFICIAL_VISCOSITY =                      \
      deneb::ArtificialViscosity::GetArtificialViscosity(name)
#define DENEB_ARTIFICIAL_VISCOSITY_FINALIZE() DENEB_ARTIFICIAL_VISCOSITY.reset()

namespace avocado {
class Communicate;
}

namespace deneb {
class ArtificialViscosity {
 protected:
  std::vector<double> artificial_viscosity_;

 public:
  static std::shared_ptr<ArtificialViscosity> GetArtificialViscosity(
      const std::string& name);

  ArtificialViscosity(){};
  virtual ~ArtificialViscosity(){};

  virtual void BuildData(void) = 0;
  virtual void ComputeArtificialViscosity(const double* solution) = 0;

  virtual double GetArtificialViscosityValue(const int icell,
                                             const int ipoint) const = 0;
  const std::vector<double>& GetArtificialViscosity() const {
    return artificial_viscosity_;
  };
};
extern std::shared_ptr<ArtificialViscosity> DENEB_ARTIFICIAL_VISCOSITY_NAME;

class NoArtificialViscosity : public ArtificialViscosity {
 public:
  NoArtificialViscosity();
  virtual ~NoArtificialViscosity(){};

  virtual void BuildData(void);
  virtual void ComputeArtificialViscosity(const double* solution){};

  virtual double GetArtificialViscosityValue(const int icell,
                                             const int ipoint) const {
    return 0.0;
  };
};

// ArtificialViscosity = Peclet, kappa
class LaplacianP0 : public ArtificialViscosity {
 protected:
  int target_state_;
  int num_bases_m1_;  // num basis of P(n-1)
  double Peclet_;
  double kappa_;
  double S0_;
  double dLmax_;

  std::shared_ptr<avocado::Communicate> communicate_;

 public:
  LaplacianP0();
  virtual ~LaplacianP0(){};

  virtual void BuildData(void);
  virtual void ComputeArtificialViscosity(const double* solution);

  virtual double GetArtificialViscosityValue(const int icell,
                                             const int ipoint) const {
    if (icell >= 0)
      return artificial_viscosity_[icell];
    else
      return 0.0;
  };

 protected:
  double SmoothnessIndicator(const double* solution);

  double MaxArtificialViscosity(const double* solution,
                                const double cell_volumes,
                                const double cell_basis_value);
};
}  // namespace deneb
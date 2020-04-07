#pragma once

#include "GeneralUserObject.h"

/**
 * Computes the partitioning of the interphase heat exchange coefficients
 */
class HeatExchangeCoefficientPartitioning : public GeneralUserObject
{
public:
  HeatExchangeCoefficientPartitioning(const InputParameters & parameters);

  virtual void execute(){};
  virtual void initialize(){};
  virtual void finalize(){};

  virtual Real getPartition(Real alpha_liquid, Real dalpha_liquid_dt) const;
  virtual Real getPartitionDer(Real alpha_liquid, Real dalpha_liquid_dt, Real area) const;

protected:
  /// Upper cut-off limit
  Real _lower;
  /// Lower cut-off limit
  Real _upper;

public:
  static InputParameters validParams();
};

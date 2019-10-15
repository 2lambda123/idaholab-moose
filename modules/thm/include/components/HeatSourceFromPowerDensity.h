#pragma once

#include "HeatSourceBase.h"

class HeatSourceFromPowerDensity;

template <>
InputParameters validParams<HeatSourceFromPowerDensity>();

/**
 * Heat source from power density
 */
class HeatSourceFromPowerDensity : public HeatSourceBase
{
public:
  HeatSourceFromPowerDensity(const InputParameters & parameters);

  virtual void addMooseObjects() override;

protected:
  /// The name of the power density variable (typically an aux variable)
  const VariableName _power_density_name;
};

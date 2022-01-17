#pragma once

#include "OneDNodalBC.h"

/**
 * Mass flow rate (m_dot) and temperature (T) BC
 */
class OneDMomentumMassFlowRateTemperatureBC : public OneDNodalBC
{
public:
  OneDMomentumMassFlowRateTemperatureBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

  const Real & _m_dot;

public:
  static InputParameters validParams();
};

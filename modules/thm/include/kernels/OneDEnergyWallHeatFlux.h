#ifndef ONEDENERGYWALLHEATFLUX_H
#define ONEDENERGYWALLHEATFLUX_H

#include "Kernel.h"

class OneDEnergyWallHeatFlux;

template <>
InputParameters validParams<OneDEnergyWallHeatFlux>();

class OneDEnergyWallHeatFlux : public Kernel
{
public:
  OneDEnergyWallHeatFlux(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  /// Heat flux function
  const VariableValue & _q_wall;
  /// Heat flux perimeter
  const VariableValue & _P_hf;
};

#endif // ONEDENERGYWALLHEATFLUX_H

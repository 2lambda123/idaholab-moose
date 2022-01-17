#pragma once

#include "HeatConductionTimeDerivative.h"
#include "RZSymmetry.h"

/**
 * Time derivative kernel used by heat conduction equation in arbitrary RZ symmetry
 */
class HeatConductionTimeDerivativeRZ : public HeatConductionTimeDerivative, public RZSymmetry
{
public:
  HeatConductionTimeDerivativeRZ(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

public:
  static InputParameters validParams();
};

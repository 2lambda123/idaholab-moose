#pragma once

#include "ODETimeDerivative.h"

/**
 * Time derivative multiplied by a coefficient for ODEs
 */
class ODECoefTimeDerivative : public ODETimeDerivative
{
public:
  ODECoefTimeDerivative(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  /// Coefficient that the time derivative terms is multiplied with
  const Real & _coef;

public:
  static InputParameters validParams();
};

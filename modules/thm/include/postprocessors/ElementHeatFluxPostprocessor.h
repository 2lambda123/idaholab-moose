#pragma once

#include "ElementIntegralPostprocessor.h"

/**
 * Postprocessor to compute total heat flux going into the fluid
 *
 * This is used for debugging.
 */
class ElementHeatFluxPostprocessor : public ElementIntegralPostprocessor
{
public:
  ElementHeatFluxPostprocessor(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral();

  /// Wall temperature
  const MaterialProperty<Real> & _T_wall;
  const VariableValue & _Tfluid;
  /// convective heat transfer coefficient, W/m^2-K
  const MaterialProperty<Real> & _Hw;
  /// Heat flux perimeter
  const VariableValue & _P_hf;

public:
  static InputParameters validParams();
};

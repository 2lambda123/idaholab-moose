#pragma once

#include "OneDIntegratedBC.h"
#include "DerivativeMaterialInterfaceTHM.h"

/**
 *
 */
class OneDMassHRhoUBC : public DerivativeMaterialInterfaceTHM<OneDIntegratedBC>
{
public:
  OneDMassHRhoUBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  /// Volume fraction
  const MaterialProperty<Real> & _alpha;
  const MaterialProperty<Real> * const _dalpha_dbeta;
  /// Specified momentum
  const Real & _rhou;
  /// Cross-sectional area
  const VariableValue & _area;

  unsigned int _beta_var_num;

public:
  static InputParameters validParams();
};

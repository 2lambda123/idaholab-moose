#ifndef ONEDENERGYGRAVITY_H
#define ONEDENERGYGRAVITY_H

#include "Kernel.h"
#include "DerivativeMaterialInterfaceRelap.h"

class OneDEnergyGravity;

template <>
InputParameters validParams<OneDEnergyGravity>();

/**
 * Computes gravity term for the energy equation
 */
class OneDEnergyGravity : public DerivativeMaterialInterfaceRelap<Kernel>
{
public:
  OneDEnergyGravity(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  const bool _has_beta;

  const VariableValue & _A;

  const MaterialProperty<Real> & _alpha;
  const MaterialProperty<Real> * const _dalpha_dbeta;

  const MaterialProperty<Real> & _rho;
  const MaterialProperty<Real> * const _drho_dbeta;
  const MaterialProperty<Real> & _drho_darhoA;

  const MaterialProperty<Real> & _vel;
  const MaterialProperty<Real> & _dvel_darhoA;
  const MaterialProperty<Real> & _dvel_darhouA;

  /// x-component of gravity
  const VariableValue & _gx;

  const unsigned int _beta_var_number;
  const unsigned int _arhoA_var_number;
  const unsigned int _arhouA_var_number;
};

#endif /* ONEDENERGYGRAVITY_H */

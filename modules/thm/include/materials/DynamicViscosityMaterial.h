#pragma once

#include "Material.h"
#include "DerivativeMaterialInterfaceTHM.h"

class SinglePhaseFluidProperties;

/**
 * Computes dynamic viscosity
 */
class DynamicViscosityMaterial : public DerivativeMaterialInterfaceTHM<Material>
{
public:
  static InputParameters validParams();

  DynamicViscosityMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// Dynamic viscosity name
  const MaterialPropertyName & _mu_name;

  // Dynamic viscosity derivatives
  MaterialProperty<Real> & _mu;
  MaterialProperty<Real> * const _dmu_dbeta;
  MaterialProperty<Real> & _dmu_darhoA;
  MaterialProperty<Real> & _dmu_darhouA;
  MaterialProperty<Real> & _dmu_darhoEA;

  /// Specific volume
  const MaterialProperty<Real> & _v;
  const MaterialProperty<Real> * const _dv_dbeta;
  const MaterialProperty<Real> & _dv_darhoA;

  /// Specific internal energy
  const MaterialProperty<Real> & _e;
  const MaterialProperty<Real> & _de_darhoA;
  const MaterialProperty<Real> & _de_darhouA;
  const MaterialProperty<Real> & _de_darhoEA;

  /// Single-phase fluid properties
  const SinglePhaseFluidProperties & _fp_1phase;
};

#pragma once

#include "VolumeJunctionBase.h"

class VolumeJunction1Phase;

template <>
InputParameters validParams<VolumeJunctionBase>();

/**
 * Junction between 1-phase flow channels that has a non-zero volume
 */
class VolumeJunction1Phase : public VolumeJunctionBase
{
public:
  VolumeJunction1Phase(const InputParameters & params);

  virtual void addVariables() override;
  virtual void addMooseObjects() override;

  /// Enumeration for junction variable/equation indices
  enum VolumeJunction1PhaseIndices
  {
    RHOV_INDEX = 0,
    RHOUV_INDEX = 1,
    RHOVV_INDEX = 2,
    RHOWV_INDEX = 3,
    RHOEV_INDEX = 4
  };
  /// Number of equations for the junction
  static const unsigned int N_EQ;

protected:
  virtual void check() const override;

  /**
   * Builds user object for computing and storing the fluxes
   */
  virtual const UserObjectName buildVolumeJunctionUserObject();

  /// Scaling factor for rho*V
  const Real & _scaling_factor_rhoV;
  /// Scaling factor for rho*u*V
  const Real & _scaling_factor_rhouV;
  /// Scaling factor for rho*v*V
  const Real & _scaling_factor_rhovV;
  /// Scaling factor for rho*w*V
  const Real & _scaling_factor_rhowV;
  /// Scaling factor for rho*E*V
  const Real & _scaling_factor_rhoEV;

  /// rho*V variable name for junction
  const VariableName _rhoV_var_name;
  /// rho*u*V variable name for junction
  const VariableName _rhouV_var_name;
  /// rho*v*V variable name for junction
  const VariableName _rhovV_var_name;
  /// rho*w*V variable name for junction
  const VariableName _rhowV_var_name;
  /// rho*E*V variable name for junction
  const VariableName _rhoEV_var_name;
  /// pressure variable name for junction
  const VariableName _pressure_var_name;
  /// temperature variable name for junction
  const VariableName _temperature_var_name;
  /// velocity variable name for junction
  const VariableName _velocity_var_name;

  /// Form loss coefficient
  const Real & _K;
  /// Reference area
  const Real & _A_ref;
};

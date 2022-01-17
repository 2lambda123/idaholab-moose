#pragma once

#include "ScalarKernel.h"
#include "Assembly.h"

class VolumeJunctionBaseUserObject;

/**
 * Adds advective fluxes for the junction variables for a volume junction
 */
class VolumeJunctionAdvectionScalarKernel : public ScalarKernel
{
public:
  VolumeJunctionAdvectionScalarKernel(const InputParameters & params);

  virtual void reinit() override;
  virtual void computeResidual() override;
  virtual void computeJacobian() override;

protected:
  /// Index within local system of the equation upon which this object acts
  const unsigned int _equation_index;

  /// Volume junction user object
  const VolumeJunctionBaseUserObject & _volume_junction_uo;

public:
  static InputParameters validParams();
};

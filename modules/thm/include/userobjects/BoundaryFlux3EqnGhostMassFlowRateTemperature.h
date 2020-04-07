#pragma once

#include "BoundaryFlux3EqnGhostBase.h"

class SinglePhaseFluidProperties;

/**
 * Computes a boundary flux from a specified mass flow rate and temperature for
 * the 1-D, 1-phase, variable-area Euler equations using a ghost cell.
 */
class BoundaryFlux3EqnGhostMassFlowRateTemperature : public BoundaryFlux3EqnGhostBase
{
public:
  BoundaryFlux3EqnGhostMassFlowRateTemperature(const InputParameters & parameters);

protected:
  virtual std::vector<Real> getGhostCellSolution(const std::vector<Real> & U) const override;
  virtual DenseMatrix<Real>
  getGhostCellSolutionJacobian(const std::vector<Real> & U) const override;

  /// Specified mass flow rate
  const Real & _rhouA;

  /// Specified temperature
  const Real & _T;
  /// Reversible flag
  const bool & _reversible;

  /// Fluid properties object
  const SinglePhaseFluidProperties & _fp;

public:
  static InputParameters validParams();
};

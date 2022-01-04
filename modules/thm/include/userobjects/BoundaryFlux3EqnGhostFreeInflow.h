#pragma once

#include "BoundaryFlux3EqnGhostBase.h"
#include "SinglePhaseFluidProperties.h"

/**
 * Free inflow boundary conditions from a ghost cell for the 1-D, 1-phase, variable-area Euler
 * equations
 */
class BoundaryFlux3EqnGhostFreeInflow : public BoundaryFlux3EqnGhostBase
{
public:
  BoundaryFlux3EqnGhostFreeInflow(const InputParameters & parameters);

protected:
  virtual std::vector<Real> getGhostCellSolution(const std::vector<Real> & U1) const override;
  virtual DenseMatrix<Real>
  getGhostCellSolutionJacobian(const std::vector<Real> & U1) const override;

  /// farstream density
  const Real _rho_inf;
  /// farstream velocity
  const Real _vel_inf;
  /// farstream pressure
  const Real _p_inf;

  /// fluid properties object
  const SinglePhaseFluidProperties & _fp;

public:
  static InputParameters validParams();
};

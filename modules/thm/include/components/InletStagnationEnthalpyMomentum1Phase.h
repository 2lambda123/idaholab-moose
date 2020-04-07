#pragma once

#include "FlowBoundary.h"

/**
 * Boundary condition with prescribed stagnation enthalpy and momentum for 1-phase flow channels
 */
class InletStagnationEnthalpyMomentum1Phase : public FlowBoundary
{
public:
  InletStagnationEnthalpyMomentum1Phase(const InputParameters & params);

  virtual void addMooseObjects() override;

protected:
  virtual void check() const override;

  /// True to allow the flow to reverse, otherwise false
  bool _reversible;

public:
  static InputParameters validParams();
};

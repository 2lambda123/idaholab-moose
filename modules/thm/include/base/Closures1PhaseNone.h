#pragma once

#include "Closures1PhaseBase.h"

class Closures1PhaseNone;

template <>
InputParameters validParams<Closures1PhaseNone>();

/**
 * Sets up no 1-phase closures
 *
 * Users have to implement the closure materials and supply them in the input file
 */
class Closures1PhaseNone : public Closures1PhaseBase
{
public:
  Closures1PhaseNone(const InputParameters & params);

  virtual void check(const FlowChannelBase & flow_channel) const override;
  virtual void check(const HeatTransferBase & heat_transfer) const override;
  virtual void addMooseObjects(const FlowChannelBase & flow_channel) override;
  virtual void addMooseObjects(const HeatTransferBase & heat_transfer) override;
};

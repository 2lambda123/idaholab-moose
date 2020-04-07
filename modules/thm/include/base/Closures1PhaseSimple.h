#pragma once

#include "Closures1PhaseBase.h"

/**
 * Simple 1-phase closures
 */
class Closures1PhaseSimple : public Closures1PhaseBase
{
public:
  Closures1PhaseSimple(const InputParameters & params);

  virtual void check(const FlowChannelBase & flow_channel) const override;
  virtual void check(const HeatTransferBase & heat_transfer) const override;
  virtual void addMooseObjects(const FlowChannelBase & flow_channel) override;
  virtual void addMooseObjects(const HeatTransferBase & heat_transfer) override;

protected:
  /**
   * Adds material to compute wall temperature from heat flux
   *
   * @param[in] flow_channel   Flow channel component
   */
  void addWallTemperatureFromHeatFluxMaterial(const FlowChannel1Phase & flow_channel) const;

public:
  static InputParameters validParams();
};

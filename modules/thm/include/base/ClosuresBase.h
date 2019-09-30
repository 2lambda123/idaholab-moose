#pragma once

#include "MooseObject.h"
#include "LoggingInterface.h"
#include "NamingInterface.h"

class ClosuresBase;
class FlowChannelBase;
class HeatTransferBase;
class THMProblem;
class Factory;

template <>
InputParameters validParams<ClosuresBase>();

/**
 * Base class for closures implementations
 *
 * The responsibilities of the closures objects depends on the flow model that
 * uses them. Examples of responsibilities will be to provide material properties
 * for friction factors and heat transfer coefficients.
 */
class ClosuresBase : public MooseObject, public LoggingInterface, public NamingInterface
{
public:
  ClosuresBase(const InputParameters & params);

  /**
   * Checks for errors
   *
   * @param[in] flow_channel   Flow channel component
   */
  virtual void check(const FlowChannelBase & flow_channel) const = 0;
  virtual void check(const HeatTransferBase & heat_transfer) const = 0;

  /**
   * Adds MOOSE objects
   *
   * @param[in] flow_channel   Flow channel component
   */
  virtual void addMooseObjects(const FlowChannelBase & flow_channel) = 0;
  virtual void addMooseObjects(const HeatTransferBase & heat_transfer) = 0;

protected:
  /**
   * Adds an arbitrary zero-value material
   *
   * @param[in] flow_channel   Flow channel component
   * @param[in] property_name   Name of the material property to create
   */
  void addZeroMaterial(const FlowChannelBase & flow_channel,
                       const std::string & property_name) const;

  /**
   * Adds a weighted average material
   *
   * @param[in] flow_channel   Flow channel component
   * @param[in] values   Values to average
   * @param[in] weights   Weights for each value
   * @param[in] property_name   Name of material property to create
   */
  void addWeightedAverageMaterial(const FlowChannelBase & flow_channel,
                                  const std::vector<MaterialPropertyName> & values,
                                  const std::vector<VariableName> & weights,
                                  const MaterialPropertyName & property_name) const;

  /**
   * Adds a material for wall temperature from an aux variable
   *
   * @param[in] flow_channel   Flow channel component
   */
  void addWallTemperatureFromAuxMaterial(const FlowChannelBase & flow_channel) const;

  /// Simulation
  THMProblem & _sim;

  /// Factory associated with the MooseApp
  Factory & _factory;
};

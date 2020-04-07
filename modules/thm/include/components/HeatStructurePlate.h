#pragma once

#include "HeatStructureBase.h"
#include "HeatConductionModel.h"

/**
 * Component to model plate heat structure
 */
class HeatStructurePlate : public HeatStructureBase
{
public:
  HeatStructurePlate(const InputParameters & params);

  virtual void addMooseObjects() override;
  virtual Real getUnitPerimeter(const SideType & side) const override;

protected:
  /// plate fuel depth
  const Real & _depth;

public:
  static InputParameters validParams();
};

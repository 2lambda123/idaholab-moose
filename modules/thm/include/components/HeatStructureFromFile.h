#pragma once

#include "HeatStructureBase.h"
#include "HeatConductionModel.h"

/**
 * Heat structure component that can load the mesh from an ExodusII file
 */
class HeatStructureFromFile : public HeatStructureBase
{
public:
  HeatStructureFromFile(const InputParameters & params);

  virtual void buildMesh() override;
  virtual void addVariables() override;
  virtual void addMooseObjects() override;

  FunctionName getInitialT() const;

  virtual Real getUnitPerimeter(const HeatStructureSideType & side) const override;

protected:
  virtual std::shared_ptr<HeatConductionModel> buildModel() override;
  virtual void init() override;
  virtual void check() const override;
  virtual bool usingSecondOrderMesh() const override;

  /// The name of the ExodusII file to load the mesh from
  const FileName & _file_name;

public:
  static InputParameters validParams();
};

#pragma once

#include "StabilizationSettings.h"

class SUPG;

template <>
InputParameters validParams<SUPG>();

class SUPG : public StabilizationSettings
{
public:
  SUPG(const InputParameters & parameters);

  virtual void addVariables(FlowModel & fm, const SubdomainName & subdomain_name) const;
  virtual void initMooseObjects(FlowModel & fm);
  virtual void addMooseObjects(FlowModel & fm, InputParameters & pars) const;

public:
  static const std::string DELTA;
  static const std::string RESIDUAL;
  static const std::string MATRIX;
  static const std::string COLUMNS;
  static const std::string DUDX;
  static const std::string DADU;
};

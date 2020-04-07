#pragma once

#include "THMControl.h"

class Function;

class TimeFunctionControl : public THMControl
{
public:
  TimeFunctionControl(const InputParameters & parameters);

  virtual void execute();

protected:
  const std::string & _component_name;
  const std::string & _param_name;
  MooseObjectParameterName _ctrl_param_name;
  const Function & _function;

public:
  static InputParameters validParams();
};

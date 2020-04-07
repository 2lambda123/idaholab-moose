#pragma once

#include "TestAction.h"

/**
 * Base class for adding common actions for Jacobian tests
 */
class JacobianTestAction : public TestAction
{
public:
  JacobianTestAction(InputParameters params);

protected:
  virtual void addPreconditioner() override;

  /// Finite differencing parameter
  const std::string _snes_test_err;

public:
  static InputParameters validParams();
};

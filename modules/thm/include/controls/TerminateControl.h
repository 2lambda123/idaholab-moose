#pragma once

#include "THMControl.h"

/**
 * This control block will terminate a run if its input indicates so.
 */
class TerminateControl : public THMControl
{
public:
  TerminateControl(const InputParameters & parameters);

  virtual void execute();

protected:
  /// Flag to throw an error if the terminate condition is met
  const bool _throw_error;

  /// Message to use if termination occurs
  const std::string & _termination_message;

  /// The control data that indicates if the simulation should be terminated
  const bool & _terminate;

public:
  static InputParameters validParams();
};

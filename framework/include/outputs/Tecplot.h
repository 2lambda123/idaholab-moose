/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#ifndef TECPLOT_H
#define TECPLOT_H

// MOOSE includes
#include "OversampleOutputter.h"

// Forward declarations
class Tecplot;

template<>
InputParameters validParams<Tecplot>();

/**
 * Class for output data to the TecplotII format
 */
class Tecplot :
  public OversampleOutputter
{
public:

  /**
   * Class constructor
   */
  Tecplot(const std::string & name, InputParameters);

protected:

  /**
   * Overload the OutputBase::output method, this is required for Tecplot
   * output due to the method utilized for outputing single/global parameters
   */
  virtual void output();

  /**
   * Returns the current filename, this method handles adding the timestep suffix
   * @return A string containing the current filename to be written
   */
  std::string filename();

  //@{
  /**
   * Individual component output is not supported for Tecplot
   */
  virtual void outputNodalVariables();
  virtual void outputElementalVariables();
  virtual void outputPostprocessors();
  virtual void outputScalarVariables();
  //@}

private:

  /// Flag for binary output
  bool _binary;
};

#endif /* TECPLOT_H */

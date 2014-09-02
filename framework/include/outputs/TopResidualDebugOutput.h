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

#ifndef TOPRESIDUALEBUGOUTPUT_H
#define TOPRESIDUALEBUGOUTPUT_H

// MOOSE includes
#include "PetscOutput.h"
#include "FEProblem.h"

// Forward declerations
class TopResidualDebugOutput;

template<>
InputParameters validParams<TopResidualDebugOutput>();


/**
 * A structure for storing data related to top residuals
 *  @see TopResidualDebugOutput::printTopResiduals()
 */
struct TopResidualDebugOutputTopResidualData
{
  unsigned int _var;
  unsigned int _nd;
  Real _residual;

  TopResidualDebugOutputTopResidualData() { _var = 0; _nd = 0; _residual = 0.; }

  TopResidualDebugOutputTopResidualData(unsigned int var, unsigned int nd, Real residual)
  {
    _var = var;
    _nd = nd;
    _residual = residual;
  }
};

/**
 * A class for producing various debug related outputs
 *
 * This class may be used from inside the [Outputs] block or via the [Debug] block (preferred)
 */
class TopResidualDebugOutput : public PetscOutput
{
public:

  /**
   * Class constructor
   * @param name Output object name
   * @param parameters Object input parameters
   */
  TopResidualDebugOutput(const std::string & name, InputParameters & parameters);

  /**
   * Class destructor
   */
  virtual ~TopResidualDebugOutput();

protected:

  /**
   * Perform the debugging output
   */
  virtual void output();

  /**
   * Prints the n top residuals for the variables in the system
   * @param residual A reference to the residual vector
   * @param n The number of residuals to print
   */
  void printTopResiduals(const NumericVector<Number> & residual, unsigned int n);

  /**
   * Method for sorting the residuals data from TopResidualDebugOutputTopResidualData structs
   * @see printTopResiduals
   */
  static bool sortTopResidualData(TopResidualDebugOutputTopResidualData i, TopResidualDebugOutputTopResidualData j) { return (fabs(i._residual) > fabs(j._residual)); }

  //@{
  /**
   * Individual component output is not supported for TopResidualDebugOutput
   */
  std::string filename();
  virtual void outputNodalVariables();
  virtual void outputElementalVariables();
  virtual void outputPostprocessors();
  virtual void outputVectorPostprocessors();
  virtual void outputScalarVariables();
  //@}

  /// Number of residuals to display
  unsigned int _num_residuals;

  /// Reference to libMesh system
  TransientNonlinearImplicitSystem & _sys;

};

#endif // TOPRESIDUALDEBUGOUTPUT_H

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

#ifndef FUNCTIONDT_H
#define FUNCTIONDT_H

#include "TimeStepper.h"
#include "LinearInterpolation.h"

class FunctionDT;

template<>
InputParameters validParams<FunctionDT>();

class FunctionDT : public TimeStepper
{
public:
  FunctionDT(const std::string & name, InputParameters parameters);

  virtual void init();

  virtual void preExecute();

  virtual void acceptStep();
  virtual void rejectStep();

protected:
  virtual Real computeInitialDT();
  virtual Real computeDT();
  virtual bool constrainStep( Real & dt );

private:
  void removeOldKnots();

  std::set<Real> _time_t;
  /// Piecewise linear definition of time stepping
  LinearInterpolation _time_ipol;
  Real _growth_factor;
  Real _min_dt;
};

#endif /* FUNCTIONDT_H_ */

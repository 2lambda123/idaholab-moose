#ifndef ADAPTIVESTEPPER_H
#define ADAPTIVESTEPPER_H

#include "TimeStepper.h"

#include "LinearInterpolation.h"

class Function;
class PiecewiseLinearFile;

/**
 *
 */
class AdaptiveDT : public TimeStepper
{
public:
  AdaptiveDT(const std::string & name, InputParameters parameters);
  virtual ~AdaptiveDT();

  virtual void init();
  virtual void preExecute();

  virtual void rejectStep();

  virtual void acceptStep();

protected:
  virtual Real computeInitialDT();
  virtual Real computeDT();
  virtual bool constrainStep(Real &dt);
  virtual Real computeFailedDT();
  void computeAdaptiveDT(Real & dt, bool allowToGrow=true, bool allowToShrink=true);
  void computeInterpolationDT(Real & dt);
  void limitDTByFunction(Real & limitedDT);

  Real & _dt_old;

  const Real _input_dt;                       ///< The dt from the input file.
  bool _tfunc_last_step;
  bool _sync_last_step;

  int _optimal_iterations;
  int _iteration_window;
  const int _linear_iteration_ratio;
  bool _adaptive_timestepping;

  Function * _timestep_limiting_function;
  PiecewiseLinearFile * _piecewise_linear_timestep_limiting_function;
  std::vector<Real> _times;
  Real _max_function_change;
  bool _force_step_every_function_point;

  std::set<Real> _tfunc_times;

  LinearInterpolation _time_ipol; ///< Piecewise linear definition of time stepping
  const bool _use_time_ipol;      ///< true if we want to use piecewise-defined time stepping
  const Real _growth_factor;
  const Real _cutback_factor;

  unsigned int & _nl_its;  /// Number of nonlinear iterations in previous solve
  unsigned int & _l_its;   /// Number of linear iterations in previous solve
  bool & _cutback_occurred;

};

template<>
InputParameters validParams<AdaptiveDT>();

#endif /* ADAPTIVESTEPPER_H */

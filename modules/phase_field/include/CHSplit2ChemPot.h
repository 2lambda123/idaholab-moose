#ifndef CHSPLIT2CHEMPOT_H
#define CHSPLIT2CHEMPOT_H

#include "Kernel.h"


//Forward Declarations
class CHSplit2ChemPot;

template<>
InputParameters validParams<CHSplit2ChemPot>();

/**
 * This file calculates the chemical potential in the split Cahn-Hilliard problem.  To create a new phase field model, just
 * override computeDFDC
 */
class CHSplit2ChemPot : public Kernel
{
public:

  CHSplit2ChemPot(const std::string & name, InputParameters parameters);
  
protected:
  
  enum PFFunctionType
  {
    Residual,
    Jacobian
  };
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);
  virtual Real computeDFDC(PFFunctionType type);
  
  unsigned int _w_var;
  VariableValue & _w;
  VariableGradient & _grad_w;
  std::string _kappa_name;
  MaterialProperty<Real> & _kappa;

private:
};
#endif //CHSPLIT2CHEMPOT_H

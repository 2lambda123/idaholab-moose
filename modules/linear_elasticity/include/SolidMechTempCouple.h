#ifndef SOLIDMECHTEMPCOUPLE_H
#define SOLIDMECHTEMPCOUPLE_H

#include "SolidMech.h"

//libMesh includes
#include "libmesh/tensor_value.h"

//Forward Declarations
class SolidMechTempCouple;

template<>
InputParameters validParams<SolidMechTempCouple>();

class SolidMechTempCouple : public SolidMech
{
public:

  SolidMechTempCouple(const std::string & name, InputParameters parameters);
  
  virtual void subdomainSetup();

  void recomputeCouplingConstants();
  
protected:
  unsigned int _temp_var;
  
  MaterialProperty<Real> & _thermal_strain;
  MaterialProperty<Real> & _alpha;

  Real _c4;
};
 

#endif //SOLIDMECHTEMPCOUPLE_H

#ifndef PLC_LSH_H
#define PLC_LSH_H

#include "MaterialModel.h"

// Forward declarations
class PLC_LSH;

template<>
InputParameters validParams<PLC_LSH>();

/**
 * Power-law creep material
 * edot = A(sigma)**n * exp(-Q/(RT))
 */

class PLC_LSH : public MaterialModel
{
public:
  PLC_LSH( const std::string & name,
                 InputParameters parameters );

protected:

  Real _coefficient;
  Real _exponent;
  Real _activation_energy;
  Real _gas_constant;
  Real _tau;
  Real _tolerance;
  Real _yield_stress;
  Real _hardening_constant;

  unsigned int _max_its;
  bool _output_iteration_info;

  MaterialProperty<SymmTensor> & _creep_strain;
  MaterialProperty<SymmTensor> & _creep_strain_old;

  MaterialProperty<SymmTensor> & _plastic_strain;
  MaterialProperty<SymmTensor> & _plastic_strain_old;

  MaterialProperty<Real> & _hardening_variable;
  MaterialProperty<Real> & _hardening_variable_old;

  MaterialProperty<Real> & _del_p;


  /// Compute the stress (sigma += deltaSigma)
  virtual void computeStress();

  void computeCreep(SymmTensor & creep_strain_increment,
                    SymmTensor & stress_new );
  void computeLSH( const SymmTensor & creep_strain_increment,
                   SymmTensor & elastic_strain_increment,
                   SymmTensor & stress_new );

private:

};

#endif //PLC_LSHMATERIAL_H

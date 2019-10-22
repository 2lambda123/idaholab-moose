#pragma once

#include "MooseObjectUnitTest.h"
#include "NumericalFlux3EqnBase.h"
#include "IdealGasFluidProperties.h"

/**
 * Base class for testing rDG objects for the 3-equation model.
 */
class Test3EqnRDGObjectBase : public MooseObjectUnitTest
{
public:
  Test3EqnRDGObjectBase()
    : MooseObjectUnitTest("THMTestApp"),

      _fp_name("fp"),
      _fp(getFluidPropertiesObject()),

      _A(2.0),
      _nLR_dot_d(1.0)
  {
  }

protected:
  /**
   * Builds and gets the fluid properties user object
   */
  const SinglePhaseFluidProperties & getFluidPropertiesObject()
  {
    const std::string class_name = "IdealGasFluidProperties";
    InputParameters params = _factory.getValidParams(class_name);
    params.set<Real>("gamma") = 1.4;
    params.set<Real>("molar_mass") = 11.640243719999999;
    _fe_problem->addUserObject(class_name, _fp_name, params);
    return _fe_problem->getUserObjectTempl<SinglePhaseFluidProperties>(_fp_name);
  }

  /**
   * Computes the conservative solution from the primitive solution
   *
   * @param[in] W   Primitive solution vector: {p, T, vel}
   */
  std::vector<Real> computeConservativeSolution(const std::vector<Real> & W) const;

  /**
   * Creates the flux object to be tested
   */
  virtual const NumericalFlux3EqnBase * createFluxObject() = 0;

  /**
   * Runs the tests
   */
  virtual void test() = 0;

  /// Fluid properties user object name
  const UserObjectName _fp_name;
  /// Fluid properties user object
  const SinglePhaseFluidProperties & _fp;

  /// Cross-sectional area
  const Real _A;
  /// Dot product of normal from L to R with direction vector
  const Real _nLR_dot_d;

  /// Flux object to be tested
  const NumericalFlux3EqnBase * _flux;
};

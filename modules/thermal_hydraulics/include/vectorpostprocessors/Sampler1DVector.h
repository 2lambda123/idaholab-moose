#pragma once

#include "Sampler1DBase.h"

/**
 * This class samples a component of a vector material property in a 1-D mesh.
 */
class Sampler1DVector : public Sampler1DBase<std::vector<Real>>
{
public:
  Sampler1DVector(const InputParameters & parameters);

  virtual Real getScalarFromProperty(const std::vector<Real> & property,
                                     const Point & curr_point) override;

protected:
  /// index of the component of the vector-valued material property
  const unsigned int _index;

public:
  static InputParameters validParams();
};

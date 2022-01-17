#pragma once

#include "NodalBC.h"

/**
 * Base class for nodal boundary conditions for 1D problems in 3D space
 */
class OneDNodalBC : public NodalBC
{
public:
  OneDNodalBC(const InputParameters & parameters);

protected:
  /// Component of outward normal along 1-D direction
  const Real _normal;

public:
  static InputParameters validParams();
};

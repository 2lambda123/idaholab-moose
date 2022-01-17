[JacobianTest1PhaseRDG]
  add_bc = true
  boundary_flux = boundary_flux
  ic_option = constant
  snes_test_err = 1e-8
[]

[UserObjects]
  [numerical_flux]
    type = NumericalFlux3EqnCentered
    fluid_properties = fluid_properties
    execute_on = 'linear nonlinear'
  []
  [boundary_flux]
    type = BoundaryFlux3EqnGhostStagnationPressureTemperature
    p0 = 1
    T0 = 2
    fluid_properties = fluid_properties
    numerical_flux = numerical_flux
    normal = 1
    execute_on = 'linear nonlinear'
  []
[]

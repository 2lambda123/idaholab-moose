#
# Testing the ability to discretize the Pipe by dividing it into
# subsections
#

[GlobalParams]
  initial_p = 1e5
  initial_T = 300
  initial_vel = 0

  closures = simple
[]

[FluidProperties]
  [./eos]
    type = StiffenedGasFluidProperties
    gamma = 2.35
    cv = 1816.0
    q = -1.167e6
    p_inf = 1.0e9
    q_prime = 0
  [../]
[]

[Components]
  [./pipe]
    type = FlowChannel1Phase
    position = '0 0 0'
    orientation = '1 0 0'
    axial_region_names = 'r1 r2'
    length = '1 1'
    n_elems = '1 2'
    A = 1
    f = 0
    fp = eos
  [../]

  [./inlet]
    type = SolidWall
    input = 'pipe:in'
  [../]

  [./outlet]
    type = SolidWall
    input = 'pipe:out'
  [../]
[]

[Problem]
  solve = false
[]

[Executioner]
  type = Transient
  num_steps = 0
[]

[Outputs]
  [./out]
    type = Exodus
    show = 'A'
  [../]
[]

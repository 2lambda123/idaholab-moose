[GlobalParams]
  gravity_vector = '0 0 0'

  initial_T = 444.447
  initial_p = 7e6
  initial_vel = 0

  scaling_factor_1phase = '1 1 1e-5'

  closures = simple

  rdg_slope_reconstruction = minmod
[]

[FluidProperties]
  [./fp]
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
    fp = fp
    # geometry
    position = '0 0 0'
    orientation = '1 0 0'
    A   = 1.0000000000e-04
    D_h  = 1.1283791671e-02
    f = 0.1
    length = 1
    n_elems = 100
  [../]

  [./inlet]
    type = InletMassFlowRateTemperature1Phase
    input = 'pipe:in'
    m_dot = 0.18
    T     = 444.447
  [../]

  [./outlet]
    type = Outlet1Phase
    input = 'pipe:out'
    p = 7e6
  [../]
[]

[Preconditioning]
  [./SMP_PJFNK]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  scheme = 'bdf2'

  dt = 0.1
  start_time = 0.0
  end_time = 5.5

  solve_type = 'PJFNK'
  line_search = 'basic'
  nl_rel_tol = 0
  nl_abs_tol = 1e-7
  nl_max_its = 20

  l_tol = 1e-3
  l_max_its = 100

  abort_on_solve_fail = true

  [./Quadrature]
    type = GAUSS
    order = SECOND
  [../]
[]

[Outputs]
  file_base = 'phy.massflowrate_3eqn'

  [./exodus]
    type = Exodus
    show = 'rhouA T p'
  [../]
[]

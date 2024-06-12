
# ############################### Tricrystal - GAMMA MODEL - GAMMAanisoGAUSS ##########################

[Mesh]

  [Base_mesh]
  #  Build a square mesh (nm)
    type = GeneratedMeshGenerator
    dim = 2
    nx = 80
    ny = 40
    nz = 1
    xmin = 0
    xmax = 800
    ymin = 0
    ymax = 400
    zmin = 0
    zmax = 1
    elem_type = QUAD4
  []

[]


[GlobalParams]
  op_num = 3 # Number of order parameters use
  var_name_base = gr # Base name of grains
  order = CONSTANT
  family = MONOMIAL
[]


[Variables]
  # Variable block, where all variables in the simulation are declared

  [PolycrystalVariables]
    order = FIRST
    family = LAGRANGE
  []

[]


[UserObjects]

  [grain_tracker]
    type = GrainTracker
    variable = 'gr0 gr1 gr2'
    threshold = 0.2
    connecting_threshold = 0.08
    compute_var_to_feature_map = true
  []

[]


[ICs]

  [InitialCondition_gr0]
     type = TricrystalTripleJunctionIC
     op_index = 1
     variable = gr0
     theta1 = 140
     junction = '165 200 0'
  []

  [InitialCondition_gr1]
     type = TricrystalTripleJunctionIC
     op_index = 2
     variable = gr1
     theta2 =  106.5
     junction = '165 200 0'
  []

  [InitialCondition_gr2]
     type = TricrystalTripleJunctionIC
     op_index = 3
     variable = gr2
     theta1 = 140
     theta2 =  106.5
     junction = '165 200 0'
  []

[]


[BCs]
  #  Default Neumann

  # [Periodic]
  #   [top_bottom]
  #     auto_direction = 'y z' # Makes problem periodic in the x and y directions
  #   []
  # []

[]


[AuxVariables]

  [bnds]
    order = FIRST
    family = LAGRANGE
  []

  [unique_grains]
  []

  [var_indices]
  []

  [ghost_regions]
  []

  [proc]
  []

  [bounds_dummy0]
    order = FIRST
    family = LAGRANGE
  []

  [bounds_dummy1]
  order = FIRST
  family = LAGRANGE
  []

  [bounds_dummy2]
    order = FIRST
    family = LAGRANGE
  []

  [free_energy]
    order = CONSTANT
    family = MONOMIAL
  []

[]


[Bounds]

  [d_upper_bound0]
    type = ConstantBounds
    variable = bounds_dummy0
    bounded_variable = gr0
    bound_type = upper
    bound_value = 1.0
  []

  [d_lower_bound0]
    type = ConstantBounds
    variable = bounds_dummy0
    bounded_variable = gr0
    bound_type = lower
    bound_value = 0.0
  []

  [d_upper_bound1]
    type = ConstantBounds
    variable = bounds_dummy1
    bounded_variable = gr1
    bound_type = upper
    bound_value = 1.0
  []

  [d_lower_bound1]
    type = ConstantBounds
    variable = bounds_dummy1
    bounded_variable = gr1
    bound_type = lower
    bound_value = 0.0
  []

  [d_upper_bound2]
    type = ConstantBounds
    variable = bounds_dummy2
    bounded_variable = gr2
    bound_type = upper
    bound_value = 1.0
  []

  [d_lower_bound2]
    type = ConstantBounds
    variable = bounds_dummy2
    bounded_variable = gr2
    bound_type = lower
    bound_value = 0.0
  []

[]


[Kernels]

  [GAMMAmodelKernel]
  []

[]


[AuxKernels]

  [bnds_aux]
    type = BndsCalcAux
    variable = bnds
    execute_on = timestep_end
  []

  [unique_grains]
    type = FeatureFloodCountAux
    variable = unique_grains
    flood_counter = grain_tracker
    field_display = UNIQUE_REGION
    execute_on = 'initial timestep_end'
  []

  [var_indices]
    type = FeatureFloodCountAux
    variable = var_indices
    flood_counter = grain_tracker
    field_display = VARIABLE_COLORING
    execute_on = 'initial timestep_end'
  []

  [ghosted_entities]
    type = FeatureFloodCountAux
    variable = ghost_regions
    flood_counter = grain_tracker
    field_display = GHOSTED_ENTITIES
    execute_on = 'initial timestep_end'
  []

  [proc]
    type = ProcessorIDAux
    variable = proc
    execute_on = 'initial timestep_end'
  []

  [energy_density]
    #  Total free energy for the bicrystal
    type = TotalFreeEnergy
    variable = free_energy
    kappa_names = "kappa_op kappa_op kappa_op"
    interfacial_vars = "gr0 gr1 gr2"
  []

[]


[Materials]

  [5DGaussian]
    # Material properties adding anisotropy to gamma;
    # L is isotropic and computed from a given mobility.

    # Minima Library file: Normalized_Axis_X  Normalized_Axis_Y Normalized_Axis_Z  Angle(Radian) GB_Normal_X  GB_Normal_Y  GB_Normal_Z  Minimun_GB_energy(J/m^2)
    # Grains orienttation quaternions file: qw qx qy qz


    type = GAMMAanisoGAUSS
    outputs = exodus

    grain_tracker = grain_tracker

    Ggamma_name = Ggamma
    Ggamma_EN_name = Ggamma_EN
    dGgammadx_name = dGgammadx
    dGgammadx_EN_name = dGgammadx_EN
    dGgammady_name = dGgammady
    dGgammady_EN_name = dGgammady_EN
    dGgammadz_name = dGgammadz
    dGgammadz_EN_name = dGgammadz_EN

    dGgammadxplus_name = dGgammadxplus
    dGgammadxplus_EN_name = dGgammadxplus_EN
    dGgammadyplus_name = dGgammadyplus
    dGgammadyplus_EN_name = dGgammadyplus_EN
    dGgammadzplus_name = dGgammadzplus
    dGgammadzplus_EN_name = dGgammadzplus_EN

    gamma_name = gamma
    gamma_EN_name = gamma_asymm

    S_switch_name = switch

    dgammadx_name = dgammadx
    dgammadx_EN_name = dgammadx_EN
    dgammady_name = dgammady
    dgammady_EN_name = dgammady_EN
    dgammadz_name = dgammadz
    dgammadz_EN_name = dgammadz_EN

    dgammadxplus_name = dgammadxplus
    dgammadxplus_EN_name = dgammadxplus_EN
    dgammadyplus_name = dgammadyplus
    dgammadyplus_EN_name = dgammadyplus_EN
    dgammadzplus_name = dgammadzplus
    dgammadzplus_EN_name = dgammadzplus_EN

    m_name = m
    m_EN_name = mu

    kappa_name = kappa

    kappa_EN_name = kappa_op

    L_name = L_AD
    L_EN_name = L

    MGBVALUE_name = MGBVALUE

    lgb_name = lgb
    lgb_EN_name = lgb_EN

    sigma_name = sigma
    sigma_EN_name = sigma_EN
    sigmaORIUNIT_name = sigmaORIUNIT
    sigmaORIUNIT_EN_name = sigmaORIUNIT_EN

    qwg_name = qwg
    qxg_name = qxg
    qyg_name = qyg
    qzg_name = qzg

    Ggammabar_name = Ggammabar
    Ggammamin2grains_name = Ggammamin2grains

    TotGauss_name = TotGauss

    VwlibR_name = VwlibR
    VxlibR_name = VxlibR
    VylibR_name = VylibR
    VzlibR_name = VzlibR
    Vx_name = Vx
    Vy_name = Vy
    Vz_name = Vz

    Library_file_name = "Minima_Library_TRI_GAMMAGAUSS" # Name given to the minima libray file

    Quaternion_file_name = "3grains_Quaternions" # Name given to the grains orientation quaternions file

    gammaBASE = 1.500000 # Base value of gamma from which gaussian is substracted or added
    GgammaBASE =  0.471404 # Value of g(gamma) corresponding to gammaBASE
    f0gammaBASE =  0.124961 # Value of function f corresponding to gammaBASE

    sigmaBASE = 1.5 # (J/m^2); Base value of energy from  which gaussian is substracted or added

    alphaswitch = 3 # Range of variation for switch 1

    betaswitch = 3 # Range of variation for switch 2

    libnum = 3 # Number of lines in minima libray file

    op_num = 3 # Number of lines in the grains orientation quaternions file

    lgbBASE_minimum = 20 # (nm);  Grain boundary width at Base energy

    sharpness = 30 # Width of gaussian

    amplitudeScale = 1 # MAINTAIN AT 1; Value to sacle the units.

    Gaussian_Tolerance = 1e-30 # Used to control gaussian effect; 1e-30 is enough

    ADDGaussian = false # 'true' to add to base value; 'false' to substract from base value

    BoundaryNormal = 5 # If = 0,the boundary normal in the simulation reference frame is computed using the orientations quaternions, else it is taken directly from the minima library file

    # Commpute L from Temperature-dependant mobility
    Mob = 2.5e-6 # (m^4/J.s)
    Q = 0.23 # (eV)
    T = 300 # (k)

    # Alternate Commpute L from directly given mobility
    #GBMobility =  1 # (m^4/(J.s))
  []

  [free_energy]
  # Free energy defined for thre grains variables 'gr0 gr1 gr2'
    type = DerivativeParsedMaterial
    property_name = F
    material_property_names = 'mu gamma_asymm'
    constant_names =       'phi  omega '
    constant_expressions = ' 1    1    '
    coupled_variables = 'gr0 gr1 gr2'
    expression = 'mu * ( (phi * (((gr0^4)/4) + ((gr1^4)/4) + ((gr2^4)/4)) )  - (omega * (((gr0^2)/2) + ((gr1^2)/2) + ((gr2^2)/2)))
                +  ( gamma_asymm * ( ((gr0^2)*(gr1^2))+((gr0^2)*(gr2^2))+((gr1^2)*(gr2^2)) ) ) + (1.0/4.0) )'
  []

[]


[Postprocessors]

  [dt]
    # Outputs the current time step
    type = TimestepSize
  []

  [total_energy]
    # Total free energy (eV)
    type = ElementIntegralVariablePostprocessor
    variable = free_energy
    execute_on = 'initial timestep_end'
  []

[]


[Preconditioning]

  [SMP]
    type = SMP
    full = true
  []

[]


[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart  -snes_type'
  petsc_options_value = ' hypre    boomeramg      31                  vinewtonrsls'

  nl_abs_tol = 1e-5
  nl_rel_tol = 1e-5
  l_max_its = 30

  end_time = 400000 # (ns)

  [TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 6
    iteration_window = 2
    dt = 1 # (ns)
    growth_factor = 1.1
    cutback_factor = 0.75
  []

[]

[Debug]
  show_var_residual_norms = true
[]

[Outputs]
  exodus = true
[]

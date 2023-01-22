[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
  temperature = temp
  order = SECOND
  family = LAGRANGE
[]

[Problem]
  type = ReferenceResidualProblem
  extra_tag_vectors = 'ref'
  reference_vector = 'ref'
[]

[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 1
  ny = 1
  nz = 1
  xmin = 0.0
  xmax = 4.0
  ymin = 0.0
  ymax = 0.5
  zmin = 0.0
  zmax = 1.0
  elem_type = HEX27
[]

[Variables]
  [temp]
    initial_condition = 600.0
  []
[]

[Kernels]
  [heat]
    type = HeatConduction
    variable = temp
  []
  [heat_ie]
    type = HeatConductionTimeDerivative
    variable = temp
  []
[]

[Modules/TensorMechanics/Master]
  [all]
    strain = FINITE
    add_variables = true
    generate_output = 'stress_xx stress_xy stress_yy strain_xx strain_xy strain_yy creep_strain_xx creep_strain_yy creep_strain_zz'
    extra_vector_tags = 'ref'
  []
[]

[Functions]
  [loading_func]
    type = PiecewiseLinear
    x = '0.  1.   5.'
    y = '0. 50.0 100.0'
  []
[]

[BCs]
  [free_end_moment]
    type = Pressure
    variable = disp_y
    component = 0
    boundary = right
    factor = 1
    function = loading_func
  []
  [FixedCenterLineX]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0.0
  []
  [FixedCenterLineY]
    type = DirichletBC
    variable = disp_y
    boundary = left
    value = 0.0
  []
  [FixedCenterLineZ]
    type = DirichletBC
    variable = disp_z
    boundary = left
    value = 0.0
  []
  [temp_fix]
    type = DirichletBC
    variable = temp
    boundary = 'left right bottom top front back'
    value = 600.0
  []
[]

[Materials]
  [elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 1e6
    poissons_ratio = 0.31
  []
  [stress]
    type = ComputeMultipleInelasticStress
    inelastic_models = 'isoplas powerlawcrp'
  []
  [isoplas]
    type = IsotropicPlasticityStressUpdate
    yield_stress = 2.5e3
    hardening_constant = 1.85e5
  []
  [powerlawcrp]
    type = PowerLawCreepStressUpdate
#    coefficient = 3.125e-14
#    n_exponent = 5.0
    coefficient = 5.0e-15
    n_exponent = 3.0
    m_exponent = 0.0
    activation_energy = 0.0
  []
  [thermal]
    type = HeatConductionMaterial
    specific_heat = 1.0
    thermal_conductivity = 100.
  []
  [density]
    type = Density
    density = 1.0
  []
[]

[Preconditioning]
  [pc]
    type = SMP
    full = True
  []
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu     superlu_dist'

  l_max_its = 20
  l_tol = 1e-3

  nl_max_its = 50
  nl_rel_tol = 1e-4
  nl_abs_tol = 1e-8

  start_time = 0.0
  dt = 1
  end_time = 5
[]

[Postprocessors]
  [num_lin_it]
    type = NumLinearIterations
  []
  [num_nonlin_it]
    type = NumNonlinearIterations
  []
  [tot_lin_it]
    type = CumulativeValuePostprocessor
    postprocessor = num_lin_it
  []
  [tot_nonlin_it]
    type = CumulativeValuePostprocessor
    postprocessor = num_nonlin_it
  []
  [alive_time]
    type = PerfGraphData
    section_name = Root
    data_type = TOTAL
  []
  [max_beam_deflection]
    type = NodalExtremeValue
    variable = disp_y
    boundary = 'right'
  []
[]

[Outputs]
  exodus = true
  csv = true
[]

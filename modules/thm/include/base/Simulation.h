#pragma once

#include "THMApp.h"
#include "FlowModel.h"
#include "Logger.h"
#include "ControlData.h"
#include "LoggingInterface.h"
#include "NamingInterface.h"

class ActionWarehouse;
class Component;
class THMMesh;
class THMProblem;

/**
 * Main class for simulation (the driver of the simulation)
 */
class Simulation : public LoggingInterface, public NamingInterface
{
public:
  Simulation(FEProblemBase & fe_problem, const InputParameters & params);
  virtual ~Simulation();

  /**
   * Gets the FE type for the flow in this simulation
   */
  const FEType & getFlowFEType() const { return _flow_fe_type; }

  /**
   * Sets up quadrature rules
   */
  virtual void setupQuadrature();

  /**
   * Initialize this simulation
   */
  virtual void initSimulation();

  /**
   * Initialize this simulation's components
   */
  virtual void initComponents();

  /**
   * Identifies the component loops
   */
  void identifyLoops();

  /**
   * Prints the component loops
   */
  void printComponentLoops() const;

  /**
   * Run the simulation
   */
  virtual void run();

  /**
   * Add a component into this simulation
   * @param type Type (the registered class name) of the component
   * @param name Name of the component
   * @param params Input parameters
   */
  void addComponent(const std::string & type, const std::string & name, InputParameters params);

  /**
   * Find out if simulation has a component with the given name
   * @param name The name of the component
   * @return true if the components exists, otherwise false
   */
  bool hasComponent(const std::string & name) const;

  /**
   * Find out if simulation has a component with the given name and specified type
   * @tparam T the type of the component we are requesting
   * @param name The name of the component
   * @return true if the component exists and has specified type, otherwise false
   */
  template <typename T>
  bool hasComponentOfType(const std::string & name) const;

  /**
   * Get component by its name
   * @tparam T the type of the component we are requesting
   * @param name The name of the component
   * @return Pointer to the component if found, otherwise throws and error
   */
  template <typename T>
  const T & getComponentByName(const std::string & name) const;

  /**
   * Called by a component to announce a variable
   * @param nl True is nonlinear variable is being added
   * @param name The name of the variable
   * @param type Type of the variable
   * @param subdomain_id Subdomain of the variable
   * @param scaling_factor Scaling factor for the variable
   */
  void addSimVariable(bool nl, const VariableName & name, FEType type, Real scaling_factor = 1.);
  void addSimVariable(bool nl,
                      const VariableName & name,
                      FEType type,
                      const std::vector<SubdomainName> & subdomain_names,
                      Real scaling_factor = 1.);

  void addConstantIC(const VariableName & var_name,
                     Real value,
                     const std::vector<SubdomainName> & block_names);
  void addFunctionIC(const VariableName & var_name,
                     const std::string & func_name,
                     const std::vector<SubdomainName> & block_names);
  void addConstantScalarIC(const VariableName & var_name, Real value);
  void addComponentScalarIC(const VariableName & var_name, const std::vector<Real> & value);

  void addSimInitialCondition(const std::string & type,
                              const std::string & name,
                              InputParameters params);

  /**
   * Add a control
   * @param type Type (registered name) of the control
   * @param name Name of the control
   * @param params Input parameters
   */
  void addControl(const std::string & type, const std::string & name, InputParameters params);

  void addFileOutputter(const std::string & name);
  void addScreenOutputter(const std::string & name);

  /**
   * Gets the vector of output names corresponding to a 1-word key string
   *
   * @param[in] key  string key that corresponds to an output names vector
   * @returns output names vector corresponding to key
   */
  std::vector<OutputName> getOutputsVector(const std::string & key) const;

  /**
   * Create mesh for this simulation
   */
  void buildMesh();

  /**
   * Add variables involved in this simulation
   */
  void addVariables();

  /**
   * Add components based physics
   */
  void addComponentPhysics();

  /**
   * Perform mesh setup actions such as setting up the coordinate system(s) and
   * creating ghosted elements.
   */
  void setupMesh();

  /**
   * Get the THMApp
   */
  THMApp & getApp() { return _app; }

  /**
   * Check the integrity of the simulation
   */
  virtual void integrityCheck() const;

  /**
   * Advance all of the state holding vectors / datastructures so that we can move to the next
   * timestep.
   */
  virtual void advanceState();

  /**
   * Check the integrity of the control data
   */
  virtual void controlDataIntegrityCheck();

  /**
   * Check integrity of coupling matrix used by the preconditioner
   */
  virtual void couplingMatrixIntegrityCheck() const;

  /**
   * Query if control data with name 'name' exists
   *
   * @param name The unique name of the control data
   * @return true if control data 'name' exists, false otherwise
   */
  template <typename T>
  bool hasControlData(const std::string & name)
  {
    if (_control_data.find(name) == _control_data.end())
      return false;
    else
      return dynamic_cast<ControlData<T> *>(_control_data[name]) != NULL;
  }

  /**
   * Get control data of type T and name 'name', if it does not exist it will be created
   *
   * @param name The unique name of the control data
   * @param name The control object that declared this data
   * @return Pointer to the control data of type T
   */
  template <typename T>
  ControlData<T> * getControlData(const std::string & name)
  {
    ControlData<T> * data = nullptr;
    if (_control_data.find(name) == _control_data.end())
    {
      data = new ControlData<T>(name);
      _control_data[name] = data;
    }
    else
      data = dynamic_cast<ControlData<T> *>(_control_data[name]);

    return data;
  }

  /**
   * Declare control data of type T and name 'name', if it does not exist it will be created
   *
   * @param name The unique name of the control data
   * @return Pointer to the control data of type T
   */
  template <typename T>
  ControlData<T> * declareControlData(const std::string & name, THMControl * ctrl)
  {
    ControlData<T> * data = getControlData<T>(name);
    if (!data->getDeclared())
    {
      // Mark the data for error checking
      data->setDeclared();
      data->setControl(ctrl);
    }
    else
      logError("Trying to declare '", name, "', but it was already declared.");

    return data;
  }

  /**
   * Gets the flag indicating whether an implicit time integration scheme is being used
   */
  const bool & getImplicitTimeIntegrationFlag() { return _implicit_time_integration; }

  /**
   * Are initial conditions specified from a file
   *
   * @return true if initial conditions are specified from a file
   */
  bool hasInitialConditionsFromFile();

  Logger & log() { return _log; }

  /**
   * Enable Jacobian checking
   *
   * @param state True for Jacobian checking, otherwise false
   */
  void setCheckJacobian(bool state) { _check_jacobian = state; }

  /**
   * Hint how to augment sparsity pattern between two elements.
   *
   * The augmentation will be symmetric
   */
  virtual void augmentSparsity(const dof_id_type & elem_id1, const dof_id_type & elem_id2);

protected:
  struct VariableInfo
  {
    bool _nl; ///< true if the variable is non-linear
    FEType _type;
    std::set<SubdomainName> _subdomain;
    Real _scaling_factor;
  };
  THMMesh & _mesh;

  /// Pointer to FEProblem representing this simulation
  FEProblemBase & _fe_problem;

  /// The application this is associated with
  THMApp & _app;

  /// The Factory associated with the MooseApp
  Factory & _factory;

  /// List of components in this simulation
  std::vector<std::shared_ptr<Component>> _components;
  /// Map of components by their names
  std::map<std::string, std::shared_ptr<Component>> _comp_by_name;
  /// Map of component name to component loop name
  std::map<std::string, std::string> _component_name_to_loop_name;
  /// Map of loop name to model type
  std::map<std::string, THM::FlowModelID> _loop_name_to_model_id;

  /// variables for this simulation (name and info about the var)
  std::map<VariableName, VariableInfo> _vars;

  struct ICInfo
  {
    std::string _type;
    InputParameters _params;

    ICInfo() : _params(emptyInputParameters()) {}
    ICInfo(const std::string & type, const InputParameters & params) : _type(type), _params(params)
    {
    }
  };
  std::map<std::string, ICInfo> _ics;

  /// "Global" of this simulation
  const InputParameters & _pars;

  /// finite element type for the flow in the simulation
  FEType _flow_fe_type;

  /**
   * Setup equations to be solved in this simulation
   */
  void setupEquations();

  /**
   * Setup reading initial conditions from a specified file, see 'initial_from_file' and
   * 'initial_from_file_timestep' parameters
   */
  void setupInitialConditionsFromFile();

  void setupInitialConditions();

  /**
   * Sets the coordinate system for each subdomain
   */
  void setupCoordinateSystem();

  /**
   * Setup ctirical heat flux table user object
   */
  void setupCriticalHeatFluxTable();

  std::vector<OutputName> _outputters_all;
  std::vector<OutputName> _outputters_file;
  std::vector<OutputName> _outputters_screen;

  /// Control data created in the control logic system
  std::map<std::string, ControlDataValue *> _control_data;

  /// true if using implicit time integration scheme
  bool _implicit_time_integration;

  Logger _log;

  /// True if checking jacobian
  bool _check_jacobian;

  /// Additional sparsity pattern that needs to be added into the Jacobian matrix
  std::map<dof_id_type, std::vector<dof_id_type>> _sparsity_elem_augmentation;

public:
  Real _zero;
};

template <typename T>
bool
Simulation::hasComponentOfType(const std::string & name) const
{
  auto it = _comp_by_name.find(name);
  if (it != _comp_by_name.end())
    return dynamic_cast<T *>((it->second).get()) != nullptr;
  else
    return false;
}

template <typename T>
const T &
Simulation::getComponentByName(const std::string & name) const
{
  auto it = _comp_by_name.find(name);
  if (it != _comp_by_name.end())
    return *dynamic_cast<T *>((it->second).get());
  else
    mooseError("Component '",
               name,
               "' does not exist in the simulation. Use hasComponent or "
               "checkComponnetByName before calling getComponent.");
}

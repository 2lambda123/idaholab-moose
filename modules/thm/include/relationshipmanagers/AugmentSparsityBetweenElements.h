#pragma once

#include "RelationshipManager.h"

using libMesh::processor_id_type;

class AugmentSparsityBetweenElements : public RelationshipManager
{
public:
  AugmentSparsityBetweenElements(const InputParameters &);

  /**
   * This function must be overriden by application codes to add
   * required elements from (range_begin, range_end) to the
   * coupled_elements map.
   */
  virtual void operator()(const MeshBase::const_element_iterator & range_begin,
                          const MeshBase::const_element_iterator & range_end,
                          processor_id_type p,
                          map_type & coupled_elements) override;

  /**
   * According to the base class docs, "We call mesh_reinit() whenever
   * the relevant Mesh has changed, but before remote elements on a
   * distributed mesh are deleted."
   */
  virtual void mesh_reinit() override;

  /**
   * Update the cached _lower_to_upper map whenever our Mesh has been
   * redistributed.  We'll be lazy and just recalculate from scratch.
   */
  virtual void redistribute() override { this->mesh_reinit(); }

  std::string getInfo() const override;

  virtual bool operator==(const RelationshipManager & other) const override;

protected:
  virtual void internalInit() override;

  const std::map<dof_id_type, std::vector<dof_id_type>> & _elem_map;

public:
  static InputParameters validParams();
};

#ifndef SCREAM_GRIDS_MANAGER_HPP
#define SCREAM_GRIDS_MANAGER_HPP

#include "share/grid/abstract_grid.hpp"
#include "share/grid/remap/abstract_remapper.hpp"
#include "share/grid/remap/identity_remapper.hpp"

#include "ekat/util/ekat_factory.hpp"
#include "ekat/util/ekat_string_utils.hpp"
#include "ekat/ekat_parameter_list.hpp"
#include "ekat/ekat_assert.hpp"
#include "ekat/mpi/ekat_comm.hpp"

#include <map>
#include <set>
#include <memory>

namespace scream
{

class GridsManager
{
public:
  using grid_type         = AbstractGrid;
  using grid_ptr_type     = std::shared_ptr<const grid_type>;
  using grid_repo_type    = std::map<std::string, grid_ptr_type>;
  using remapper_type     = AbstractRemapper;
  using remapper_ptr_type = std::shared_ptr<remapper_type>;

  GridsManager () = default;
  virtual ~GridsManager () = default;

  virtual std::string name () const = 0;

  grid_ptr_type get_grid (const std::string& name) const;

  grid_ptr_type get_reference_grid () const {
    return get_grid(get_reference_grid_name());
  }

  // Check if the given grid has been built
  bool has_grid (const std::string& grid_name) const;

  virtual void build_grids () = 0;

  remapper_ptr_type
  create_remapper (const grid_ptr_type& from_grid,
                   const grid_ptr_type& to_grid) const;

  remapper_ptr_type
  create_remapper (const std::string& from_grid,
                   const std::string& to_grid) const {
    return create_remapper(get_grid(from_grid),get_grid(to_grid));
  }

  remapper_ptr_type
  create_remapper_from_ref_grid(const grid_ptr_type& grid) const {
    return create_remapper(get_reference_grid(),grid);
  }

  remapper_ptr_type
  create_remapper_to_ref_grid(const grid_ptr_type& grid) const {
    return create_remapper(grid,get_reference_grid());
  }

  const grid_repo_type& get_repo () const { return m_grids; }

protected:

  virtual std::string get_reference_grid_name  () const = 0;

  virtual remapper_ptr_type
  do_create_remapper (const grid_ptr_type from_grid,
                      const grid_ptr_type to_grid) const = 0;

  std::string print_available_grids () const;


  grid_repo_type            m_grids;
};

// A short name for the factory for grid managers
using GridsManagerFactory 
    = ekat::Factory<GridsManager,
                    ekat::CaseInsensitiveString,
                    std::shared_ptr<GridsManager>,
                    const ekat::Comm&,const ekat::ParameterList&>;

} // namespace scream

#endif // SCREAM_GRIDS_MANAGER_HPP

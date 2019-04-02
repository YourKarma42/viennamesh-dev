#include "viennameshpp/plugin.hpp"
#include "cgal_mesh.hpp"
#include "cgal_mesh_simplification.hpp"
#include "cgal_automatic_mesh_simplification.hpp"
#include "cgal_curve.hpp"
#include "cgal_simplify_curve.hpp"
#include "cgal_remove_short_edges.hpp"
#include "cgal_compare_curvature_calculations.hpp"
#include "cgal_data_creator.hpp"




viennamesh_error viennamesh_plugin_init(viennamesh_context context)
{
  viennamesh::register_data_type<viennamesh::cgal::polyhedron_surface_mesh>(context); //CGAL's triangulated surface mesh data type

  viennamesh::register_conversion<viennagrid_mesh, viennamesh::cgal::polyhedron_surface_mesh>(context); //viennagrid datastructure (DS) --> CGAL DS
  viennamesh::register_conversion<viennamesh::cgal::polyhedron_surface_mesh, viennagrid_mesh>(context); // CGAL DS --> viennagrid DS

  viennamesh::register_algorithm<viennamesh::cgal::cgal_mesh_simplification>(context); //mesh simplification with user defined policies and parameters
  viennamesh::register_algorithm<viennamesh::cgal::cgal_automatic_mesh_simplification>(context); //mesh simplifiaction with automatically chosen policies and parameters
  
  viennamesh::register_algorithm<viennamesh::cgal_curve>(context);
  //viennamesh::register_algorithm<viennamesh::cgal::cgal_statistic_2>(context); 

  viennamesh::register_algorithm<viennamesh::cgal::cgal_simplify_curve>(context);


  viennamesh::register_algorithm<viennamesh::cgal::cgal_remove_short_edges>(context);

  viennamesh::register_algorithm<viennamesh::cgal::cgal_compare_curvature_calculations>(context);

  viennamesh::register_algorithm<viennamesh::cgal::cgal_data_creator>(context);

  return VIENNAMESH_SUCCESS;
}

int viennamesh_version()
{
  return VIENNAMESH_VERSION;
}

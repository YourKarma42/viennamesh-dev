#include "viennameshpp/plugin.hpp"
#include "vtk_mesh.hpp"
#include "vtk_decimate_pro.hpp"
#include "vtk_quadric_clustering.hpp"
#include "vtk_quadric_decimation.hpp"
#include "vtk_mesh_quality.hpp"
#include "vtk_clean.hpp"
#include "vtk_print_mesh_data.hpp"
#include "vtk_mesh_distance.hpp"
#include "vtk_calculate_mesh_quality_metrics.hpp"

viennamesh_error viennamesh_plugin_init(viennamesh_context context)
{
    viennamesh::register_data_type<viennamesh::vtk::mesh>(context);

    viennamesh::register_conversion<viennagrid_mesh, viennamesh::vtk::mesh>(context);
    viennamesh::register_conversion<viennamesh::vtk::mesh, viennagrid_mesh>(context);

    viennamesh::register_algorithm<viennamesh::vtk::decimate_pro>(context);
    viennamesh::register_algorithm<viennamesh::vtk::quadric_clustering>(context);
    viennamesh::register_algorithm<viennamesh::vtk::quadric_decimation>(context);
    viennamesh::register_algorithm<viennamesh::vtk::mesh_quality>(context);
    viennamesh::register_algorithm<viennamesh::vtk::clean>(context);
    viennamesh::register_algorithm<viennamesh::vtk::mesh_distance>(context);

    viennamesh::register_algorithm<viennamesh::vtk::print_Mesh_data>(context);

    viennamesh::register_algorithm<viennamesh::vtk::calculate_mesh_quality_metrics>(context);


    return VIENNAMESH_SUCCESS;
}

int viennamesh_version()
{
    return VIENNAMESH_VERSION;
}

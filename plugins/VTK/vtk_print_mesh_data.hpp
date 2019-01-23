#ifndef VIENNAMESH_VTK_PRINT_MESH_DATA_HPP
#define VIENNAMESH_VTK_PRINT_MESH_DATA_HPP

#include "viennameshpp/plugin.hpp"

namespace viennamesh
{
    namespace vtk
    {
        class print_Mesh_data :
                public plugin_algorithm
        {
            public:
                print_Mesh_data();

                static std::string name() { return "vtk_print_mesh_data"; }
                bool run(viennamesh::algorithm_handle&);
        };
    }
}

#endif //VIENNAMESH_VTK_PRINT_MESH_DATA_HPP
#ifndef VIENNAMESH_VTK_MESH_DISTANCE_HPP
#define VIENNAMESH_VTK_MESH_DISTANCE_HPP

#include "viennameshpp/plugin.hpp"

namespace viennamesh
{
    namespace vtk
    {
        class mesh_distance :
                public plugin_algorithm
        {
            public:
                mesh_distance();

                static std::string name() { return "vtk_mesh_distance"; }
                bool run(viennamesh::algorithm_handle&);
        };
    }
}

#endif //VIENNAMESH_VTK_VTK_MESH_DISTANCE_HPP
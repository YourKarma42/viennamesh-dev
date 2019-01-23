#ifndef VIENNAMESH_VTK_CLEAN_HPP
#define VIENNAMESH_VTK_CLEAN_HPP

#include "viennameshpp/plugin.hpp"

namespace viennamesh
{
    namespace vtk
    {
        class clean :
                public plugin_algorithm
        {
            public:
                clean();

                static std::string name() { return "vtk_clean"; }
                bool run(viennamesh::algorithm_handle&);
        };
    }
}

#endif //VIENNAMESH_VTK_CLEAN_HPP
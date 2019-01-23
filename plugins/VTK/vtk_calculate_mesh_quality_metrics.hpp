#ifndef VIENNAMESH_CALCULATE_MESH_QUALITY_METRICS_HPP
#define VIENNAMESH_CALCULATE_MESH_QUALITY_METRICS_HPP

#include "viennameshpp/plugin.hpp"

namespace viennamesh
{
    namespace vtk
    {
        class calculate_mesh_quality_metrics :
                public plugin_algorithm
        {
            public:
                calculate_mesh_quality_metrics();

                static std::string name() { return "vtk_calculate_mesh_quality_metrics"; }
                bool run(viennamesh::algorithm_handle&);
        };
    }
}

#endif //VIENNAMESH_CALCULATE_MESH_QUALITY_METRICS_HPP
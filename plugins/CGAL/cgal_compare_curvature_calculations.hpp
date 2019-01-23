#ifndef VIENNAMESH_ALGORITHM_CGAL_COMPARE_CURVATURE_CALCULATIONS_HPP
#define VIENNAMESH_ALGORITHM_CGAL_COMPARE_CURVATURE_CALCULATIONS_HPP

/* ============================================================================
   Copyright (c) 2011-2016, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.

                            -----------------
                ViennaMesh - The Vienna Meshing Framework
                            -----------------

                    http://viennamesh.sourceforge.net/

   License:         MIT (X11), see file LICENSE in the base directory
=============================================================================== */

/*
 *   This algorithm depends on the mesh comparison capabilities of the statistics plugin. Thus, statistics plugin must be enabled!
*/

#include "viennameshpp/plugin.hpp"

namespace viennamesh
{
    namespace cgal
    {
        class cgal_compare_curvature_calculations : public plugin_algorithm
        {
        public:
            cgal_compare_curvature_calculations();

            static std::string name();
            bool run(viennamesh::algorithm_handle &);
        };
    }

}



#endif

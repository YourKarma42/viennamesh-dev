#ifndef VIENNAMESH_ALGORITHM_CGAL_FIRST_TEST_HPP
#define VIENNAMESH_ALGORITHM_CGAL_FIRST_TEST_HPP

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


#include "viennameshpp/plugin.hpp"





namespace viennamesh
{



        class cgal_curve : public plugin_algorithm
        {
        public:
            cgal_curve();

            static std::string name();

            bool run(viennamesh::algorithm_handle &);

        };


}



#endif
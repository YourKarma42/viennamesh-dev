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

/*
 *   testclass
*/




#include "viennameshpp/plugin.hpp"

/*#include <limits>
#include <set>
#include <cmath>*/




//libigl includes
#include <igl/bfs_orient.h>

#include <igl/hausdorff.h>
#include <igl/bounding_box_diagonal.h>
#include <igl/point_mesh_squared_distance.h>

#include <igl/gaussian_curvature.h>
#include <igl/barycentric_coordinates.h>

#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/principal_curvature.h>






namespace viennamesh
{



        class cgal_first_test : public plugin_algorithm
        {
        public:
            cgal_first_test();

            static std::string name();

            bool run(viennamesh::algorithm_handle &);

            void curvature(mesh_handle &);
        };


}



#endif
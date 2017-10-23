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


#include "cgal_curve.hpp"


namespace viennamesh
{

        cgal_curve::cgal_curve(){}

        std::string cgal_curve::name()
        {
           
            return "cgal_curve";
        }





        bool cgal_curve::run(viennamesh::algorithm_handle &){
            


 
        }




}
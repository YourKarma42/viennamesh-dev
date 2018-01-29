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

#include "cgal_mesh.hpp"

#include "cgal_simplify_curve.hpp"

#include<chrono>

#include <boost/algorithm/string.hpp>  // sto_lower(string str)
#include <typeinfo>   // operator typeid

// Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>

// Stop-condition policies
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>

// Cost policies
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_cost.h>

// Placement policies
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_placement.h>




//my test placement policy
#include <CGAL/Curvature_policies/cgal_curve_placement.h>

//my test cost policy
#include <CGAL/Curvature_policies/cgal_curve_cost.h>

//my test stop prediction
#include <CGAL/Curvature_policies/cgal_curve_stop.h>

//hash table to store curvatures
#include <unordered_map>

//my hash function for unordered map

#include <CGAL/Curvature_policies/cgal_hash_points.hpp>



//Feature preservation
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Constrained_placement.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Unique_hash_map.h>




namespace SMS = CGAL::Surface_mesh_simplification;


namespace viennamesh
{
    namespace cgal
    {

        cgal_simplify_curve::cgal_simplify_curve(){}

        std::string cgal_simplify_curve::name()
        {
           
            return "cgal_simplify_curve";
        }


         // constrained edges are set via constrained_edge_collapse
        template <class STOP, class COST, class PLACEMENT>
        int smart_edge_collapse (cgal::polyhedron_surface_mesh & surface_mesh, STOP stop, COST cost, PLACEMENT placement )
        {
            //print_selected_parameters (stop, cost, placement);

            info(5) << "cost:        " << std::endl;
            return SMS::edge_collapse
                   (surface_mesh
                    ,stop
                    ,CGAL::parameters::vertex_index_map(get(CGAL::vertex_external_index,surface_mesh))
                    .halfedge_index_map  (get(CGAL::halfedge_external_index, surface_mesh))
                    .get_cost (cost)
                    .get_placement(placement)
                   );
        }



        bool cgal_simplify_curve::run(viennamesh::algorithm_handle &){

            data_handle<viennagrid_numeric> input_ratio = get_input<viennagrid_numeric>("ratio");

            data_handle<viennamesh_string> policy = get_required_input<viennamesh_string>("policy");          

            //_________________________Testing_Area_______________________________________________________
            // clear costs file

            std::ofstream myfile;
            myfile.open ("costs.txt");
            if (myfile.is_open()){  
        
              myfile << "";
              myfile.close();
            }

            //_________________________END_Testing_Area_______________________________________________________
            
            //Get mesh data (REQUIRED)
            data_handle<cgal::polyhedron_surface_mesh> input_mesh = get_required_input<cgal::polyhedron_surface_mesh>("mesh");


            //Create output mesh handle
            data_handle<cgal::polyhedron_surface_mesh> output_mesh = make_data<cgal::polyhedron_surface_mesh>();

            //Declare internal reference to output mesh
            cgal::polyhedron_surface_mesh & my_mesh = const_cast<cgal::polyhedron_surface_mesh&> (output_mesh());
            my_mesh = input_mesh();


            
            //Calculate or import curvature

            //Lindstrom Turk policies

            viennagrid_numeric volume_weight = 0.1;
            viennagrid_numeric boundary_weight = 0.5;
            viennagrid_numeric shape_weight = 0.4;

            viennagrid_numeric ratio_of_edges = input_ratio();

            
            
            
            int removed_edges=0;

            std::vector<double> time_for_curves;


            auto start = std::chrono::high_resolution_clock::now();

            auto finish = std::chrono::high_resolution_clock::now();


            if(policy() == "lt"){

                //Lindstrom Turk policys

                SMS::LindstromTurk_cost<cgal::polyhedron_surface_mesh> cost(SMS::LindstromTurk_params(volume_weight,boundary_weight,shape_weight));

                SMS::Count_ratio_stop_predicate<cgal::polyhedron_surface_mesh> stop( ratio_of_edges );


                SMS::LindstromTurk_placement<cgal::polyhedron_surface_mesh> placement(SMS::LindstromTurk_params(volume_weight,boundary_weight,shape_weight));

                start = std::chrono::high_resolution_clock::now();

                removed_edges = smart_edge_collapse(my_mesh, stop, cost, placement);

                finish = std::chrono::high_resolution_clock::now();

            }else if(policy() == "m"){

                //Testing area my policies



                //create Hashtable to store curvatures

                std::unordered_map<Point_3, std::pair<double, double>, Point_3_Hash, Point_3_Equal> curvatures;



                SMS::Curvature_placement<cgal::polyhedron_surface_mesh> placement_test; //atm midpoint

                SMS::LindstromTurk_placement<cgal::polyhedron_surface_mesh> placement(SMS::LindstromTurk_params(volume_weight,boundary_weight,shape_weight));

                SMS::Curvature_cost<cgal::polyhedron_surface_mesh> cost_test(curvatures, time_for_curves); //atm my curve cost

                SMS::Curvature_stop_predicate<cgal::polyhedron_surface_mesh> stop_test( ratio_of_edges );

                start = std::chrono::high_resolution_clock::now();

                removed_edges = smart_edge_collapse(my_mesh, stop_test, cost_test, placement_test);

                finish = std::chrono::high_resolution_clock::now();


                //outputs for development

                info(1) << Hashtable_Statistics(curvatures) << std::endl;

                info(1) << curvature_output(curvatures) << std::endl;


            }



            int i = 0;
            /*for(mesh_t::Vertex_iterator at=mesh.vertices_begin(),end=mesh.vertices_end();at!=end;++at){

                if(i%300 == 0){
                    //at->vertex_id
                }


                i++;

            }*/

            long t=0;
            long max =0;
            for ( auto &i : time_for_curves ) {
                if(i > max)
                    max = i;
                t += i;
            }

            info(1)  << "----------------++++++++++++++++++ Time output ++++++++++++++++++----------------" << std::endl;

            info(1) << "runtime : " << std::chrono::duration_cast<std::chrono::seconds>(finish-start).count() << " s" << std::endl 
                    << "          " << std::chrono::duration_cast<std::chrono::microseconds>(finish-start).count() << " micro s" <<  std::endl
                    << "curves  : " << t<< " micro s" <<  std::endl
                    << "Max time: " << max << " micro s" << std::endl
                    << time_for_curves.size() << std::endl;

            info(1) << "--------------------------------------------------------------------------------" << std::endl;

            info(1) << "removed edges:" << removed_edges << std::endl;

            info(1) << "saving the mesh..." << std::endl;

            set_output("mesh", my_mesh);

            return true;
 
        }

       
    }
}
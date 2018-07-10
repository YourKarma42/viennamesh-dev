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

#include<fstream>
#include<chrono>

#include <boost/algorithm/string.hpp>  // sto_lower(string str)
#include <typeinfo>   // operator typeid

//#define CGAL_SURFACE_SIMPLIFICATION_ENABLE_TRACE 4

void Surface_simplification_external_trace( std::string s )
{
   static std::ofstream out("log.txt");
    // std::ofstream out;
  //out.open ("log.txt");
  

   out << s << std::endl ;

     //out.close();
} 

#include "cgal_mesh.hpp"

#include "cgal_remove_short_edges.hpp"



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




//hash table to store curvatures
#include <unordered_map>


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

        cgal_remove_short_edges::cgal_remove_short_edges(){}

        std::string cgal_remove_short_edges::name()
        {
           
            return "cgal_remove_short_edges";
        }


         // constrained edges are set via constrained_edge_collapse
        template <class STOP, class COST, class PLACEMENT>
        int smart_edge_collapse (cgal::polyhedron_surface_mesh & surface_mesh, STOP stop, COST cost, PLACEMENT placement )
        {

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



        bool cgal_remove_short_edges::run(viennamesh::algorithm_handle &){

            data_handle<viennagrid_numeric> min_edge_length = get_input<viennagrid_numeric>("min_edge_size");      


            //Get mesh data (REQUIRED)
            data_handle<cgal::polyhedron_surface_mesh> input_mesh = get_required_input<cgal::polyhedron_surface_mesh>("mesh");


            //Create output mesh handle
            data_handle<cgal::polyhedron_surface_mesh> output_mesh = make_data<cgal::polyhedron_surface_mesh>();

            //Declare internal reference to output mesh
            cgal::polyhedron_surface_mesh & my_mesh = const_cast<cgal::polyhedron_surface_mesh&> (output_mesh());
            my_mesh = input_mesh();

            std::vector<double> length;

           
            info(5)<< "stuff_started" <<std::endl;

            //info(5)<<min_edge_length()*min_edge_length()<<std::endl;

           // length.push_back(CGAL::squared_distance(at->vertex()->point(), at->opposite()->vertex()->point()));

           int prev_v = my_mesh.size_of_vertices();
           int prev_he = my_mesh.size_of_halfedges();
           int prev_f = my_mesh.size_of_facets();

           info(5)<< "Mesh is valid: "<< my_mesh.is_valid() << std::endl;
           info(5)<< "Vertices vorher: "<< prev_v <<std::endl;
           info(5)<< "Edges vorher: "<< prev_he <<std::endl;
           info(5)<< "Facests vorher: "<< prev_f<<std::endl;
 

            for ( auto e = my_mesh.edges_begin(); e != my_mesh.edges_end(); ++e){

                length.push_back(CGAL::squared_distance(e->prev()->vertex()->point(),e->vertex()->point()));

            }

            double t = min_edge_length()*min_edge_length();

            std::sort(length.begin(), length.end());

            info(5)<< length.back()  <<std::endl;

            for(int i =0; i < 8; i++){
                info(5)<< length[i]  <<std::endl;
            }

            cgal::polyhedron_surface_mesh::Vertex_handle v;

            cgal::polyhedron_surface_mesh::Facet_handle f1, f2;
            

            //find the short edge
            for ( auto e = my_mesh.edges_begin(); e != my_mesh.edges_end(); ++e){

                //dosnt remove Facets

                if (CGAL::squared_distance(e->prev()->vertex()->point(),e->vertex()->point()) == length[0]){

                    f1 = e->face();

                    f2 = e->opposite()->face();

                    v = (my_mesh.join_vertex(e))->vertex();
                    break;
                }

            }

          

            //clear old faces from the Mesh (creates holes)
            cgal::polyhedron_surface_mesh::Vertex::Halfedge_around_vertex_circulator at=v->vertex_begin(), end = at;
            do{  
                

                if(at->opposite()->face() == f1){
                    
                    my_mesh.erase_facet(at->opposite()->opposite());
                    my_mesh.erase_facet(at->opposite());

                }

                if(at->opposite()->face() == f2){
                    
                    my_mesh.erase_facet(at->opposite()->opposite());
                    my_mesh.erase_facet(at->opposite());
                    
                }
                           
                at++;
            }while(at != end);


            //fill the holes
            at=v->vertex_begin();
            end = at;
            do{ 
                if(at->opposite()->is_border()){
                    my_mesh.fill_hole(at->opposite());

                }
                         
                at++;
            }while(at != end);

            at=v->vertex_begin();
            end = at;
            do{ 
              /*  info(5)<< "point origin: "<< at->opposite()->opposite()->vertex()->point() << std::endl;
                info(5)<< "point: "<< at->opposite()->vertex()->point() << std::endl;

                info(5)<< "bordde1: "<< at->opposite()->opposite()->is_border() << std::endl;
                info(5)<< "border2: "<< at->opposite()->is_border() << std::endl;
*/
                         
                at++;
            }while(at != end);


           info(5)<< "Mesh is valid: "<< my_mesh.is_valid() << std::endl;
           info(5)<< "Vertices delted: "<< prev_v - my_mesh.size_of_vertices() <<std::endl;
           info(5)<< "Edges deleted: "<< prev_he - my_mesh.size_of_halfedges() <<std::endl;
           info(5)<< "Facests deleted: "<< prev_f - my_mesh.size_of_facets()<<std::endl;




            info(5)<< "stuff finished" <<std::endl;



            set_output("mesh", my_mesh);


            return true;

        }

    }
}
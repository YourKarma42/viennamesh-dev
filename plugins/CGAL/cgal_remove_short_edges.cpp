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


#include "cgal_mesh.hpp"

#include "cgal_remove_short_edges.hpp"


//hash table to store curvatures
#include <unordered_map>


//Feature preservation
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Constrained_placement.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Unique_hash_map.h>

//██╗    ██╗ █████╗ ██████╗ ███╗   ██╗██╗███╗   ██╗ ██████╗                           
//██║    ██║██╔══██╗██╔══██╗████╗  ██║██║████╗  ██║██╔════╝                           
//██║ █╗ ██║███████║██████╔╝██╔██╗ ██║██║██╔██╗ ██║██║  ███╗                          
//██║███╗██║██╔══██║██╔══██╗██║╚██╗██║██║██║╚██╗██║██║   ██║                          
//╚███╔███╔╝██║  ██║██║  ██║██║ ╚████║██║██║ ╚████║╚██████╔╝                          
// ╚══╝╚══╝ ╚═╝  ╚═╝╚═╝  ╚═╝╚═╝  ╚═══╝╚═╝╚═╝  ╚═══╝ ╚═════╝    

/*
WARNING this code is WITCHCRAFT

Idea:
- scan the mesh and find all short edges 

- edges are considred short by a user supplyed squard length (edge_length < length*length)

- each short edge is collapsed into a vertex with the CGAL function  Polyhedron_3.join_vertex(he)

- this step produces facets that are not fully connectet to the mesh

- delete non connected facets (creates holes)

- fill the holes

This code work for the supplyed meshes exp_conf.vtu and exp_noox_surf_conf.str.vtu.
It needs to be adepted for other meshes!!
*/
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

            std::vector<cgal::polyhedron_surface_mesh::Vertex_handle> new_vertices;

           
            info(5)<< "stuff_started" <<std::endl;


           int prev_v = my_mesh.size_of_vertices();
           int prev_he = my_mesh.size_of_halfedges();
           int prev_f = my_mesh.size_of_facets();

           info(5)<< "Mesh is valid: "<< my_mesh.is_valid() << std::endl;
           info(5)<< "Vertices vorher: "<< prev_v <<std::endl;
           info(5)<< "Edges vorher: "<< prev_he <<std::endl;
           info(5)<< "Facests vorher: "<< prev_f<<std::endl;


 
            //search short edges 
            for(cgal::polyhedron_surface_mesh::Halfedge_iterator he = my_mesh.halfedges_begin(), end = my_mesh.halfedges_end(); he !=end; ++he){

                length.push_back(CGAL::squared_distance(he->prev()->vertex()->point(),he->vertex()->point()));

            }


            std::sort(length.begin(), length.end());


            for(int i =0; i <50; i++){
                info(5)<< length[i]  << std::endl;
            }

        

            int removed = 0;

            while(length[0] < min_edge_length() && removed <= 5){

                cgal::polyhedron_surface_mesh::Vertex_handle v;

                cgal::polyhedron_surface_mesh::Halfedge_handle h;

                cgal::polyhedron_surface_mesh::Facet_handle f1, f2;
                
                //find short edges and collapse them
                for(cgal::polyhedron_surface_mesh::Halfedge_iterator he = my_mesh.halfedges_begin(), end = my_mesh.halfedges_end(); he !=end; ++he){

                    //dosnt remove Facets

                    if (CGAL::squared_distance(he->prev()->vertex()->point(),he->vertex()->point()) == length[0] && removed<=3){

                        f1 = he->face();

                        f2 = he->opposite()->face();
                     
                        v = (my_mesh.join_vertex(he))->vertex();
                        
                        break;
                    }

                    if (CGAL::squared_distance(he->prev()->vertex()->point(),he->vertex()->point()) == length[4] && removed == 4){

                        f1 = he->face();

                        f2 = he->opposite()->face();

                        v = (my_mesh.join_vertex(he))->vertex();
                       
                        break;
                    }

                    if (CGAL::squared_distance(he->prev()->vertex()->point(),he->vertex()->point()) == length[3] && removed == 5){

                        f1 = he->face();

                        f2 = he->opposite()->face();

                        v = (my_mesh.join_vertex(he))->vertex();

                        
                        break;
                    }


                }

            


                //clear old faces from the Mesh (creates holes)
                //ERROR OCCURES HERE
                //some halfedge (of different length) are strangely linked with otehrs and produce errors depending on wich halfedges facet is choosen
                cgal::polyhedron_surface_mesh::Vertex::Halfedge_around_vertex_circulator at=v->vertex_begin(), end = at;
                do{  
                       

                    if(at->opposite()->face() == f1){
                        

                        my_mesh.erase_facet(at);


                        my_mesh.erase_facet(at->opposite());

                    }

                    if(at->opposite()->face() == f2){
                    
                        my_mesh.erase_facet(at);

                        my_mesh.erase_facet(at->opposite());
                                              
                    }

                    new_vertices.push_back(v);
                            
                    at++;
                }while(at != end);

                               



                length.clear();
                //search short edges
                for(cgal::polyhedron_surface_mesh::Halfedge_iterator he = my_mesh.halfedges_begin(), end = my_mesh.halfedges_end(); he !=end; ++he){

                    length.push_back(CGAL::squared_distance(he->prev()->vertex()->point(),he->vertex()->point()));

                }
                std::sort(length.begin(), length.end());

                info(5)<< "next 2 Half-Edges: " << length[0] << " and " << length[1]<< std::endl;

                for(int i =0; i < 10; i++){
                    info(5)<< length[i]  << std::endl;
                }


                removed++;

            }

            info(5)<< "filling holes..." << std::endl;

            //fill the holes

            for(auto v: new_vertices){
                
                cgal::polyhedron_surface_mesh::Vertex::Halfedge_around_vertex_circulator at=v->vertex_begin(), end = at;
                do{ 
                    if(at->opposite()->is_border()){
                        my_mesh.fill_hole(at->opposite());

                        if(at->is_border()){                            
                            info(5)<< "BAD EDGE nebeneinander" << std::endl;
                        }

                    }

                    if(at->is_border()){
                        my_mesh.fill_hole(at);
                        info(5)<< "BAD EDGE" << std::endl;

                    }
                            
                    at++;
                }while(at != end);
            }


            std::cout << std::endl;


           info(5)<< "Edges removed: "<< removed << std::endl;
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



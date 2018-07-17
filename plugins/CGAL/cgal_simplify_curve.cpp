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

//curvature testfunction
#include <CGAL/curve_calc_dev/curve_test_1.hpp>

//curvature testfunction
#include <CGAL/analizing/cgal_mesh_analytics.hpp>


#include "CGAL/Curvature_policies/boundary/cgal_curve_boundary_cost.h"
#include "CGAL/Curvature_policies/boundary/cgal_curve_boundary_placement.h"
#include "CGAL/Curvature_policies/boundary/cgal_curve_boundary_stop.h"

#include "CGAL/Curvature_policies/flat/cgal_curve_flat_cost.hpp"
#include "CGAL/Curvature_policies/flat/cgal_curve_flat_placement.hpp"
#include "CGAL/Curvature_policies/flat/cgal_curve_flat_stop.hpp"


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

            int removed_edges_this_step=0;

            std::vector<double> time_for_curves;


            auto start = std::chrono::high_resolution_clock::now();

            auto finish = std::chrono::high_resolution_clock::now();


            start = std::chrono::high_resolution_clock::now();

            cgal_mesh_analytics<cgal::polyhedron_surface_mesh> analytics(my_mesh);

            finish = std::chrono::high_resolution_clock::now();

            info(1) << "runtime analytics: " << std::chrono::duration_cast<std::chrono::seconds>(finish-start).count() << " s" << std::endl;

            //cgal_mesh_analytics<cgal::polyhedron_surface_mesh> test(1.0,2.0,3.0);




            //create Hashtable to store curvatures
            //move into if statment
            std::unordered_map<Point_3, viennamesh::cgal::vertex_statistics, Point_3_Hash, Point_3_Equal> curvatures;

           


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





                int edge_edges=0;

                int ecke=0;

                //------------------------------------step one edge------------------------------

                SMS::Curvature_boundary_placement<cgal::polyhedron_surface_mesh> placement_test(analytics); //atm midpoint
           
                SMS::Curvature_boundary_cost<cgal::polyhedron_surface_mesh> cost_test(analytics); //atm my curve cost

                SMS::Curvature_boundary_stop<cgal::polyhedron_surface_mesh> stop_test(analytics, ratio_of_edges);

                start = std::chrono::high_resolution_clock::now();

                removed_edges_this_step = smart_edge_collapse(my_mesh, stop_test, cost_test, placement_test);

                finish = std::chrono::high_resolution_clock::now();

                removed_edges += removed_edges_this_step;

                info(1) << "removed boundary edges:" << removed_edges_this_step << std::endl;
                info(1) << "runtime boundary edges: " << std::chrono::duration_cast<std::chrono::seconds>(finish-start).count() << " s" << std::endl; 

                //------------------------------- Step 2 flat plane------------------------------
                SMS::Curvature_flat_placement<cgal::polyhedron_surface_mesh> placement_flat(analytics); //atm midpoint

                SMS::Curvature_flat_cost<cgal::polyhedron_surface_mesh> cost_flat(analytics); //atm my curve cost

                SMS::Curvature_flat_stop<cgal::polyhedron_surface_mesh> stop_flat(0.316);

                start = std::chrono::high_resolution_clock::now();

                removed_edges_this_step = smart_edge_collapse(my_mesh, stop_flat, cost_flat, placement_flat);

                finish = std::chrono::high_resolution_clock::now();

                info(1) << "removed flat edges:" << removed_edges_this_step << std::endl;
                info(1) << "runtime flat edges: " << std::chrono::duration_cast<std::chrono::seconds>(finish-start).count() << " s" << std::endl;

                removed_edges += removed_edges_this_step;



                //outputs for development

                //info(1) << Hashtable_Statistics(curvatures) << std::endl;

               // info(1) << curvature_output(curvatures) << std::endl;

               // info(1) << "!!!!!!!!!! Eck Kanten: " << edge_edges << std::endl;

               // info(1) << "!!!!!!!!!! ECKEN placement: " << ecke << std::endl;


                /*std::ofstream myfile;
                myfile.open ("costs.txt", std::ios_base::app);
                if (myfile.is_open()){

                    for ( auto &it : curvatures ){

                        myfile << (it.second).first << " " << (it.second).second << "\n";
                    }
                    
                    }
                myfile.close();*/


            }



            int i = 0;
            /*for(mesh_t::Vertex_iterator at=mesh.vertices_begin(),end=mesh.vertices_end();at!=end;++at){

                if(i%300 == 0){
                    //at->vertex_id
                }


                i++;

            }*/

           /* long t=0;
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
*/

            info(1) << "overall removed edges:" << removed_edges << std::endl;


            //___________________________________Testing_writing data into VTU file________________________________________________
            
            if(policy() == "m" || policy() == "p"){
            
                //create Vector to store quantities
                std::vector<viennagrid_quantity_field> quantity_fields;

                //create one quantity
                viennagrid_quantity_field quantity_field=0;
                viennagrid_quantity_field_create(&quantity_field);

                //set the type of the quantity
                viennagrid_quantity_field_init(quantity_field,
                                VIENNAGRID_ELEMENT_TYPE_VERTEX ,                                        // topological dimension of the elements
                                VIENNAGRID_QUANTITY_FIELD_TYPE_NUMERIC,   // floats
                                1,                                        // one float per element
                                VIENNAGRID_QUANTITY_FIELD_STORAGE_SPARSE);

                            //create one quantity
                viennagrid_quantity_field quantity_field1=0;
                viennagrid_quantity_field_create(&quantity_field1);

                //set the type of the quantity
                viennagrid_quantity_field_init(quantity_field1,
                                VIENNAGRID_ELEMENT_TYPE_VERTEX ,                                        // topological dimension of the elements
                                VIENNAGRID_QUANTITY_FIELD_TYPE_NUMERIC,   // floats
                                1,                                        // one float per element
                                VIENNAGRID_QUANTITY_FIELD_STORAGE_SPARSE);



                //create one quantity
                viennagrid_quantity_field quantity_field2=0;
                viennagrid_quantity_field_create(&quantity_field2);

                //set the type of the quantity
                viennagrid_quantity_field_init(quantity_field2,
                                VIENNAGRID_ELEMENT_TYPE_VERTEX ,                                        // topological dimension of the elements
                                VIENNAGRID_QUANTITY_FIELD_TYPE_NUMERIC,   // floats
                                1,                                        // one float per element
                                VIENNAGRID_QUANTITY_FIELD_STORAGE_SPARSE);
                

                //create one quantity
                viennagrid_quantity_field quantity_field3=0;
                viennagrid_quantity_field_create(&quantity_field3);

                //set the type of the quantity
                viennagrid_quantity_field_init(quantity_field3,
                                VIENNAGRID_ELEMENT_TYPE_VERTEX ,                                        // topological dimension of the elements
                                VIENNAGRID_QUANTITY_FIELD_TYPE_NUMERIC,   // floats
                                1,                                        // one float per element
                                VIENNAGRID_QUANTITY_FIELD_STORAGE_SPARSE);


                //create one quantity
                viennagrid_quantity_field quantity_field4=0;
                viennagrid_quantity_field_create(&quantity_field4);

                //set the type of the quantity
                viennagrid_quantity_field_init(quantity_field4,
                                VIENNAGRID_ELEMENT_TYPE_VERTEX ,                                        // topological dimension of the elements
                                VIENNAGRID_QUANTITY_FIELD_TYPE_NUMERIC,   // floats
                                1,                                        // one float per element
                                VIENNAGRID_QUANTITY_FIELD_STORAGE_SPARSE);


                                                
                info(1)  << "----------------++++++++++++++++++ cgal_curve ++++++++++++++++++----------------" << std::endl;


#pragma region RECALCULATE principle curvatures and mean curvature write VTU
               /* long id=0;

                start = std::chrono::high_resolution_clock::now();

                for(cgal::polyhedron_surface_mesh::Vertex_iterator at=my_mesh.vertices_begin(),end=my_mesh.vertices_end();at!=end;++at)
                {
                    //add data into the quantitiy

                        double blub[2];

                        double a = principal_curvatures_cgal(*at,my_mesh, blub);

                        a=fabs((blub[0]+blub[1])/2.0);
                    
                        viennagrid_quantity_field_value_set(quantity_field,id,&a);

                        double b = blub[0]*blub[1];

                        viennagrid_quantity_field_value_set(quantity_field1,id,&b);
  
                    id++;
                }

                finish = std::chrono::high_resolution_clock::now();

                info(1)  << std::chrono::duration_cast<std::chrono::microseconds>(finish-start).count() << " micro s" << std::endl;

                // add current quantity to the vector
                viennagrid_quantity_field_name_set(quantity_field,"cgal_mean_curve");
                quantity_fields.push_back(quantity_field);

                viennagrid_quantity_field_name_set(quantity_field1,"cgal_gauss_curve");
                quantity_fields.push_back(quantity_field1);

                info(1)  << "----------------++++++++++++++++++ my_curve ++++++++++++++++++----------------" << std::endl;

                id=0;

                start = std::chrono::high_resolution_clock::now();

                for(cgal::polyhedron_surface_mesh::Vertex_iterator at=my_mesh.vertices_begin(),end=my_mesh.vertices_end();at!=end;++at)
                {

                        double blub[2];

                        double a = mean_curvature_my(*at);

                       // double a = 20;

                        if(std::isnan(a)){
                            info(1) << id << std::endl;
                            info(1) << (*at).point() << std::endl;
                            info(1) << my_curvature_test_print(*at) << std::endl;
                            a=150.0;
                        }
                        if(std::isinf(a)){
                           // info(1) << id << std::endl;
                            //info(1) << (*at).point() << std::endl;
                            
                            a=150.0;
                        }

                        //a=(blub[0]+blub[1])/2.0;

                        double b= gauss_curvatures_my(*at);


                        if(std::isnan(b)){

                            info(1) << my_curvature_test_print(*at) << std::endl;

                            info(1) << (*at).point() << std::endl;
                            //info(1) << id << std::endl;
                            b=100.0;
                        }
  
                    
                        viennagrid_quantity_field_value_set(quantity_field2,id,&a);

                        viennagrid_quantity_field_value_set(quantity_field3,id,&b);
  
                    id++;
                }

                finish = std::chrono::high_resolution_clock::now();

                info(1) << std::chrono::duration_cast<std::chrono::microseconds>(finish-start).count() << " micro s" << std::endl;

                

                // add current quantity to the vector

                viennagrid_quantity_field_name_set(quantity_field2,"my_mean_curve");
                quantity_fields.push_back(quantity_field2);

                viennagrid_quantity_field_name_set(quantity_field3,"my_gauss_curve");
                quantity_fields.push_back(quantity_field3);

                //export quantities
                data_handle<viennagrid_quantity_field> quantity_field_handle=viennamesh::plugin_algorithm::make_data(to_cpp(quantity_fields[0]));  
                quantity_field_handle.push_back(to_cpp(quantity_fields[1]));
                quantity_field_handle.push_back(to_cpp(quantity_fields[2]));    
                quantity_field_handle.push_back(to_cpp(quantity_fields[3]));        
                set_output("quantities",quantity_field_handle);*/

#pragma endregion

            

#pragma region write principle curvatures and mean curvature in VTU
                
                long id=0;

                for(cgal::polyhedron_surface_mesh::Vertex_iterator at=my_mesh.vertices_begin(),end=my_mesh.vertices_end();at!=end;++at)
                {

                    double a = analytics.get_mean(at);

                    double b = analytics.get_gauss(at);

                    double c = analytics.get_p1(at);

                    double d = analytics.get_p2(at);

                    double e=analytics.get_color(at);

                    /*if(analytics.get_color(at) >= 0){
                        if(analytics.get_color(at) == 1){
                            e = 50.0;
                        }else if(analytics.get_color(at) == 3) { 
                            e = 100.0;
                        }

                    }*/

                        viennagrid_quantity_field_value_set(quantity_field,id,&(c));

                        viennagrid_quantity_field_value_set(quantity_field1,id,&(d));

                        viennagrid_quantity_field_value_set(quantity_field2,id,
                            &(a));

                        viennagrid_quantity_field_value_set(quantity_field3,id,
                            &(b));

                        viennagrid_quantity_field_value_set(quantity_field4,id,
                            &(e));


                    id++;
                }

                /*for(cgal::polyhedron_surface_mesh::Point_iterator at=my_mesh.points_begin(),end=my_mesh.points_end();at!=end;++at)
                {



                    //add data into the quantitiy

                    if(curvatures.find(*at) == curvatures.end()){

                        info(1) << (*at) << std::endl;

                    }else{
                    
                        viennagrid_quantity_field_value_set(quantity_field,id,&((curvatures.find(*at)->second).curvature_1));

                        viennagrid_quantity_field_value_set(quantity_field1,id,&((curvatures.find(*at)->second).curvature_2));

                        double blub = fabs((curvatures.find(*at)->second).curvature_1)+fabs((curvatures.find(*at)->second).curvature_2);

                        viennagrid_quantity_field_value_set(quantity_field2,id,
                            &(blub));
                    }
                    
                    id++;
                }*/



                // add current quantity to the vector
                viennagrid_quantity_field_name_set(quantity_field,"PC1");
                quantity_fields.push_back(quantity_field);

                viennagrid_quantity_field_name_set(quantity_field1,"PC2");
                quantity_fields.push_back(quantity_field1);

                viennagrid_quantity_field_name_set(quantity_field2,"mean");
                quantity_fields.push_back(quantity_field2);

                viennagrid_quantity_field_name_set(quantity_field3,"gauss");
                quantity_fields.push_back(quantity_field3);

                viennagrid_quantity_field_name_set(quantity_field4,"tracing_test");
                quantity_fields.push_back(quantity_field4);

                

                //export quantities
                data_handle<viennagrid_quantity_field> quantity_field_handle=viennamesh::plugin_algorithm::make_data(to_cpp(quantity_fields[0]));  
                quantity_field_handle.push_back(to_cpp(quantity_fields[1]));    
                quantity_field_handle.push_back(to_cpp(quantity_fields[2])); 
                quantity_field_handle.push_back(to_cpp(quantity_fields[3]));
                quantity_field_handle.push_back(to_cpp(quantity_fields[4]));       
                set_output("quantities",quantity_field_handle);

#pragma endregion

            }


        

            //___________________________________Testing_writing data into VTU file________________________________________________


            //first new curvature tests

            int h = 0;
            for(cgal::polyhedron_surface_mesh::Vertex_iterator at=my_mesh.vertices_begin(),end=my_mesh.vertices_end();at!=end;++at)
            {

               
               /* if(fabs((*at).point().x() - 1.0) <= 0.00001 &&
                   fabs((*at).point().y() - 0.0) <= 0.00001 &&
                   fabs((*at).point().z() - 0.5) <= 0.00001) {   
                    
                    info(1) << "&& Schlecht &&" <<std::endl;
                    info(1) << (*at).point() <<std::endl; 
                    info(1) << mean_test(*at) <<std::endl; 
                    info(1) << my_curvature_test_print(*at) <<std::endl;
                    
                }

               
                if(fabs((*at).point().x() - 1.0) <= 0.00001 &&
                   fabs((*at).point().y() - 0.0) <= 0.00001 &&
                   fabs((*at).point().z() - 0.0) <= 0.00001) {   
                       
                    info(1) << "&& Gut &&" <<std::endl; 
                    info(1) << (*at).point() <<std::endl; 
                    info(1) << mean_test(*at) <<std::endl;  
                    info(1) << my_curvature_test_print(*at) <<std::endl;

                }*/

                /*if (h==5){
                    info(1) << my_curvature_test_print(*at) <<std::endl;
                    info(1) << (*at).point() <<std::endl;
                }*/


                h++;
    
            }        



          

            info(1) << "exporting the mesh..." << std::endl;


            set_output("mesh", my_mesh);


            return true;
 
        }


        /*data_handle<viennagrid_quantity_field> create_quanteties_test(){

            return NULL;
        }*/




       
    }
}
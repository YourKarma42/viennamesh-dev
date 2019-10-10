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




#include "CGAL/Curvature_policies/cgal_curve_cost.h"
#include "CGAL/Curvature_policies/cgal_flat_cost_p.h"


//hash table to store curvatures
#include <unordered_map>

//my hash function for unordered map
//#include <CGAL/Curvature_policies/cgal_hash_points.hpp>

//curvature testfunction
//#include <CGAL/curve_calc_dev/curve_test_1.hpp>

//curvature testfunction
#include <CGAL/analizing/cgal_mesh_analytics.hpp>


#include "CGAL/Curvature_policies/boundary/cgal_curve_boundary_cost.h"
#include "CGAL/Curvature_policies/boundary/cgal_curve_boundary_placement.h"
#include "CGAL/Curvature_policies/boundary/cgal_curve_boundary_stop.h"

#include "CGAL/Curvature_policies/flat/cgal_curve_flat_cost.hpp"
#include "CGAL/Curvature_policies/flat/cgal_curve_flat_placement.hpp"
#include "CGAL/Curvature_policies/flat/cgal_curve_flat_stop.hpp"

#include "CGAL/Curvature_policies/features/cgal_curve_features_cost.hpp"
#include "CGAL/Curvature_policies/features/cgal_curve_features_placement.hpp"
#include "CGAL/Curvature_policies/features/cgal_curve_features_stop.hpp"

#include "CGAL/Curvature_policies/curved/cgal_curve_curved_cost.hpp"
#include "CGAL/Curvature_policies/curved/cgal_curve_curved_placement.hpp"
#include "CGAL/Curvature_policies/curved/cgal_curve_curved_stop.hpp"


#include "CGAL/Curvature_policies/cgal_curve_stop.h"


#include <CGAL/Plane_3.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/intersections.h>

//Feature preservation
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Constrained_placement.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Unique_hash_map.h>



typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef K::Segment_2 Segment_2;
typedef K::Line_2 Line_2;

typedef K::Line_3 Line_3;
//typedef K::Plane_3 Plane_3;
typedef K::Intersect_3 Intersect_3;
typedef K::Intersect_2 Intersect_2;


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

            data_handle<viennagrid_numeric> anz_zones = get_required_input<viennagrid_numeric>("anz_zones");

            data_handle<viennagrid_numeric> zone_size_mult = get_required_input<viennagrid_numeric>("zone_multiplicator");

            std::list<double> zones;

            for(int i = 0; i< anz_zones(); i++){

                std::string zone = "zone" + std::to_string(i);

                data_handle<viennagrid_numeric> z = get_required_input<viennagrid_numeric>(zone);
                zones.push_back(z());

            }

            for(out: zones){
                std::cout << out << std::endl;
            }

           

   

            
            //Get mesh data (REQUIRED)
            data_handle<cgal::polyhedron_surface_mesh> input_mesh = get_required_input<cgal::polyhedron_surface_mesh>("mesh");


            //Create output mesh handle
            data_handle<cgal::polyhedron_surface_mesh> output_mesh = make_data<cgal::polyhedron_surface_mesh>();

            //Declare internal reference to output mesh
            cgal::polyhedron_surface_mesh & my_mesh = const_cast<cgal::polyhedron_surface_mesh&> (output_mesh());
            my_mesh = input_mesh();


            
            //Calculate or import curvature

            //Lindstrom Turk policies

            viennagrid_numeric volume_weight = 0.4;
            viennagrid_numeric boundary_weight = 0.3;
            viennagrid_numeric shape_weight = 0.3;

            viennagrid_numeric ratio_of_edges = input_ratio();

            
            
            
            int removed_edges=0;

            int removed_edges_this_step=0;

            std::vector<double> time_for_curves;

            info(1) << "Running analytics" << std::endl;


            auto start = std::chrono::high_resolution_clock::now();

            auto finish = std::chrono::high_resolution_clock::now();


            start = std::chrono::high_resolution_clock::now();

            cgal_mesh_analytics<cgal::polyhedron_surface_mesh> analytics(my_mesh);

            finish = std::chrono::high_resolution_clock::now();

            info(1) << "runtime analytics: " << std::chrono::duration_cast<std::chrono::microseconds>(finish-start).count() << " mus" << std::endl;

            double runtime_analytics = std::chrono::duration_cast<std::chrono::microseconds>(finish-start).count();


            double total_runtime = runtime_analytics;


            if(policy() == "lt"){

                //Lindstrom Turk policys


                //test

                SMS::Curvature_curved_cost<cgal::polyhedron_surface_mesh> 
                    cost_test(analytics, SMS::LindstromTurk_params(volume_weight,boundary_weight,shape_weight));


                total_runtime = 0.0;

                SMS::LindstromTurk_cost<cgal::polyhedron_surface_mesh> cost(SMS::LindstromTurk_params(volume_weight,boundary_weight,shape_weight));

                SMS::Count_ratio_stop_predicate<cgal::polyhedron_surface_mesh> stop( ratio_of_edges );

                //edge length
                SMS::Curvature_stop_predicate<cgal::polyhedron_surface_mesh> stop_2( ratio_of_edges);
                


                SMS::LindstromTurk_placement<cgal::polyhedron_surface_mesh> placement(SMS::LindstromTurk_params(volume_weight,boundary_weight,shape_weight));

                start = std::chrono::high_resolution_clock::now();

                removed_edges = smart_edge_collapse(my_mesh, stop_2, cost, placement);

                finish = std::chrono::high_resolution_clock::now();

                total_runtime = total_runtime + std::chrono::duration_cast<std::chrono::microseconds>(finish-start).count();



            }else if(policy() == "m"){

                //Testing area my policies





                int edge_edges=0;

                int ecke=0;

 

//------------------------------- Step 1 curved area------------------------------

                viennagrid_numeric curve_volume_weight = 0.4;
                viennagrid_numeric curve_boundary_weight = 0.3;
                viennagrid_numeric curve_shape_weight = 0.3;


                //SMS::Curvature_curved_placement<cgal::polyhedron_surface_mesh> placement_curved(analytics); //atm midpoint

               // SMS::Curvature_curved_cost<cgal::polyhedron_surface_mesh> 
                //     cost_curved(analytics, SMS::LindstromTurk_params(volume_weight,boundary_weight,shape_weight)); //atm my curve cost

                SMS::LindstromTurk_placement<cgal::polyhedron_surface_mesh> 
                    lt_placement(SMS::LindstromTurk_params(curve_volume_weight,curve_boundary_weight,curve_shape_weight));

                SMS::Curvature_curved_cost<cgal::polyhedron_surface_mesh> 
                    cost_curved(analytics, SMS::LindstromTurk_params(curve_volume_weight,curve_boundary_weight,curve_shape_weight));

                SMS::Curvature_curved_stop<cgal::polyhedron_surface_mesh> stop_curved(ratio_of_edges, analytics);
                //SMS::Count_ratio_stop_predicate<cgal::polyhedron_surface_mesh> stop( 0.2 );

                start = std::chrono::high_resolution_clock::now();



               removed_edges_this_step = smart_edge_collapse(my_mesh, /*stop*/ stop_curved, cost_curved, lt_placement /*placement_curved*/);

                finish = std::chrono::high_resolution_clock::now();

                info(1) << "removed curved edges:" << removed_edges_this_step << std::endl;
                info(1) << "runtime curved edges: " << std::chrono::duration_cast<std::chrono::microseconds>(finish-start).count() << " mus" << std::endl;

                total_runtime = total_runtime + std::chrono::duration_cast<std::chrono::microseconds>(finish-start).count();

                
                
               
                removed_edges += removed_edges_this_step;


                info(1) << "AVERAGE EDGE LENGHT CURVED:" << analytics.average_length_curved_edges() << std::endl;


               


 

 //------------------------------------step 2 transition area------------------------------

                //analytics.extend_transition_area();

                start = std::chrono::high_resolution_clock::now();

                analytics.calculate_curvatures();

                analytics.new_extrect_features();

                finish = std::chrono::high_resolution_clock::now();

               // analytics.calculate_transition_areas();



                viennagrid_numeric transition_volume_weight = 0.0;
                viennagrid_numeric transition_boundary_weight = 0.2;
                viennagrid_numeric transition_shape_weight = 0.8;

                SMS::Curvature_features_placement<cgal::polyhedron_surface_mesh> 
                  placement_transition(analytics, SMS::LindstromTurk_params(transition_volume_weight,transition_boundary_weight,transition_shape_weight));

                total_runtime = total_runtime + std::chrono::duration_cast<std::chrono::microseconds>(finish-start).count();



                double area_start=0.0;

                start = std::chrono::high_resolution_clock::now();

                /*for(auto s: zones){


                    analytics.calculate_transition_areas();

                    SMS::Curvature_features_cost<cgal::polyhedron_surface_mesh> cost_transition2(area_start * zone_size_mult(),
                    analytics, SMS::LindstromTurk_params(transition_volume_weight,transition_boundary_weight,transition_shape_weight)); 
      

                    SMS::Curvature_features_stop<cgal::polyhedron_surface_mesh> stop_transition2(analytics, s);

                    removed_edges_this_step = smart_edge_collapse(my_mesh, stop_transition2, cost_transition2, placement_transition);

                    //finish = std::chrono::high_resolution_clock::now();

                    removed_edges += removed_edges_this_step;
                   

                    area_start = area_start + s*zone_size_mult();
                    info(1) << "Area Start: " << area_start << std::endl;

                }*/

                double s =  ratio_of_edges + zone_size_mult();

                while(removed_edges_this_step != 0){


                    analytics.calculate_transition_areas();

                    SMS::Curvature_features_cost<cgal::polyhedron_surface_mesh> cost_transition2(area_start,
                    analytics, SMS::LindstromTurk_params(transition_volume_weight,transition_boundary_weight,transition_shape_weight)); 
      

                    SMS::Curvature_features_stop<cgal::polyhedron_surface_mesh> stop_transition2(analytics, s);

                    removed_edges_this_step = smart_edge_collapse(my_mesh, stop_transition2, cost_transition2, placement_transition);

                    //finish = std::chrono::high_resolution_clock::now();

                    removed_edges += removed_edges_this_step;
                   

                    area_start = area_start + s;
                    
                    s += zone_size_mult();

                    info(1) << "Area Start: " << area_start << std::endl;

                }

                finish = std::chrono::high_resolution_clock::now();

                info(1) << "removed transition edges:" << removed_edges_this_step << std::endl;
                info(1) << "runtime transition edges: " << std::chrono::duration_cast<std::chrono::microseconds>(finish-start).count() << " mus" << std::endl; 
                

                total_runtime = total_runtime + std::chrono::duration_cast<std::chrono::microseconds>(finish-start).count();


                //------------------------------- Simplification end -----------------------------

        
                analytics.calculate_metrics();




                info(1) << "NUMBER OF EDGES: " << my_mesh.size_of_halfedges() << std::endl;

            } 

            int i = 0;


            info(1) << "overall removed edges: " << removed_edges << std::endl;

            std::string output_file = "times_compare.csv";

            std::ofstream csv_my;
            csv_my.open(output_file.c_str(),  std::ios::app);

            csv_my << std::endl << total_runtime;

            csv_my.close();

            info(1) << "overall runtime: " << total_runtime <<  " mus" << std::endl;




            //___________________________________Testing_writing data into VTU file________________________________________________
            
            if(policy() == "m" || policy() == "p" || policy() == "lt"){


                //TODO: rethink
                cgal_mesh_analytics<cgal::polyhedron_surface_mesh> analytics2(my_mesh);
            
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

                //create one quantity
                viennagrid_quantity_field quantity_field5=0;
                viennagrid_quantity_field_create(&quantity_field5);

                //set the type of the quantity
                viennagrid_quantity_field_init(quantity_field5,
                                VIENNAGRID_ELEMENT_TYPE_VERTEX ,                                        // topological dimension of the elements
                                VIENNAGRID_QUANTITY_FIELD_TYPE_NUMERIC,   // floats
                                1,                                        // one float per element
                                VIENNAGRID_QUANTITY_FIELD_STORAGE_SPARSE);
 
                //create one quantity
                viennagrid_quantity_field quantity_field6=0;
                viennagrid_quantity_field_create(&quantity_field6);

                //set the type of the quantity
                viennagrid_quantity_field_init(quantity_field6,
                                VIENNAGRID_ELEMENT_TYPE_VERTEX ,                                        // topological dimension of the elements
                                VIENNAGRID_QUANTITY_FIELD_TYPE_NUMERIC,   // floats
                                1,                                        // one float per element
                                VIENNAGRID_QUANTITY_FIELD_STORAGE_SPARSE);

                //create one quantity
                viennagrid_quantity_field quantity_field7=0;
                viennagrid_quantity_field_create(&quantity_field7);

                //set the type of the quantity
                viennagrid_quantity_field_init(quantity_field7,
                                VIENNAGRID_ELEMENT_TYPE_VERTEX ,                                        // topological dimension of the elements
                                VIENNAGRID_QUANTITY_FIELD_TYPE_NUMERIC,   // floats
                                3,                                        // one float per element
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

                    double a = analytics2.get_mean(at);

                    double b = analytics2.get_gauss(at);

                    double c = analytics2.get_p1(at);

                    double d = analytics2.get_p2(at);

                    //double e = analytics.get_transition_area(at);

                    double e = analytics2.get_area(at);

                    if(isnan(e)){
                        std::cout << "blublbublublublublub" << std::endl;
                        e=0.0;

                    }
                        
                    

                    double f = analytics2.get_feature(at);

                    double g = analytics2.get_curved_area(at);


                    Vector_3 gg = analytics2.get_normal_vec(at);


                    double test[3] = {gg.x(),gg.y(),gg.z()};


                    /*double test[3] = {0.0,0.0,0.0};

                    if(isnan(gg.x())){
                        test[0]=100.0;
                        std::cout << gg.x() << "  " << gg.y() << "  " << gg.z() << std::endl;
                    }else{
                        test[0]=(double)gg.x();
                    }

                    if(isnan(gg.y())){
                        test[1]=100.0;
                        std::cout << gg.x() << "  " << gg.y() << "  " << gg.z() << std::endl;
                    }else{
                        test[1]=(double)gg.y();
                    }

                    if(isnan(gg.z())){
                        test[2]=100.0;
                        std::cout << gg.x() << "  " << gg.y() << "  " << gg.z() << std::endl;
                    }else{
                        test[2]=(double)gg.z();
                    }*/



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

                        viennagrid_quantity_field_value_set(quantity_field5,id,
                            &(f));

                        viennagrid_quantity_field_value_set(quantity_field6,id,
                            &(g));

                        viennagrid_quantity_field_value_set(quantity_field7,id,&(test));


                    id++;
                }   


                // add current quantity to the vector
                viennagrid_quantity_field_name_set(quantity_field,"PC1");
                quantity_fields.push_back(quantity_field);

                viennagrid_quantity_field_name_set(quantity_field1,"PC2");
                quantity_fields.push_back(quantity_field1);

                viennagrid_quantity_field_name_set(quantity_field2,"mean");
                quantity_fields.push_back(quantity_field2);

                viennagrid_quantity_field_name_set(quantity_field3,"gauss");
                quantity_fields.push_back(quantity_field3);

                viennagrid_quantity_field_name_set(quantity_field4,"area_around_vertex");
                quantity_fields.push_back(quantity_field4);

                viennagrid_quantity_field_name_set(quantity_field5,"features");
                quantity_fields.push_back(quantity_field5);

                viennagrid_quantity_field_name_set(quantity_field6,"curved area");
                quantity_fields.push_back(quantity_field6);

                viennagrid_quantity_field_name_set(quantity_field7,"normal_vector");
                quantity_fields.push_back(quantity_field7);

                

                //export quantities
                data_handle<viennagrid_quantity_field> quantity_field_handle=viennamesh::plugin_algorithm::make_data(to_cpp(quantity_fields[0]));  
                quantity_field_handle.push_back(to_cpp(quantity_fields[1]));    
                quantity_field_handle.push_back(to_cpp(quantity_fields[2])); 
                quantity_field_handle.push_back(to_cpp(quantity_fields[3]));
                quantity_field_handle.push_back(to_cpp(quantity_fields[4]));
                quantity_field_handle.push_back(to_cpp(quantity_fields[5]));  
                quantity_field_handle.push_back(to_cpp(quantity_fields[6]));
                quantity_field_handle.push_back(to_cpp(quantity_fields[7]));          
                set_output("quantities",quantity_field_handle);

#pragma endregion

            }


        

            //___________________________________Testing_writing data into VTU file________________________________________________


            //first new curvature tests

            int h = 0;
            for(cgal::polyhedron_surface_mesh::Vertex_iterator at=my_mesh.vertices_begin(),end=my_mesh.vertices_end();at!=end;++at)
            {

               
                if(fabs((*at).point().x() +0.53125) <= 0.00001 &&
                   fabs((*at).point().y() - 0.5) <= 0.00001 &&
                   fabs((*at).point().z() + 0.00625) <= 0.00001) {   
                    
                    
                    cgal::polyhedron_surface_mesh::Vertex::Halfedge_around_vertex_circulator at1=(*at).vertex_begin(), end1 = at1;

                    std::vector<Point_3> points;

                    std::vector<Point_3> points2;

                    std::cout << "GAGAGAGAGAGAGGAGA";

                    Plane_3 test_plane1, test_plane2;

                    do{
                        if(analytics.get_feature(at1->opposite()->vertex()) == -1){
                            points.push_back(at->point());
                        }

                        if (points.size() == 2){

                            test_plane1 = Plane_3(points[0], points[1], (*at).point());

                            break;
                        }

                        
                        at1++;
                    }while(at1 != end1);

                    at1=(*at).vertex_begin();
                    end1 = at1;

                    do{
                        if(analytics.get_feature(at1->opposite()->vertex()) == 35){
                            points2.push_back(at->point());
                        }

                        if (points2.size() == 2){

                            test_plane2 = Plane_3(points2[0], points2[1], (*at).point());

                            break;
                        }

                        
                        at1++;
                    }while(at1 != end1);
                    
                }

                h++;
    
            }        
    

            info(1) << "exporting the mesh..." << std::endl;


            set_output("mesh", my_mesh);


            return true;
 
        }
      
    }
}
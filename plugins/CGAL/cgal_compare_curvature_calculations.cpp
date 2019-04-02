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

#include "cgal_compare_curvature_calculations.hpp"

#include "CGAL/cgal_curvature_functions.hpp"

//#include "curve_calc_dev/curvature_calculator.hpp"

#include "curvature_calculator.hpp"

#include<chrono>

#include <boost/algorithm/string.hpp>  // sto_lower(string str)
#include <typeinfo>   // operator typeid



//hash table to store curvatures
#include <unordered_map>

#include <fstream>







namespace viennamesh
{
    namespace cgal
    {

        cgal_compare_curvature_calculations::cgal_compare_curvature_calculations(){}

        std::string cgal_compare_curvature_calculations::name()
        {
           
            return "cgal_compare_curvature_calculations";
        }

        bool cgal_compare_curvature_calculations::run(viennamesh::algorithm_handle &){

            data_handle<viennagrid_numeric> number_of_loops = get_input<viennagrid_numeric>("number_of_loops");

            data_handle<viennagrid_numeric> number_of_points = get_input<viennagrid_numeric>("number_of_points");

            data_handle<bool> create_outputmesh = get_input<bool>("outputmesh");
 
            
            //Get mesh data (REQUIRED)
            data_handle<cgal::polyhedron_surface_mesh> input_mesh = get_required_input<cgal::polyhedron_surface_mesh>("mesh");


            //Declare internal reference to output mesh
            cgal::polyhedron_surface_mesh & my_mesh = const_cast<cgal::polyhedron_surface_mesh&> (input_mesh());
            my_mesh = input_mesh();


            auto start = std::chrono::high_resolution_clock::now();

            auto finish = std::chrono::high_resolution_clock::now();

            double curvatures[2];

            double summed_up_runtimes_my=0.0;

            double summed_up_runtimes_cgal=0.0;

            info(1) << "Start my calculation"  << std::endl; 

            for(int i = 0; i < number_of_loops(); i++){


                //calculate all curvatures my
                start = std::chrono::high_resolution_clock::now();

                for(cgal::polyhedron_surface_mesh::Vertex_iterator at=my_mesh.vertices_begin(),end=my_mesh.vertices_end();at!=end;++at)
                {

                    cgal::polyhedron_surface_mesh::Vertex_handle t = at;
                    Curves c = calc_curvatures(*at);
                    

                }

                finish = std::chrono::high_resolution_clock::now();

                summed_up_runtimes_my += std::chrono::duration_cast<std::chrono::microseconds>(finish-start).count();

                
            }

            info(1) << "start CGAL calculation"  << std::endl; 

            for(int i = 0; i < number_of_loops(); i++){

                //calculate all curvatures cgal

                start = std::chrono::high_resolution_clock::now();

            char cwd[PATH_MAX];
            std::string output_filename_my = getcwd(cwd, sizeof(cwd));
            std::string output_filename_cgal = getcwd(cwd, sizeof(cwd));
            output_filename_my += "/distance.csv";
            output_filename_cgal += "/times_cgal.csv";

            std::ofstream csv_my, csv_cgal;
            csv_my.open(output_filename_my.c_str(),  std::ios::app);

            std::stringstream out1, out2, out3, out4, out5;

            bool first = true;

 

                for(cgal::polyhedron_surface_mesh::Vertex_iterator at=my_mesh.vertices_begin(),end=my_mesh.vertices_end();at!=end;++at)
                {
                    if(!first){
                        out1 << ", ";
                        out2 << ", ";
                        out3 << ", ";
                        out4 << ", "; 
                        out5 << ", ";   
                    }

                    double test_c1[2];

                    
                    principal_curvatures_cgal_time_test(*at, my_mesh, test_c1, 6);

                    /*double c = (test_c1[0]+test_c1[1])/2.0;
                    
                    out1 << std::to_string(c);

                    
                    principal_curvatures_cgal_time_test(*at, my_mesh, test_c1, 18);

                    c = (test_c1[0]+test_c1[1])/2.0;
                    
                    out2 << std::to_string(c);

                    
                    principal_curvatures_cgal_time_test(*at, my_mesh, test_c1, 35);

                    c = (test_c1[0]+test_c1[1])/2.0;
                    
                    out3 << std::to_string(c);

                    
                    principal_curvatures_cgal_time_test(*at, my_mesh, test_c1, 40);

                    c = (test_c1[0]+test_c1[1])/2.0;
                    
                    out4 << std::to_string(c);

                    Curves cu = calc_curvatures(*at);

                    out5 << std::to_string(cu.mean);







                    first = false;*/

                }

                csv_my << out1.str() << std::endl;

                csv_my << out2.str() << std::endl;

                csv_my << out3.str() << std::endl;

                csv_my << out4.str() << std::endl;

                csv_my << out5.str();

                csv_my.close();

                finish = std::chrono::high_resolution_clock::now();

                summed_up_runtimes_cgal += std::chrono::duration_cast<std::chrono::microseconds>(finish-start).count();

                


            }
            info(1) << "number of loops:   " << number_of_loops()  << std::endl; 
            info(1) << "runtime my curvature:   " << summed_up_runtimes_my << " mus" << std::endl; 
            info(1) << "runtime cgal curvature " << summed_up_runtimes_cgal << " mus" << std::endl; 




            /*std::ofstream csv_my, csv_cgal;
            csv_my.open(output_filename_my.c_str(),  std::ios::app);
            csv_cgal.open(output_filename_cgal.c_str(),  std::ios::app);

            csv_my << summed_up_runtimes_my/(double)number_of_loops() << ", "; 
            csv_cgal << summed_up_runtimes_cgal/(double)number_of_loops() << ", ";*/

            //close csv file
            //csv_my.close();
            //csv_cgal.close();



    
            if(create_outputmesh() == true){

                info(1) << "creating output mesh..." << std::endl; 


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
                                1,                                        // one float per element
                                VIENNAGRID_QUANTITY_FIELD_STORAGE_SPARSE);


                long id=0;

                for(cgal::polyhedron_surface_mesh::Vertex_iterator at=my_mesh.vertices_begin(),end=my_mesh.vertices_end();at!=end;++at)
                    {

                            Curves c = calc_curvatures(*at);

                            principal_curvatures_cgal_time_test(*at, my_mesh, curvatures,number_of_points());


                            viennagrid_quantity_field_value_set(quantity_field,id,&(c.p1));

                            viennagrid_quantity_field_value_set(quantity_field1,id,&(c.p2));

                            viennagrid_quantity_field_value_set(quantity_field2,id,&(c.mean));

                            viennagrid_quantity_field_value_set(quantity_field3,id,&(c.gauss));

                            viennagrid_quantity_field_value_set(quantity_field4,id,&(curvatures[0]));

                            viennagrid_quantity_field_value_set(quantity_field5,id,&(curvatures[1]));

                            double a = (double)((curvatures[0]+curvatures[1])/2.0);

                            viennagrid_quantity_field_value_set(quantity_field6,id,&a);

                            double b = curvatures[0]*curvatures[1];

                            viennagrid_quantity_field_value_set(quantity_field7,id,&b);


                        id++;
                    }


                    // add current quantity to the vector
                    viennagrid_quantity_field_name_set(quantity_field,"my PC1");
                    quantity_fields.push_back(quantity_field);

                    viennagrid_quantity_field_name_set(quantity_field1,"my PC2");
                    quantity_fields.push_back(quantity_field1);

                    viennagrid_quantity_field_name_set(quantity_field2,"my mean");
                    quantity_fields.push_back(quantity_field2);

                    viennagrid_quantity_field_name_set(quantity_field3,"my gauss");
                    quantity_fields.push_back(quantity_field3);

                    viennagrid_quantity_field_name_set(quantity_field4,"cgal PC1");
                    quantity_fields.push_back(quantity_field4);

                    viennagrid_quantity_field_name_set(quantity_field5,"cgal PC2");
                    quantity_fields.push_back(quantity_field5);

                    viennagrid_quantity_field_name_set(quantity_field6,"cgal mean");
                    quantity_fields.push_back(quantity_field6);
                    
                    viennagrid_quantity_field_name_set(quantity_field7,"cgal gauss");
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
                    set_output("mesh", my_mesh);
                }

      

            return true;
 
        }
       
    }
}
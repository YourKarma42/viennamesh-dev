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

#include "curvature_calculator.hpp"


namespace viennamesh
{

        cgal_curve::cgal_curve(){}

        std::string cgal_curve::name()
        {
           
            return "cgal_curve";
        }

        typedef typename cgal::polyhedron_surface_mesh::Vertex_handle Vertex_handle;

        std::unordered_map<Vertex_handle, cgal::Curves> _curvatures_scaled;

        std::unordered_map<Vertex_handle, cgal::Curves> _curvatures_original;

        //TODO: remove
        

        bool cgal_curve::run(viennamesh::algorithm_handle &){

                data_handle<cgal::polyhedron_surface_mesh> input_mesh = get_required_input<cgal::polyhedron_surface_mesh>("mesh");

                //Create output mesh handle
                data_handle<cgal::polyhedron_surface_mesh> output_mesh = make_data<cgal::polyhedron_surface_mesh>();

                //Declare internal reference to output mesh
                cgal::polyhedron_surface_mesh & my_mesh = const_cast<cgal::polyhedron_surface_mesh&> (output_mesh());

                //viennamesh::cgal::scale_test =1.0f;

                my_mesh = input_mesh();

                for(cgal::polyhedron_surface_mesh::Vertex_iterator at=my_mesh.vertices_begin(),end=my_mesh.vertices_end();at!=end;++at)
                {
                    _curvatures_original[at] = cgal::calc_curvatures(*at);

                }

                data_handle<cgal::polyhedron_surface_mesh> input_mesh2 = get_required_input<cgal::polyhedron_surface_mesh>("mesh");

                data_handle<cgal::polyhedron_surface_mesh> output_mesh2 = make_data<cgal::polyhedron_surface_mesh>();

                //Declare internal reference to output mesh
                cgal::polyhedron_surface_mesh & my_mesh2 = const_cast<cgal::polyhedron_surface_mesh&> (output_mesh2());

                //viennamesh::cgal::scale_test =1.0f;

                my_mesh2 = input_mesh2();


                for(cgal::polyhedron_surface_mesh::Vertex_iterator at=my_mesh2.vertices_begin(),end=my_mesh2.vertices_end();at!=end;++at)
                {
                    _curvatures_scaled[at] = cgal::calc_curvatures(*at);

                }

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
                

                long id=0;
                cgal::polyhedron_surface_mesh::Vertex_iterator itorig   = my_mesh.vertices_begin();
                cgal::polyhedron_surface_mesh::Vertex_iterator itscaled = my_mesh2.vertices_begin();           


                for( ;itorig != my_mesh.vertices_end(); )
                {

                    double a = _curvatures_original[itorig].mean;

                    double b = _curvatures_scaled[itscaled].mean;
                  

                    double c = std::abs(a)-std::abs(b*10000);

      

                        viennagrid_quantity_field_value_set(quantity_field,id,&(a));

                        viennagrid_quantity_field_value_set(quantity_field1,id,&(b));

                        viennagrid_quantity_field_value_set(quantity_field2,id,&(c));


                    id++;
                    itorig++;
                    itscaled++;
                }   




        // add current quantity to the vector
        viennagrid_quantity_field_name_set(quantity_field,"orig");
        quantity_fields.push_back(quantity_field);

        viennagrid_quantity_field_name_set(quantity_field1,"scaled");
        quantity_fields.push_back(quantity_field1);

        viennagrid_quantity_field_name_set(quantity_field2,"diff");
        quantity_fields.push_back(quantity_field2);


                

        //export quantities
        data_handle<viennagrid_quantity_field> quantity_field_handle=viennamesh::plugin_algorithm::make_data(to_cpp(quantity_fields[0]));  
        quantity_field_handle.push_back(to_cpp(quantity_fields[1]));    
        quantity_field_handle.push_back(to_cpp(quantity_fields[2]));  

        set_output("quantities",quantity_field_handle);
                   
        set_output("mesh", my_mesh);

        }


}
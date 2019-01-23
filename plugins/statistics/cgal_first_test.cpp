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


#include "cgal_first_test.hpp"

#include "libigl_convert.hpp"


//libigl includes



namespace viennamesh
{

        cgal_first_test::cgal_first_test(){}

        std::string cgal_first_test::name()
        {
           
            return "cgal_first_test";
        }





        bool cgal_first_test::run(viennamesh::algorithm_handle &){
            
            //Get mesh data
            
            mesh_handle input_mesh = get_required_input<mesh_handle>("mesh");

            mesh_handle corsed_mesh = get_required_input<mesh_handle>("corsed_mesh");

            data_handle<viennagrid_numeric> numb = get_required_input<viennagrid_numeric>("num");

            data_handle<viennamesh_string> input = get_required_input<viennamesh_string>("test");

            //viennamesh::info(1) << input_mesh_cgal()[0] << std::endl;

            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Vertices;
            Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> Facets;

            convert_to_igl_mesh(input_mesh(), Vertices, Facets);

            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> VerticesS;
            Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> FacetsS;

            convert_to_igl_mesh(corsed_mesh(), VerticesS, FacetsS);

            

            viennamesh::info(1) << "Not Simplyfied:"  << std::endl;

            viennamesh::info(1) << "rows: " << Vertices.rows() << std::endl;

            viennamesh::info(1) << "cols: " << Vertices.cols() << std::endl;

            viennamesh::info(1) << "Simplyfied:"  << std::endl;

            viennamesh::info(1) << "rows: " << VerticesS.rows() << std::endl;

            viennamesh::info(1) << "cols: " << VerticesS.cols() << std::endl;



            Eigen::Matrix<double, Eigen::Dynamic, 1> curvature;

            igl::gaussian_curvature(Vertices, Facets, curvature);

            Eigen::Matrix<double, Eigen::Dynamic, 1> curvatureS;

            igl::gaussian_curvature(VerticesS, FacetsS, curvatureS);


            for(int i =0; i < Vertices.rows(); i++){
                for(int j=0; j< VerticesS.rows(); j++){

                    if(Vertices(i,0)==VerticesS(j,0) && Vertices(i,1)==VerticesS(j,1) && Vertices(i,2)==VerticesS(j,2)){
                        viennamesh::info(1) << "curve:" << curvature(i, 0) << std::endl;
                        viennamesh::info(1) << "curveS:" << curvatureS(j, 0) << std::endl;
                    }


                }
            }


            //viennamesh::info(1) << Vertices(5,0) << std::endl;

            //viennamesh::info(1) << input() << std::endl;

            //viennamesh::info(1) << numb() << std::endl;

            //viennamesh::info(1) << "echt gangstaaaaa cout" << std::endl;


 
        }

        void cgal_first_test::curvature(mesh_handle & mesh){

            ///mesh_handle mesh = get_required_input<mesh_handle>("mesh");

           

            //Comparison metrics are implemented using libigl, which uses matrices provided by Eigen library.

            
        }





}
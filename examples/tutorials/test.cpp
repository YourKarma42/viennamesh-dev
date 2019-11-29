#include "viennameshpp/core.hpp"


//testclass

int main(int argc, char **argv, char **argv2)
{
    
    if (argc > 3)
    {
        std::cout << "1 argument required: 1 meshe" << std::endl;
        return 0;
    }

    std::string meshPath(argv[1]);

    //std::string meshPath2(argv[2]);

    //squared edge length


    viennamesh::context_handle context;

    viennamesh::algorithm_handle mesh_reader = context.make_algorithm("mesh_reader");
    mesh_reader.set_input( "filename", meshPath );
    mesh_reader.run();


    viennamesh::algorithm_handle test = context.make_algorithm("cgal_curve");
    test.set_default_source(mesh_reader);
    test.run();


    std::string file_name = std::string(basename(argv[1]));

    viennamesh::algorithm_handle mesh_writer = context.make_algorithm("mesh_writer");
    mesh_writer.set_default_source(test);
    mesh_writer.set_input( "filename", file_name + "_compare.vtu" );
    mesh_writer.run();



//IMPORTANT QUALITY METICS AND STUFF
/*
    viennamesh::algorithm_handle test = context.make_algorithm("vtk_mesh_distance");
    test.set_default_source(mesh_reader);
    test.set_input("mesh1", meshPath);
    test.set_input("mesh2", meshPath2);
    //test.set_input("number_of_loops", number_of_loops);
    test.run();

    viennamesh::algorithm_handle test1 = context.make_algorithm("vtk_calculate_mesh_quality_metrics");
    test1.set_default_source(mesh_reader);
    test1.set_input("mesh_path", meshPath2);
    //test.set_input("number_of_loops", number_of_loops);
    test1.run();


*/
    




    //viennamesh reader
   /* viennamesh::algorithm_handle mesh_reader = context.make_algorithm("mesh_reader");
    mesh_reader.set_input( "filename", meshPath );
    mesh_reader.run();

    viennamesh::algorithm_handle test = context.make_algorithm("cgal_compare_curvature_calculations");

    test.set_default_source(mesh_reader);
    test.set_input("number_of_loops", number_of_loops);
    test.run();

*/


   /* viennamesh::algorithm_handle test = context.make_algorithm("remove_short_edges");

    test.set_default_source(mesh_reader);
    test.set_input("min_edge_size", sqrt_edge_length);
    test.run();

     std::string file_name = std::string(basename(argv[1]));

    //the resulting coarsened mesh is written to a vtu file.
    viennamesh::algorithm_handle mesh_writer = context.make_algorithm("mesh_writer");
    mesh_writer.set_default_source(test);
    mesh_writer.set_input( "filename", file_name + "_short_removed.vtu" );
    mesh_writer.run();*/

}
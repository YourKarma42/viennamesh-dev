#include "viennameshpp/core.hpp"


int main(int argc, char **argv, char **argv2)
{
    
    if (argc > 5 || argc < 3)
    {
        std::cout << "3 argument required: 1 triangle mesh. Number of reptitions, cgal number of points " << std::endl;
        std::cout << "1 additonal argument if a output mesh should be created " << std::endl;
        return 0;
    }

    std::string meshPath(argv[1]);

    int number_of_loops = atoi(argv[2]);

    int number_of_points = atoi(argv[3]);

    bool output_mesh = false;

    if(argc == 5){
        output_mesh = true;
    }

    //squared edge length


    viennamesh::context_handle context;

    viennamesh::algorithm_handle mesh_reader = context.make_algorithm("mesh_reader");
    mesh_reader.set_input( "filename", meshPath );
    mesh_reader.run();

    viennamesh::algorithm_handle comperator = context.make_algorithm("cgal_compare_curvature_calculations");
    comperator.set_default_source(mesh_reader);
    comperator.set_input("number_of_loops", number_of_loops);
    comperator.set_input("number_of_points", number_of_points);
    comperator.set_input("outputmesh", output_mesh);
    comperator.run();

    std::cout << "number of points: " << number_of_points << std::endl;

    if(output_mesh){
        viennamesh::algorithm_handle mesh_writer = context.make_algorithm("mesh_writer");
        mesh_writer.set_default_source(comperator);
        mesh_writer.set_input( "filename", "comparison.vtu" );
        mesh_writer.run();
    }

}
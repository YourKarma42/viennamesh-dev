#include "viennameshpp/core.hpp"



//testclass

int main(int argc, char **argv, char **argv2)
{
    
    if (argc != 2)
    {
        std::cout << "1 argument required: Path to mesh Victory process mesh (Volume surface?)" << std::endl;
        return 0;
    }

    std::string meshPath(argv[1]);


    viennamesh::context_handle context;

    viennamesh::algorithm_handle mesh_reader = context.make_algorithm("mesh_reader");
    mesh_reader.set_input( "filename", meshPath );
    mesh_reader.run();

    viennamesh::algorithm_handle mesh_writer = context.make_algorithm("mesh_writer");
    mesh_writer.set_default_source(mesh_reader);
    mesh_writer.set_input( "filename", "output_test.vtu" );
    mesh_writer.run();

/*
    viennamesh::algorithm_handle distance_calculator = context.make_algorithm("vtk_mesh_distance");
    distance_calculator.set_default_source(mesh_reader);
    distance_calculator.set_input("mesh1", meshPath);
    distance_calculator.set_input("mesh2", meshPath2);
    distance_calculator.run();

*/


}
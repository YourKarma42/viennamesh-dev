#include "viennameshpp/core.hpp"


//testclass

int main(int argc, char **argv, char **argv2)
{
    
    if (argc != 3)
    {
        std::cout << "2 argument required: 2 Meshes. Mesh 1 not_simplified. Mesh 2 simplified" << std::endl;
        return 0;
    }

    std::string meshPath(argv[1]);

    std::string meshPath2(argv[2]);



    viennamesh::context_handle context;

    viennamesh::algorithm_handle mesh_reader = context.make_algorithm("mesh_reader");
    mesh_reader.set_input( "filename", meshPath );
    mesh_reader.run();


    viennamesh::algorithm_handle distance_calculator = context.make_algorithm("vtk_mesh_distance");
    distance_calculator.set_default_source(mesh_reader);
    distance_calculator.set_input("mesh1", meshPath);
    distance_calculator.set_input("mesh2", meshPath2);
    distance_calculator.run();




}
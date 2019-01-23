#include "viennameshpp/core.hpp"


//testclass

int main(int argc, char **argv, char **argv2)
{
    
    if (argc > 2)
    {
        std::cout << "1 argument required: 1 triangle mesh to calculate quality metrics of" << std::endl;
        return 0;
    }

    std::string meshPath(argv[1]);

    //squared edge length


    viennamesh::context_handle context;

    viennamesh::algorithm_handle mesh_reader = context.make_algorithm("mesh_reader");
    mesh_reader.set_input( "filename", meshPath );
    mesh_reader.run();

    viennamesh::algorithm_handle quality_calculator = context.make_algorithm("vtk_calculate_mesh_quality_metrics");
    quality_calculator.set_default_source(mesh_reader);
    quality_calculator.set_input("mesh_path", meshPath);
    //test.set_input("number_of_loops", number_of_loops);
    quality_calculator.run();


}
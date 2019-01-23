#include "viennameshpp/core.hpp"


//tests for curvature calculation in cgal

int main(int argc, char **argv)
{
    
    if (argc < 2)
    {
        std::cout << "at least 1 argument required: 1 .vtu file" << std::endl;
        return 0;
    }


    std::string vtkPath(argv[1]);

    viennamesh::context_handle context;


    viennamesh::algorithm_handle printer = context.make_algorithm("vtk_print_mesh_data");

    printer.set_input( "file_name", vtkPath );
    

    printer.run();

}
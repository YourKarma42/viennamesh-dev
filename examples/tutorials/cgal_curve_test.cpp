#include "viennameshpp/core.hpp"


//tests for curvature calculation in cgal

int main(int argc, char **argv, char **argv2)
{
    
    if (argc != 2)
    {
        std::cout << "1 argument required: 1 meshes to compare curvature" << std::endl;
        return 0;
    }

    std::string meshPath(argv[1]);

    viennamesh::context_handle context;



        //viennamesh reader
    viennamesh::algorithm_handle mesh_reader = context.make_algorithm("mesh_reader");
    mesh_reader.set_input( "filename", meshPath );
    mesh_reader.run();

    viennamesh::algorithm_handle test = context.make_algorithm("cgal_curve");

    test.set_default_source(mesh_reader);

    test.run();

}
#include "viennameshpp/core.hpp"


//tests for curvature calculation in cgal

int main(int argc, char **argv)
{
    
    if (argc < 2)
    {
        std::cout << "at least 1 argument required: 1 .vtu file to clean 2 tolerance for point merging" << std::endl;
        return 0;
    }

    double tolerance = 0.0;

    if (argc == 3){
        tolerance = atof(argv[2]);
    }

    std::string vtkPath(argv[1]);

    viennamesh::context_handle context;

    //viennamesh reader
    viennamesh::algorithm_handle mesh_reader = context.make_algorithm("mesh_reader");
    mesh_reader.set_input( "filename", vtkPath );
    mesh_reader.run();


    viennamesh::algorithm_handle cleaner = context.make_algorithm("vtk_clean");

    cleaner.set_default_source(mesh_reader);

    cleaner.set_input( "tolerance", tolerance );
    

    cleaner.run();


    std::string file_name = vtkPath;

    file_name = file_name + "_cleaned.vtu";


    viennamesh::algorithm_handle mesh_writer = context.make_algorithm("mesh_writer");
    mesh_writer.set_default_source(cleaner);
    mesh_writer.set_input( "filename", file_name );
    mesh_writer.run();
}
#include "viennameshpp/core.hpp"

#include <libgen.h>


//testclass

int main(int argc, char **argv)
{
    
    if (argc < 2)
    {
        std::cout << "1 argument required: 1 meshes to compare curvature" << std::endl;
        return 0;
    }

    std::string meshPath(argv[1]);

    double ratio = 0.5;

    std::string policy("m");
    

    if (argc > 2){
        ratio = atof(argv[2]);
    }
        

    if (argc > 3){
        std::string new_pol(argv[3]);
        if(new_pol == "m" || new_pol == "lt" || new_pol == "p"){
            policy = new_pol;
        }else{
            std::cout << "Wrong argument in 3: must be lt or m or p" << std::endl;
            std::cout << "Using standard m" << std::endl;
        }
    }


    if (argc > 4)
    {
        std::cout << "2 optional parameters: ratio as double | lindstom turk (lt) or my stuff (m)" << std::endl;
        return 0;
    }

    std::cout << "choosen ratio is: " << ratio << std::endl;

    std::cout << "choosen policy: " << policy << std::endl;


    viennamesh::context_handle context;


    //viennamesh reader
    viennamesh::algorithm_handle mesh_reader = context.make_algorithm("mesh_reader");
    mesh_reader.set_input( "filename", meshPath );
    mesh_reader.run();


    viennamesh::algorithm_handle test = context.make_algorithm("cgal_simplify_curve");

    test.set_default_source(mesh_reader);

    test.set_input("ratio", ratio);
    test.set_input("policy", policy);

    test.run();

    std::string file_name = std::string(basename(argv[1]));


    //the resulting coarsened mesh is written to a vtu file.
    viennamesh::algorithm_handle mesh_writer = context.make_algorithm("mesh_writer");
    mesh_writer.set_default_source(test);
    mesh_writer.set_input( "filename", file_name + "_" + std::to_string(ratio) + "_"+ policy + "_" + "_simplified.vtu" );
    mesh_writer.run();

}
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
        if(new_pol == "m" || new_pol == "lt" || new_pol == "p" || new_pol == "mN1" || new_pol == "mN2" || new_pol == "mC1" || new_pol == "fcp"){
            policy = new_pol;
        }else{
            std::cout << "Wrong argument in 3: must be lt or m or p" << std::endl;
            std::cout << "Using standard m" << std::endl;
        }
    }

    double transition_el=0;
    double flat_el=0;
    
    if (argc > 4){

      transition_el = atof(argv[4]);
      flat_el = atof(argv[5]);
    }

    if (argc > 6)
    {
        std::cout << "4 optional parameters: ratio as double | lindstom turk (lt) or my stuff (m)" << std::endl;
        std::cout << "transition area edge length as double | flat area edge length as double" << std::endl;
        return 0;
    }

    std::cout << "choosen ratio is: " << ratio << std::endl;

    std::cout << "choosen policy: " << policy << std::endl;

    std::cout << "choosen transition are edge length: " << transition_el << std::endl;

    std::cout << "choosen flat area edge lenth: " << flat_el << std::endl;

    viennamesh::context_handle context;


    //viennamesh reader
    viennamesh::algorithm_handle mesh_reader = context.make_algorithm("mesh_reader");
    mesh_reader.set_input( "filename", meshPath );
    mesh_reader.run();


    viennamesh::algorithm_handle test = context.make_algorithm("cgal_simplify_curve");

    test.set_default_source(mesh_reader);

    test.set_input("ratio", ratio);
    test.set_input("policy", policy);
    test.set_input("transition_el", transition_el);
    test.set_input("flat_el", flat_el); 

    test.run();

//komisch nochmal überlegen
    std::string file_name = std::string(basename(argv[1])).substr(0,std::string(basename(argv[1])).find_last_of("/"));



    //the resulting coarsened mesh is written to a vtu file.
    viennamesh::algorithm_handle mesh_writer = context.make_algorithm("mesh_writer");
    mesh_writer.set_default_source(test);
    mesh_writer.set_input( "filename", file_name + "_" + std::to_string(ratio) + "_"+ policy + "_" + "_simplified.vtu" );
    mesh_writer.run();

}
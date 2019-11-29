#include "viennameshpp/core.hpp"

#include <libgen.h>

#include <list>

#include <fstream>


//testclass

int main(int argc, char **argv)
{
    
    if (argc < 2)
    {
        std::cout << "at least 2 argument required:" << std::endl;
        std::cout << "for option p: 'path to mesh' 'p'" << std::endl;
        std::cout << "for option lt: 'path to mesh' 'lt' 'number of remaining triangles'" << std::endl;
        std::cout << "for option m: 'path to mesh' 'm' 'start length' 'end length' 'step size' 'flat boundary'" << std::endl;
        return 0;
    }

    std::string meshPath(argv[1]);

    std::string policy("m");

    if (argc > 2){
        std::string new_pol(argv[2]);
        if(new_pol == "m" || new_pol == "lt" || new_pol == "p"){
            policy = new_pol;
        }else{
            std::cout << "Wrong argument in 3: must be lt or m or p" << std::endl;
            return 0;
        }
              
    }
        

       
    double arg1 = 0.0;

    if (argc > 3){
        arg1 = atof(argv[3]);
    }

    double end_length = 0.0;
    double step_size = 0.0;
    double flat_boundary = 0.0;

    if (argc > 4){
        end_length = atof(argv[4]);
        step_size = atof(argv[5]);
        flat_boundary = atof(argv[6]);
    }


    std::cout << "choosen policy       : " << policy << std::endl;
    std::cout << "choosen start length : " << arg1 << std::endl;
    std::cout << "choosen end length   : " << end_length << std::endl;
    std::cout << "choosen step length  : " << step_size << std::endl;
    std::cout << "choosen flat boundary: " << flat_boundary << std::endl;

    viennamesh::context_handle context;

    //viennamesh reader
    viennamesh::algorithm_handle mesh_reader = context.make_algorithm("mesh_reader");
    mesh_reader.set_input( "filename", meshPath );
    mesh_reader.run();


    viennamesh::algorithm_handle simplification = context.make_algorithm("cgal_simplify_curve");

    simplification.set_default_source(mesh_reader);


    simplification.set_input("ratio", arg1);
    simplification.set_input("policy", policy);
    simplification.set_input("end_length", end_length);
    simplification.set_input("step_size", step_size);
    simplification.set_input("flat_boudary", flat_boundary);
    

    simplification.run();

    
//komisch nochmal überlegen
    std::stringstream file_name_tmp;
    file_name_tmp << std::string(basename(argv[1])).substr(0,std::string(basename(argv[1])).find_last_of("/"));

    std::string segment;
    std::vector<std::string> seglist;

    while(std::getline(file_name_tmp, segment, '_'))
    {
        seglist.push_back(segment);
    }

    std::stringstream file_name;
    file_name << seglist[0] << "_" << seglist[1];

    std::cout << file_name.str() << std::endl;




    //the resulting coarsened mesh is written to a vtu file.
    viennamesh::algorithm_handle mesh_writer = context.make_algorithm("mesh_writer");
    mesh_writer.set_default_source(simplification);
    mesh_writer.set_input( "filename", file_name.str() + "_" + policy + "_simpl.vtu" );
    mesh_writer.run();

}

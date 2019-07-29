#include "viennameshpp/core.hpp"

#include <libgen.h>

#include <list>

#include <fstream>


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

    int anz_zones = 0;

    double zone_size_mult = 0.0;

    std::list<double>  zones;
    
    if (argc > 4){

      anz_zones = atoi(argv[4]);

      zone_size_mult = atof(argv[5]);

      for(int i = 0; i<anz_zones; i++){
          zones.push_back(atof(argv[6+i]));
      }

    

    }

    /*if (argc > 6)
    {
        std::cout << "4 optional parameters: ratio as double | lindstom turk (lt) or my stuff (m)" << std::endl;
        std::cout << "transition area edge length as double | flat area edge length as double" << std::endl;
        return 0;
    }*/



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
    test.set_input("zone_multiplicator", zone_size_mult);

    test.set_input("anz_zones", anz_zones);

    int i = 0;

    for(auto z: zones){

        std::string zone = "zone" + std::to_string(i);

        test.set_input(zone, z);

        i++;
          
    }

    

    test.run();

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

    std::string output_filename;
    output_filename += file_name.str();
    output_filename += "_params.txt";

    //create csv file handle
    std::ofstream output;
    output.open(output_filename.c_str());

        output << policy << std::endl;
        output << ratio << std::endl;

        if(policy == "m"){

            output << "Zone mult: "<< zone_size_mult << std::endl;

            output << "Zone sizes:" << std::endl;
            for(auto n: zones){
                output << n << std::endl;
            }
        }

    output.close();




    //the resulting coarsened mesh is written to a vtu file.
    viennamesh::algorithm_handle mesh_writer = context.make_algorithm("mesh_writer");
    mesh_writer.set_default_source(test);
    mesh_writer.set_input( "filename", file_name.str() + "_" + policy + "_simpl.vtu" );
    mesh_writer.run();

}

#include "vtk_print_mesh_data.hpp"
#include "vtk_mesh.hpp"

#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <vtkGenericDataObjectReader.h>


#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkCleanPolyData.h>

#include <vtkUnstructuredGrid.h>
#include <vtkCellArray.h>

#include <set>


namespace viennamesh {

    namespace vtk {

        print_Mesh_data::print_Mesh_data() {}

        bool print_Mesh_data::run(viennamesh::algorithm_handle &) {

            data_handle<viennamesh_string> file_name = get_required_input<viennamesh_string>("file_name");


            //read all the data from the file



            vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
            vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
            reader->SetFileName(file_name().c_str());
            reader->Update();

            vtkUnstructuredGrid* ugrid = reader->GetOutput();


            info(5) << "Number of Cells: " << ugrid->GetNumberOfCells() << std::endl;


            std::vector<std::set<int>> faces_around_vertex(ugrid->GetNumberOfPoints());


            int numCells = ugrid->GetNumberOfCells();

            int numPoints = ugrid->GetNumberOfPoints();

            //numPoints

            for(int i = 0; i < numPoints ; i++){

                if(i % 100 == 0){

                    info(5) << "Point Number " << i << std::endl;
                }

                for (int j = 0; j < numCells; j++){


                    //ugrid->GetCell(i)->GetNumberOfPoints();

                    //info(5) << "Points " << ugrid->GetCell(i)->GetNumberOfPoints() << std::endl;

                    vtkIdList* cell_points = ugrid->GetCell(j)->GetPointIds();

                    if(cell_points->GetId(0) == i){
                        faces_around_vertex[i].insert(j);

                    }

                    if(cell_points->GetId(1) == i){
                        faces_around_vertex[i].insert(j);

                    }

                    if(cell_points->GetId(2) == i){
                        faces_around_vertex[i].insert(j);

                    }





                    //info(5) << "Points " << test->GetId(0) << std::endl;



    

                }
            }

            std::ofstream myfile;
            myfile.open ("Triangles_around_each_Vertex.txt");

            if (myfile.is_open()){  

                for(int i = 0; i < faces_around_vertex.size(); i++){
                    myfile << i << " ";
                    for (auto face_id : faces_around_vertex[i]){
                        myfile << face_id << " ";
                    }

                    myfile << std::endl;
                }     
              
              myfile.close();
            }

            std::ofstream myfile2;
            myfile2.open ("Vertex_IDs_Points.txt");

            if (myfile2.is_open()){  

                 for (int j = 0; j < numPoints; j++){

                    myfile2 << j << " ";
                     

                    double p[3];

                        
                    ugrid->GetPoint(j,p);


                    myfile2 << p[0] << " "<< p[1] << " "<< p[2] << " ";
             

                    myfile2 << std::endl;

                }
              
              myfile2.close();
            }
        

            std::ofstream myfile3;
            myfile3.open ("Triangle_IDs_Vertices.txt");

            if (myfile3.is_open()){  

                 for (int j = 0; j < numCells; j++){

                    myfile3 << j << " ";

                    vtkPoints* cell_points = ugrid->GetCell(j)->GetPoints();

                    for(int k = 0; k<3; k++){                   

                        double p[3];

                            
                        cell_points->GetPoint(k,p);


                        myfile3 << p[0] << " "<< p[1] << " "<< p[2] << " ";
                    }

                    myfile3 << std::endl;

                }
              
              myfile3.close();
            }



 
            return true;

        }

    }

}
#include "vtk_mesh_distance.hpp"
#include "vtk_mesh.hpp"

#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPolyDataReader.h>

#include <vtkPolyDataMapper.h>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkCleanPolyData.h>

#include <vtkDistancePolyDataFilter.h>

#include <vtkUnstructuredGrid.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>

#include <vtkGeometryFilter.h>
#include <vtkAppendFilter.h>
#include <vtkPolyDataReader.h>
#include <vtkDoubleArray.h>



#include <vtkPointData.h>



namespace viennamesh {

    namespace vtk {

        mesh_distance::mesh_distance() {}

        bool mesh_distance::run(viennamesh::algorithm_handle &) {

            data_handle<vtk::mesh> is_mesh = get_required_input<vtk::mesh>("mesh");

            data_handle<viennamesh_string> meshPath1 = get_required_input<viennamesh_string>("mesh1"); 

            data_handle<viennamesh_string> meshPath2 = get_required_input<viennamesh_string>("mesh2");  

            //std::string a = meshPath1();

            std::string filename = meshPath1();
            //std::string BaseOutputPath = meshPath2().substr(0,meshPath2().find_last_of("/"));
            std::string BaseOutputPath = "";
        
            //store all quality values in a vector
            std::vector<double> qualities;
            std::vector<double> tmp_data;
            
            //create output csv-file for storing the mesh quality data
            //construct filename for csv and vtu file
            size_t found = (filename).find_last_of("/");
            size_t find_vtu = (filename).find_last_of(".");

            std::string output_filename_csv = BaseOutputPath;
            output_filename_csv += filename.substr(found+1, find_vtu-found-1);;
            output_filename_csv += "_distance.csv";

            std::string output_filename_vtu = "Distances.vtu";

            //std::string output_filename_vtu = BaseOutputPath;
            //output_filename_vtu += filename.substr(found+1, find_vtu-found-1);
            //output_filename_vtu += "_distance.vtu";

            //create csv file handle
            ofstream csv;
            csv.open(output_filename_csv.c_str());



            vtkSmartPointer<vtkPolyData> input1;
            vtkSmartPointer<vtkPolyData> input2;

            vtkSmartPointer<vtkXMLUnstructuredGridReader> reader1 =
            vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
            reader1->SetFileName(meshPath1().c_str());
            reader1->Update();

            vtkSmartPointer<vtkGeometryFilter> geometryFilter = 
            vtkSmartPointer<vtkGeometryFilter>::New();
            #if VTK_MAJOR_VERSION <= 5
                geometryFilter->SetInput(reader1->GetOutput());
            #else
                geometryFilter->SetInputData(reader1->GetOutput());
            #endif
            geometryFilter->Update(); 




            input1 = geometryFilter->GetOutput();

            vtkSmartPointer<vtkXMLUnstructuredGridReader> reader2 =
            vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
            reader2->SetFileName(meshPath2().c_str());
            reader2->Update();

            vtkSmartPointer<vtkGeometryFilter> geometryFilter2 = 
            vtkSmartPointer<vtkGeometryFilter>::New();
            #if VTK_MAJOR_VERSION <= 5
                geometryFilter2->SetInput(reader2->GetOutput());
            #else
                geometryFilter2->SetInputData(reader2->GetOutput());
            #endif
            geometryFilter2->Update(); 

            input2 = geometryFilter2->GetOutput();



            vtkSmartPointer<vtkCleanPolyData> clean1 =
                vtkSmartPointer<vtkCleanPolyData>::New();
            #if VTK_MAJOR_VERSION <= 5
            clean1->SetInputConnection( input1->GetProducerPort());
            #else
            clean1->SetInputData( input1);
            #endif

            vtkSmartPointer<vtkCleanPolyData> clean2 =
                vtkSmartPointer<vtkCleanPolyData>::New();
            #if VTK_MAJOR_VERSION <= 5
            clean2->SetInputConnection( input2->GetProducerPort());
            #else
            clean2->SetInputData( input2);
            #endif


            std::cout << "number of points Mesh 1: " << (reader1->GetOutput()->GetNumberOfPoints()) << std::endl;


            std::cout << "number of points Mesh 2: "  << (reader2->GetOutput()->GetNumberOfPoints()) << std::endl;

            vtkSmartPointer<vtkDistancePolyDataFilter> distanceFilter =
            vtkSmartPointer<vtkDistancePolyDataFilter>::New();

            distanceFilter->SetInputConnection( 0, clean1->GetOutputPort() );
            distanceFilter->SetInputConnection( 1, clean2->GetOutputPort() );
            std::cout << "calculating distance..." << std::endl;
            distanceFilter->Update();


            vtkDataSet* distancesToMesh = distanceFilter->GetOutput();
            vtkSmartPointer<vtkDoubleArray> qualityArray = vtkDoubleArray::SafeDownCast(distancesToMesh->GetPointData()->GetArray("Distance"));
            qualityArray->SetName("distance"); 
            reader1->GetOutput()->GetPointData()->AddArray(qualityArray);

            /*tmp_data.resize(reader1->GetNumberOfCells());    
            for (size_t i = 0; i < reader1->GetNumberOfCells(); ++i)
            {
                tmp_data[i] = *qualityArray->GetTuple(i);
            }

            qualities.push_back(tmp_data);
            tmp_data.clear();*/


            std::cout << (double)(distanceFilter->GetOutput()->GetNumberOfPoints()) << std::endl;

            //double test = vtkDoubleArray::SafeDownCast(distanceFilter->GetOutput()->GetPointData()->GetScalars())->GetValue(0);


            double sum_dist = 0.0;
            double sum_dist_abs = 0.0;

            for(int i = 0; i <= distanceFilter->GetOutput()->GetNumberOfPoints(); i++){
               sum_dist = sum_dist + vtkDoubleArray::SafeDownCast(distanceFilter->GetOutput()->GetPointData()->GetScalars())->GetValue(i);
               sum_dist_abs = sum_dist_abs + 
                        fabs(vtkDoubleArray::SafeDownCast(distanceFilter->GetOutput()->GetPointData()->GetScalars())->GetValue(i));

                qualities.push_back(vtkDoubleArray::SafeDownCast(distanceFilter->GetOutput()->GetPointData()->GetScalars())->GetValue(i));

               //probably write data into file to make statistic
            }

            std::cout << "Points in distance filter: " << distanceFilter->GetOutput()->GetNumberOfPoints() << std::endl;


            std::cout << "summed up distance: " << sum_dist << std::endl;

            std::cout << "summed up abs distance: " << sum_dist_abs << std::endl;

            std::cout << "avarage abs distance: " << sum_dist_abs/(double)(distanceFilter->GetOutput()->GetNumberOfPoints()) << std::endl;



            std::cout << "max dist: " << (double)(distanceFilter->GetOutput()->GetPointData()->GetScalars()->GetRange()[0]) << std::endl;


            std::cout << "min dist: " <<(double)(distanceFilter->GetOutput()->GetPointData()->GetScalars()->GetRange()[1]) << std::endl;







            std::cout << "writing didstances into .csv file" << std::endl;

            csv << "distance" << endl;

            std::cout << "Writing csv file " << output_filename_csv << std::endl;

           for(size_t vid = 0; vid < reader1->GetNumberOfCells(); ++vid)
            {
                csv << qualities[vid];


                //csv << ",";
                    
        
                csv << std::endl;
            }

            //close csv file
            csv.close();

            // Write file
            std::cout << "Writing vtu file ";
            std::cout << output_filename_vtu << std::endl;
            
            vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
            writer->SetFileName(output_filename_vtu.c_str());
            //writer->SetInputData(mq->GetOutput());
            writer->SetInputData(reader1->GetOutput());
            writer->SetDataModeToAscii();
            writer->Write();








 
            return true;

        }

    }

}


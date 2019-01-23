#include "vtk_calculate_mesh_quality_metrics.hpp"
#include "vtk_mesh.hpp"

#include <vtkSmartPointer.h>

#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkMeshQuality.h>
#include <vtkFieldData.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <vtkThreshold.h>
#include <vtkDoubleArray.h>
#include <vtkCellData.h>




namespace viennamesh {

    namespace vtk {

        calculate_mesh_quality_metrics::calculate_mesh_quality_metrics() {}

        bool calculate_mesh_quality_metrics::run(viennamesh::algorithm_handle &) {

            data_handle<viennamesh_string> meshpath = get_required_input<viennamesh_string>("mesh_path"); 

            std::string filename = meshpath();
            std::string BaseOutputPath = meshpath().substr(0,meshpath().find_last_of("/"));
        
            //store all quality values in a vector
            std::vector<std::vector<double>> qualities;
            std::vector<double> tmp_data;
            
            //create output csv-file for storing the mesh quality data
            //construct filename for csv and vtu file
            size_t found = (filename).find_last_of("/");
            size_t find_vtu = (filename).find_last_of(".");

            std::string output_filename_csv = BaseOutputPath;
            output_filename_csv += filename.substr(found+1, find_vtu-found-1);;
            output_filename_csv += "_quality.csv";

            std::string output_filename_vtu = BaseOutputPath;
            output_filename_vtu += filename.substr(found+1, find_vtu-found-1);
            output_filename_vtu += "_quality.vtu";

            //create csv file handle
            ofstream csv;
            csv.open(output_filename_csv.c_str());

            //read all the data from the file
            vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
                vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
            reader->SetFileName(filename.c_str());
            reader->Update();

            vtkSmartPointer<vtkMeshQuality> mq = 
            vtkSmartPointer<vtkMeshQuality>::New();

            mq->SetInputConnection( reader->GetOutputPort() );
            mq->SaveCellQualityOn();

            //Edge Ratio
            mq->SetTriangleQualityMeasureToEdgeRatio();
            mq->Update();

            vtkDataSet* qualityMesh = mq->GetOutput();
            vtkSmartPointer<vtkDoubleArray> qualityArray = vtkDoubleArray::SafeDownCast(qualityMesh->GetCellData()->GetArray("Quality"));
            qualityArray->SetName("Edge Ratio"); //CHANGE FOR TRIS
            reader->GetOutput()->GetCellData()->AddArray(qualityArray);

            tmp_data.resize(reader->GetNumberOfCells());    
            for (size_t i = 0; i < reader->GetNumberOfCells(); ++i)
            {
                tmp_data[i] = *qualityArray->GetTuple(i);
            }

            qualities.push_back(tmp_data);
            tmp_data.clear();

            //end Edge Ratio

            //Aspect Ratio
            mq->SetTriangleQualityMeasureToAspectRatio();
            mq->Update();

            qualityMesh = mq->GetOutput();
            qualityArray = vtkDoubleArray::SafeDownCast(qualityMesh->GetCellData()->GetArray("Quality"));
            qualityArray->SetName("Aspect Ratio");
            reader->GetOutput()->GetCellData()->AddArray(qualityArray);

            tmp_data.resize(reader->GetNumberOfCells());    
            for (size_t i = 0; i < reader->GetNumberOfCells(); ++i)
            {
                tmp_data[i] = *qualityArray->GetTuple(i);
            }

            qualities.push_back(tmp_data);
            tmp_data.clear();
            //End of Aspect Ratio

            //Radius Ratio
            mq->SetTriangleQualityMeasureToRadiusRatio();
            mq->Update();

            qualityMesh = mq->GetOutput();
            qualityArray = vtkDoubleArray::SafeDownCast(qualityMesh->GetCellData()->GetArray("Quality"));
            qualityArray->SetName("Radius Ratio");
            reader->GetOutput()->GetCellData()->AddArray(qualityArray);

            tmp_data.resize(reader->GetNumberOfCells());    
            for (size_t i = 0; i < reader->GetNumberOfCells(); ++i)
            {
                tmp_data[i] = *qualityArray->GetTuple(i);
            }

            qualities.push_back(tmp_data);
            tmp_data.clear();
            //End of Radius Ratio

            //Frobenius Norm
            mq->SetTriangleQualityMeasureToAspectFrobenius();
            mq->Update();
            

            qualityMesh = mq->GetOutput();
            qualityArray = vtkDoubleArray::SafeDownCast(qualityMesh->GetCellData()->GetArray("Quality"));
            qualityArray->SetName("Frobenius Norm");
            reader->GetOutput()->GetCellData()->AddArray(qualityArray);

            tmp_data.resize(reader->GetNumberOfCells());    
            for (size_t i = 0; i < reader->GetNumberOfCells(); ++i)
            {
                tmp_data[i] = *qualityArray->GetTuple(i);
            }

            qualities.push_back(tmp_data);
            tmp_data.clear();
            //End of Frobenius Norm

            //Minimal Dihedral Angle
            mq->SetTriangleQualityMeasureToMinAngle();
            mq->Update();

            qualityMesh = mq->GetOutput();
            qualityArray = vtkDoubleArray::SafeDownCast(qualityMesh->GetCellData()->GetArray("Quality"));
            qualityArray->SetName("Minimal Dihedral Angle");
            reader->GetOutput()->GetCellData()->AddArray(qualityArray);

            tmp_data.resize(reader->GetNumberOfCells());    
            for (size_t i = 0; i < reader->GetNumberOfCells(); ++i)
            {
                tmp_data[i] = *qualityArray->GetTuple(i);
            }

            qualities.push_back(tmp_data);
            tmp_data.clear();
            //End of Minimal Dihedral Angle



           
            //Condition
            mq->SetTriangleQualityMeasureToCondition();
            mq->Update();

            qualityMesh = mq->GetOutput();
            qualityArray = vtkDoubleArray::SafeDownCast(qualityMesh->GetCellData()->GetArray("Quality"));
            qualityArray->SetName("Condition");
            reader->GetOutput()->GetCellData()->AddArray(qualityArray);

            tmp_data.resize(reader->GetNumberOfCells());    
            for (size_t i = 0; i < reader->GetNumberOfCells(); ++i)
            {
                tmp_data[i] = *qualityArray->GetTuple(i);
            }

            qualities.push_back(tmp_data);
            tmp_data.clear();
            //End of Condition

            //Scaled Jacobian
            mq->SetTriangleQualityMeasureToScaledJacobian();
            mq->Update();

            qualityMesh = mq->GetOutput();
            qualityArray = vtkDoubleArray::SafeDownCast(qualityMesh->GetCellData()->GetArray("Quality"));
            qualityArray->SetName("Scaled Jacobian");
            reader->GetOutput()->GetCellData()->AddArray(qualityArray);

            tmp_data.resize(reader->GetNumberOfCells());    
            for (size_t i = 0; i < reader->GetNumberOfCells(); ++i)
            {
                tmp_data[i] = *qualityArray->GetTuple(i);
            }

            qualities.push_back(tmp_data);
            tmp_data.clear();
            //End of Scaled Jacobian

            //Shape
            mq->SetTriangleQualityMeasureToShape();
            mq->Update();
            

            qualityMesh = mq->GetOutput();
            qualityArray = vtkDoubleArray::SafeDownCast(qualityMesh->GetCellData()->GetArray("Quality"));
            qualityArray->SetName("Shape");
            reader->GetOutput()->GetCellData()->AddArray(qualityArray);

            tmp_data.resize(reader->GetNumberOfCells());    
            for (size_t i = 0; i < reader->GetNumberOfCells(); ++i)
            {
                tmp_data[i] = *qualityArray->GetTuple(i);
            }

            qualities.push_back(tmp_data);
            tmp_data.clear();
            //End of Shape

            //Relative Size Squared
            mq->SetTriangleQualityMeasureToRelativeSizeSquared();
            mq->Update();

            qualityMesh = mq->GetOutput();
            qualityArray = vtkDoubleArray::SafeDownCast(qualityMesh->GetCellData()->GetArray("Quality"));
            qualityArray->SetName("Relative Size Squared");
            reader->GetOutput()->GetCellData()->AddArray(qualityArray);

            tmp_data.resize(reader->GetNumberOfCells());    
            for (size_t i = 0; i < reader->GetNumberOfCells(); ++i)
            {
                tmp_data[i] = *qualityArray->GetTuple(i);
            }

            qualities.push_back(tmp_data);
            tmp_data.clear();
            //End of Relative Size Squared

            //Shape And Size
            mq->SetTriangleQualityMeasureToShapeAndSize();
            mq->Update();

            qualityMesh = mq->GetOutput();
            qualityArray = vtkDoubleArray::SafeDownCast(qualityMesh->GetCellData()->GetArray("Quality"));
            qualityArray->SetName("Shape And Size");
            reader->GetOutput()->GetCellData()->AddArray(qualityArray);

            tmp_data.resize(reader->GetNumberOfCells());    
            for (size_t i = 0; i < reader->GetNumberOfCells(); ++i)
            {
                tmp_data[i] = *qualityArray->GetTuple(i);
            }

            qualities.push_back(tmp_data);
            tmp_data.clear();
            //End of Shape And Size

            //Distortion
            mq->SetTriangleQualityMeasureToDistortion();
            mq->Update();

            qualityMesh = mq->GetOutput();
            qualityArray = vtkDoubleArray::SafeDownCast(qualityMesh->GetCellData()->GetArray("Quality"));
            qualityArray->SetName("Distortion");
            reader->GetOutput()->GetCellData()->AddArray(qualityArray);

            tmp_data.resize(reader->GetNumberOfCells());    
            for (size_t i = 0; i < reader->GetNumberOfCells(); ++i)
            {
                tmp_data[i] = *qualityArray->GetTuple(i);
            }

            qualities.push_back(tmp_data);
            tmp_data.clear();
            //End of Distortion


            //write output

            std::cout << "Computed " << qualities.size() << " different quality metrics" << std::endl;

            csv << "edge_ratio, aspect_ratio, radius_ratio, frobenius_norm, min_dihedral_angle, condition, jacobian, scaled_jacobian, shape, rel_size_squared, shape_and_size, distortion" << endl;

            std::cout << "Writing csv file " << output_filename_csv << std::endl;

            for(size_t vid = 0; vid < reader->GetNumberOfCells(); ++vid)
            {
                for (size_t mid = 0; mid < qualities.size(); ++mid)
                {
                    csv << qualities[mid][vid];

                    if (mid != (qualities.size()-1))
                    {
                        csv << ",";
                    }
                }
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
            writer->SetInputData(reader->GetOutput());
            writer->SetDataModeToAscii();
            writer->Write();






 
            return true;

        }

    }

}


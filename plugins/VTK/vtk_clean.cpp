#include "vtk_clean.hpp"
#include "vtk_mesh.hpp"

#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>


#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkCleanPolyData.h>

#include <vtkUnstructuredGrid.h>
#include <vtkCellArray.h>

#include <vtkGeometryFilter.h>
#include <vtkAppendFilter.h>


namespace viennamesh {

    namespace vtk {

        clean::clean() {}

        bool clean::run(viennamesh::algorithm_handle &) {

            data_handle<viennagrid_numeric> tolerance = get_input<viennagrid_numeric>("tolerance");

            data_handle<vtk::mesh> input_mesh = get_required_input<vtk::mesh>("mesh");


            vtkSmartPointer<vtkCleanPolyData> testclean = 
            vtkSmartPointer<vtkCleanPolyData>::New();

            testclean->SetInputData(input_mesh().GetMesh());

            testclean->PointMergingOn();

            testclean->SetAbsoluteTolerance(tolerance());
            
            testclean->Update();

            data_handle<vtk::mesh> output_mesh = make_data<vtk::mesh>();

            output_mesh().GetMesh()->DeepCopy((vtkDataObject*)testclean->GetOutput());

  
            set_output("mesh", output_mesh());






 
            return true;

        }

    }

}
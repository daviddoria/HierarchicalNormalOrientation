#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkSphereSource.h>
#include <vtkPolyDataNormals.h>
#include <vtkPointData.h>
#include <vtkCellData.h>

#include "vtkSuperPoints.h"


int main(int argc, char* argv[])
{
  if(argc != 4)
    {
    std::cerr << "Required arguments: InputFile.vtp OutputSuperPointLabels.vtp OuputSuperPointCenters.vtp" << std::endl;
    return EXIT_FAILURE;
    }
    
  std::string inputFileName = argv[1];
  std::string labelsFileName = argv[2];
  std::string centersFileName = argv[3];
  
  std::cout << "Reading " << inputFileName << std::endl;
  
  vtkSmartPointer<vtkXMLPolyDataReader> reader =
      vtkSmartPointer<vtkXMLPolyDataReader>::New();
  reader->SetFileName(inputFileName.c_str());
  reader->Update();

  vtkSmartPointer<vtkSuperPoints> superPoints =
    vtkSmartPointer<vtkSuperPoints>::New();
  superPoints->SetInputConnection(reader->GetOutputPort());
  superPoints->SetAutomaticRadiusRatio(5);
  superPoints->Update();

  std::cout << "Writing " << labelsFileName << std::endl;
  {
  vtkSmartPointer<vtkXMLPolyDataWriter> writer =
    vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName(labelsFileName.c_str());
  writer->SetInputConnection(superPoints->GetOutputPort(0));
  writer->Write();
  }

  std::cout << "Writing " << centersFileName << std::endl;
  {
  vtkSmartPointer<vtkXMLPolyDataWriter> writer =
    vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName(centersFileName.c_str());
  writer->SetInputConnection(superPoints->GetOutputPort(1));
  writer->Write();
  }
  return EXIT_SUCCESS;
}

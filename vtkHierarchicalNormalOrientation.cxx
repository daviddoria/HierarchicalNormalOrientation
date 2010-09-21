#include "vtkHierarchicalNormalOrientation.h"

#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkInformationVector.h"
#include "vtkInformation.h"
#include "vtkDataObject.h"
#include "vtkSmartPointer.h"

vtkHierarchicalNormalOrientation::vtkHierarchicalNormalOrientation()
{

}



int vtkHierarchicalNormalOrientation::RequestData(vtkInformation *vtkNotUsed(request),
                                             vtkInformationVector **inputVector,
                                             vtkInformationVector *outputVector)
{
  // Get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // Get the input and ouptut
  vtkPolyData* input = vtkPolyData::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkPolyData *output = vtkPolyData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));


  return 1;
}


//----------------------------------------------------------------------------
void vtkHierarchicalNormalOrientation::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
  //os << indent << "NNRadius: (" << this->NNRadius << ")\n";

}

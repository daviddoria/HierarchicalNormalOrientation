#ifndef __vtkHierarchicalNormalOrientation_h
#define __vtkHierarchicalNormalOrientation_h

#include "vtkPolyDataAlgorithm.h"

class vtkHierarchicalNormalOrientation : public vtkPolyDataAlgorithm
{
public:
  vtkTypeMacro(vtkHierarchicalNormalOrientation,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  static vtkHierarchicalNormalOrientation *New();

protected:
  vtkHierarchicalNormalOrientation();
  ~vtkHierarchicalNormalOrientation(){}

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
  vtkHierarchicalNormalOrientation(const vtkHierarchicalNormalOrientation&);  // Not implemented.
  void operator=(const vtkHierarchicalNormalOrientation&);  // Not implemented.

};

#endif
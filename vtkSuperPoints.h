#ifndef __vtkSuperPoints_h
#define __vtkSuperPoints_h

#include "vtkPolyDataAlgorithm.h"

class vtkSuperPoints : public vtkPolyDataAlgorithm
{
public:
  vtkTypeMacro(vtkSuperPoints,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  static vtkSuperPoints *New();

  vtkSetMacro(NNRadius, double);
  vtkSetMacro(AutomaticRadius, int);
  vtkSetMacro(AutomaticRadiusRatio, int);

protected:
  vtkSuperPoints();
  ~vtkSuperPoints(){}

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:

  class SuperPoint
  {
    public:
      std::vector<int> Points;
      int Center; // ID of the center point

      int Count(){return this->Points.size();}
      void AddPoint(int point){this->Points.push_back(point);}
      int GetPoint(int id){return this->Points[id];}

  };

  vtkSuperPoints(const vtkSuperPoints&);  // Not implemented.
  void operator=(const vtkSuperPoints&);  // Not implemented.

  double NNRadius;
  int AutomaticRadius;
  int AutomaticRadiusRatio;
  int TooSmall;
  void LabelSuperPoints(vtkPolyData* polydata);
  void ColorSuperPoints(vtkPolyData* polydata);

  void ComputeAutomaticRadius();

  void AverageIntensities(vtkPolyData* superPoints);
  void AverageNormals(vtkPolyData* superPoints);
  void AverageColors(vtkPolyData* superPoints);

  void ValidateInput(vtkPolyData*);

  std::vector<int> Used;

  int CountValidSuperPoints();
  std::vector<SuperPoint> SuperPoints;

  vtkPolyData* Input;

  void CreateSuperPointCenters(vtkPolyData*);
};

#endif
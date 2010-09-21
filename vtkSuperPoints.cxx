#include "vtkSuperPoints.h"

#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkInformationVector.h"
#include "vtkInformation.h"
#include "vtkDataObject.h"
#include "vtkSmartPointer.h"

#include <vtkMath.h>
#include <vtkLookupTable.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnsignedCharArray.h>
#include <vtkIdFilter.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkPolyData.h>
#include <vtkKdTreePointLocator.h>
#include <vtkIdList.h>
#include <vtkSelectionNode.h>
#include <vtkSelection.h>
#include <vtkExtractSelection.h>
#include <vtkIdTypeArray.h>
#include <vtkGraphWriter.h>
#include <vtkExtractSelectedGraph.h>
#include <vtkMutableUndirectedGraph.h>

#include <vtkFullyConnectedGraphFilter.h>
#include <vtkGraphBFSIterator.h>
#include <vtkUnstructuredGridToGraph.h>
#include <vtkNearestNeighborGraphFilter.h>
#include <vtkConnectGraph.h>
#include <vtkGraphVertexDataConditionalIterator.h>

#include <numeric>
#include <vector>
#include <algorithm>
#include <cmath>

// For testing only
#include <vtkGraphToPolyData.h>
#include <vtkXMLPolyDataWriter.h>

#define UNASSIGNED -1
#define TOOSMALL -2

vtkStandardNewMacro(vtkSuperPoints);

vtkSuperPoints::vtkSuperPoints()
{
  this->SetNumberOfOutputPorts(2);

  this->NNRadius = 1.0;
  this->AutomaticRadius = 1;
  this->AutomaticRadiusRatio = 10;
  this->TooSmall = 10;
}


class vtkCustomGraphIterator : public vtkGraphBFSIterator
{
public:
  static vtkCustomGraphIterator* New();
  vtkTypeMacro(vtkCustomGraphIterator, vtkGraphBFSIterator);

  vtkCustomGraphIterator(){}
  ~vtkCustomGraphIterator(){}

  protected:
    virtual void AddVertex(vtkIdType);

    vtkCustomGraphIterator(const vtkCustomGraphIterator &);  // Not implemented.
    void operator=(const vtkCustomGraphIterator &);        // Not implemented.

};
vtkStandardNewMacro(vtkCustomGraphIterator);

// Weight normals, colors, and intensities by normal quality
void vtkCustomGraphIterator::AddVertex(vtkIdType val)
{
  // Hard "use normals" threshold.
  // If normal error is above 6.6e-06, weight .5(color) + .5(intensity)
  // If normal error is low, use .33 weight for color, intensity, and normal

  // We will refer to things as "XYZSimilarity". These are scaled from 0 (not similar) to 1 (exactly the same)
  vtkFloatArray* normalError = vtkFloatArray::SafeDownCast(this->Graph->GetVertexData()->GetArray("NormalError"));
  vtkFloatArray* normals = vtkFloatArray::SafeDownCast(this->Graph->GetVertexData()->GetNormals());
  vtkUnsignedCharArray* colors = vtkUnsignedCharArray::SafeDownCast(this->Graph->GetVertexData()->GetArray("Colors"));
  vtkDoubleArray* intensities = vtkDoubleArray::SafeDownCast(this->Graph->GetVertexData()->GetArray("Intensities"));

  if(!normalError)
    {
    std::cout << "There is no NormalError array!" << std::endl;
    exit(-1);
    }
  if(!normals)
    {
    std::cout << "There is no normals array!" << std::endl;
    exit(-1);
    }
  if(!colors)
    {
    std::cout << "There is no Colors array!" << std::endl;
    exit(-1);
    }
  if(!intensities)
    {
    std::cout << "There is no Intensities array!" << std::endl;
    exit(-1);
    }

  for(vtkIdType i = 0; i < this->Graph->GetNumberOfVertices(); i++)
    {
    // Color
    unsigned char c1[3];
    colors->GetTupleValue(this->GetVertex(), c1);

    unsigned char c2[3];
    colors->GetTupleValue(val, c2);

    double colorDistance = sqrt(pow(c1[0] - c2[0],2) + pow(c1[1] - c2[1],2) + pow(c1[2] - c2[2],2)); //low if colors are similar
    double maxColorDistance = sqrt(pow(255,2) + pow(255,2) + pow(255,2));
    double colorSimilarity = 1.0 - (colorDistance/maxColorDistance);

    // Normal
    double n1[3];
    normals->GetTuple(this->GetVertex(), n1);
    vtkMath::Normalize(n1);
    //std::cout << "Point normal " << i << ": " << n[0] << " " << n[1] << " " << n[2] << std::endl;

    double n2[3];
    normals->GetTuple(val, n2);
    vtkMath::Normalize(n2);

    double angle = acos(vtkMath::Dot(n1,n2));

    double degAngle = fabs(vtkMath::DegreesFromRadians(angle));

    if(degAngle > 90)
      {
      degAngle = 180 - degAngle;
      }

    double normalSimilarity = 1.0 - (degAngle / 90.);

    // Intensity
    double intensitySimilarity = 1.0 - fabs(intensities->GetValue(this->GetVertex()) - intensities->GetValue(val));

    // Weight according to normalError
    double similarity = 0.0;
    double normalErrorThreshold = 6.6e-06;

    // If both normals are trustable, include the normals
    if((normalError->GetValue(this->GetVertex()) < normalErrorThreshold) && (normalError->GetValue(val) < normalErrorThreshold))
      {
      similarity = .33 * colorSimilarity + .33 * normalSimilarity + .33 * intensitySimilarity;
      }
    else // Ignore the normals
      {
      similarity = .5 * colorSimilarity + .5 * intensitySimilarity;
      }

    if(similarity > .9)
      {
      this->QueuePushBack(val);
      }

    } // end for
}

int vtkSuperPoints::RequestData(vtkInformation *vtkNotUsed(request),
                                             vtkInformationVector **inputVector,
                                             vtkInformationVector *outputVector)
{
  // Get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo0 = outputVector->GetInformationObject(0);
  vtkInformation *outInfo1 = outputVector->GetInformationObject(1);

  // Get the input and ouptut
  //vtkPolyData *input
  this->Input = vtkPolyData::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));

  ValidateInput(this->Input);

  vtkPolyData *outputLabeled = vtkPolyData::SafeDownCast(
    outInfo0->Get(vtkDataObject::DATA_OBJECT()));

  vtkPolyData *outputCenters = vtkPolyData::SafeDownCast(
    outInfo1->Get(vtkDataObject::DATA_OBJECT()));

  if(this->AutomaticRadius)
    {
    ComputeAutomaticRadius();
    }

  // Build a graph on all of the input points
  vtkSmartPointer<vtkIdFilter> idFilter =
    vtkSmartPointer<vtkIdFilter>::New();
  idFilter->SetInputConnection(this->Input->GetProducerPort());
  idFilter->SetIdsArrayName("OriginalIds");
  idFilter->Update();

  vtkSmartPointer<vtkNearestNeighborGraphFilter> nnFilter =
    vtkSmartPointer<vtkNearestNeighborGraphFilter>::New();
  nnFilter->SetInputConnection(idFilter->GetOutputPort());
  nnFilter->Update();

  // Connect the NN graph
  vtkSmartPointer<vtkConnectGraph> connectFilter =
    vtkSmartPointer<vtkConnectGraph>::New();
  connectFilter->SetInputConnection(nnFilter->GetOutputPort());
  connectFilter->Update();

  // Create a kd tree
  vtkSmartPointer<vtkKdTreePointLocator> kdTree =
    vtkSmartPointer<vtkKdTreePointLocator>::New();
  kdTree->SetDataSet(this->Input);
  kdTree->BuildLocator();

  this->Used.assign(this->Input->GetNumberOfPoints(), 0);

  for(vtkIdType pointID = 0; pointID < this->Input->GetNumberOfPoints(); ++pointID)
    {
    if(pointID % 1000 == 0)
      {
      std::cout << "point " << pointID << " out of " << this->Input->GetNumberOfPoints() << std::endl;
      }

    if(this->Used[pointID] != 0)
      {
      continue;
      }

    // Find all the points around the query point
    vtkSmartPointer<vtkIdList> neighbors =
      vtkSmartPointer<vtkIdList>::New();
    kdTree->FindPointsWithinRadius(this->NNRadius, pointID, neighbors, true);

    // Extract the points around the query point
    vtkSmartPointer<vtkIdTypeArray> ids =
      vtkSmartPointer<vtkIdTypeArray>::New();
    ids->SetNumberOfComponents(1);

    // Only keep points that have not already been used
    for(vtkIdType i = 0; i < neighbors->GetNumberOfIds(); i++)
      {
      if(this->Used[neighbors->GetId(i)] == 0)
        {
        ids->InsertNextValue(neighbors->GetId(i));
        }
      }

    if(ids->GetNumberOfTuples() == 0) // All points in this neighborhood already belong to a superpixel
      {
      SuperPoint superpoint;
      superpoint.Center = pointID;
      superpoint.AddPoint(pointID);
      continue;
      }

    vtkSmartPointer<vtkSelectionNode> selectionNode =
      vtkSmartPointer<vtkSelectionNode>::New();
    selectionNode->SetFieldType(vtkSelectionNode::VERTEX);
    selectionNode->SetContentType(vtkSelectionNode::INDICES);
    selectionNode->SetSelectionList(ids);

    vtkSmartPointer<vtkSelection> selection =
      vtkSmartPointer<vtkSelection>::New();
    selection->AddNode(selectionNode);

    vtkSmartPointer<vtkExtractSelectedGraph> extractSelection =
      vtkSmartPointer<vtkExtractSelectedGraph>::New();
    extractSelection->SetInput(0, connectFilter->GetOutput());
    extractSelection->SetInput(1, selection);
    extractSelection->Update();

    //std::cout << "There are " << extractSelection->GetOutput()->GetNumberOfVertices() << " extracted points." << std::endl;

    vtkDataSetAttributes* vertexData = extractSelection->GetOutput()->GetVertexData();

    vtkDataArray* normalsAfter = vtkDataArray::SafeDownCast(vertexData->GetNormals());

    vtkIdTypeArray* originalIds =
      vtkIdTypeArray::SafeDownCast(vertexData->GetArray("OriginalIds"));

    vtkSmartPointer<vtkIdList> idList =
      vtkSmartPointer<vtkIdList>::New();
    for(vtkIdType i = 0; i < originalIds->GetNumberOfTuples(); i++)
      {
      idList->InsertNextId(originalIds->GetValue(i));
      }
    vtkIdType query = idList->IsId(pointID);
    /*
    vtkIdType query = -1;
    // Map from original id to current id in selected graph
    for(vtkIdType j = 0; j < originalIds->GetNumberOfTuples(); j++)
      {
      if(originalIds->GetValue(j) == pointID)
        {
        query = j;
        //std::cout << "Original Id " << pointID << " is now ID " << j << " in extracted points" << std::endl;
        break;
        }
      }
    */

    if(query == -1)
      {
      std::cerr << "Original Id " << pointID << " was not found in extracted points - this should never happen!" << std::endl;
      exit(-1);
      }

    vtkSmartPointer<vtkCustomGraphIterator> graphIterator =
      vtkSmartPointer<vtkCustomGraphIterator>::New();
    graphIterator->SetGraph(extractSelection->GetOutput());
    graphIterator->SetStartVertex(query);

    vtkSmartPointer<vtkMutableUndirectedGraph> connectedRegion =
      vtkSmartPointer<vtkMutableUndirectedGraph>::New();
    graphIterator->GetSelectedRegion(connectedRegion);
    int numberOfConnectedPoints = connectedRegion->GetNumberOfVertices();

    vtkIdTypeArray* connectedIds =
      vtkIdTypeArray::SafeDownCast(connectedRegion->GetVertexData()->GetArray("OriginalIds"));

    SuperPoint superpoint;
    superpoint.Center = pointID;
    for(vtkIdType connectedId = 0; connectedId < connectedRegion->GetNumberOfVertices(); connectedId++)
      {
      superpoint.AddPoint(connectedIds->GetValue(connectedId));
      this->Used[connectedIds->GetValue(connectedId)] = 1;
      }

    this->SuperPoints.push_back(superpoint);

    } //end for

  outputLabeled->ShallowCopy(this->Input);
  LabelSuperPoints(outputLabeled);
  ColorSuperPoints(outputLabeled);

  CreateSuperPointCenters(outputCenters);

  return 1;
}


//----------------------------------------------------------------------------
void vtkSuperPoints::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
  os << indent << "NNRadius: (" << this->NNRadius << ")\n";
  os << indent << "AutomaticRadius: (" << this->AutomaticRadius << ")\n";
  os << indent << "AutomaticRadiusRatio: (" << this->AutomaticRadiusRatio << ")\n";
}

void vtkSuperPoints::ColorSuperPoints(vtkPolyData* polydata)
{
  // Randomly color the super points. Color superpoints that are
  // too small black.

  vtkSmartPointer<vtkLookupTable> lut =
    vtkSmartPointer<vtkLookupTable>::New();
  lut->SetNumberOfTableValues(this->SuperPoints.size());

  vtkSmartPointer<vtkUnsignedCharArray> superColors =
    vtkSmartPointer<vtkUnsignedCharArray>::New();
  superColors->SetNumberOfComponents(3);
  superColors->SetName("SuperColors");
  superColors->SetNumberOfTuples(this->Input->GetNumberOfPoints());

  vtkMath::RandomSeed(time(NULL));


  for(unsigned int i = 0; i < this->SuperPoints.size(); i++)
    {
    unsigned char randomColor[3];
    randomColor[0] = static_cast<unsigned char>(255. * vtkMath::Random(0.0,1.0));
    randomColor[1] = static_cast<unsigned char>(255. * vtkMath::Random(0.0,1.0));
    randomColor[2] = static_cast<unsigned char>(255. * vtkMath::Random(0.0,1.0));

    SuperPoint superpoint = this->SuperPoints[i];

    if(superpoint.Count() > this->TooSmall)
      {
      for(unsigned int point = 0; point < superpoint.Count(); point++)
        {
        superColors->SetTupleValue(superpoint.GetPoint(point), randomColor);
        }
      }
    else
      {
      for(unsigned int point = 0; point < superpoint.Count(); point++)
        {
        unsigned char black[3] = {0,0,0};
        superColors->SetTupleValue(superpoint.GetPoint(point), black);
        }
      }
    }

  polydata->GetPointData()->AddArray(superColors);
}

void vtkSuperPoints::LabelSuperPoints(vtkPolyData* polydata)
{
  vtkSmartPointer<vtkIntArray> labelArray =
    vtkSmartPointer<vtkIntArray>::New();
  labelArray->SetNumberOfComponents(1);
  labelArray->SetName("SuperPointLabels");
  labelArray->SetNumberOfTuples(this->Input->GetNumberOfPoints());

  for(unsigned int superpointIndex = 0; superpointIndex < this->SuperPoints.size(); superpointIndex++)
    {
    SuperPoint superpoint = this->SuperPoints[superpointIndex];
    for(unsigned int point = 0; point < superpoint.Count(); point++)
      {
      labelArray->SetValue(superpoint.GetPoint(point), superpointIndex);
      }
    }

  polydata->GetPointData()->PassData(this->Input->GetPointData());

}

void vtkSuperPoints::AverageNormals(vtkPolyData* superPointCenters)
{
  // Average normals
  vtkDataArray* inputNormals =
    vtkDataArray::SafeDownCast(this->Input->GetPointData()->GetNormals());

  //vtkSmartPointer<vtkDoubleArray> averageNormals =
    //vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkFloatArray> averageNormals =
    vtkSmartPointer<vtkFloatArray>::New();
  averageNormals->SetNumberOfComponents(3);
  averageNormals->SetName("Normals");

  for(vtkIdType superID = 0; superID < this->SuperPoints.size(); superID++)
    {
    SuperPoint superpoint = this->SuperPoints[superID];

    if(superpoint.Count() < this->TooSmall)
      {
      continue;
      }

    double x=0; double y=0; double z=0;

    for(vtkIdType point = 0; point < superpoint.Count(); point++)
      {
      double n[3];
      inputNormals->GetTuple(superpoint.GetPoint(point), n);
      x += n[0];
      y += n[1];
      z += n[2];
      }

    x /= static_cast<float>(superpoint.Count());
    y /= static_cast<float>(superpoint.Count());
    z /= static_cast<float>(superpoint.Count());

    float superNormal[3];
    superNormal[0] = x;
    superNormal[1] = y;
    superNormal[2] = z;
    averageNormals->InsertNextTupleValue(superNormal);
    }
  std::cout << "There are " << averageNormals->GetNumberOfTuples() << " average normals." << std::endl;

  superPointCenters->GetPointData()->SetNormals(averageNormals);
}

void vtkSuperPoints::AverageColors(vtkPolyData* superPointCenters)
{
  // Average colors
  vtkUnsignedCharArray* inputColors =
    vtkUnsignedCharArray::SafeDownCast(this->Input->GetPointData()->GetArray("Colors"));

  vtkSmartPointer<vtkUnsignedCharArray> averageColors =
    vtkSmartPointer<vtkUnsignedCharArray>::New();
  averageColors->SetNumberOfComponents(3);
  averageColors->SetName("Colors");

  for(vtkIdType superID = 0; superID < this->SuperPoints.size(); superID++)
    {
    SuperPoint superpoint = this->SuperPoints[superID];

    if(superpoint.Count() < this->TooSmall)
      {
      continue;
      }

    double r=0; double g=0; double b=0;

    for(unsigned int point = 0; point < superpoint.Count(); point++)
      {
      unsigned char c[3];
      inputColors->GetTupleValue(superpoint.GetPoint(point), c);
      r += c[0];
      g += c[1];
      b += c[2];
      }

    r /= static_cast<double>(superpoint.Count());
    g /= static_cast<double>(superpoint.Count());
    b /= static_cast<double>(superpoint.Count());

    unsigned char averageColor[3];
    averageColor[0] = static_cast<unsigned char>(r);
    averageColor[1] = static_cast<unsigned char>(g);
    averageColor[2] = static_cast<unsigned char>(b);
    averageColors->InsertNextTupleValue(averageColor);
    }
  superPointCenters->GetPointData()->AddArray(averageColors);

  std::cout << "There are " << averageColors->GetNumberOfTuples() << " average colors." << std::endl;
}

void vtkSuperPoints::AverageIntensities(vtkPolyData* superPointCenters)
{

  // Average intensities
  //vtkDataArray* inputIntensities =
    //vtkDataArray::SafeDownCast(input->GetPointData()->GetArray("Intensities"));
  vtkDoubleArray* inputIntensities =
    vtkDoubleArray::SafeDownCast(this->Input->GetPointData()->GetArray("Intensities"));

  vtkSmartPointer<vtkDoubleArray> averageIntensities =
    vtkSmartPointer<vtkDoubleArray>::New();
  averageIntensities->SetNumberOfComponents(1);
  averageIntensities->SetName("Intensities");

  for(vtkIdType superID = 0; superID < this->SuperPoints.size(); superID++)
    {
    SuperPoint superpoint = this->SuperPoints[superID];

    if(superpoint.Count() < this->TooSmall)
      {
      continue;
      }

    double intensities = 0.0;

    for(unsigned int point = 0; point < superpoint.Count(); point++)
      {
      double intensity = inputIntensities->GetValue(superpoint.GetPoint(point));
      intensities += intensity;
      }

    intensities /= static_cast<double>(superpoint.Count());
    averageIntensities->InsertNextValue(intensities);
    }

  superPointCenters->GetPointData()->AddArray(averageIntensities);
  std::cout << "There are " << averageIntensities->GetNumberOfTuples() << " average intensities." << std::endl;
}

void vtkSuperPoints::ComputeAutomaticRadius()
{
    // Decide on a reasonable RBNN radius
    double bounds[6];
    this->Input->GetBounds(bounds);

    double delx = bounds[1] - bounds[0];
    double dely = bounds[3] - bounds[2];
    double delz = bounds[5] - bounds[4];
    std::cout << "delx: " << delx << " dely: " << dely << " delz: " << delz << std::endl;
    double minDim = std::min(delx, std::min(dely,delz));

    this->NNRadius = minDim / static_cast<double>(this->AutomaticRadiusRatio);
    std::cout << "Automatic radius: " << this->NNRadius << std::endl;
}

void vtkSuperPoints::ValidateInput(vtkPolyData* input)
{
  /*
  vtkDataArray* normals = vtkDataArray::SafeDownCast(input->GetPointData()->GetArray("Normals"));
  vtkDataArray* colors = vtkDataArray::SafeDownCast(input->GetPointData()->GetArray("Colors"));
  vtkDataArray* intensities = vtkDataArray::SafeDownCast(input->GetPointData()->GetArray("Intensities"));
  */
  vtkFloatArray* normals = vtkFloatArray::SafeDownCast(input->GetPointData()->GetArray("Normals"));
  vtkUnsignedCharArray* colors = vtkUnsignedCharArray::SafeDownCast(input->GetPointData()->GetArray("Colors"));
  vtkDoubleArray* intensities = vtkDoubleArray::SafeDownCast(input->GetPointData()->GetArray("Intensities"));
  if(!input)
    {
    std::cerr << "Input is invalid!" << std::endl;
    exit(-1);
    }
  if(!normals)
    {
    std::cerr << "No input normals!" << std::endl;
    exit(-1);
    }
  if(!colors)
    {
    std::cerr << "No input colors!" << std::endl;
    exit(-1);
    }
  if(!intensities)
    {
    std::cerr << "No input intensities!" << std::endl;
    exit(-1);
    }
}

int vtkSuperPoints::CountValidSuperPoints()
{
  int count = 0;
  for(unsigned int i = 0; i < this->SuperPoints.size(); i++)
    {
    if(this->SuperPoints[i].Count() > this->TooSmall)
      {
      count++;
      }
    }
  return count;
}

void vtkSuperPoints::CreateSuperPointCenters(vtkPolyData* outputCenters)
{
  std::cout << "There are " << this->SuperPoints.size() << " super points." << std::endl;

  vtkSmartPointer<vtkPoints> points =
    vtkSmartPointer<vtkPoints>::New();
  for(unsigned int i = 0; i < this->SuperPoints.size(); i++)
    {
    if(this->SuperPoints[i].Count() < this->TooSmall)
      {
      continue;
      }
    double p[3];
    this->Input->GetPoint(this->SuperPoints[i].Center, p);
    points->InsertNextPoint(p);
    }

  vtkSmartPointer<vtkPolyData> polyData =
    vtkSmartPointer<vtkPolyData>::New();
  polyData->SetPoints(points);

  vtkSmartPointer<vtkVertexGlyphFilter> glyphFilter =
    vtkSmartPointer<vtkVertexGlyphFilter>::New();
  glyphFilter->SetInputConnection(polyData->GetProducerPort());
  glyphFilter->Update();

  outputCenters->ShallowCopy(glyphFilter->GetOutput());

  AverageIntensities(outputCenters);
  AverageNormals(outputCenters);
  AverageColors(outputCenters);
}
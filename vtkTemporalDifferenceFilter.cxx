/*=========================================================================

Program:   Visualization Toolkit
Module:    $RCSfile: vtkTemporalDifferenceFilter.cxx,v $

Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
All rights reserved.
See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/*

The Custom particle interpolator currently only works with the modified
time pipeline which is a work in progress in the pv-meshless branch of paraview.

This filter will not produce any output when compiled against the standard paraview

*/

#include "vtkTemporalDifferenceFilter.h"

#include "vtkCellData.h"
#include "vtkCompositeDataIterator.h"
#include "vtkDataSet.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkIdTypeArray.h"
#include "vtkIntArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPointSet.h"
#include "vtkSmartPointer.h"
#include "vtkCompositeDataPipeline.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkCellArray.h"
#include "vtkPolyData.h"
#include "vtkUnstructuredGrid.h"
#include "vtkDataArraySelection.h"

#include <vector>
#include <algorithm>
#include <functional>
#include <cmath>
//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkTemporalDifferenceFilter);
//----------------------------------------------------------------------------
vtkTemporalDifferenceFilter::vtkTemporalDifferenceFilter()
{
  this->ComputeDerivative         = 0;
  this->ComputeSpatialDifferences = 1;
  this->ComputeMagnitudes         = 0;
  this->ArrayNamePrefix           = NULL;
  this->ExplicitTimeStepCaching   = 0;
  this->LastTimeStep              = NULL;
  this->LastUpdateTime            = 0.0;
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
  this->PointDataArraySelection   = vtkSmartPointer<vtkDataArraySelection>::New();
}
//----------------------------------------------------------------------------
vtkTemporalDifferenceFilter::~vtkTemporalDifferenceFilter()
{
  this->SetArrayNamePrefix(NULL);
  this->LastTimeStep = NULL;
}
//----------------------------------------------------------------------------
void vtkTemporalDifferenceFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
  os << indent << "ComputeDerivative: " << this->ComputeDerivative << "\n";
}
//----------------------------------------------------------------------------
int vtkTemporalDifferenceFilter::FillInputPortInformation(int port, vtkInformation* info)
{
  if (port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
#ifdef JB_NEW_TIME
    info->Set(vtkCompositeDataPipeline::TEMPORAL_DATA_INPUT(), 1);
#endif
  }
  return 1;
}
//----------------------------------------------------------------------------
int vtkTemporalDifferenceFilter::FillOutputPortInformation(int vtkNotUsed(port), vtkInformation* info)
{
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPointSet");
  return 1;
}
//----------------------------------------------------------------------------
int vtkTemporalDifferenceFilter::RequestDataObject( vtkInformation*,
                                                   vtkInformationVector** inputVector ,
                                                   vtkInformationVector* outputVector)
{
  if (this->GetNumberOfInputPorts() == 0 || this->GetNumberOfOutputPorts() == 0)
  {
    return 1;
  }

  vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
  if (!inInfo)
  {
    return 0;
  }
  vtkDataObject *input = inInfo->Get(vtkDataObject::DATA_OBJECT());

  if (input)
  {
    // for each output
    for(int i=0; i < this->GetNumberOfOutputPorts(); ++i)
    {
      vtkInformation* info = outputVector->GetInformationObject(i);
      vtkDataObject *output = info->Get(vtkDataObject::DATA_OBJECT());

      if (!output || !output->IsA(input->GetClassName()))
      {
        vtkDataObject* newOutput = input->NewInstance();
        info->Set(vtkDataObject::DATA_OBJECT(), newOutput);
        newOutput->FastDelete();
      }
    }
    return 1;
  }
  return 0;
}
//----------------------------------------------------------------------------
int vtkTemporalDifferenceFilter::RequestUpdateExtent (
  vtkInformation * vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);

  // find the required input time steps and request them
  if (outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP()))
  {
    // get the update times
    double upTime =  outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());

    // get the available input times
    double *inTimes =
      inInfo->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
    int numInTimes =
      inInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());

    // only if the input is not continuous should we do anything
    if (inTimes)
    {
      //compute request times f
      double inUpTimes[2];
      int numInUpTimes(0);

      // for each requested time mark the required input times
      numInUpTimes= 0;
      // below the range
      if (upTime <= inTimes[0])
      {
        inUpTimes[numInUpTimes++] = inTimes[0];
      }
      // above the range?
      else if (upTime >= inTimes[numInTimes-1])
      {
        inUpTimes[numInUpTimes++] = inTimes[numInTimes-1];
      }
      // in the middle
      else
      {
        int i = 0;
        while (upTime > inTimes[i])
        {
          ++i;
        }
        inUpTimes[numInUpTimes++] = inTimes[i-1];
        inUpTimes[numInUpTimes++] = inTimes[i];
      }

      inInfo->Set(vtkMultiTimeStepAlgorithm::UPDATE_TIME_STEPS(), inUpTimes,numInUpTimes);
    }
  }
  return 1;
}
//----------------------------------------------------------------------------
// Change the information
int vtkTemporalDifferenceFilter::RequestInformation (
  vtkInformation * vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);

  //
  // find time on input
  //
  int     numTimes = 0;
  double *inTimes  = NULL;
  double  outRange[2];

  if (inInfo->Has(vtkStreamingDemandDrivenPipeline::TIME_STEPS()))
  {
    inTimes =
      inInfo->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
    numTimes =
      inInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());

    outRange[0] = inTimes[0];
    outRange[1] = inTimes[numTimes-1];
    outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(),
      outRange,2);
  }

  // If the time source is a traditional file reader or other, then it declares times in information
  // but if it is a live object such as a simulation, then we must explicityl cache the last time step
  // in order to compute differences
  if (numTimes>1) 
  {
    this->ExplicitTimeStepCaching = 0;
    this->LastTimeStep = NULL;
  }
  else {
    this->ExplicitTimeStepCaching = 1;
  }

  vtkDataSet *inData = vtkDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPointData *pd = inData->GetPointData();

  //
  // when information is regenerated, the arrays might change (bah!)
  //
  if (pd->GetNumberOfArrays()>0) {
    vtkSmartPointer<vtkDataArraySelection> TempArraySelection = vtkSmartPointer<vtkDataArraySelection>::New();
    for (int i=0; i<pd->GetNumberOfArrays(); i++) {
      const char *name = pd->GetArray(i)->GetName();
      TempArraySelection->AddArray(name);
      if (this->PointDataArraySelection->ArrayIsEnabled(name)) {
        TempArraySelection->EnableArray(name);
      }
      else {
        TempArraySelection->DisableArray(name);
      }
    }
    this->PointDataArraySelection->CopySelections(TempArraySelection);
  }
  return 1;
}
//----------------------------------------------------------------------------
class vtkTemporalDifferenceToleranceCheck: public std::binary_function<double, double, bool>
{
public:
  vtkTemporalDifferenceToleranceCheck(double tol) { this->tolerance = tol; }
  double tolerance;
  //
  result_type operator()(first_argument_type a, second_argument_type b) const
  {
    bool result = (fabs(a-b)<=(this->tolerance));
    return (result_type)result;
  }
};
//----------------------------------------------------------------------------
int vtkTemporalDifferenceFilter::ComputeInputUpdateExtent (
  vtkInformation * vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // if the input is not a time source, we can't do much other than cache time steps as they arrive
  if (this->ExplicitTimeStepCaching) {
    return 1;
  }

  // otherwise ... get the info objects
  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);

  // find the required input time steps and request them
  if (outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP()))
  {
    double requestedTimeValue = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
    double timestep = std::find_if(
      this->TimeStepValues.begin(), this->TimeStepValues.end(),
      std::bind2nd( vtkTemporalDifferenceToleranceCheck( this->TimeStepTolerance ), requestedTimeValue ))
      - this->TimeStepValues.begin();

    // We need this time step and the previous one ...
    if (timestep>0) {
      inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP(), this->TimeStepValues[timestep-1]);
    }
    else if (this->TimeStepValues.size()>0) {
      inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP(), this->TimeStepValues[timestep]);
    }
  }
  return 1;
}
//----------------------------------------------------------------------------
// templated difference function
template <class T>
void vtkTemporalDifferencedVdt(vtkTemporalDifferenceFilter *tdf,
                               double dt,
                               vtkDataArray *output,
                               vtkDataArray **arrays,
                               vtkIdType numComp,                                    
                               T *)
{
  T *outData = static_cast<T*>(output->GetVoidPointer(0));
  T *inData0 = static_cast<T*>(arrays[0]->GetVoidPointer(0));
  T *inData1 = static_cast<T*>(arrays[1]->GetVoidPointer(0));
  //
  vtkIdType N = arrays[0]->GetNumberOfTuples();
  for (vtkIdType t=0; t<N; ++t)
  {
    T *x0 = &inData0[t*numComp];
    T *x1 = &inData1[t*numComp];
    if (!tdf->GetComputeMagnitudes()) {
        for (int c=0; c<numComp; ++c)
          *outData++ = static_cast<T>((x1[c]-x0[c])/dt);
      }
    else {
      for (int c=0; c<numComp; ++c)
        *outData++ = static_cast<T>(std::abs((x1[c]-x0[c])/dt));
    }
  }
  output->SetNumberOfTuples(N);
}
//----------------------------------------------------------------------------
// templated zero function
template <class T>
void vtkTemporalDifferenceZero(vtkDataArray *output,
                               vtkIdType nC, vtkIdType nT,                                    
                               T *)
{
  T *outData = static_cast<T*>(output->GetVoidPointer(0));
  //
  for (vtkIdType t=0; t<nT*nC; ++t)
  {
    *outData++ = static_cast<T>(0);
  }
}
//----------------------------------------------------------------------------
vtkDataArray *vtkTDFNewArray(vtkDataArray *da, vtkIdType Nc, vtkIdType Nt, const char *prefix)
{
  //
  // Create the array
  //
  vtkAbstractArray *aa = da->CreateArray(da->GetDataType());
  vtkDataArray *output = vtkDataArray::SafeDownCast(aa);
  //
  // initialize 
  //
  output->SetNumberOfComponents(Nc);
  output->SetNumberOfTuples(Nt);
  std::string newname = std::string(prefix) + da->GetName();
  output->SetName(newname.c_str());
  return output;
}  
//----------------------------------------------------------------------------
vtkDataArray *ZeroDataArray(vtkDataArray *da, vtkIdType Nc, vtkIdType Nt, const char *prefix)
{
  //
  // Create the output array
  //
  vtkDataArray *output = vtkTDFNewArray(da, Nc, Nt, prefix);

  // now do the interpolation
  switch (da->GetDataType())
  {
    vtkTemplateMacro(vtkTemporalDifferenceZero
      (output, Nc, Nt, static_cast<VTK_TT *>(0)));
  default:
    return 0;
  }

  return output;
}
//----------------------------------------------------------------------------
vtkDataArray *vtkTemporalDifferenceFilter::DifferenceDataArray(double dt, vtkDataArray **arrays, vtkIdType Nt)
{
  //
  // Create the output array
  //
  int Nc = arrays[0]->GetNumberOfComponents();
  vtkDataArray *output = vtkTDFNewArray(arrays[0], Nc, Nt, this->ArrayNamePrefix);

  // now do the interpolation
  switch (arrays[0]->GetDataType())
  {
    vtkTemplateMacro(vtkTemporalDifferencedVdt
      (this, dt, output, arrays, Nc, static_cast<VTK_TT *>(0)));
  default:
    vtkErrorMacro(<< "Execute: Unknown ScalarType");
    return 0;
  }

  return output;
}
//----------------------------------------------------------------------------
vtkDataSet *vtkTemporalDifferenceFilter::DifferenceDataSet(vtkDataSet *in1, vtkDataSet *in2, double dt)
{
  vtkDataSet *input[2] = { in1, in2 };

  vtkDataSet *output = input[0]->NewInstance();
  output->CopyStructure(input[0]);
  output->GetPointData()->ShallowCopy(input[0]->GetPointData());
  output->GetCellData()->ShallowCopy(input[0]->GetCellData());

  //
  // Compute spatial difference if the dataset is a vtkPointSet
  //
  std::vector<vtkDataArray*> arrays;
  vtkDataArray *outarray;
  //
  if (this->ComputeSpatialDifferences) {
    vtkPointSet *inPointSet[2] = { vtkPointSet::SafeDownCast(input[0]), vtkPointSet::SafeDownCast(input[1]) };
    for (int i=0; i<2; ++i) {
      if (inPointSet[i]) {
        arrays.push_back(inPointSet[i]->GetPoints()->GetData());
      }
    }
    if (arrays.size()>1) {
      outarray = this->DifferenceDataArray(dt, &arrays[0], arrays[0]->GetNumberOfTuples());
    }
    else {
      outarray = ZeroDataArray(arrays[0], 
        arrays[0]->GetNumberOfComponents(), arrays[0]->GetNumberOfTuples(), this->ArrayNamePrefix);
    }
    output->GetPointData()->AddArray(outarray);
    outarray->FastDelete();
  }
  //
  // Loop over all pointdata 
  //
  for (int s=0; s<input[0]->GetPointData()->GetNumberOfArrays(); ++s) {
    arrays.clear();
    //
    // On some data, the scalar arrays are consistent but ordered
    // differently on each time step, so we will fetch them by name if
    // possible.
    //
    vtkDataArray *dataarray = input[0]->GetPointData()->GetArray(s);
    char *scalarname = dataarray->GetName();
    if (this->GetPointArrayStatus(scalarname)) {
      arrays.push_back(dataarray);
      if (input[1] && this->GetPointArrayStatus(scalarname)) {
        vtkDataArray *dataarray = input[1]->GetPointData()->GetArray(scalarname);
        arrays.push_back(dataarray);
      }
      if (arrays.size()>1) {
        outarray = this->DifferenceDataArray(dt, &arrays[0], arrays[0]->GetNumberOfTuples());
      }
      else {
        outarray = ZeroDataArray(arrays[0], 
          arrays[0]->GetNumberOfComponents(), arrays[0]->GetNumberOfTuples(), this->ArrayNamePrefix);
      }
      output->GetPointData()->AddArray(outarray);
      outarray->FastDelete();
    }
  }

  //
  // Interpolate celldata if present
  //
  for (int s=0; s<input[0]->GetCellData()->GetNumberOfArrays(); ++s) 
  {
    arrays.clear();
    char *scalarname = NULL;
    for (int i=0; i<2; ++i) 
    {
      //
      // On some data, the scalar arrays are consistent but ordered
      // differently on each time step, so we will fetch them by name if
      // possible.
      //
      if (i==0 || (scalarname==NULL)) 
      {
        vtkDataArray *dataarray = input[i]->GetCellData()->GetArray(s);
        scalarname = dataarray->GetName();
        arrays.push_back(dataarray);
      }
      else if (input[i])
      {
        vtkDataArray *dataarray = 
          input[i]->GetCellData()->GetArray(scalarname);
        arrays.push_back(dataarray);
      }
    }
    if (arrays.size()>1) 
    {
      outarray = this->DifferenceDataArray(dt, &arrays[0], arrays[0]->GetNumberOfTuples());
    }
    else 
    {
      outarray = ZeroDataArray(arrays[0], 
        arrays[0]->GetNumberOfComponents(), arrays[0]->GetNumberOfTuples(), this->ArrayNamePrefix);
    }
    output->GetPointData()->AddArray(outarray);
    outarray->FastDelete();
  }

  if (in1->GetInformation()->Has(vtkDataObject::DATA_GEOMETRY_UNMODIFIED()) &&
    in2->GetInformation()->Has(vtkDataObject::DATA_GEOMETRY_UNMODIFIED()))
  {
    output->GetInformation()->Set(vtkDataObject::DATA_GEOMETRY_UNMODIFIED(),1);
  }
  return output;
}
//----------------------------------------------------------------------------
int vtkTemporalDifferenceFilter::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  vtkDataSet *outData = NULL;

  // get the requested update times
  double upTime = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());

  vtkMultiBlockDataSet *inData = vtkMultiBlockDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  int numTimeSteps  = inData->GetNumberOfBlocks();

  if (numTimeSteps==1)
  {
    vtkDataSet* data0 = vtkDataSet::SafeDownCast(inData->GetBlock(0));
    // Do a diff, with a NULL second dataset to create the correct delta_ arrays
    outData = this->DifferenceDataSet(data0, NULL, 0.0);

    if (data0->GetInformation()->Has(vtkDataObject::DATA_GEOMETRY_UNMODIFIED()))
    {
      outData->GetInformation()->Set(vtkDataObject::DATA_GEOMETRY_UNMODIFIED(),1);
    }
    outInfo->Set(vtkDataObject::DATA_OBJECT(),outData);
  }
  else
  {
    vtkDataSet* data0 = vtkDataSet::SafeDownCast(inData->GetBlock(0));
    vtkDataSet* data1 = vtkDataSet::SafeDownCast(inData->GetBlock(1));
    if (data0==NULL && data1==NULL)
    {
      vtkErrorMacro("Null data set");
      return 0;
    }
    // interpolate i-1 and i
    double t0 = data0->GetInformation()->Get(vtkDataObject::DATA_TIME_STEP());
    double t1 = data1->GetInformation()->Get(vtkDataObject::DATA_TIME_STEP());

    // compute delta T
    if (t0!=t1) {
      this->DeltaT = t1-t0;  
    }
    else {
      this->DeltaT = 1.0;  
    }

    outData = this->DifferenceDataSet(data0, data1, this->DeltaT);
    outInfo->Set(vtkDataObject::DATA_OBJECT(),outData);
    outData->FastDelete();
  }
  // stamp this new dataset with a time key
  outData->GetInformation()->Set(vtkDataObject::DATA_TIME_STEP(),upTime);
  //
  this->LastUpdateTime = upTime;
  return 1;
}
//----------------------------------------------------------------------------
int vtkTemporalDifferenceFilter::GetNumberOfPointArrays()
{
  return this->PointDataArraySelection->GetNumberOfArrays();
}
//----------------------------------------------------------------------------
const char* vtkTemporalDifferenceFilter::GetPointArrayName(int index)
{
  return this->PointDataArraySelection->GetArrayName(index);
}
//----------------------------------------------------------------------------
int vtkTemporalDifferenceFilter::GetPointArrayStatus(const char* name)
{
  return this->PointDataArraySelection->ArrayIsEnabled(name);
}
//----------------------------------------------------------------------------
void vtkTemporalDifferenceFilter::SetPointArrayStatus(const char* name, int status)
{
  if (status!=this->GetPointArrayStatus(name))
  {
    //    this->MeshParamsModifiedTime.Modified();
    if (status)
    {
      this->PointDataArraySelection->EnableArray(name);
    }
    else
    {
      this->PointDataArraySelection->DisableArray(name);
    }
    this->Modified();
  }
}

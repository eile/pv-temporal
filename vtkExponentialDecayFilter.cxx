/*=========================================================================

Program:   Visualization Toolkit
Module:    $RCSfile: vtkExponentialDecayFilter.cxx,v $

Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
All rights reserved.
See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkExponentialDecayFilter.h"

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
#include "vtkMath.h"
//
#include <vector>
#include <algorithm>
#include <functional>
#include <cmath>
//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkExponentialDecayFilter);
//----------------------------------------------------------------------------
double decaycompute(double avg, double decay, double lasttime, double thistime, double value) 
{
  double alpha = 1.0 - 1.0/exp(decay*abs(thistime-lasttime));
  avg += alpha*(value - avg);
  return avg;
}
//----------------------------------------------------------------------------
vtkExponentialDecayFilter::vtkExponentialDecayFilter()
{
  this->DecayFactor               = 100.0;
  this->HighFrequencyResponse     = 0;
  this->HighFrequencyDelta        = 20.0;
  this->ArrayNamePrefix           = NULL;
  this->LastPointData             = vtkSmartPointer<vtkPointData>::New();
  this->LastUpdateTime            = 0.0;
  this->FirstIteration            = true;
  this->OutputAbsoluteValue       = 1;
  this->ClampAndNormalizeOutput   = 1;
  this->NormalizedRange[0] = -2;
  this->NormalizedRange[1] =  2;

  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
  this->PointDataArraySelection   = vtkSmartPointer<vtkDataArraySelection>::New();
}
//----------------------------------------------------------------------------
vtkExponentialDecayFilter::~vtkExponentialDecayFilter()
{
  this->SetArrayNamePrefix(NULL);
  this->LastPointData = NULL;
}
//----------------------------------------------------------------------------
// Change the information
int vtkExponentialDecayFilter::ExecuteInformation (
  vtkInformation * vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);

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

  double upTime = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
  if (upTime<=this->LastUpdateTime && !this->FirstIteration) {
    this->FirstIteration = true;
    this->LastPointData = vtkSmartPointer<vtkPointData>::New();
  }
  return 1;
}
//----------------------------------------------------------------------------
// templated difference function
template <class T>
void vtkExponentialDecayCompute(vtkExponentialDecayFilter *tdf,
  bool first, double decay, double lasttime, double thistime,
                               vtkDataArray *output,
                               vtkDataArray **arrays,
                               vtkIdType numComp,                                    
                               T *)
{
  T *outData = static_cast<T*>(output->GetVoidPointer(0));
  T *inData0 = static_cast<T*>(arrays[0]->GetVoidPointer(0));
  T *inData1 = NULL;
  if (arrays[1]) {
    inData1 = static_cast<T*>(arrays[1]->GetVoidPointer(0));
  }
  bool    HFR = tdf->GetHighFrequencyResponse();
  double  HFD = tdf->GetHighFrequencyDelta();
  bool    ABS = tdf->GetOutputAbsoluteValue();
  bool    NOR = tdf->GetClampAndNormalizeOutput();
  double *NRM = tdf->GetNormalizedRange();
  //
  vtkIdType N = arrays[0]->GetNumberOfTuples();
  for (vtkIdType t=0; t<N; ++t)
  {
    T *value = &inData0[t*numComp];
    double vv = *value;
    double pv = 0.0;
    if (arrays[1]) {
      inData1 = static_cast<T*>(arrays[1]->GetVoidPointer(0));
      pv = inData1[t*numComp];
    }
    for (int c=0; c<numComp; ++c) {
      double temp;
      if (first || (HFR && abs(vv)>HFD)) {
        temp = ABS ? (vv) : static_cast<T>(vv);
      }
      else {
        temp = decaycompute(pv, decay, lasttime, thistime, vv);
        temp = ABS ? abs(temp) : temp;
      }
      if (NRM) {
        *outData++ = vtkMath::ClampAndNormalizeValue(temp, NRM);
      }
      else {
        *outData++ = static_cast<T>(temp);
      }
    }
  }
  output->SetNumberOfTuples(N);
}
//----------------------------------------------------------------------------
vtkDataArray *vtkEDFNewArray(vtkDataArray *da, vtkIdType Nc, vtkIdType Nt, const char *prefix)
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
vtkDataArray *vtkExponentialDecayFilter::DecayDataArray(double timevalue, vtkDataArray **arrays, vtkIdType Nt)
{
  //
  // Create the output array
  //
  int Nc = arrays[0]->GetNumberOfComponents();
  vtkDataArray *output = vtkEDFNewArray(arrays[0], Nc, Nt, this->ArrayNamePrefix);

  // now do the interpolation
  switch (arrays[0]->GetDataType())
  {
    vtkTemplateMacro(vtkExponentialDecayCompute
      (this, 
       this->FirstIteration, this->DecayFactor, this->LastUpdateTime, timevalue,
       output, arrays, Nc, static_cast<VTK_TT *>(0)));
  default:
    vtkErrorMacro(<< "Execute: Unknown ScalarType");
    return 0;
  }

  return output;
}
//----------------------------------------------------------------------------
vtkDataSet *vtkExponentialDecayFilter::DecayDataSet(vtkDataSet *in1, double timevalue)
{
  vtkDataSet *output = in1->NewInstance();
  output->CopyStructure(in1);
  output->GetPointData()->ShallowCopy(in1->GetPointData());
  output->GetCellData()->ShallowCopy(in1->GetCellData());

  //
  // Compute spatial difference if the dataset is a vtkPointSet
  //
  std::vector<vtkDataArray*> arrays;
  vtkDataArray *outarray;
  //
  //
  // Loop over all pointdata 
  //
  for (int s=0; s<in1->GetPointData()->GetNumberOfArrays(); ++s) {
    arrays.clear();
    //
    // On some data, the scalar arrays are consistent but ordered
    // differently on each time step, so we will fetch them by name if
    // possible.
    //
    vtkDataArray *dataarray = in1->GetPointData()->GetArray(s);
    char *scalarname = dataarray->GetName();
    if (this->GetPointArrayStatus(scalarname)) {
      arrays.push_back(dataarray);
      vtkDataArray *dataarray = LastPointData->GetArray((std::string(this->ArrayNamePrefix)+std::string(scalarname)).c_str());
      arrays.push_back(dataarray); // NULL is OK
      outarray = this->DecayDataArray(timevalue, &arrays[0], arrays[0]->GetNumberOfTuples());
      this->FirstIteration = false;
      output->GetPointData()->AddArray(outarray);
      outarray->FastDelete();
    }
  }

/*
  if (in1->GetInformation()->Has(vtkDataObject::DATA_GEOMETRY_UNMODIFIED()))
  {
    output->GetInformation()->Set(vtkDataObject::DATA_GEOMETRY_UNMODIFIED(),1);
  }
  */
  return output;
}
//----------------------------------------------------------------------------
int vtkExponentialDecayFilter::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  vtkDataSet *outData = NULL;

  // get the requested update times
  double upTime = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());

  vtkPointSet *data0 = vtkPointSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

  // Do a diff, with a NULL second dataset to create the correct delta_ arrays
  outData = this->DecayDataSet(data0, upTime);

/*
  if (data0->GetInformation()->Has(vtkDataObject::DATA_GEOMETRY_UNMODIFIED()))
  {
    outData->GetInformation()->Set(vtkDataObject::DATA_GEOMETRY_UNMODIFIED(),1);
  }
  */
  outInfo->Set(vtkDataObject::DATA_OBJECT(),outData);

  // stamp this new dataset with a time key
  outData->GetInformation()->Set(vtkDataObject::DATA_TIME_STEP(),upTime);
  //
  this->LastUpdateTime = upTime;
  this->LastPointData->ShallowCopy(outData->GetPointData());
  return 1;
}
//----------------------------------------------------------------------------
int vtkExponentialDecayFilter::GetNumberOfPointArrays()
{
  return this->PointDataArraySelection->GetNumberOfArrays();
}
//----------------------------------------------------------------------------
const char* vtkExponentialDecayFilter::GetPointArrayName(int index)
{
  return this->PointDataArraySelection->GetArrayName(index);
}
//----------------------------------------------------------------------------
int vtkExponentialDecayFilter::GetPointArrayStatus(const char* name)
{
  return this->PointDataArraySelection->ArrayIsEnabled(name);
}
//----------------------------------------------------------------------------
void vtkExponentialDecayFilter::SetPointArrayStatus(const char* name, int status)
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

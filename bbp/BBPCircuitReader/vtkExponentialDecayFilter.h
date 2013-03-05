/*=========================================================================

Program:   Visualization Toolkit
Module:    $RCSfile: vtkExponentialDecayFilter.h,v $

Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
All rights reserved.
See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkExponentialDecayFilter - Create a pseudo moving average of data using an Exponential Decay Function
// .SECTION Description
// 

#ifndef __vtkExponentialDecayFilter_h
#define __vtkExponentialDecayFilter_h

#include "vtkPointSetAlgorithm.h"
#include "vtkSmartPointer.h" // required
#include <vector>     // required

class vtkDataSet;
class vtkDataArraySelection;
class vtkPointData;

class VTK_EXPORT vtkExponentialDecayFilter : public vtkPointSetAlgorithm
{
public:
  static vtkExponentialDecayFilter *New();
  vtkTypeMacro(vtkExponentialDecayFilter, vtkPointSetAlgorithm);

  // Description:
  vtkSetMacro(DecayFactor, double);
  vtkGetMacro(DecayFactor, int);

  // Description:
  vtkSetStringMacro(ArrayNamePrefix);
  vtkGetStringMacro(ArrayNamePrefix); 

  // Description:
  // Get the number of point or cell arrays available in the input.
  int GetNumberOfPointArrays();

  // Description:
  // Get the name of the point or cell array with the given index in
  // the input.
  const char* GetPointArrayName(int index);

  // Description:
  // Get/Set whether the point or cell array with the given name is to
  // be read.
  int GetPointArrayStatus(const char* name);
  void SetPointArrayStatus(const char* name, int status);

protected:
   vtkExponentialDecayFilter();
  ~vtkExponentialDecayFilter();

  virtual int ExecuteInformation(vtkInformation *,
    vtkInformationVector **,
    vtkInformationVector *);

  virtual int RequestData(vtkInformation *,
    vtkInformationVector **,
    vtkInformationVector *);

  // Description:
  virtual vtkDataSet *DecayDataSet(vtkDataSet *in1, double timevalue);

  // Description:
  virtual vtkDataArray *DecayDataArray(double timevalue, vtkDataArray **arrays, vtkIdType N);

  // internal variables
  bool   FirstIteration;
  double DecayFactor;
  char  *ArrayNamePrefix;

  double TimeStepTolerance;
  double LastUpdateTime;
  //BTX
  std::vector<double>           TimeStepValues;
  vtkSmartPointer<vtkPointData> LastPointData;
  //ETX

  // To allow paraview gui to enable/disable scalars
  vtkSmartPointer<vtkDataArraySelection> PointDataArraySelection;

private:
  vtkExponentialDecayFilter(const vtkExponentialDecayFilter&);  // Not implemented.
  void operator=(const vtkExponentialDecayFilter&);  // Not implemented.
};

#endif

/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkTemporalDifferenceFilter.h,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkTemporalDifferenceFilter - display displacement vectors for a mesh
// .SECTION Description
// vtkTemporalDifferenceFilter takes two successive meshes in a time series and adds
// a vector array which represents the displacement of each pt to the corresponding 
// one at the next time step. 
// The displacement attached to data at T=1 is what would be applied for T=0 to reach T=1

#ifndef __vtkTemporalDifferenceFilter_h
#define __vtkTemporalDifferenceFilter_h

#include "vtkMultiTimeStepAlgorithm.h"
#include "vtkSmartPointer.h" // required
#include <vector>     // required

class vtkDataSet;

class VTK_EXPORT vtkTemporalDifferenceFilter : public vtkMultiTimeStepAlgorithm
{
public:
  static vtkTemporalDifferenceFilter *New();
  vtkTypeMacro(vtkTemporalDifferenceFilter, vtkMultiTimeStepAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  vtkSetMacro(ComputeMagnitudes, int);
  vtkGetMacro(ComputeMagnitudes, int);
  vtkBooleanMacro(ComputeMagnitudes, int);

  // Description:
  vtkSetMacro(ComputeScalarDifferences, int);
  vtkGetMacro(ComputeScalarDifferences, int);
  vtkBooleanMacro(ComputeScalarDifferences, int);

  // Description:
  vtkSetMacro(ComputeSpatialDifferences, int);
  vtkGetMacro(ComputeSpatialDifferences, int);
  vtkBooleanMacro(ComputeSpatialDifferences, int);

  // Description:
  vtkSetMacro(ComputeDerivative, int);
  vtkGetMacro(ComputeDerivative, int);
  vtkBooleanMacro(ComputeDerivative, int);

  // Description:
  vtkSetStringMacro(ArrayNamePrefix);
  vtkGetStringMacro(ArrayNamePrefix); 

protected:
   vtkTemporalDifferenceFilter();
  ~vtkTemporalDifferenceFilter();

  virtual int FillInputPortInformation(int port, vtkInformation* info);
  virtual int FillOutputPortInformation(int vtkNotUsed(port), vtkInformation* info);

  virtual int RequestInformation(vtkInformation *,
                                 vtkInformationVector **,
                                 vtkInformationVector *);

  virtual int RequestData(vtkInformation *,
                          vtkInformationVector **,
                          vtkInformationVector *);

  virtual int ComputeInputUpdateExtent(vtkInformation *,
                          vtkInformationVector **,
                          vtkInformationVector *);

  virtual int RequestUpdateExtent(vtkInformation *,
                          vtkInformationVector **,
                          vtkInformationVector *);

  virtual int RequestDataObject(vtkInformation *,
                          vtkInformationVector **,
                          vtkInformationVector *);

  // Description:
  virtual vtkDataSet *DifferenceDataSet(vtkDataSet *in1, 
                                        vtkDataSet *in2,
                                        double dt);

  // Description:
  virtual vtkDataArray *DifferenceDataArray(double ratio, 
                                            vtkDataArray **arrays, 
                                            vtkIdType N);

  // internal variables
  int    ComputeMagnitudes;
  int    ComputeScalarDifferences;
  int    ComputeSpatialDifferences;
  int    ComputeDerivative;
  char  *ArrayNamePrefix;

  double TimeStepTolerance;
  double DeltaT;
  int    ExplicitTimeStepCaching;
  double LastUpdateTime;
  //BTX
  std::vector<double>      TimeStepValues;
  vtkSmartPointer<vtkDataSet> LastTimeStep;
  //ETX

private:
  vtkTemporalDifferenceFilter(const vtkTemporalDifferenceFilter&);  // Not implemented.
  void operator=(const vtkTemporalDifferenceFilter&);  // Not implemented.
};

#endif

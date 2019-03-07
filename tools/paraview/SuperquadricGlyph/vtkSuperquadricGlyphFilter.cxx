/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile$

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// This plugin has been developed and contributed by Daniel Queteschiner, JKU
// Linz, Austria.

#include "vtkSuperquadricGlyphFilter.h"

#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkSuperquadricSource.h"
#include "vtkCellArray.h"
#include "vtkAppendPolyData.h"
#include "vtkMath.h"
#include "vtkDoubleArray.h"
#include "vtkIdTypeArray.h"
#include "vtkQuaternion.h"
#include "vtkTransform.h"


vtkStandardNewMacro(vtkSuperquadricGlyphFilter);

vtkSuperquadricGlyphFilter::vtkSuperquadricGlyphFilter()
{
  this->SetNumberOfInputPorts(1); // Set the number of input ports used by the algorithm.
  this->ColorMode = COLOR_BY_SCALARS;
  this->ThetaResolution = 16;
  this->PhiResolution = 16;
  this->ThetaRoundness = 0.3;
  this->PhiRoundness = 0.3;
  this->HalfAxisVectorArray = NULL;
  this->HalfAxisXArray = NULL;
  this->HalfAxisYArray = NULL;
  this->HalfAxisZArray = NULL;
  this->Blockiness1Array = NULL;
  this->Blockiness2Array = NULL;
  this->Quat1Array = NULL;
  this->Quat2Array = NULL;
  this->Quat3Array = NULL;
  this->Quat4Array = NULL;
  this->RotationMatrixArray = NULL;

  // by default process active point scalars
  this->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,
                               vtkDataSetAttributes::SCALARS);
  // by default process active point vectors
  this->SetInputArrayToProcess(1,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,
                               vtkDataSetAttributes::VECTORS);
  // by default process active point normals
  this->SetInputArrayToProcess(2,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,
                               vtkDataSetAttributes::NORMALS);
  // by default process active point scalars (color scalars)
  this->SetInputArrayToProcess(3,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,
                               vtkDataSetAttributes::SCALARS);
}

vtkSuperquadricGlyphFilter::~vtkSuperquadricGlyphFilter()
{
  delete [] HalfAxisVectorArray;
  delete [] HalfAxisXArray;
  delete [] HalfAxisYArray;
  delete [] HalfAxisZArray;
  delete [] Blockiness1Array;
  delete [] Blockiness2Array;
  delete [] Quat1Array;
  delete [] Quat2Array;
  delete [] Quat3Array;
  delete [] Quat4Array;
  delete [] RotationMatrixArray;
}

int vtkSuperquadricGlyphFilter::FillInputPortInformation(
  int vtkNotUsed(port), vtkInformation* info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
  return 1;
}


int vtkSuperquadricGlyphFilter::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkDataObject  *doInput = inInfo->Get(vtkDataObject::DATA_OBJECT());
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the input and output
  vtkDataSet *input = vtkDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  if (!input)
  {
    if (doInput)
    {
      vtkErrorMacro("This filter cannot process input of type: " << doInput->GetClassName());
    }
    return 0;
  }

  vtkPointData *inputPD;
  vtkPointData *outputPD;

  vtkSuperquadricSource *superquadric;
  vtkAppendPolyData* append = NULL;
  vtkIdType numPtsIn, numGlyphPtsTotal, inPtId;
  int axis, abort = 0;
  vtkSuperquadricSource** superquadricsArray = NULL;
  vtkIdType superquadricsArrayCount = 0;
  vtkDataArray *inScalars = NULL;
  vtkDataArray *inVectors = NULL;
  vtkDataArray *inNormals = NULL;
  vtkDataArray *inCScalars = NULL;
  vtkIdTypeArray *numGlyphPts;
  double theta, phi;
  double scalex, scaley, scalez;
  double quatw, quatx, quaty, quatz;
  vtkQuaterniond quat;
  vtkTransform *trans = vtkTransform::New();

  inputPD = input->GetPointData();
  if (!inputPD)
  {
    vtkErrorMacro(<<"No data to glyph!");
    return 1;
  }

  numPtsIn = input->GetNumberOfPoints();
  if (numPtsIn < 1)
  {
    vtkErrorMacro(<<"No points to glyph!");
    return 1;
  }

  inScalars = this->GetInputArrayToProcess(0, inputVector);
  inVectors = this->GetInputArrayToProcess(1, inputVector);
  inNormals = this->GetInputArrayToProcess(2, input);
  inCScalars = this->GetInputArrayToProcess(3, input);
  if (inCScalars == NULL)
  {
    inCScalars = inScalars;
  }

  vtkDataArray *halfaxisdata = this->HalfAxisVectorArray ? input->GetPointData()->GetArray(this->HalfAxisVectorArray) : NULL;
  vtkDataArray *halfaxisxdata = this->HalfAxisXArray ? input->GetPointData()->GetArray(this->HalfAxisXArray) : NULL;
  vtkDataArray *halfaxisydata = this->HalfAxisYArray ? input->GetPointData()->GetArray(this->HalfAxisYArray) : NULL;
  vtkDataArray *halfaxiszdata = this->HalfAxisZArray ? input->GetPointData()->GetArray(this->HalfAxisZArray) : NULL;
  vtkDataArray *blockiness1data = this->Blockiness1Array ? input->GetPointData()->GetArray(this->Blockiness1Array) : NULL;
  vtkDataArray *blockiness2data = this->Blockiness2Array ? input->GetPointData()->GetArray(this->Blockiness2Array) : NULL;
  vtkDataArray *quat1data = this->Quat1Array ? input->GetPointData()->GetArray(this->Quat1Array) : NULL;
  vtkDataArray *quat2data = this->Quat2Array ? input->GetPointData()->GetArray(this->Quat2Array) : NULL;
  vtkDataArray *quat3data = this->Quat3Array ? input->GetPointData()->GetArray(this->Quat3Array) : NULL;
  vtkDataArray *quat4data = this->Quat4Array ? input->GetPointData()->GetArray(this->Quat4Array) : NULL;
  vtkDataArray *rotmatrixdata = this->RotationMatrixArray ? input->GetPointData()->GetArray(this->RotationMatrixArray) : NULL;

  append = vtkAppendPolyData::New();

  // Not every glyph necessarily has the same number of points.
  numGlyphPts = vtkIdTypeArray::New();
  numGlyphPts->SetNumberOfValues(numPtsIn);
  superquadricsArray = new vtkSuperquadricSource*[numPtsIn];

  for (inPtId = 0; inPtId < numPtsIn; ++inPtId)
  {
    if (! (inPtId % 500))
    {
      this->UpdateProgress(static_cast<double>(inPtId)/numPtsIn);
      abort = this->GetAbortExecute();
    }

    if (!abort)
    {
      superquadricsArray[inPtId] = NULL;

      {
        scalex = scaley = scalez = 1.;
        superquadric = vtkSuperquadricSource::New();

        // scale superquadric according to point data
        if (this->HalfAxisVector)
        {
          if (halfaxisdata)
          {
            double* halfscale = halfaxisdata->GetTuple3(inPtId);
            scalex = 2.*halfscale[0];
            scaley = 2.*halfscale[1];
            scalez = 2.*halfscale[2];
          }
        }
        else
        {
          if (halfaxisxdata) scalex = 2.*halfaxisxdata->GetTuple1(inPtId);
          if (halfaxisydata) scaley = 2.*halfaxisydata->GetTuple1(inPtId);
          if (halfaxiszdata) scalez = 2.*halfaxiszdata->GetTuple1(inPtId);
        }

        superquadric->SetScale(scalex, scaley, scalez);
        superquadric->SetThetaResolution(this->ThetaResolution);
        superquadric->SetPhiResolution(this->PhiResolution);
        superquadric->ToroidalOff();
        axis = 2;  // superquadrics with rotational symmetriy around the z axis

        // Set roundness per tensor, if requested
        theta = this->ThetaRoundness;
        phi = this->PhiRoundness;

        if (!this->FixedThetaPhiRoundness)
        {
          // theta,phi == 2/blockiness1, 2/blockiness2
          if (blockiness1data)
            theta = 2./blockiness1data->GetTuple1(inPtId);
          if (blockiness2data)
            phi   = 2./blockiness2data->GetTuple1(inPtId);
        }

        superquadric->SetAxisOfSymmetry(axis);
        superquadric->SetThetaRoundness(theta);
        superquadric->SetPhiRoundness(phi);
        superquadric->Update();

        numGlyphPts->SetValue(inPtId, superquadric->GetOutput()->GetNumberOfPoints());

        superquadricsArray[superquadricsArrayCount++] = superquadric;
      }
    }
  } // end for loop

  vtkIdType count;

  append->SetUserManagedInputs(1);
  append->SetNumberOfInputs(superquadricsArrayCount);

  for (count = 0; count < superquadricsArrayCount; count++)
  {
    superquadric = superquadricsArray[count];
    if (superquadric != NULL)
    {
      append->SetInputConnectionByNumber(count, superquadric->GetOutputPort());
    }
  }
  append->Update();

  for (count = 0; count < superquadricsArrayCount; count++)
  {
    superquadric = superquadricsArray[count];
    if(superquadric != NULL)
    {
      superquadric->Delete();
      superquadricsArray[count] = NULL;
      superquadric = NULL;
    }
  }
  delete [] superquadricsArray;
  superquadricsArray = NULL;

  if (append)
    output->ShallowCopy(append->GetOutput());

  numGlyphPtsTotal = 0;
  for (inPtId = 0; inPtId < numPtsIn; inPtId++)
  {
    numGlyphPtsTotal += numGlyphPts->GetValue(inPtId);
  }

  outputPD = output->GetPointData();
  // Allocate storage for output PolyData
  outputPD->CopyAllocate(inputPD, numGlyphPtsTotal);
  outputPD->CopyVectorsOff();
  outputPD->CopyNormalsOff();
  outputPD->CopyTCoordsOff();

  {
    vtkIdType ptIncr=0;
    vtkPoints * outputPoints = output->GetPoints();
    for (inPtId=0; inPtId < numPtsIn; inPtId++)
    {
      double* pos = input->GetPoint(inPtId);
      trans->PostMultiply();

      if (this->QuatOrientation)
      {
        trans->Identity();
        quatw = 1.;
        quatx = quaty = quatz = 0.;
        if (quat1data) quatw = quat1data->GetTuple1(inPtId);
        if (quat2data) quatx = quat2data->GetTuple1(inPtId);
        if (quat3data) quaty = quat3data->GetTuple1(inPtId);
        if (quat4data) quatz = quat4data->GetTuple1(inPtId);

        quat.Set(quatw,quatx,quaty,quatz);
        double axis[3] = {};
        double angle;
        angle = vtkMath::DegreesFromRadians(quat.GetRotationAngleAndAxis(axis));
        trans->RotateWXYZ(angle, axis[0], axis[1], axis[2]);
      }
      else
      {
        double *m3x3;
        double m4x4[16] = {};
        m4x4[0] = m4x4[5] = m4x4[10]= m4x4[15] = 1.;
        if (rotmatrixdata)
        {
          m3x3 = rotmatrixdata->GetTuple9(inPtId);
          m4x4[0]  = m3x3[0];
          m4x4[1]  = m3x3[1];
          m4x4[2]  = m3x3[2];
          m4x4[3]  = 0.;
          m4x4[4]  = m3x3[3];
          m4x4[5]  = m3x3[4];
          m4x4[6]  = m3x3[5];
          m4x4[7]  = 0.;
          m4x4[8]  = m3x3[6];
          m4x4[9]  = m3x3[7];
          m4x4[10] = m3x3[8];
          m4x4[11] = 0.;
          m4x4[12] = 0.;
          m4x4[13] = 0.;
          m4x4[14] = 0.;
          m4x4[15] = 1.;
        }
        trans->SetMatrix(m4x4);
      }

      trans->Translate(pos[0], pos[1], pos[2]);

      vtkIdType numSourcePts = numGlyphPts->GetValue(inPtId);
      for (vtkIdType i=0; i < numSourcePts; i++)
      {
        double pointPos[3];
        outputPoints->GetPoint(ptIncr+i, pointPos);
        trans->TransformPoint(pointPos,pointPos);
        outputPoints->SetPoint(ptIncr+i, pointPos);

        if(inputPD) outputPD->CopyData(inputPD, inPtId, ptIncr+i);
        else printf("no input PointData to copy data from\n");
      }
      ptIncr += numSourcePts;
    }
  }

  trans->Delete();
  numGlyphPts->Delete();
  append->Delete();

  return 1;
}


void vtkSuperquadricGlyphFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}


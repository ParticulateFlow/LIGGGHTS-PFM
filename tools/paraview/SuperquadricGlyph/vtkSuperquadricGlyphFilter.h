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

// .NAME vtkSuperquadricGlyphFilter - scale and orient superquadric glyph

// .SECTION Description
// vtkSuperquadricGlyphFilter generates a superquadric glyph at every point in
// the input data set. The glyphs are oriented and scaled according to the input
// data set. By default a different glyph is used per point. By setting
// FixedThetaPhiRoundness to true fixed theta and phi roundness values for all
// points can be chosen. Set theta and phi roundness to 0.0 to get rectangular
// glyphs, set them to 1.0 to get ellipsoidal glyphs, set theta roundness to 1.0
// and phi roundness to 0.0 to get cylindrical glyphs.
// vtkSuperquadricGlyphFilter operates on any type of data set. Its output is
// polygonal.

// .SECTION Thanks
// This plugin has been developed and contributed by Daniel Queteschiner, JKU
// Linz, Austria.

// .SECTION See Also
// vtkSuperquadricTensorGlyph

#ifndef vtkSuperquadricGlyphFilter_h
#define vtkSuperquadricGlyphFilter_h

#include "vtkPolyDataAlgorithm.h"

class vtkSuperquadricSource;

class VTK_EXPORT vtkSuperquadricGlyphFilter : public vtkPolyDataAlgorithm
{
public:
  vtkTypeMacro(vtkSuperquadricGlyphFilter, vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
  static vtkSuperquadricGlyphFilter *New();

  // Description:
  // Set the number of points in the longitude direction. Initial value is 16.
  vtkSetMacro(ThetaResolution,int);

  // Description:
  // Set the number of points in the latitude direction. Initial value is 16.
  vtkSetMacro(PhiResolution,int);

  // Description:
  // Set Superquadric east/west roundness.
  // Values range from 0 (rectangular) to 1 (circular) to higher orders.
  // Initial value is 1.0.
  vtkSetMacro(ThetaRoundness,double);

  // Description:
  // Set Superquadric north/south roundness.
  // Values range from 0 (rectangular) to 1 (circular) to higher orders.
  // Initial value is 1.0.
  vtkSetMacro(PhiRoundness,double);

  // Description:
  // If true, then set theta and phi roundness settings for superquadric per tensor. False by default.
  vtkGetMacro(FixedThetaPhiRoundness, int);
  vtkSetMacro(FixedThetaPhiRoundness, int);
  vtkBooleanMacro(FixedThetaPhiRoundness, int);

  enum
  {
      COLOR_BY_SCALARS
  };

  // Description:
  // Set the color mode to be used for the glyphs. The recognized values are:
  // COLOR_BY_SCALARS = 0 (default)
  vtkSetClampMacro(ColorMode, int, COLOR_BY_SCALARS, COLOR_BY_SCALARS);
  vtkGetMacro(ColorMode, int);
  void SetColorModeToScalars()
    {this->SetColorMode(COLOR_BY_SCALARS);}

  // Description:
  // If true, then set size from vector array. True by default.
  vtkGetMacro(HalfAxisVector, int);
  vtkSetMacro(HalfAxisVector, int);
  vtkBooleanMacro(HalfAxisVector, int);

  // Description:
  // Array to use to control all 3 half axis
  vtkSetStringMacro(HalfAxisVectorArray);
  vtkGetStringMacro(HalfAxisVectorArray);
  // Description:
  // Array to use to control half axis x
  vtkSetStringMacro(HalfAxisXArray);
  vtkGetStringMacro(HalfAxisXArray);
  // Description:
  // Array to use to control half axis y
  vtkSetStringMacro(HalfAxisYArray);
  vtkGetStringMacro(HalfAxisYArray);
  // Description:
  // Array to use to control half axis z
  vtkSetStringMacro(HalfAxisZArray);
  vtkGetStringMacro(HalfAxisZArray);

  // Description:
  // Array to use to control blockiness 1
  vtkSetStringMacro(Blockiness1Array);
  vtkGetStringMacro(Blockiness1Array);
  // Description:
  // Array to use to control blockiness 2
  vtkSetStringMacro(Blockiness2Array);
  vtkGetStringMacro(Blockiness2Array);

  // Description:
  // If true, then set orientation via quaternion. True by default.
  vtkGetMacro(QuatOrientation, int);
  vtkSetMacro(QuatOrientation, int);
  vtkBooleanMacro(QuatOrientation, int);
  // Description:
  // Array to use to control quat 1
  vtkSetStringMacro(Quat1Array);
  vtkGetStringMacro(Quat1Array);
  // Description:
  // Array to use to control quat 2
  vtkSetStringMacro(Quat2Array);
  vtkGetStringMacro(Quat2Array);
  // Description:
  // Array to use to control quat 3
  vtkSetStringMacro(Quat3Array);
  vtkGetStringMacro(Quat3Array);
  // Description:
  // Array to use to control quat 4
  vtkSetStringMacro(Quat4Array);
  vtkGetStringMacro(Quat4Array);
  // Description:
  // Array to use to control glyph rotation matrix
  vtkSetStringMacro(RotationMatrixArray);
  vtkGetStringMacro(RotationMatrixArray);

protected:
  vtkSuperquadricGlyphFilter();
  ~vtkSuperquadricGlyphFilter();

  /* implementation of algorithm */
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
  virtual int FillInputPortInformation(int port, vtkInformation* info);

  int ThetaResolution;
  int PhiResolution;
  double ThetaRoundness;
  double PhiRoundness;
  int FixedThetaPhiRoundness;
  int ColorMode; // The coloring mode to use for the glyphs.

  int HalfAxisVector; // use vector instead of 3 separate values
  char* HalfAxisVectorArray;
  char* HalfAxisXArray;
  char* HalfAxisYArray;
  char* HalfAxisZArray;
  char* Blockiness1Array;
  char* Blockiness2Array;
  int QuatOrientation;
  char* Quat1Array;
  char* Quat2Array;
  char* Quat3Array;
  char* Quat4Array;
  char* RotationMatrixArray;

private:
  vtkSuperquadricGlyphFilter(const vtkSuperquadricGlyphFilter&);  // Not implemented.
  void operator=(const vtkSuperquadricGlyphFilter&);  // Not implemented.
};

#endif

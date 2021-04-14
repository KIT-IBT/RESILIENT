/*! This file contains all header files. that would be used.*/

//
//  ioHeader.h
//  RESILIENT
//
//  Created by Andreas Wachter on 01.09.15.
//  Copyright (c) 2015 IBT. All rights reserved.
//

#ifndef _RESILIENT_Headers_h_
#define _RESILIENT_Headers_h_

#include <vector>
#include <list>
#include <sstream>
#include <string>
#include <iostream>
#include <map>
#include <set>
#include <cmath>
#include <stdexcept>
#include <unordered_map>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <unistd.h>

#include "DataFormat.h"
#include "Material.h"
#include "tetgen.h"
#include "ConvertTetgenio.h"
#include "Config.h"
#include "AveragingOrientation.h"
#include "Methods.h"
#include "Reader.h"
#include "Writer.h"

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkSelectEnclosedPoints.h>
#include <vtkCell.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkPointLocator.h>
#include <vtkPointSet.h>
#include <vtkPointSource.h>
#include <vtkMath.h>
#include <vtkLine.h>
#include <vtkPlane.h>
#include <vtkDataArray.h>
#include <vtkSortDataArray.h>
#include <vtkGeometryFilter.h>
#include <vtkCleanPolyData.h>
#include <vtkSelectionNode.h>
#include <vtkSelection.h>
#include <vtkExtractSelection.h>
#include <vtkMergePoints.h>
#include <vtkParametricEllipsoid.h>
#include <vtkParametricFunctionSource.h>
#include <vtkSphereSource.h>
#include <vtkAppendPolyData.h>
#include <vtkLineSource.h>
#include <vtkTubeFilter.h>
#include <vtkTriangleFilter.h>
#include <vtkTriangle.h>
#include <vtkParametricSpline.h>
#include <vtkSTLWriter.h>
#include <vtkSTLReader.h>
#include <vtkSmoothPolyDataFilter.h>

#endif  // ifndef _RESILIENT_Headers_h_

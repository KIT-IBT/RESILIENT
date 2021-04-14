/*!For converting an vtk polydata in a tetgen format and a tetgen mesh in an vtk unstrucerd grid.*/

//
//  convertTetgenio.h
//  RESILIENT
//
//  Created by Andreas Wachter on 01.09.15.
//  Copyright (c) 2015 IBT. All rights reserved.
//

#ifndef _RESILIENT_ConvertTetgenio_
#define _RESILIENT_ConvertTetgenio_

#include <stdio.h>
#include "Headers.h"

class ConvertTetgenio {
 public:
  static tetgenio convertVTPtoTetgenio(vtkSmartPointer<vtkPolyData> &polyData);
  static vtkSmartPointer<vtkUnstructuredGrid> convertTetgeniotoVTU(tetgenio &tetData);
};


#endif /* defined(_RESILIENT_convertTetgenio_) */

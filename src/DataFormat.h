/*! It stores and contains all mesh informations (the actually mesh, the input tpye, the centerpoint of the cells, double layer, the surface of the mesh and the centrepoints
 of the surface cells).*/

 //
 // DataFormat.h
 // RESILIENT
 //
 // Created by Andreas Wachter on 01.09.15.
 // Copyright (c) 2015 IBT. All rights reserved.
 //

#ifndef _RESILIENT_DataFormat_
#define _RESILIENT_DataFormat_

#include <iostream>
#include <vtkUnstructuredGrid.h>
#include <vtkPolyData.h>
#include <vector>
#include <vtkPointLocator.h>
#include <vtkPointSet.h>
#include <vtkSmartPointer.h>
#include <vtkIdList.h>


// class for new Datatype with methods to add a polydata object
class DataFormat {
public:
	// constructor
	DataFormat();
	enum PossibleInputType { undef = 0, vtu = 2, vtp = 3, vtk = 4 };

	// methods
	void setInputType(PossibleInputType);
	void setCentrePoints(vtkSmartPointer<vtkUnstructuredGrid>);
	void setDoubleLayer(bool doubleLayer);
	void setVtkData(vtkSmartPointer<vtkPointSet>);
	void setRemoveCells(std::vector<std::vector<double>> removeCells);
	void setSurfacewithNormals(vtkSmartPointer<vtkPolyData> surfacewithNormals);
	void setCentrePointsSurface(vtkSmartPointer<vtkUnstructuredGrid>);
	void setSurfaceForBridges(vtkSmartPointer<vtkPolyData> surface);
	PossibleInputType getInputType();
	vtkSmartPointer<vtkUnstructuredGrid> getCentrePoints();
	vtkSmartPointer<vtkPointSet> getVtkData();
	bool getDoubleLayer();
	std::vector<std::vector<double>> getRemoveCells();
	vtkSmartPointer<vtkPolyData> getSurfacewithNormals();
	vtkSmartPointer<vtkPolyData> getSurfaceForBridges();
	vtkSmartPointer<vtkUnstructuredGrid> getCentrePointsSurface();
	static DataFormat deepCopy(DataFormat &oldData);

private:
	PossibleInputType inputType;
	vtkSmartPointer<vtkPointSet> vtkData;
	vtkSmartPointer<vtkUnstructuredGrid> centrePoints;
	bool DoubleLayer;
	std::vector<std::vector<double>> RemoveCells;
	vtkSmartPointer<vtkPolyData> SurfacewithNormals;
	vtkSmartPointer<vtkPolyData> SurfaceForBridges;
	vtkSmartPointer<vtkUnstructuredGrid> CentrePointsSurface;
}; // class DataFormat

#endif /* defined(_RESILIENT_DataFormat_) */

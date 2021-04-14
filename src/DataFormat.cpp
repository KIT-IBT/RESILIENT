/*! \file DataFormat

   \brief It stores and contains all mesh informations (the actually mesh, the input tpye, the centerpoint of the cells, double layer, the surface of the mesh and the
	  centrepoints of the surface cells).

   \version 1.0.0

   \date Andreas Wachter 01.09.15

   \author Andreas Wachter\n
   Institute of Biomedical Engineering\n
   Kit\n
   http://www.ibt.kit.de\n

   \sa Synopsis \ref DataFormat
 */


#include "DataFormat.h"


 /*! Constructor for the fiber orientation and alation. It contains all informations (the actually mesh, the input tpye,
	 the centerpoint of the cells, double layer, the  surface of the mesh and the centrepoints of the surface cells). */
DataFormat::DataFormat() {
	inputType = PossibleInputType::undef;
	centrePoints = vtkSmartPointer<vtkUnstructuredGrid>::New();
	DoubleLayer = true;
	SurfacewithNormals = vtkSmartPointer<vtkPolyData>::New();
	CentrePointsSurface = vtkSmartPointer<vtkUnstructuredGrid>::New();
	SurfaceForBridges = vtkSmartPointer<vtkPolyData>::New();
}

// methods for class DataFormat

/*! Set function for the input type.
   \param type input type of the mesh code as PossibleInputType
 */
void DataFormat::setInputType(PossibleInputType type) {
	inputType = type;
}

/*! Get function for the input type.
   \return The type input type of the mesh code as PossibleInputType
 */
DataFormat::PossibleInputType DataFormat::getInputType() {
	return inputType;
}

/*! Set function for the centerpoints of the mesh cells.
   \param cPoints Centerpoints of the mesh cells as vtkUnstructuredGrid. The point-ids and the cell-id must be the same.
 */
void DataFormat::setCentrePoints(vtkSmartPointer<vtkUnstructuredGrid> cPoints) {
	centrePoints = cPoints;
}

/*! Get function for the centerpoints of the mesh cells.
   \return Centerpoints of the mesh cells as vtkUnstructuredGrid.
 */
vtkSmartPointer<vtkUnstructuredGrid> DataFormat::getCentrePoints() {
	return centrePoints;
}

/*! Set function for the orginal mesh.
   \param data The orginal mesh as a vtk pointSet, but it stores also polydata and unstructured grid.
 */
void DataFormat::setVtkData(vtkSmartPointer<vtkPointSet> data) {
	vtkData = data;
}

/*! Get function for the orginal mesh.
   \return The orginal mesh as a vtk pointSet.
 */
vtkSmartPointer<vtkPointSet> DataFormat::getVtkData() {
	return vtkData;
}

/*! Set function for the double layer.
   \param doubleLayer True for polydata and else False.
 */
void DataFormat::setDoubleLayer(bool doubleLayer) {
	DoubleLayer = doubleLayer;
}

/*! Get function for the the double layer.
   \return True if it is a volume mesh, else False.
 */
bool DataFormat::getDoubleLayer() {
	return DoubleLayer;
}

/*! Set function for the remove cells between the atrials.
   \param removeCells A vector that contains a vector with the coordination (x,y,z) of the removed cells between the
	  atrials.
 */
void DataFormat::setRemoveCells(std::vector<std::vector<double>> removeCells) {
	RemoveCells = removeCells;
}

/*! Get function for the remove cells between the atrials.
   \return A vector that contains a vector with the coordination (x,y,z) of the removed cells between the atrials.
 */
std::vector<std::vector<double>> DataFormat::getRemoveCells() {
	return RemoveCells;
}

/*! Set function for the surface of the mesh.
   \param surfaceWithNormals The Surface of the mesh as vtk polydata.
 */
void DataFormat::setSurfacewithNormals(vtkSmartPointer<vtkPolyData> surfaceWithNormals) {
	SurfacewithNormals = surfaceWithNormals;
}

/*! Get function for the surface of the mesh.
   \param The Surface of the mesh as vtk polydata.
 */
vtkSmartPointer<vtkPolyData> DataFormat::getSurfacewithNormals() {
	return SurfacewithNormals;
}

/*! Set function for the surface of the mesh.
   \param surfaceWithNormals The Surface of the mesh as vtk polydata.
 */
void DataFormat::setSurfaceForBridges(vtkSmartPointer<vtkPolyData> surface) {
	SurfaceForBridges = surface;
}

/*! Get function for the surface of the mesh.
   \param The Surface of the mesh as vtk polydata.
 */
vtkSmartPointer<vtkPolyData> DataFormat::getSurfaceForBridges() {
	return SurfaceForBridges;
}

/*! Set function for the centerpoints of the mesh surface cells.
   \param cPoints Centerpoints of the mesh surface cells as vtkUnstructuredGrid. The point-ids and the cell-id must be
	  the same.
 */
void DataFormat::setCentrePointsSurface(vtkSmartPointer<vtkUnstructuredGrid> cPointsSurface) {
	CentrePointsSurface = cPointsSurface;
}

/*! Set function for the centerpoints of the mesh surface cells.
   \return The centerpoints of the mesh surface cells as vtkUnstructuredGrid.*/
vtkSmartPointer<vtkUnstructuredGrid> DataFormat::getCentrePointsSurface() {
	return CentrePointsSurface;
}

/*! Make a deep copy form a Datafomat.
   \param oldData Pointer to the orginal Data.
   \return A new idetically Dataformat.
 */
DataFormat DataFormat::deepCopy(DataFormat &oldData) {
	DataFormat newData;

	vtkSmartPointer<vtkUnstructuredGrid> centerPoints = vtkSmartPointer<vtkUnstructuredGrid>::New();

	centerPoints->DeepCopy(oldData.getCentrePoints());
	newData.setCentrePoints(centerPoints);

	bool doubleLayer = oldData.getDoubleLayer();
	newData.setDoubleLayer(doubleLayer);

	DataFormat::PossibleInputType inputType = oldData.getInputType();
	newData.setInputType(inputType);

	vtkSmartPointer<vtkPointSet> vtkData;
	if (oldData.getInputType() == DataFormat::vtp) {
		vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
		polyData->DeepCopy(vtkPolyData::SafeDownCast(oldData.getVtkData()));
		vtkData = polyData;
	}
	else {
		vtkSmartPointer<vtkUnstructuredGrid> uGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
		uGrid->DeepCopy(vtkUnstructuredGrid::SafeDownCast(oldData.getVtkData()));
		vtkData = uGrid;
	}
	newData.setVtkData(vtkData);

	vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
	polyData->DeepCopy(oldData.getSurfacewithNormals());
	newData.setSurfacewithNormals(polyData);

	vtkSmartPointer<vtkPolyData> polyData2 = vtkSmartPointer<vtkPolyData>::New();
	polyData2->DeepCopy(oldData.getSurfaceForBridges());
	newData.setSurfaceForBridges(polyData2);

	vtkSmartPointer<vtkUnstructuredGrid> centerPointsSurface = vtkSmartPointer<vtkUnstructuredGrid>::New();
	centerPointsSurface->DeepCopy(oldData.getCentrePointsSurface());
	newData.setCentrePointsSurface(centerPointsSurface);

	return newData;
}  // DataFormat::deepCopy

/*!
   \page DataFormat

   \section DESCRIPTION_DataFormat DESCRIPTION
   It stores and contains all mesh informations (the actually mesh, the input tpye, the centerpoint of the cells, double layer, the surface of the mesh and the
	  centrepoints of the surface cells).

   \section SOURCE_DataFormat SOURCE

   DataFormat.cpp

   \section CHANGELOG_DataFormat CHANGELOG
   V1.0.0 - 01.09.2015 (Andreas Wachter): Starting with changelog\n
 */

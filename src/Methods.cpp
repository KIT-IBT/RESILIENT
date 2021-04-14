/*! \file Methods.cpp

 \brief This is the methods container that contains the most methods for the RESILIENT.

 \version 1.0.0

 \date Andreas Wachter 01.09.15

 \author Andreas Wachter\n
 Institute of Biomedical Engineering\n
 Kit\n
 http://www.ibt.kit.de\n

 \sa Synopsis \ref Methods
 */

#include "Methods.h"
#include "Writer.h"
#include "tetgen.h"
#include "cork.h"

#include <vector>
#include <limits>
#include <cstddef>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>

using namespace std;

// methods for Algorithm class

/*! Function for searching a path it the atrial mesh. It used a modificated dijkstra algorithm.
 \param data Pointer to the orginal mesh.
 \param startPoint Start point of the path as Vector (x,y,z).
 \param targetPoint Target point of the path as Vector (x,y,z).
 \param atrium In which atrium the path would be searched (right, left, leftEpi, leftEndo or in both atrial).
 \return VtkIdList that conatins the cell-ids of the path.
 */
vtkSmartPointer<vtkIdList> Methods::pathSearch(DataFormat &data, vector<double> startPoint, vector<double> targetPoint, string atrium) {
	vtkSmartPointer<vtkPoints> centrePoints = data.getCentrePoints()->GetPoints();
	vtkSmartPointer<vtkUnstructuredGrid> centrePointsUGrid = data.getCentrePoints();


	list<vtkIdType> tempPath;
	set<vtkIdType> arrayIDs;

	vtkSmartPointer<vtkDoubleArray> possibleCells = vtkSmartPointer<vtkDoubleArray>::New();

	possibleCells->SetNumberOfComponents(2);
	possibleCells->SetComponentName(0, "CellID");
	possibleCells->SetComponentName(1, "Distance to the starting point");

	// initialize
	vtkSmartPointer<vtkDoubleArray> distanceMatrix = vtkSmartPointer<vtkDoubleArray>::New();
	distanceMatrix->SetName("Distance");
	distanceMatrix->SetNumberOfComponents(2);
	distanceMatrix->SetNumberOfTuples(centrePoints->GetNumberOfPoints());
	distanceMatrix->SetComponentName(0, "CellID parent");
	distanceMatrix->SetComponentName(1, "Distance to the starting point");
	distanceMatrix->FillComponent(0, -1);
	distanceMatrix->FillComponent(1, std::numeric_limits<double>::max());


	vtkIdType startPointId = -1;
	vtkIdType targetPointId = -1;
	double startPointArray[3] = { 0, 0, 0 };
	double targetPointArray[3] = { 0, 0, 0 };
	if (atrium.compare("right") == 0) {
		startPointId = findClosedPointIdinMaterialInRight(data, startPoint);
		targetPointId = findClosedPointIdinMaterialInRight(data, targetPoint);
		startPointArray[0] = centrePoints->GetPoint(startPointId)[0];
		startPointArray[1] = centrePoints->GetPoint(startPointId)[1];
		startPointArray[2] = centrePoints->GetPoint(startPointId)[2];
		targetPointArray[0] = centrePoints->GetPoint(targetPointId)[0];
		targetPointArray[1] = centrePoints->GetPoint(targetPointId)[1];
		targetPointArray[2] = centrePoints->GetPoint(targetPointId)[2];
	}
	else if (atrium.compare("left") == 0) {
		startPointId = findClosedPointIdinMaterialInLeft(data, startPoint);
		targetPointId = findClosedPointIdinMaterialInLeft(data, targetPoint);
		startPointArray[0] = centrePoints->GetPoint(startPointId)[0];
		startPointArray[1] = centrePoints->GetPoint(startPointId)[1];
		startPointArray[2] = centrePoints->GetPoint(startPointId)[2];
		targetPointArray[0] = centrePoints->GetPoint(targetPointId)[0];
		targetPointArray[1] = centrePoints->GetPoint(targetPointId)[1];
		targetPointArray[2] = centrePoints->GetPoint(targetPointId)[2];
	}
	else if (atrium.compare("leftEpi") == 0) {
		startPointId = findClosedPointIdinMaterialInLeftEpi(data, startPoint);
		targetPointId = findClosedPointIdinMaterialInLeftEpi(data, targetPoint);
		startPointArray[0] = centrePoints->GetPoint(startPointId)[0];
		startPointArray[1] = centrePoints->GetPoint(startPointId)[1];
		startPointArray[2] = centrePoints->GetPoint(startPointId)[2];
		targetPointArray[0] = centrePoints->GetPoint(targetPointId)[0];
		targetPointArray[1] = centrePoints->GetPoint(targetPointId)[1];
		targetPointArray[2] = centrePoints->GetPoint(targetPointId)[2];
	}
	else if (atrium.compare("leftEndo") == 0) {
		startPointId = findClosedPointIdinMaterialInLeftEndo(data, startPoint);
		targetPointId = findClosedPointIdinMaterialInLeftEndo(data, targetPoint);
		startPointArray[0] = centrePoints->GetPoint(startPointId)[0];
		startPointArray[1] = centrePoints->GetPoint(startPointId)[1];
		startPointArray[2] = centrePoints->GetPoint(startPointId)[2];
		targetPointArray[0] = centrePoints->GetPoint(targetPointId)[0];
		targetPointArray[1] = centrePoints->GetPoint(targetPointId)[1];
		targetPointArray[2] = centrePoints->GetPoint(targetPointId)[2];
	}
	else {
		startPointArray[0] = startPoint.at(0);
		startPointArray[1] = startPoint.at(1);
		startPointArray[2] = startPoint.at(2);
		targetPointArray[0] = targetPoint.at(0);
		targetPointArray[1] = targetPoint.at(1);
		targetPointArray[2] = targetPoint.at(2);

		// cout << "Point was not test of in left or right atrium"<< endl;
		startPointId = centrePointsUGrid->FindPoint(startPointArray);
		targetPointId = centrePointsUGrid->FindPoint(targetPointArray);
	}

	if ((startPointId < 0) || (targetPointId < 0)) {
		cout << "Point was not found" << endl;
		writeIntermediateData(data, "error");
		exit(1);
	}

	if (startPointId == targetPointId) {
		vtkSmartPointer<vtkIdList> path = vtkSmartPointer<vtkIdList>::New();
		path->InsertNextId(startPointId);

		return path;
	}
    
	// starting point insert

	vtkIdType currentCellId = startPointId;

	double distanceMatrixArray[2];
	distanceMatrixArray[0] = 0;
	distanceMatrixArray[1] = 0;
	distanceMatrix->InsertTuple(currentCellId, distanceMatrixArray);

	double possibleCellsArray[2];
	possibleCellsArray[0] = currentCellId;
	possibleCellsArray[1] = 0;
	possibleCells->InsertTuple(0, possibleCellsArray);
	arrayIDs.insert(currentCellId);


	vtkSmartPointer<vtkDoubleArray> distanceMap = vtkSmartPointer<vtkDoubleArray>::New();
	distanceMap->SetName("DistanceMap");
	distanceMap->SetNumberOfComponents(1);
	distanceMap->SetNumberOfTuples(centrePoints->GetNumberOfPoints());
	distanceMap->FillComponent(0, 0);

	int count = 0;
	bool exit = false;

	// forward search (distance map)
	while (!exit) {
		vtkIdType parentID = possibleCells->vtkDataArray::GetComponent(0, 0);
		distanceMap->SetComponent(parentID, 0, count);
		count++;

		possibleCells->RemoveFirstTuple();

		arrayIDs.erase(parentID);

		double parentDistance = distanceMatrix->GetComponent(parentID, 1);

		vtkSmartPointer<vtkIdList> cellPointIds = vtkSmartPointer<vtkIdList>::New();
		vtkSmartPointer<vtkIdList> currentNeighborCellIds = vtkSmartPointer<vtkIdList>::New();

		data.getVtkData()->GetCellPoints(parentID, cellPointIds);

		for (vtkIdType i = 0; i < cellPointIds->GetNumberOfIds(); i++) {
			vtkSmartPointer<vtkIdList> PointIdList = vtkSmartPointer<vtkIdList>::New();
			PointIdList->InsertNextId(cellPointIds->GetId(i));

			vtkSmartPointer<vtkIdList> tempNeighborCellIds = vtkSmartPointer<vtkIdList>::New();

			data.getVtkData()->GetCellNeighbors(parentID, PointIdList, tempNeighborCellIds);

			for (vtkIdType j = 0; j < tempNeighborCellIds->GetNumberOfIds(); j++) {
				int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(tempNeighborCellIds->GetId(j), 0);
				if (atrium.compare("right") == 0) {
					if (Material::isInRight(material)) {
						currentNeighborCellIds->InsertUniqueId(tempNeighborCellIds->GetId(j));
					}
				}
				else if (atrium.compare("left") == 0) {
					if (Material::isInLeft(material)) {
						currentNeighborCellIds->InsertUniqueId(tempNeighborCellIds->GetId(j));
					}
				}
				else if (atrium.compare("leftEpi") == 0) {
					if (Material::isInLeftEpi(material)) {
						currentNeighborCellIds->InsertUniqueId(tempNeighborCellIds->GetId(j));
					}
				}
				else if (atrium.compare("leftEndo") == 0) {
					if (Material::isInLeftEndo(material)) {
						currentNeighborCellIds->InsertUniqueId(tempNeighborCellIds->GetId(j));
					}
				}
				else if (Material::isInLeft(material) || Material::isInRight(material)) {
					currentNeighborCellIds->InsertUniqueId(tempNeighborCellIds->GetId(j));
				}
			}
		}


		for (vtkIdType i = 0; i < currentNeighborCellIds->GetNumberOfIds(); i++) {
			vtkIdType cellId = currentNeighborCellIds->GetId(i);
			double parentPoint[3] = { centrePoints->GetPoint(parentID)[0], centrePoints->GetPoint(parentID)[1], centrePoints->GetPoint(parentID)[2] };
			double cellPoint[3] = { centrePoints->GetPoint(cellId)[0], centrePoints->GetPoint(cellId)[1], centrePoints->GetPoint(cellId)[2] };

			double distanceCellpointToParentpoint = sqrt(abs(vtkMath::Distance2BetweenPoints(parentPoint, cellPoint)));
			double distanceCellpointToLine = sqrt(abs(vtkLine::DistanceToLine(cellPoint, startPointArray, targetPointArray))); // negative square possible
			// from DistanceToLine

			if ((parentDistance + distanceCellpointToParentpoint + distanceCellpointToLine) <
				(distanceMatrix->vtkDataArray::GetComponent(cellId, 1))) {
				distanceMatrixArray[0] = parentID;
				distanceMatrixArray[1] = parentDistance + distanceCellpointToParentpoint + distanceCellpointToLine;
				distanceMatrix->InsertTuple(cellId, distanceMatrixArray);

				if (arrayIDs.find(cellId) == arrayIDs.end()) {
					possibleCellsArray[0] = cellId;
					possibleCellsArray[1] = parentDistance + distanceCellpointToParentpoint + distanceCellpointToLine;
					possibleCells->InsertNextTuple(possibleCellsArray);
					arrayIDs.insert(cellId);
				}
				else {
					vtkIdType posInArray = Methods::isIDInsideArrayLinear(possibleCells, cellId);

					if (possibleCells->GetComponent(posInArray, 1) > parentDistance + distanceCellpointToParentpoint + distanceCellpointToLine) {
						possibleCells->SetComponent(posInArray, 1, parentDistance + distanceCellpointToParentpoint + distanceCellpointToLine);
					}
				}
			}
			else if ((parentDistance + distanceCellpointToParentpoint + distanceCellpointToLine) == (distanceMatrix->vtkDataArray::GetComponent(cellId, 1))) {
				double newparentTargetDistance = sqrt(abs(vtkLine::DistanceToLine(centrePoints->GetPoint(parentID), startPointArray, targetPointArray)));
				double oldparentTargetDistance = sqrt(abs(vtkLine::DistanceToLine(centrePoints->GetPoint(distanceMatrix->vtkDataArray::GetComponent(cellId, 0)), startPointArray, targetPointArray)));

				if (newparentTargetDistance < oldparentTargetDistance) {
					distanceMatrixArray[0] = parentID;
					distanceMatrixArray[1] = parentDistance + distanceCellpointToParentpoint + distanceCellpointToLine;
					distanceMatrix->InsertTuple(cellId, distanceMatrixArray);

					if (arrayIDs.find(cellId) == arrayIDs.end()) {
						possibleCellsArray[0] = cellId;
						possibleCellsArray[1] = parentDistance + distanceCellpointToParentpoint + distanceCellpointToLine;
						possibleCells->InsertNextTuple(possibleCellsArray);
						arrayIDs.insert(cellId);
					}
					else {
						vtkIdType posInArray = Methods::isIDInsideArrayLinear(possibleCells, cellId);
						if (possibleCells->GetComponent(posInArray, 1) >
							parentDistance + distanceCellpointToParentpoint + distanceCellpointToLine) {
							possibleCells->SetComponent(posInArray, 1, parentDistance + distanceCellpointToParentpoint + distanceCellpointToLine);
						}
					}
				}
			}
		}

		if (parentID == targetPointId)
			exit = true;

		if (possibleCells->GetNumberOfTuples() == 0) {
			cout << "Target point can't be arrived" << endl;

			if (ConfigGeneral::debug) {
				data.getVtkData()->GetCellData()->AddArray(distanceMap);
				Methods::writeIntermediateData(data, "debug");
				data.getVtkData()->GetCellData()->vtkFieldData::RemoveArray("DistanceMap");
			}

			vtkSmartPointer<vtkIdList> path0 = vtkSmartPointer<vtkIdList>::New();
			return path0;
		}
		vtkSortDataArray::SortArrayByComponent(possibleCells, 1);
	}

	vtkIdType pointId = targetPointId;


	tempPath.push_front(pointId);

	// backward search
	do {
		pointId = distanceMatrix->vtkDataArray::GetComponent(pointId, 0);
		tempPath.push_front(pointId);
	} while (pointId != startPointId);


	long size = tempPath.size();

	vtkSmartPointer<vtkIdList> path = vtkSmartPointer<vtkIdList>::New();
	for (long i = 0; i < size; i++) {
		path->InsertNextId(tempPath.front());
		tempPath.pop_front();
	}


	return path;
} // Methods::pathSearch

/*! Function for searching a path closed to another path. It used a modificated dijkstra algorithm.
 \param data Pointer to the orginal mesh.
 \param startPoint Start point of the path as Vector (x,y,z).
 \param targetPoint Target point of the path as Vector (x,y,z).
 \param corrPath The corresponding path for the path search.
 \param atrium In which atrium the path would be searched (right, left, leftEpi, leftEndo or in both atrial).
 \return VtkIdList that conatins the cell-ids of the path.
 */
vtkSmartPointer<vtkIdList> Methods::pathSearchNearToPath(DataFormat &data, vector<double> startPoint, vector<double> targetPoint, vtkSmartPointer<vtkIdList> corrPath, string atrium) {
	vtkSmartPointer<vtkPoints> centrePoints = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkUnstructuredGrid> centrePointsUGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

	if ((atrium.compare("surfaceEndo") == 0) || (atrium.compare("surfaceEpi") == 0)) {
		centrePoints = data.getCentrePointsSurface()->GetPoints();
		centrePointsUGrid = data.getCentrePointsSurface();
	}
	else {
		centrePoints = data.getCentrePoints()->GetPoints();
		centrePointsUGrid = data.getCentrePoints();
	}

	list<vtkIdType> tempPath;
	set<vtkIdType> arrayIDs;

	vtkSmartPointer<vtkDoubleArray> possibleCells = vtkSmartPointer<vtkDoubleArray>::New();
	possibleCells->SetNumberOfComponents(2);
	possibleCells->SetComponentName(0, "CellID");
	possibleCells->SetComponentName(1, "Distance to the starting point");

	// initialize
	vtkSmartPointer<vtkDoubleArray> distanceMatrix = vtkSmartPointer<vtkDoubleArray>::New();
	distanceMatrix->SetName("Distance");
	distanceMatrix->SetNumberOfComponents(2);
	distanceMatrix->SetNumberOfTuples(centrePoints->GetNumberOfPoints());
	distanceMatrix->SetComponentName(0, "CellID parent");
	distanceMatrix->SetComponentName(1, "Distance to the starting point");
	distanceMatrix->FillComponent(0, -1);
	distanceMatrix->FillComponent(1, std::numeric_limits<double>::max());

	vtkIdType startPointId = -1;
	vtkIdType targetPointId = -1;
	double startPointArray[3] = { 0, 0, 0 };
	double targetPointArray[3] = { 0, 0, 0 };
	if (atrium.compare("right") == 0) {
		startPointId = findClosedPointIdinMaterialInRight(data, startPoint);
		targetPointId = findClosedPointIdinMaterialInRight(data, targetPoint);
		startPointArray[0] = centrePoints->GetPoint(startPointId)[0];
		startPointArray[1] = centrePoints->GetPoint(startPointId)[1];
		startPointArray[2] = centrePoints->GetPoint(startPointId)[2];
		targetPointArray[0] = centrePoints->GetPoint(targetPointId)[0];
		targetPointArray[1] = centrePoints->GetPoint(targetPointId)[1];
		targetPointArray[2] = centrePoints->GetPoint(targetPointId)[2];
	}
	else if (atrium.compare("left") == 0) {
		startPointId = findClosedPointIdinMaterialInLeft(data, startPoint);
		targetPointId = findClosedPointIdinMaterialInLeft(data, targetPoint);
		startPointArray[0] = centrePoints->GetPoint(startPointId)[0];
		startPointArray[1] = centrePoints->GetPoint(startPointId)[1];
		startPointArray[2] = centrePoints->GetPoint(startPointId)[2];
		targetPointArray[0] = centrePoints->GetPoint(targetPointId)[0];
		targetPointArray[1] = centrePoints->GetPoint(targetPointId)[1];
		targetPointArray[2] = centrePoints->GetPoint(targetPointId)[2];
	}
	else if (atrium.compare("leftEpi") == 0) {
		startPointId = findClosedPointIdinMaterialInLeftEpi(data, startPoint);
		targetPointId = findClosedPointIdinMaterialInLeftEpi(data, targetPoint);
		startPointArray[0] = centrePoints->GetPoint(startPointId)[0];
		startPointArray[1] = centrePoints->GetPoint(startPointId)[1];
		startPointArray[2] = centrePoints->GetPoint(startPointId)[2];
		targetPointArray[0] = centrePoints->GetPoint(targetPointId)[0];
		targetPointArray[1] = centrePoints->GetPoint(targetPointId)[1];
		targetPointArray[2] = centrePoints->GetPoint(targetPointId)[2];
	}
	else if (atrium.compare("leftEndo") == 0) {
		startPointId = findClosedPointIdinMaterialInLeftEndo(data, startPoint);
		targetPointId = findClosedPointIdinMaterialInLeftEndo(data, targetPoint);
		startPointArray[0] = centrePoints->GetPoint(startPointId)[0];
		startPointArray[1] = centrePoints->GetPoint(startPointId)[1];
		startPointArray[2] = centrePoints->GetPoint(startPointId)[2];
		targetPointArray[0] = centrePoints->GetPoint(targetPointId)[0];
		targetPointArray[1] = centrePoints->GetPoint(targetPointId)[1];
		targetPointArray[2] = centrePoints->GetPoint(targetPointId)[2];
	}
	else if (atrium.compare("surfaceEndo") == 0) {
		DataFormat surface;
		surface.setVtkData(data.getSurfacewithNormals());
		surface.setCentrePoints(data.getCentrePointsSurface());

		startPointId = findClosedPointIdinMaterialInEndo(surface, startPoint);
		targetPointId = findClosedPointIdinMaterialInEndo(surface, targetPoint);
		startPointArray[0] = centrePoints->GetPoint(startPointId)[0];
		startPointArray[1] = centrePoints->GetPoint(startPointId)[1];
		startPointArray[2] = centrePoints->GetPoint(startPointId)[2];
		targetPointArray[0] = centrePoints->GetPoint(targetPointId)[0];
		targetPointArray[1] = centrePoints->GetPoint(targetPointId)[1];
		targetPointArray[2] = centrePoints->GetPoint(targetPointId)[2];
	}
	else if (atrium.compare("surfaceEpi") == 0) {
		DataFormat surface;
		surface.setVtkData(data.getSurfacewithNormals());
		surface.setCentrePoints(data.getCentrePointsSurface());

		startPointId = findClosedPointIdinMaterialInEpi(surface, startPoint);
		targetPointId = findClosedPointIdinMaterialInEpi(surface, targetPoint);
		startPointArray[0] = centrePoints->GetPoint(startPointId)[0];
		startPointArray[1] = centrePoints->GetPoint(startPointId)[1];
		startPointArray[2] = centrePoints->GetPoint(startPointId)[2];
		targetPointArray[0] = centrePoints->GetPoint(targetPointId)[0];
		targetPointArray[1] = centrePoints->GetPoint(targetPointId)[1];
		targetPointArray[2] = centrePoints->GetPoint(targetPointId)[2];
	}
	else {
		startPointArray[0] = startPoint.at(0);
		startPointArray[1] = startPoint.at(1);
		startPointArray[2] = startPoint.at(2);
		targetPointArray[0] = targetPoint.at(0);
		targetPointArray[1] = targetPoint.at(1);
		targetPointArray[2] = targetPoint.at(2);

		// cout << "Point was not test of in left or right atrium"<< endl;
		startPointId = centrePointsUGrid->FindPoint(startPointArray);
		targetPointId = centrePointsUGrid->FindPoint(targetPointArray);
	}


	if ((startPointId < 0) || (targetPointId < 0)) {
		cout << "Point was not found" << endl;
		writeIntermediateData(data, "error");
		exit(1);
	}

	if (startPointId == targetPointId) {
		vtkSmartPointer<vtkIdList> path = vtkSmartPointer<vtkIdList>::New();
		path->InsertNextId(startPointId);

		return path;
	}

	vtkSmartPointer<vtkUnstructuredGrid> corresPathList = vtkSmartPointer<vtkUnstructuredGrid>::New();
	vtkSmartPointer<vtkPoints> pointsCorresPath = vtkSmartPointer<vtkPoints>::New();

	for (vtkIdType i = 0; i < corrPath->GetNumberOfIds(); i++) {
		double point[3] = { data.getCentrePoints()->GetPoint(corrPath->GetId(i))[0], data.getCentrePoints()->GetPoint(corrPath->GetId(i))[1], data.getCentrePoints()->GetPoint(corrPath->GetId(i))[2] };
		pointsCorresPath->InsertNextPoint(point);
	}
	corresPathList->SetPoints(pointsCorresPath);

	// starting point insert

	vtkIdType currentCellId = startPointId;

	double distanceMatrixArray[2];
	distanceMatrixArray[0] = 0;
	distanceMatrixArray[1] = 0;
	distanceMatrix->InsertTuple(currentCellId, distanceMatrixArray);

	double possibleCellsArray[2];
	possibleCellsArray[0] = currentCellId;
	possibleCellsArray[1] = 0;
	possibleCells->InsertTuple(0, possibleCellsArray);
	arrayIDs.insert(currentCellId);


	vtkSmartPointer<vtkDoubleArray> distanceMap = vtkSmartPointer<vtkDoubleArray>::New();
	distanceMap->SetName("DistanceMap");
	distanceMap->SetNumberOfComponents(1);
	distanceMap->SetNumberOfTuples(centrePoints->GetNumberOfPoints());
	distanceMap->FillComponent(0, 0);

	int count = 0;
	bool exit = false;
	while (!exit) {
		vtkIdType parentID = possibleCells->vtkDataArray::GetComponent(0, 0);
		distanceMap->SetComponent(parentID, 0, count);
		count++;
		possibleCells->RemoveFirstTuple();

		arrayIDs.erase(parentID);

		double parentDistance = distanceMatrix->GetComponent(parentID, 1);

		vtkSmartPointer<vtkIdList> cellPointIds = vtkSmartPointer<vtkIdList>::New();
		vtkSmartPointer<vtkIdList> currentNeighborCellIds = vtkSmartPointer<vtkIdList>::New();


		if ((atrium.compare("surfaceEndo") == 0) || (atrium.compare("surfaceEpi") == 0)) {
			data.getSurfacewithNormals()->GetCellPoints(parentID, cellPointIds);
		}
		else {
			data.getVtkData()->GetCellPoints(parentID, cellPointIds);
		}


		for (vtkIdType i = 0; i < cellPointIds->GetNumberOfIds(); i++) {
			vtkSmartPointer<vtkIdList> PointIdList = vtkSmartPointer<vtkIdList>::New();
			PointIdList->InsertNextId(cellPointIds->GetId(i));

			vtkSmartPointer<vtkIdList> tempNeighborCellIds = vtkSmartPointer<vtkIdList>::New();

			if ((atrium.compare("surfaceEndo") == 0) || (atrium.compare("surfaceEpi") == 0)) {
				data.getSurfacewithNormals()->GetCellNeighbors(parentID, PointIdList, tempNeighborCellIds);
			}
			else {
				data.getVtkData()->GetCellNeighbors(parentID, PointIdList, tempNeighborCellIds);
			}

			for (vtkIdType j = 0; j < tempNeighborCellIds->GetNumberOfIds(); j++) {
				int material = -1;

				if ((atrium.compare("surfaceEndo") == 0) || (atrium.compare("surfaceEpi") == 0)) {
					material = data.getSurfacewithNormals()->GetCellData()->GetArray("Material")->GetComponent(tempNeighborCellIds->GetId(j), 0);
				}
				else {
					material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(tempNeighborCellIds->GetId(j), 0);
				}

				if (atrium.compare("right") == 0) {
					if (Material::isInRight(material)) {
						currentNeighborCellIds->InsertUniqueId(tempNeighborCellIds->GetId(j));
					}
				}
				else if (atrium.compare("left") == 0) {
					if (Material::isInLeft(material)) {
						currentNeighborCellIds->InsertUniqueId(tempNeighborCellIds->GetId(j));
					}
				}
				else if (atrium.compare("leftEpi") == 0) {
					if (Material::isInLeftEpi(material)) {
						currentNeighborCellIds->InsertUniqueId(tempNeighborCellIds->GetId(j));
					}
				}
				else if (atrium.compare("leftEndo") == 0) {
					if (Material::isInEndo(material)) {
						currentNeighborCellIds->InsertUniqueId(tempNeighborCellIds->GetId(j));
					}
				}
				else if (atrium.compare("surfaceEndo") == 0) {
					if (Material::isInEndo(material)) {
						currentNeighborCellIds->InsertUniqueId(tempNeighborCellIds->GetId(j));
					}
				}
				else if (atrium.compare("surfaceEpi") == 0) {
					if (Material::isInEpi(material)) {
						currentNeighborCellIds->InsertUniqueId(tempNeighborCellIds->GetId(j));
					}
				}
				else if (Material::isInLeft(material) || Material::isInRight(material)) {
					currentNeighborCellIds->InsertUniqueId(tempNeighborCellIds->GetId(j));
				}
			}
		}


		for (vtkIdType i = 0; i < currentNeighborCellIds->GetNumberOfIds(); i++) {
			vtkIdType cellId = currentNeighborCellIds->GetId(i);
			double parentPoint[3] = { centrePoints->GetPoint(parentID)[0], centrePoints->GetPoint(parentID)[1], centrePoints->GetPoint(parentID)[2] };
			double cellPoint[3] = { centrePoints->GetPoint(cellId)[0], centrePoints->GetPoint(cellId)[1], centrePoints->GetPoint(cellId)[2] };

			double distanceCellpointToParentpoint = sqrt(abs(vtkMath::Distance2BetweenPoints(parentPoint, cellPoint)));

			vtkIdType pathPointId = corresPathList->FindPoint(cellPoint);
			double pathPoint[3] = { corresPathList->GetPoint(pathPointId)[0], corresPathList->GetPoint(pathPointId)[1], corresPathList->GetPoint(pathPointId)[2] };
			double distanceCelltoPath = sqrt(abs(vtkMath::Distance2BetweenPoints(cellPoint, pathPoint)));

			if ((parentDistance + distanceCellpointToParentpoint + distanceCelltoPath) <
				(distanceMatrix->vtkDataArray::GetComponent(cellId, 1))) {
				distanceMatrixArray[0] = parentID;
				distanceMatrixArray[1] = parentDistance + distanceCellpointToParentpoint + distanceCelltoPath;
				distanceMatrix->InsertTuple(cellId, distanceMatrixArray);

				if (arrayIDs.find(cellId) == arrayIDs.end()) {
					possibleCellsArray[0] = cellId;
					possibleCellsArray[1] = parentDistance + distanceCellpointToParentpoint + distanceCelltoPath;
					possibleCells->InsertNextTuple(possibleCellsArray);
					arrayIDs.insert(cellId);
				}
				else {
					vtkIdType posInArray = Methods::isIDInsideArrayLinear(possibleCells, cellId);

					if (possibleCells->GetComponent(posInArray, 1) > parentDistance + distanceCellpointToParentpoint + distanceCelltoPath) {
						possibleCells->SetComponent(posInArray, 1, parentDistance + distanceCellpointToParentpoint + distanceCelltoPath);
					}
				}
			}
			else if ((parentDistance + distanceCellpointToParentpoint + distanceCelltoPath) == (distanceMatrix->vtkDataArray::GetComponent(cellId, 1))) {
				double newparentTargetDistance = sqrt(abs(vtkLine::DistanceToLine(centrePoints->GetPoint(parentID), startPointArray, targetPointArray)));
				double oldparentTargetDistance = sqrt(abs(vtkLine::DistanceToLine(centrePoints->GetPoint(distanceMatrix->vtkDataArray::GetComponent(cellId, 0)), startPointArray, targetPointArray)));

				if (newparentTargetDistance < oldparentTargetDistance) {
					distanceMatrixArray[0] = parentID;
					distanceMatrixArray[1] = parentDistance + distanceCellpointToParentpoint + distanceCelltoPath;
					distanceMatrix->InsertTuple(cellId, distanceMatrixArray);

					if (arrayIDs.find(cellId) == arrayIDs.end()) {
						possibleCellsArray[0] = cellId;
						possibleCellsArray[1] = parentDistance + distanceCellpointToParentpoint + distanceCelltoPath;
						possibleCells->InsertNextTuple(possibleCellsArray);
						arrayIDs.insert(cellId);
					}
					else {
						vtkIdType posInArray = Methods::isIDInsideArrayLinear(possibleCells, cellId);
						if (possibleCells->GetComponent(posInArray, 1) > parentDistance + distanceCellpointToParentpoint + distanceCelltoPath) {
							possibleCells->SetComponent(posInArray, 1, parentDistance + distanceCellpointToParentpoint + distanceCelltoPath);
						}
					}
				}
			}
		}


		if (parentID == targetPointId)
			exit = true;

		if (possibleCells->GetNumberOfTuples() == 0) {
			cout << "Target point can't be arrived" << endl;

			if (ConfigGeneral::debug) {
				data.getVtkData()->GetCellData()->AddArray(distanceMap);
				Methods::writeIntermediateData(data, "debug");
			}

			vtkSmartPointer<vtkIdList> path0 = vtkSmartPointer<vtkIdList>::New();
			return path0;
		}
		vtkSortDataArray::SortArrayByComponent(possibleCells, 1);
	}

	vtkIdType pointId = targetPointId;


	tempPath.push_front(pointId);

	do {
		pointId = distanceMatrix->vtkDataArray::GetComponent(pointId, 0);
		tempPath.push_front(pointId);
	} while (pointId != startPointId);

	long size = tempPath.size();

	vtkSmartPointer<vtkIdList> path = vtkSmartPointer<vtkIdList>::New();
	for (long i = 0; i < size; i++) {
		if ((atrium.compare("surfaceEndo") == 0) || (atrium.compare("surfaceEpi") == 0)) {
			vtkIdType tempId = centrePointsUGrid->FindPoint(data.getCentrePointsSurface()->GetPoint(tempPath.front())[0], data.getCentrePointsSurface()->GetPoint(tempPath.front())[1], data.getCentrePointsSurface()->GetPoint(tempPath.front())[2]);
			path->InsertNextId(tempId);
		}
		else {
			path->InsertNextId(tempPath.front());
		}
		tempPath.pop_front();
	}


	return path;
} // Methods::pathSearchNearToPath

/*! Function for searching a path in one tissue class. It used a modificated dijkstra algorithm.
 \param data Pointer to the orginal mesh.
 \param startPoint Start point of the path as Vector (x,y,z).
 \param targetPoint Target point of the path as Vector (x,y,z).
 \param inMaterial1 The tissue clas in which the parth would be searched.
 \param atrium In which atrium the path would be searched (right, left, leftEpi, leftEndo or in both atrial).
 \return VtkIdList that conatins the cell-ids of the path.
 */

vtkSmartPointer<vtkIdList> Methods::mapPathFromSurfaceToVolume(DataFormat &data, vtkSmartPointer<vtkIdList> path) {
	vtkSmartPointer<vtkIdList> returnPath = vtkSmartPointer<vtkIdList>::New();

	for (vtkIdType i = 0; i < path->GetNumberOfIds(); i++) {
		vtkIdType tempId = data.getCentrePoints()->FindPoint(data.getCentrePointsSurface()->GetPoint(path->GetId(i)));
		returnPath->InsertNextId(tempId);
	}

	return returnPath;
}

// Methods::mapPathFromSurfaceToVolume

/*! Function for searching a path corresponding path in the volumatric mesh.
 \param data Pointer to the orginal mesh.
 \param path conatins the cell-ids of the path on the surface.
 \return VtkIdList that conatins the cell-ids of the path in the volumatric mesh.
 */

vtkSmartPointer<vtkIdList> Methods::pathSearch(DataFormat &data, vector<double> startPoint, vector<double> targetPoint, Material::Mat inMaterial1, string atrium) {
	return pathSearch(data, startPoint, targetPoint, inMaterial1, inMaterial1, inMaterial1, inMaterial1, inMaterial1, inMaterial1, atrium);
}

/*! Function for searching a path in two tissue class. It used a modificated dijkstra algorithm.
 \param data Pointer to the orginal mesh.
 \param startPoint Start point of the path as Vector (x,y,z).
 \param targetPoint Target point of the path as Vector (x,y,z).
 \param inMaterial1 The first tissue clas in which the parth would be searched.
 \param inMaterial2 The second tissue clas in which the parth would be searched.
 \param atrium In which atrium the path would be searched (right, left, leftEpi, leftEndo or in both atrial).
 \return VtkIdList that conatins the cell-ids of the path.
 */
vtkSmartPointer<vtkIdList> Methods::pathSearch(DataFormat &data, vector<double> startPoint, vector<double> targetPoint, Material::Mat inMaterial1, Material::Mat inMaterial2, string atrium) {
	return pathSearch(data, startPoint, targetPoint, inMaterial1, inMaterial2, inMaterial1, inMaterial1, inMaterial1, inMaterial1, atrium);
}

/*! Function for searching a path in three tissue class. It used a modificated dijkstra algorithm.
 \param data Pointer to the orginal mesh.
 \param startPoint Start point of the path as Vector (x,y,z).
 \param targetPoint Target point of the path as Vector (x,y,z).
 \param inMaterial1 The first tissue clas in which the parth would be searched.
 \param inMaterial2 The second tissue clas in which the parth would be searched.
 \param inMaterial3 The third tissue clas in which the parth would be searched.
 \param atrium In which atrium the path would be searched (right, left, leftEpi, leftEndo or in both atrial).
 \return VtkIdList that conatins the cell-ids of the path.
 */
vtkSmartPointer<vtkIdList> Methods::pathSearch(DataFormat &data, vector<double> startPoint, vector<double> targetPoint, Material::Mat inMaterial1, Material::Mat inMaterial2, Material::Mat inMaterial3, string atrium) {
	return pathSearch(data, startPoint, targetPoint, inMaterial1, inMaterial2, inMaterial3, inMaterial1, inMaterial1, inMaterial1, atrium);
}

/*! Function for searching a path in four tissue class. It used a modificated dijkstra algorithm.
 \param data Pointer to the orginal mesh.
 \param startPoint Start point of the path as Vector (x,y,z).
 \param targetPoint Target point of the path as Vector (x,y,z).
 \param inMaterial1 The first tissue clas in which the parth would be searched.
 \param inMaterial2 The second tissue clas in which the parth would be searched.
 \param inMaterial3 The third tissue clas in which the parth would be searched.
 \param inMaterial4 The fourth tissue clas in which the parth would be searched.
 \param atrium In which atrium the path would be searched (right, left, leftEpi, leftEndo or in both atrial).
 \return VtkIdList that conatins the cell-ids of the path.
 */
vtkSmartPointer<vtkIdList> Methods::pathSearch(DataFormat &data, vector<double> startPoint, vector<double> targetPoint, Material::Mat inMaterial1, Material::Mat inMaterial2, Material::Mat inMaterial3, Material::Mat inMaterial4, string atrium) {
	return pathSearch(data, startPoint, targetPoint, inMaterial1, inMaterial2, inMaterial3, inMaterial4, inMaterial1, inMaterial1, atrium);
}

/*! Function for searching a path in five tissue class. It used a modificated dijkstra algorithm.
 \param data Pointer to the orginal mesh.
 \param startPoint Start point of the path as Vector (x,y,z).
 \param targetPoint Target point of the path as Vector (x,y,z).
 \param inMaterial1 The first tissue clas in which the parth would be searched.
 \param inMaterial2 The second tissue clas in which the parth would be searched.
 \param inMaterial3 The third tissue clas in which the parth would be searched.
 \param inMaterial4 The fourth tissue clas in which the parth would be searched.
 \param inMaterial5 The fifth tissue clas in which the parth would be searched.
 \param atrium In which atrium the path would be searched (right, left, leftEpi, leftEndo or in both atrial).
 \return VtkIdList that conatins the cell-ids of the path.
 */
vtkSmartPointer<vtkIdList> Methods::pathSearch(DataFormat &data, vector<double> startPoint, vector<double> targetPoint, Material::Mat inMaterial1, Material::Mat inMaterial2, Material::Mat inMaterial3, Material::Mat inMaterial4, Material::Mat inMaterial5, string atrium) {
	return pathSearch(data, startPoint, targetPoint, inMaterial1, inMaterial2, inMaterial3, inMaterial4, inMaterial5, inMaterial1, atrium);
}

/*! Function for searching a path in six tissue class. It used a modificated dijkstra algorithm.
 \param data Pointer to the orginal mesh.
 \param startPoint Start point of the path as Vector (x,y,z).
 \param targetPoint Target point of the path as Vector (x,y,z).
 \param inMaterial1 The first tissue clas in which the parth would be searched.
 \param inMaterial2 The second tissue clas in which the parth would be searched.
 \param inMaterial3 The third tissue clas in which the parth would be searched.
 \param inMaterial4 The fourth tissue clas in which the parth would be searched.
 \param inMaterial5 The fifth tissue clas in which the parth would be searched.
 \param inMaterial6 The sixth tissue clas in which the parth would be searched.
 \param atrium In which atrium the path would be searched (right, left, leftEpi, leftEndo or in both atrial).
 \return VtkIdList that conatins the cell-ids of the path.
 */
vtkSmartPointer<vtkIdList> Methods::pathSearch(DataFormat &data, vector<double> startPoint, vector<double> targetPoint, Material::Mat inMaterial1, Material::Mat inMaterial2, Material::Mat inMaterial3, Material::Mat inMaterial4, Material::Mat inMaterial5, Material::Mat inMaterial6, string atrium) {
	vtkSmartPointer<vtkPoints> centrePoints = data.getCentrePoints()->GetPoints();
	vtkSmartPointer<vtkUnstructuredGrid> centrePointsUGrid = data.getCentrePoints();

	list<vtkIdType> tempPath;
	vtkSmartPointer<vtkIdList> path = vtkSmartPointer<vtkIdList>::New();

	set<vtkIdType> arrayIDs;

	vtkSmartPointer<vtkDoubleArray> possibleCells = vtkSmartPointer<vtkDoubleArray>::New();

	possibleCells->SetNumberOfComponents(2);
	possibleCells->SetComponentName(0, "CellID");
	possibleCells->SetComponentName(1, "Distance to the starting point");

	// initialize
	vtkSmartPointer<vtkDoubleArray> distanceMatrix = vtkSmartPointer<vtkDoubleArray>::New();
	distanceMatrix->SetName("Distance");
	distanceMatrix->SetNumberOfComponents(2);
	distanceMatrix->SetNumberOfTuples(centrePoints->GetNumberOfPoints());
	distanceMatrix->SetComponentName(0, "CellID parent");
	distanceMatrix->SetComponentName(1, "Distance to the starting point");
	distanceMatrix->FillComponent(0, -1);
	distanceMatrix->FillComponent(1, std::numeric_limits<double>::max());

	vtkIdType startPointId = -1;
	vtkIdType targetPointId = -1;
	double startPointArray[3] = { 0, 0, 0 };
	double targetPointArray[3] = { 0, 0, 0 };
	if (atrium.compare("right") == 0) {
		startPointId = findClosedPointIdinMaterialInRight(data, startPoint);
		targetPointId = findClosedPointIdinMaterialInRight(data, targetPoint);
		startPointArray[0] = centrePoints->GetPoint(startPointId)[0];
		startPointArray[1] = centrePoints->GetPoint(startPointId)[1];
		startPointArray[2] = centrePoints->GetPoint(startPointId)[2];
		targetPointArray[0] = centrePoints->GetPoint(targetPointId)[0];
		targetPointArray[1] = centrePoints->GetPoint(targetPointId)[1];
		targetPointArray[2] = centrePoints->GetPoint(targetPointId)[2];
	}
	else if (atrium.compare("left") == 0) {
		startPointId = findClosedPointIdinMaterialInLeft(data, startPoint);
		targetPointId = findClosedPointIdinMaterialInLeft(data, targetPoint);
		startPointArray[0] = centrePoints->GetPoint(startPointId)[0];
		startPointArray[1] = centrePoints->GetPoint(startPointId)[1];
		startPointArray[2] = centrePoints->GetPoint(startPointId)[2];
		targetPointArray[0] = centrePoints->GetPoint(targetPointId)[0];
		targetPointArray[1] = centrePoints->GetPoint(targetPointId)[1];
		targetPointArray[2] = centrePoints->GetPoint(targetPointId)[2];
	}
	else {
		startPointArray[0] = startPoint.at(0);
		startPointArray[1] = startPoint.at(1);
		startPointArray[2] = startPoint.at(2);
		targetPointArray[0] = targetPoint.at(0);
		targetPointArray[1] = targetPoint.at(1);
		targetPointArray[2] = targetPoint.at(2);
		cout << "Point was not test if it in left or right atrium" << endl;
		startPointId = centrePointsUGrid->FindPoint(startPointArray);
		targetPointId = centrePointsUGrid->FindPoint(targetPointArray);
	}

	if ((startPointId < 0) || (targetPointId < 0)) {
		cout << "Point was not found" << endl;
		if (ConfigGeneral::debug) {
			writeIntermediateData(data, "error");
		}
		exit(1);
	}

	if (startPointId == targetPointId) {
		path->Reset();
		return path;
	}

	// starting point insert
	vtkIdType currentCellId = startPointId;

	double distanceMatrixArray[2];
	distanceMatrixArray[0] = 0;
	distanceMatrixArray[1] = 0;
	distanceMatrix->InsertTuple(currentCellId, distanceMatrixArray);

	double possibleCellsArray[2];
	possibleCellsArray[0] = currentCellId;
	possibleCellsArray[1] = 0;
	possibleCells->InsertTuple(0, possibleCellsArray);
	arrayIDs.insert(currentCellId);

	vtkSmartPointer<vtkDoubleArray> distanceMap = vtkSmartPointer<vtkDoubleArray>::New();
	distanceMap->SetName("DistanceMap");
	distanceMap->SetNumberOfComponents(1);
	distanceMap->SetNumberOfTuples(centrePoints->GetNumberOfPoints());
	distanceMap->FillComponent(0, 0);

	int count = 0;
	bool exit = false;
	while (!exit) {
		vtkIdType parentID = possibleCells->vtkDataArray::GetComponent(0, 0);
		distanceMap->SetComponent(parentID, 0, count);
		count++;

		possibleCells->RemoveFirstTuple();

		arrayIDs.erase(parentID);

		double parentDistance = distanceMatrix->GetComponent(parentID, 1);

		vtkSmartPointer<vtkIdList> cellPointIds = vtkSmartPointer<vtkIdList>::New();
		vtkSmartPointer<vtkIdList> currentNeighborCellIds = vtkSmartPointer<vtkIdList>::New();

		data.getVtkData()->GetCellPoints(parentID, cellPointIds);

		for (vtkIdType i = 0; i < cellPointIds->GetNumberOfIds(); i++) {
			vtkSmartPointer<vtkIdList> PointIdList = vtkSmartPointer<vtkIdList>::New();
			PointIdList->InsertNextId(cellPointIds->GetId(i));

			vtkSmartPointer<vtkIdList> tempNeighborCellIds = vtkSmartPointer<vtkIdList>::New();

			data.getVtkData()->GetCellNeighbors(parentID, PointIdList, tempNeighborCellIds);

			for (vtkIdType j = 0; j < tempNeighborCellIds->GetNumberOfIds(); j++) {
				if (!data.getDoubleLayer()) {
					int material = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(tempNeighborCellIds->GetId(j), 0);
					if ((material == inMaterial1) || (material == inMaterial2) || (material == inMaterial3) ||
						(material == inMaterial4) || (material == inMaterial5) || (material == inMaterial6)) {
						currentNeighborCellIds->InsertUniqueId(tempNeighborCellIds->GetId(j));
					}
				}

				int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(tempNeighborCellIds->GetId(j), 0);
				if ((material == inMaterial1) || (material == inMaterial2) || (material == inMaterial3) ||
					(material == inMaterial4) || (material == inMaterial5) || (material == inMaterial6)) {
					currentNeighborCellIds->InsertUniqueId(tempNeighborCellIds->GetId(j));
				}
			}
		}


		for (vtkIdType i = 0; i < currentNeighborCellIds->GetNumberOfIds(); i++) {
			vtkIdType cellId = currentNeighborCellIds->GetId(i);
			double parentPoint[3] = { centrePoints->GetPoint(parentID)[0], centrePoints->GetPoint(parentID)[1], centrePoints->GetPoint(parentID)[2] };
			double cellPoint[3] = { centrePoints->GetPoint(cellId)[0], centrePoints->GetPoint(cellId)[1], centrePoints->GetPoint(cellId)[2] };

			double distanceCellpointToParentpoint = sqrt(abs(vtkMath::Distance2BetweenPoints(parentPoint, cellPoint)));
			double distanceCellpointToLine = sqrt(abs(vtkLine::DistanceToLine(cellPoint, startPointArray, targetPointArray))); // negative square possible
			// from DistanceToLine

			if ((parentDistance + distanceCellpointToParentpoint + distanceCellpointToLine) <
				(distanceMatrix->vtkDataArray::GetComponent(cellId, 1))) {
				distanceMatrixArray[0] = parentID;
				distanceMatrixArray[1] = parentDistance + distanceCellpointToParentpoint + distanceCellpointToLine;
				distanceMatrix->InsertTuple(cellId, distanceMatrixArray);

				if (arrayIDs.find(cellId) == arrayIDs.end()) {
					possibleCellsArray[0] = cellId;
					possibleCellsArray[1] = parentDistance + distanceCellpointToParentpoint + distanceCellpointToLine;
					possibleCells->InsertNextTuple(possibleCellsArray);
					arrayIDs.insert(cellId);
				}
				else {
					vtkIdType posInArray = Methods::isIDInsideArrayLinear(possibleCells, cellId);

					if (possibleCells->GetComponent(posInArray, 1) > parentDistance + distanceCellpointToParentpoint + distanceCellpointToLine) {
						possibleCells->SetComponent(posInArray, 1, parentDistance + distanceCellpointToParentpoint + distanceCellpointToLine);
					}
				}
			}
			else if ((parentDistance + distanceCellpointToParentpoint + distanceCellpointToLine) == (distanceMatrix->vtkDataArray::GetComponent(cellId, 1))) {
				double newparentTargetDistance = sqrt(abs(vtkLine::DistanceToLine(centrePoints->GetPoint(parentID), startPointArray, targetPointArray)));
				double oldparentTargetDistance = sqrt(abs(vtkLine::DistanceToLine(centrePoints->GetPoint(distanceMatrix->vtkDataArray::GetComponent(cellId, 0)), startPointArray, targetPointArray)));

				if (newparentTargetDistance < oldparentTargetDistance) {
					distanceMatrixArray[0] = parentID;
					distanceMatrixArray[1] = parentDistance + distanceCellpointToParentpoint + distanceCellpointToLine;
					distanceMatrix->InsertTuple(cellId, distanceMatrixArray);

					if (arrayIDs.find(cellId) == arrayIDs.end()) {
						possibleCellsArray[0] = cellId;
						possibleCellsArray[1] = parentDistance + distanceCellpointToParentpoint + distanceCellpointToLine;
						possibleCells->InsertNextTuple(possibleCellsArray);
						arrayIDs.insert(cellId);
					}
					else {
						vtkIdType posInArray = Methods::isIDInsideArrayLinear(possibleCells, cellId);
						if (possibleCells->GetComponent(posInArray, 1) >
							parentDistance + distanceCellpointToParentpoint + distanceCellpointToLine) {
							possibleCells->SetComponent(posInArray, 1, parentDistance + distanceCellpointToParentpoint + distanceCellpointToLine);
						}
					}
				}
			}
		}

		if (parentID == targetPointId) {
			exit = true;
		}
		else if (possibleCells->GetNumberOfTuples() == 0) {
			cout << "Target point can't be arrived" << endl;
			if (ConfigGeneral::debug) {
				data.getVtkData()->GetCellData()->AddArray(distanceMap);
				Methods::writeIntermediateData(data, "pathSearchDebug");
				data.getVtkData()->GetCellData()->vtkFieldData::RemoveArray("DistanceMap");
			}
			path->Reset();
			return path;
		}
		vtkSortDataArray::SortArrayByComponent(possibleCells, 1);
	}

	vtkIdType pointId = targetPointId;


	tempPath.push_front(pointId);

	do {
		pointId = distanceMatrix->vtkDataArray::GetComponent(pointId, 0);
		tempPath.push_front(pointId);
	} while (pointId != startPointId);


	long size = tempPath.size();

	for (long i = 0; i < size; i++) {
		path->InsertNextId(tempPath.front());
		tempPath.pop_front();
	}


	return path;
} // Methods::pathSearch

/*! Function for searching a path over a plane, that define the start, target and plane point. It used a modificated
 dijkstra algorithm.
 \param data Pointer to the orginal mesh.
 \param startPoint Start point of the path as Vector (x,y,z).
 \param targetPoint Target point of the path as Vector (x,y,z).
 \param planePoint Plane point for the plane as Vector (x,y,z).
 \return VtkIdList that conatins the cell-ids of the path.
 */
vtkSmartPointer<vtkIdList> Methods::pathSearchOverPlane(DataFormat &data, vector<double> startPoint, vector<double> targetPoint, vector<double> planePoint) {
	vtkSmartPointer<vtkPoints> centrePoints = data.getCentrePoints()->GetPoints();
	vtkSmartPointer<vtkUnstructuredGrid> centrePointsUGrid = data.getCentrePoints();
	double startPointArray[3] = { startPoint.at(0), startPoint.at(1), startPoint.at(2) };
	double targetPointArray[3] = { targetPoint.at(0), targetPoint.at(1), targetPoint.at(2) };
	double planePointArray[3] = { planePoint.at(0), planePoint.at(1), planePoint.at(2) };
	double normalVectorPlane[3];
	double vectorTargetStartpoint[3];
	double vectorPlaneStartpoint[3];

	vtkMath::Subtract(targetPointArray, startPointArray, vectorTargetStartpoint);
	vtkMath::Subtract(planePointArray, startPointArray, vectorPlaneStartpoint);
	vtkMath::Cross(vectorTargetStartpoint, vectorPlaneStartpoint, normalVectorPlane);
	vtkMath::Normalize(normalVectorPlane);


	vtkSmartPointer<vtkDoubleArray> distanceMatrix = vtkSmartPointer<vtkDoubleArray>::New();

	list<vtkIdType> tempPath;
	vtkSmartPointer<vtkIdList> path = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkDoubleArray> possibleCells = vtkSmartPointer<vtkDoubleArray>::New();
	possibleCells->SetNumberOfComponents(2);
	possibleCells->SetComponentName(0, "CellID");
	possibleCells->SetComponentName(1, "Distance to the starting point");

	set<vtkIdType> arrayIDs;

	// initialize
	distanceMatrix->SetName("Distance");
	distanceMatrix->SetNumberOfComponents(2);
	distanceMatrix->SetComponentName(0, "Cell parent");
	distanceMatrix->SetComponentName(1, "Distance to the starting point");
	distanceMatrix->SetNumberOfTuples(centrePoints->GetNumberOfPoints());
	distanceMatrix->FillComponent(0, -1);
	distanceMatrix->FillComponent(1, std::numeric_limits<double>::max());

	vtkIdType startPointId = centrePointsUGrid->FindPoint(startPoint.at(0), startPoint.at(1), startPoint.at(2));
	vtkIdType targetPointId = centrePointsUGrid->FindPoint(targetPoint.at(0), targetPoint.at(1), targetPoint.at(2));


	if ((startPointId < 0) || (targetPointId < 0)) {
		cout << "Point was not found" << endl;
		writeIntermediateData(data, "error");
		exit(1);
	}

	// starting point insert

	vtkIdType currentCellId = startPointId;

	double distanceMatrixArray[2];
	distanceMatrixArray[0] = 0;
	distanceMatrixArray[1] = 0;
	distanceMatrix->InsertTuple(currentCellId, distanceMatrixArray);

	double possibleCellsArray[2];
	possibleCellsArray[0] = currentCellId;
	possibleCellsArray[1] = 0;
	possibleCells->InsertTuple(0, possibleCellsArray);
	arrayIDs.insert(currentCellId);

	vtkSmartPointer<vtkDoubleArray> distanceMap = vtkSmartPointer<vtkDoubleArray>::New();
	distanceMap->SetName("DistanceMap");
	distanceMap->SetNumberOfComponents(1);
	distanceMap->SetNumberOfTuples(centrePoints->GetNumberOfPoints());
	distanceMap->FillComponent(0, 0);

	int count = 0;
	bool exit = false;
	while (!exit) {
		vtkIdType parentID = possibleCells->vtkDataArray::GetComponent(0, 0);
		distanceMap->SetComponent(parentID, 0, count);
		count++;

		possibleCells->RemoveFirstTuple();
		arrayIDs.erase(parentID);

		double parentDistance = distanceMatrix->GetComponent(parentID, 1);

		vtkSmartPointer<vtkIdList> cellPointIds = vtkSmartPointer<vtkIdList>::New();
		vtkSmartPointer<vtkIdList> currentNeighborCellIds = vtkSmartPointer<vtkIdList>::New();

		data.getVtkData()->GetCellPoints(parentID, cellPointIds);

		for (vtkIdType i = 0; i < cellPointIds->GetNumberOfIds(); i++) {
			vtkSmartPointer<vtkIdList> PointIdList = vtkSmartPointer<vtkIdList>::New();
			PointIdList->InsertNextId(cellPointIds->GetId(i));

			vtkSmartPointer<vtkIdList> tempNeighborCellIds = vtkSmartPointer<vtkIdList>::New();

			data.getVtkData()->GetCellNeighbors(parentID, PointIdList, tempNeighborCellIds);

			for (vtkIdType j = 0; j < tempNeighborCellIds->GetNumberOfIds(); j++) {
				currentNeighborCellIds->InsertUniqueId(tempNeighborCellIds->GetId(j));
			}
		}


		for (vtkIdType i = 0; i < currentNeighborCellIds->GetNumberOfIds(); i++) {
			vtkIdType cellId = currentNeighborCellIds->GetId(i);

			double parentPoint[3] = { centrePoints->GetPoint(parentID)[0], centrePoints->GetPoint(parentID)[1], centrePoints->GetPoint(parentID)[2] };
			double cellPoint[3] = { centrePoints->GetPoint(cellId)[0], centrePoints->GetPoint(cellId)[1], centrePoints->GetPoint(cellId)[2] };

			double distanceCellpointToParentpoint = sqrt(abs(vtkMath::Distance2BetweenPoints(parentPoint, cellPoint)));
			double distanceCellpointToPlane = sqrt(abs(vtkPlane::DistanceToPlane(cellPoint, normalVectorPlane, startPointArray)));

			if ((parentDistance + distanceCellpointToParentpoint + distanceCellpointToPlane) <
				(distanceMatrix->GetTuple(cellId)[1])) {
				double distanceMatrixArray[2];
				distanceMatrixArray[0] = parentID;
				distanceMatrixArray[1] = parentDistance + distanceCellpointToParentpoint + distanceCellpointToPlane;
				distanceMatrix->InsertTuple(cellId, distanceMatrixArray);

				if (arrayIDs.find(cellId) == arrayIDs.end()) {
					possibleCellsArray[0] = cellId;
					possibleCellsArray[1] = parentDistance + distanceCellpointToParentpoint + distanceCellpointToPlane;
					possibleCells->InsertNextTuple(possibleCellsArray);
					arrayIDs.insert(cellId);
				}
				else {
					vtkIdType posInArray = Methods::isIDInsideArrayLinear(possibleCells, cellId);
					if (possibleCells->GetComponent(posInArray, 1) > parentDistance + distanceCellpointToParentpoint + distanceCellpointToPlane) {
						possibleCells->SetComponent(posInArray, 1, parentDistance + distanceCellpointToParentpoint + distanceCellpointToPlane);
					}
				}
			}
			else if ((parentDistance + distanceCellpointToPlane) == ((distanceMatrix->GetTuple(cellId)[1]))) {
				double newparentTargetDistance = sqrt(abs(vtkPlane::DistanceToPlane(centrePoints->GetPoint(parentID), normalVectorPlane, targetPointArray)));
				double oldparentTargetDistance = sqrt(abs(vtkPlane::DistanceToPlane(centrePoints->GetPoint(distanceMatrix->GetTuple(cellId)[0]), normalVectorPlane, targetPointArray)));

				if (newparentTargetDistance < oldparentTargetDistance) {
					double distanceMatrixArray[2];
					distanceMatrixArray[0] = parentID;
					distanceMatrixArray[1] = parentDistance + distanceCellpointToParentpoint + distanceCellpointToPlane;
					distanceMatrix->InsertTuple(cellId, distanceMatrixArray);

					if (arrayIDs.find(cellId) == arrayIDs.end()) {
						possibleCellsArray[0] = cellId;
						possibleCellsArray[1] = parentDistance + distanceCellpointToParentpoint + distanceCellpointToPlane;
						possibleCells->InsertNextTuple(possibleCellsArray);
						arrayIDs.insert(cellId);
					}
					else {
						vtkIdType posInArray = Methods::isIDInsideArrayLinear(possibleCells, cellId);
						if (possibleCells->GetComponent(posInArray, 1) >
							parentDistance + distanceCellpointToParentpoint + distanceCellpointToPlane) {
							possibleCells->SetComponent(posInArray, 1, parentDistance + distanceCellpointToParentpoint + distanceCellpointToPlane);
						}
					}
				}
			}
		}

		if (parentID == targetPointId) {
			exit = true;
		}
		else if (possibleCells->GetNumberOfTuples() == 0) {
			cout << "Target point can't be arrived" << endl;
			data.getVtkData()->GetCellData()->AddArray(distanceMap);
			Methods::writeIntermediateData(data, "pathSearchDebug");
			data.getVtkData()->GetCellData()->vtkFieldData::RemoveArray("DistanceMap");
		}

		vtkSortDataArray::SortArrayByComponent(possibleCells, 1);
	}
	vtkIdType pointId = targetPointId;
	tempPath.push_front(pointId);

	do {
		pointId = distanceMatrix->vtkDataArray::GetComponent(pointId, 0);
		tempPath.push_front(pointId);
	} while (pointId != startPointId);

	long size = tempPath.size();

	for (long i = 0; i < size; i++) {
		path->InsertNextId(tempPath.front());
		tempPath.pop_front();
	}

	return path;
} // Methods::pathSearchOverPlane

/*! Initialization the DataForamt for the fieber orientation tool with all needed cell arrays.
 \param data Pointer to the orginal mesh.
 */

void Methods::init(DataFormat &data) {
	vtkSmartPointer<vtkPoints> centrePoints = data.getCentrePoints()->GetPoints();

	if (!(data.getVtkData()->GetCellData()->HasArray("RawMaterial"))) {
		vtkSmartPointer<vtkIntArray> rowMat = vtkSmartPointer<vtkIntArray>::New();
		rowMat->DeepCopy(data.getVtkData()->GetCellData()->GetArray("Material"));
		rowMat->SetName("RawMaterial");
		data.getVtkData()->GetCellData()->AddArray(rowMat);
	}

	vtkSmartPointer<vtkDoubleArray> theta = vtkSmartPointer<vtkDoubleArray>::New();
	theta->SetName("Theta");
	theta->SetNumberOfComponents(1);
	theta->SetNumberOfTuples(centrePoints->GetNumberOfPoints());
	theta->FillComponent(0, 0);

	vtkSmartPointer<vtkDoubleArray> phi = vtkSmartPointer<vtkDoubleArray>::New();
	phi->SetName("Phi");
	phi->SetNumberOfComponents(1);
	phi->SetNumberOfTuples(centrePoints->GetNumberOfPoints());
	phi->FillComponent(0, 0);

	vtkSmartPointer<vtkDoubleArray> thetaPhi = vtkSmartPointer<vtkDoubleArray>::New();
	thetaPhi->SetName("ThetaPhi");
	thetaPhi->SetNumberOfComponents(2);
	thetaPhi->SetComponentName(0, "Theta");
	thetaPhi->SetComponentName(1, "Phi");
	thetaPhi->SetNumberOfTuples(centrePoints->GetNumberOfPoints());
	thetaPhi->FillComponent(0, 0);
	thetaPhi->FillComponent(1, 0);

	vtkSmartPointer<vtkDoubleArray> diffVec = vtkSmartPointer<vtkDoubleArray>::New();
	diffVec->SetName("DifferenceVector");
	diffVec->SetNumberOfComponents(3);
	diffVec->SetComponentName(0, "Vector_comp_x");
	diffVec->SetComponentName(1, "Vector_comp_y");
	diffVec->SetComponentName(2, "Vector_comp_z");
	diffVec->SetNumberOfTuples(centrePoints->GetNumberOfPoints());
	diffVec->FillComponent(0, 0);
	diffVec->FillComponent(1, 0);
	diffVec->FillComponent(2, 0);

	vtkSmartPointer<vtkDoubleArray> baseDiffVec = vtkSmartPointer<vtkDoubleArray>::New();
	baseDiffVec->SetName("BasedPathDifferenceVector");
	baseDiffVec->SetNumberOfComponents(3);
	baseDiffVec->SetComponentName(0, "Vector_comp_x");
	baseDiffVec->SetComponentName(1, "Vector_comp_y");
	baseDiffVec->SetComponentName(2, "Vector_comp_z");
	baseDiffVec->SetNumberOfTuples(centrePoints->GetNumberOfPoints());
	baseDiffVec->FillComponent(0, 0);
	baseDiffVec->FillComponent(1, 0);
	baseDiffVec->FillComponent(2, 0);

	vtkSmartPointer<vtkDoubleArray> basePath = vtkSmartPointer<vtkDoubleArray>::New();
	basePath->SetName("basePath");
	basePath->SetNumberOfComponents(1);
	basePath->SetNumberOfTuples(centrePoints->GetNumberOfPoints());
	basePath->FillComponent(0, 0);

	vtkSmartPointer<vtkDoubleArray> status = vtkSmartPointer<vtkDoubleArray>::New();
	status->SetName("Status");
	status->SetNumberOfComponents(1);
	status->SetNumberOfTuples(centrePoints->GetNumberOfPoints());
	status->FillComponent(0, 0);



	data.getVtkData()->GetCellData()->AddArray(theta);
	data.getVtkData()->GetCellData()->AddArray(phi);
	data.getVtkData()->GetCellData()->AddArray(thetaPhi);
	data.getVtkData()->GetCellData()->AddArray(diffVec);
	data.getVtkData()->GetCellData()->AddArray(baseDiffVec);
	data.getVtkData()->GetCellData()->AddArray(basePath);
	data.getVtkData()->GetCellData()->AddArray(status);


	// if it is a surface mesh
	if (!data.getDoubleLayer()) {
		vtkSmartPointer<vtkDoubleArray> materialEndo = vtkSmartPointer<vtkDoubleArray>::New();
		materialEndo->SetName("MaterialEndo");
		materialEndo->SetNumberOfComponents(1);
		materialEndo->SetNumberOfTuples(centrePoints->GetNumberOfPoints());

		for (vtkIdType i = 0; i < data.getVtkData()->GetNumberOfCells(); i++) {
			if (data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(i, 0) == Material::Vorhof_links) {
				materialEndo->vtkDataArray::SetComponent(i, 0, Material::Vorhof_links_Endo);
			}
			else {
				materialEndo->vtkDataArray::SetComponent(i, 0, data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(i, 0));
			}
		}

		vtkSmartPointer<vtkDoubleArray> thetaEndo = vtkSmartPointer<vtkDoubleArray>::New();
		thetaEndo->SetName("ThetaEndo");
		thetaEndo->SetNumberOfComponents(1);
		thetaEndo->SetNumberOfTuples(centrePoints->GetNumberOfPoints());
		thetaEndo->FillComponent(0, 0);

		vtkSmartPointer<vtkDoubleArray> phiEndo = vtkSmartPointer<vtkDoubleArray>::New();
		phiEndo->SetName("PhiEndo");
		phiEndo->SetNumberOfComponents(1);
		phiEndo->SetNumberOfTuples(centrePoints->GetNumberOfPoints());
		phiEndo->FillComponent(0, 0);

		vtkSmartPointer<vtkDoubleArray> thetaPhiEndo = vtkSmartPointer<vtkDoubleArray>::New();
		thetaPhiEndo->SetName("ThetaPhiEndo");
		thetaPhiEndo->SetNumberOfComponents(2);
		thetaPhiEndo->SetComponentName(0, "ThetaEndo");
		thetaPhiEndo->SetComponentName(1, "PhiEndo");
		thetaPhiEndo->SetNumberOfTuples(centrePoints->GetNumberOfPoints());
		thetaPhiEndo->FillComponent(0, 0);
		thetaPhiEndo->FillComponent(1, 0);

		vtkSmartPointer<vtkDoubleArray> diffVecEndo = vtkSmartPointer<vtkDoubleArray>::New();
		diffVecEndo->SetName("DifferenceVectorEndo");
		diffVecEndo->SetNumberOfComponents(3);
		diffVecEndo->SetNumberOfTuples(centrePoints->GetNumberOfPoints());
		diffVecEndo->FillComponent(0, 0);
		diffVecEndo->FillComponent(1, 0);
		diffVecEndo->FillComponent(2, 0);

		vtkSmartPointer<vtkDoubleArray> basePathEndo = vtkSmartPointer<vtkDoubleArray>::New();
		basePathEndo->SetName("basePathEndo");
		basePathEndo->SetNumberOfComponents(1);
		basePathEndo->SetNumberOfTuples(centrePoints->GetNumberOfPoints());
		basePathEndo->FillComponent(0, 0);

		vtkSmartPointer<vtkDoubleArray> statusEndo = vtkSmartPointer<vtkDoubleArray>::New();
		statusEndo->SetName("StatusEndo");
		statusEndo->SetNumberOfComponents(1);
		statusEndo->SetNumberOfTuples(centrePoints->GetNumberOfPoints());
		statusEndo->FillComponent(0, 0);

		data.getVtkData()->GetCellData()->AddArray(materialEndo);
		data.getVtkData()->GetCellData()->AddArray(thetaEndo);
		data.getVtkData()->GetCellData()->AddArray(phiEndo);
		data.getVtkData()->GetCellData()->AddArray(thetaPhiEndo);
		data.getVtkData()->GetCellData()->AddArray(diffVecEndo);
		data.getVtkData()->GetCellData()->AddArray(basePathEndo);
		data.getVtkData()->GetCellData()->AddArray(statusEndo);

		if (!(data.getVtkData()->GetCellData()->HasArray("RawMaterialEndo"))) {
			vtkSmartPointer<vtkIntArray> rowMatEndo = vtkSmartPointer<vtkIntArray>::New();
			rowMatEndo->SetName("RawMaterialEndo");
			rowMatEndo->SetNumberOfComponents(1);
			rowMatEndo->SetNumberOfTuples(centrePoints->GetNumberOfPoints());

			for (vtkIdType i = 0; i < data.getVtkData()->GetNumberOfCells(); i++) {
				if (data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(i, 0) == Material::Vorhof_links) {
					rowMatEndo->vtkDataArray::SetComponent(i, 0, Material::Vorhof_links_Endo);
				}
				else if (data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(i, 0) == Material::Vorhof_rechts) {
					rowMatEndo->vtkDataArray::SetComponent(i, 0, Material::subendo_right_Atrial_Cardium);
				}
				else {
					rowMatEndo->vtkDataArray::SetComponent(i, 0, data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(i, 0));
				}
			}

			data.getVtkData()->GetCellData()->AddArray(rowMatEndo);
		}
	}
} // Methods::init

/*! Mark the fiber along the path in the mesh by subtracting the cell centerpoints.
 \param data Pointer to the orginal mesh.
 \param path Founded path between to points in the mesh.
 */

void Methods::pathMarker(DataFormat &data, vtkSmartPointer<vtkIdList> path) {
	vtkSmartPointer<vtkPoints> centrePoints = data.getCentrePoints()->GetPoints();

	vtkIdType currentPointId = -1;
	double tempPoint[3];

	for (long i = 0; i < path->GetNumberOfIds(); i++) {
		if (i < path->GetNumberOfIds() - 1) {
			currentPointId = path->GetId(i);
			double currentPoint[3] = { centrePoints->GetPoint(path->GetId(i))[0], centrePoints->GetPoint(path->GetId(i))[1], centrePoints->GetPoint(path->GetId(i))
			[2] };
			double nextPoint[3] = { centrePoints->GetPoint(path->GetId(i + 1))[0], centrePoints->GetPoint(path->GetId(i + 1))[1], centrePoints->GetPoint(path->GetId(i + 1))[2] };

			vtkMath::Subtract(nextPoint, currentPoint, tempPoint);
			vtkMath::Normalize(tempPoint);
		}
		else if (i == path->GetNumberOfIds() - 1) {
			currentPointId = path->GetId(i);
			double currentPoint[3] = { centrePoints->GetPoint(path->GetId(i))[0], centrePoints->GetPoint(path->GetId(i))[1], centrePoints->GetPoint(path->GetId(i))
			[2] };
			double lastPoint[3] = { centrePoints->GetPoint(path->GetId(i - 1))[0], centrePoints->GetPoint(path->GetId(i - 1))[1], centrePoints->GetPoint(path->GetId(i - 1))[2] };

			vtkMath::Subtract(currentPoint, lastPoint, tempPoint);
			vtkMath::Normalize(tempPoint);
		}
		if (!data.getDoubleLayer()) {
			data.getVtkData()->GetCellData()->GetArray("StatusEndo")->SetComponent(currentPointId, 0, 1);
			data.getVtkData()->GetCellData()->GetArray("basePathEndo")->SetComponent(currentPointId, 0, 1);
			data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->InsertTuple(currentPointId, tempPoint);
		}
		data.getVtkData()->GetCellData()->GetArray("Status")->SetComponent(currentPointId, 0, 1);
		data.getVtkData()->GetCellData()->GetArray("basePath")->SetComponent(currentPointId, 0, 1);
		data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->InsertTuple(currentPointId, tempPoint);

		data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->InsertTuple(currentPointId, tempPoint);

	}
} // Methods::pathMarker

/*! Mark the fiber along the path in one tissue class in the mesh by subtracting the cell center points.
 \param data Pointer to the orginal mesh.
 \param path Founded path between to points in the mesh.
 \param inMaterial Tissue class in that the fiber orientation would be marked.
 */

void Methods::pathMarker(DataFormat &data, vtkSmartPointer<vtkIdList> path, Material::Mat inMaterial) {
	pathMarker(data, path, inMaterial, inMaterial, inMaterial);
}

/*! Mark the fiber along the path in two tissue class in the mesh by subtracting the cell center points.
 \param data Pointer to the orginal mesh.
 \param path Founded path between to points in the mesh.
 \param inMaterial1 First tissue class in that the fiber orientation would be marked.
 \param inMaterial2 Second tissue class in that the fiber orientation would be marked.
 */
void Methods::pathMarker(DataFormat &data, vtkSmartPointer<vtkIdList> path, Material::Mat inMaterial1, Material::Mat inMaterial2) {
	pathMarker(data, path, inMaterial1, inMaterial2, inMaterial1);
}

/*! Mark the fiber along the path in three tissue class in the mesh by subtracting the cell center points.
 \param data Pointer to the orginal mesh.
 \param path Founded path between to points in the mesh.
 \param inMaterial1 First tissue class in that the fiber orientation would be marked.
 \param inMaterial2 Second tissue class in that the fiber orientation would be marked.
 \param inMaterial3 Third tissue class in that the fiber orientation would be marked.
 */
void Methods::pathMarker(DataFormat &data, vtkSmartPointer<vtkIdList> path, Material::Mat inMaterial1, Material::Mat inMaterial2, Material::Mat inMaterial3) {
	vtkSmartPointer<vtkPoints> centrePoints = data.getCentrePoints()->GetPoints();

	vtkIdType currentPointId = -1;
	double tempPoint[3];

	for (long i = 0; i < path->GetNumberOfIds(); i++) {
		if (i < path->GetNumberOfIds() - 1) {
			currentPointId = path->GetId(i);
			double currentPoint[3] = { centrePoints->GetPoint(path->GetId(i))[0], centrePoints->GetPoint(path->GetId(i))[1], centrePoints->GetPoint(path->GetId(i))[2] };
			double nextPoint[3] = { centrePoints->GetPoint(path->GetId(i + 1))[0], centrePoints->GetPoint(path->GetId(i + 1))[1], centrePoints->GetPoint(path->GetId(i + 1))[2] };

			vtkMath::Subtract(nextPoint, currentPoint, tempPoint);
			vtkMath::Normalize(tempPoint);
		}
		else if (i == path->GetNumberOfIds() - 1) {
			currentPointId = path->GetId(i);
			double currentPoint[3] = { centrePoints->GetPoint(path->GetId(i))[0], centrePoints->GetPoint(path->GetId(i))[1], centrePoints->GetPoint(path->GetId(i))[2] };
			double lastPoint[3] = { centrePoints->GetPoint(path->GetId(i - 1))[0], centrePoints->GetPoint(path->GetId(i - 1))[1], centrePoints->GetPoint(path->GetId(i - 1))[2] };

			vtkMath::Subtract(currentPoint, lastPoint, tempPoint);
			vtkMath::Normalize(tempPoint);
		}

		if (!data.getDoubleLayer()) {
			int material = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(currentPointId, 0);

			if ((data.getVtkData()->GetCellData()->GetArray("StatusEndo")->GetComponent(currentPointId, 0) == 0) &&
				((material == inMaterial1) || (material == inMaterial2) || (material == inMaterial3))) {
				data.getVtkData()->GetCellData()->GetArray("StatusEndo")->SetComponent(currentPointId, 0, 1);
				data.getVtkData()->GetCellData()->GetArray("basePathEndo")->SetComponent(currentPointId, 0, 1);
				data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->InsertTuple(currentPointId, tempPoint);
			}
		}

		int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(currentPointId, 0);
		if ((data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(currentPointId, 0) == 0) &&
			((material == inMaterial1) || (material == inMaterial2) || (material == inMaterial3))) {
			data.getVtkData()->GetCellData()->GetArray("Status")->SetComponent(currentPointId, 0, 1);
			data.getVtkData()->GetCellData()->GetArray("basePath")->SetComponent(currentPointId, 0, 1);
			data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->InsertTuple(currentPointId, tempPoint);
		}
		data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->InsertTuple(currentPointId, tempPoint);
	}
} // Methods::pathMarker

/*! Smooth the fiber along the path by averaging 5 fiber values before and 5 after the actually cell along the path.
 \param data Pointer to the orginal mesh.
 \param path Founded path between to points in the mesh.
 */

void Methods::smoothingPath(DataFormat &data, vtkSmartPointer<vtkIdList> path) {
	vtkSmartPointer<vtkDoubleArray> tempdiffVec = vtkDoubleArray::SafeDownCast(data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector"));

	for (vtkIdType i = 0; i < path->GetNumberOfIds(); i++) {
		double point[3] = { 0, 0, 0 };

		for (vtkIdType j = i - 5; j < i + 6; j++) {
			if ((j >= 0) && (j < path->GetNumberOfIds())) {
				point[0] += tempdiffVec->vtkDataArray::GetComponent(path->GetId(j), 0);
				point[1] += tempdiffVec->vtkDataArray::GetComponent(path->GetId(j), 1);
				point[2] += tempdiffVec->vtkDataArray::GetComponent(path->GetId(j), 2);
			}
		}

		vtkMath::Normalize(point);

		if (!data.getDoubleLayer()) {
			if (data.getVtkData()->GetCellData()->GetArray("StatusEndo")->GetComponent(path->GetId(i), 0) == 1) {
				data.getVtkData()->GetCellData()->GetArray("StatusEndo")->SetComponent(path->GetId(i), 0, 2);
				data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->InsertTuple(path->GetId(i), point);
			}
		}

		if (data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(path->GetId(i), 0) == 1) {
			data.getVtkData()->GetCellData()->GetArray("Status")->SetComponent(path->GetId(i), 0, 2);
			data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->InsertTuple(path->GetId(i), point);
		}
		data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->InsertTuple(path->GetId(i), point);
	}
} // Methods::smoothingPath

/*! Mark the fiber along the ring path in one tissue class in the mesh by subtracting the cell center points. A ring
 path is a closed loop.
 \param data Pointer to the orginal mesh.
 \param path Founded path between to points in the mesh.
 \param inMaterial Tissue class in that the fiber orientation would be marked.
 */
void Methods::pathMarkerRing(DataFormat &data, vtkSmartPointer<vtkIdList> path, Material::Mat inMaterial) {
	pathMarkerRing(data, path, inMaterial, inMaterial);
}

/*! Mark the fiber along the ring path in two tissue class in the mesh by subtracting the cell center points. A ring
 path is a closed loop.
 \param data Pointer to the orginal mesh.
 \param path Founded path between to points in the mesh.
 \param inMaterial1 First tissue class in that the fiber orientation would be marked.
 \param inMaterial2 Second tissue class in that the fiber orientation would be marked.
 */
void Methods::pathMarkerRing(DataFormat &data, vtkSmartPointer<vtkIdList> path, Material::Mat inMaterial1, Material::Mat inMaterial2) {
	vtkSmartPointer<vtkPoints> centrePoints = data.getCentrePoints()->GetPoints();

	vtkIdType currentPointId = -1;
	double tempPoint[3];

	for (long i = 0; i < path->GetNumberOfIds(); i++) {
		if (i < path->GetNumberOfIds() - 1) {
			currentPointId = path->GetId(i);
			double currentPoint[3] = { centrePoints->GetPoint(path->GetId(i))[0], centrePoints->GetPoint(path->GetId(i))[1], centrePoints->GetPoint(path->GetId(i))
			[2] };
			double nextPoint[3] = { centrePoints->GetPoint(path->GetId(i + 1))[0], centrePoints->GetPoint(path->GetId(i + 1))[1], centrePoints->GetPoint(path->GetId(i + 1))[2] };

			vtkMath::Subtract(nextPoint, currentPoint, tempPoint);
			vtkMath::Normalize(tempPoint);
		}
		else if (i == path->GetNumberOfIds() - 1) {
			currentPointId = path->GetId(i);
			double currentPoint[3] = { centrePoints->GetPoint(path->GetId(i))[0], centrePoints->GetPoint(path->GetId(i))[1], centrePoints->GetPoint(path->GetId(i))
			[2] };
			double nextPoint[3] = { centrePoints->GetPoint(path->GetId(0))[0], centrePoints->GetPoint(path->GetId(0))[1], centrePoints->GetPoint(path->GetId(0))
			[2] };

			vtkMath::Subtract(nextPoint, currentPoint, tempPoint);
			vtkMath::Normalize(tempPoint);
		}


		if (!data.getDoubleLayer()) {
			int material = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(currentPointId, 0);

			if ((data.getVtkData()->GetCellData()->GetArray("StatusEndo")->GetComponent(currentPointId, 0) == 0) &&
				((material == inMaterial1) || (material == inMaterial2))) {
				data.getVtkData()->GetCellData()->GetArray("StatusEndo")->SetComponent(currentPointId, 0, 1);
				data.getVtkData()->GetCellData()->GetArray("basePathEndo")->SetComponent(currentPointId, 0, 1);
				data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->InsertTuple(currentPointId, tempPoint);
			}
		}

		int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(currentPointId, 0);
		if ((data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(currentPointId, 0) == 0) &&
			((material == inMaterial1) || (material == inMaterial2))) {
			data.getVtkData()->GetCellData()->GetArray("Status")->SetComponent(currentPointId, 0, 1);
			data.getVtkData()->GetCellData()->GetArray("basePath")->SetComponent(currentPointId, 0, 1);
			data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->InsertTuple(currentPointId, tempPoint);
		}

		data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->InsertTuple(currentPointId, tempPoint);
	}
} // Methods::pathMarkerRing

/*! Smooth the fiber along the ring path by averaging 5 fiber values before and 5 after the actually cell along the
 path. A Ring path ist a closed loop.
 \param data Pointer to the orginal mesh.
 \param path Founded path between to points in the mesh.
 */
void Methods::smoothingRingPath(DataFormat &data, vtkSmartPointer<vtkIdList> path) {
	vtkSmartPointer<vtkDoubleArray> tempdiffVec = vtkDoubleArray::SafeDownCast(data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector"));

	for (vtkIdType i = 0; i < path->GetNumberOfIds(); i++) {
		double point[3] = { 0, 0, 0 };

		for (vtkIdType j = i - 5; j < i + 5; j++) {
			vtkIdType k = (j + path->GetNumberOfIds()) % path->GetNumberOfIds();
			point[0] += tempdiffVec->vtkDataArray::GetComponent(path->GetId(k), 0);
			point[1] += tempdiffVec->vtkDataArray::GetComponent(path->GetId(k), 1);
			point[2] += tempdiffVec->vtkDataArray::GetComponent(path->GetId(k), 2);
		}

		vtkMath::Normalize(point);

		if (!data.getDoubleLayer()) {
			if (data.getVtkData()->GetCellData()->GetArray("StatusEndo")->GetComponent(path->GetId(i), 0) == 1) {
				data.getVtkData()->GetCellData()->GetArray("StatusEndo")->SetComponent(path->GetId(i), 0, 2);
				data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->InsertTuple(path->GetId(i), point);
			}
		}

		if (data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(path->GetId(i), 0) == 1) {
			data.getVtkData()->GetCellData()->GetArray("Status")->SetComponent(path->GetId(i), 0, 2);
			data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->InsertTuple(path->GetId(i), point);
		}
		data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->InsertTuple(path->GetId(i), point);
	}
} // Methods::smoothingRingPath

/*! Clacuated the isolated boder layer for the atrial appendages.
 \param data Pointer to the orginal mesh.
 \param path Founded path at based of the appandage in the mesh.
 \param point Highest point oft the ovire of theatrial appandage.
 \retrun A vtkIdList with the cell-ids of the border layer.
 */
vtkSmartPointer<vtkIdList> Methods::boundaryLayerAppend(DataFormat &data, vtkSmartPointer<vtkIdList> path, std::vector<double> point) {
	vtkSmartPointer<vtkPoints> centrePoints = data.getCentrePoints()->GetPoints();

	vtkSmartPointer<vtkIdList> NeighborCellIdList = path;
	double topPoint[3] = { point.at(0), point.at(1), point.at(2) };

	vtkSmartPointer<vtkDoubleArray> tempArray = vtkSmartPointer<vtkDoubleArray>::New();

	tempArray->SetNumberOfComponents(3);
	tempArray->SetNumberOfTuples(centrePoints->GetNumberOfPoints());
	tempArray->FillComponent(0, 0);
	tempArray->FillComponent(1, 0);
	tempArray->FillComponent(2, 0);

	// init tempArray
	for (vtkIdType i = 0; i < path->GetNumberOfIds(); i++) {
		vtkIdType currentId = path->GetId(i);
		double *tempPoint = data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(currentId);
		tempArray->InsertTuple(currentId, tempPoint);
	}

	vector<double> centroid = getCentroid(data, path);
	double centrePointAppendages[3] = { centroid.at(0), centroid.at(1), centroid.at(2) };


	set<vtkIdType> viewedCells;

	for (vtkIdType i = 0; i < path->GetNumberOfIds(); i++) {
		viewedCells.insert(path->GetId(i));
	}

	set<vtkIdType> changedCells;

	for (vtkIdType j = 0; j < NeighborCellIdList->GetNumberOfIds(); j++) {
		vtkIdType parentID = NeighborCellIdList->GetId(j);
		vtkIdType preParentId = NeighborCellIdList->GetId((j + (NeighborCellIdList->GetNumberOfIds()) - 1) % (NeighborCellIdList->GetNumberOfIds()));
		vtkIdType postParentId = NeighborCellIdList->GetId((j + (NeighborCellIdList->GetNumberOfIds()) + 1) % (NeighborCellIdList->GetNumberOfIds()));

		double parentPoint[3] = { 0, 0, 0 };
		parentPoint[0] = centrePoints->GetPoint(parentID)[0];
		parentPoint[1] = centrePoints->GetPoint(parentID)[1];
		parentPoint[2] = centrePoints->GetPoint(parentID)[2];

		double preParentPoint[3] = { 0, 0, 0 };
		preParentPoint[0] = centrePoints->GetPoint(preParentId)[0];
		preParentPoint[1] = centrePoints->GetPoint(preParentId)[1];
		preParentPoint[2] = centrePoints->GetPoint(preParentId)[2];

		double postParentPoint[3] = { 0, 0, 0 };
		postParentPoint[0] = centrePoints->GetPoint(postParentId)[0];
		postParentPoint[1] = centrePoints->GetPoint(postParentId)[1];
		postParentPoint[2] = centrePoints->GetPoint(postParentId)[2];


		// define the Plane
		double vectorParentCentrePointApp[3];
		double vectorParentPreparentPoint[3];
		double normalVectorPlanePre[3];

		vtkMath::Subtract(preParentPoint, parentPoint, vectorParentPreparentPoint);
		vtkMath::Subtract(centrePointAppendages, parentPoint, vectorParentCentrePointApp);
		vtkMath::Cross(vectorParentPreparentPoint, vectorParentCentrePointApp, normalVectorPlanePre);
		vtkMath::Normalize(normalVectorPlanePre);
		vtkMath::Normalize(vectorParentPreparentPoint);

		vector<double> vecParentPreparentPoint(vectorParentPreparentPoint, vectorParentPreparentPoint + sizeof(vectorParentPreparentPoint) /
			sizeof(double));

		double vectorParentPostparentPoint[3];
		double normalVectorPlanePost[3];

		vtkMath::Subtract(postParentPoint, parentPoint, vectorParentPostparentPoint);
		vtkMath::Cross(vectorParentPostparentPoint, vectorParentCentrePointApp, normalVectorPlanePost);
		vtkMath::Normalize(normalVectorPlanePost);
		vtkMath::MultiplyScalar(normalVectorPlanePost, -1);
		vtkMath::Normalize(vectorParentPostparentPoint);
		vector<double> vecParentPostparentPoint(vectorParentPostparentPoint, vectorParentPostparentPoint + sizeof(vectorParentPostparentPoint) /
			sizeof(double));

		// left and right Plane
		double vectorPreparentTopPoint[3];
		double vectorPreparentCentrePointApp[3];
		double normalVectorLeftPlane[3];

		vtkMath::Subtract(centrePointAppendages, preParentPoint, vectorPreparentCentrePointApp);
		vtkMath::Subtract(preParentPoint, topPoint, vectorPreparentTopPoint);
		vtkMath::Cross(vectorPreparentCentrePointApp, vectorPreparentTopPoint, normalVectorLeftPlane);
		vtkMath::Normalize(normalVectorLeftPlane);

		double vectorPostparentTopPoint[3];
		double vectorPostparentCentrePointApp[3];
		double normalVectorRightPlane[3];

		vtkMath::Subtract(centrePointAppendages, postParentPoint, vectorPostparentCentrePointApp);
		vtkMath::Subtract(postParentPoint, topPoint, vectorPostparentTopPoint);
		vtkMath::Cross(vectorPostparentCentrePointApp, vectorPostparentTopPoint, normalVectorRightPlane);
		vtkMath::Normalize(normalVectorRightPlane);

		double angleBetweenPlane = Methods::angleBetween2Vectors(vecParentPreparentPoint, vecParentPostparentPoint);

		vtkSmartPointer<vtkIdList> cellPointIds = vtkSmartPointer<vtkIdList>::New();
		vtkSmartPointer<vtkIdList> currentNeighborCellIds = vtkSmartPointer<vtkIdList>::New();

		data.getVtkData()->GetCellPoints(parentID, cellPointIds);

		for (vtkIdType i = 0; i < cellPointIds->GetNumberOfIds(); i++) {
			vtkSmartPointer<vtkIdList> parentCellPointIdList = vtkSmartPointer<vtkIdList>::New();
			parentCellPointIdList->InsertNextId(cellPointIds->GetId(i));

			vtkSmartPointer<vtkIdList> tempNeighborCellIds = vtkSmartPointer<vtkIdList>::New();

			data.getVtkData()->GetCellNeighbors(parentID, parentCellPointIdList, tempNeighborCellIds);

			for (vtkIdType j = 0; j < tempNeighborCellIds->GetNumberOfIds(); j++) {
				currentNeighborCellIds->InsertUniqueId(tempNeighborCellIds->GetId(j));
			} // all neighbour cells
		}

		// calc the highest and the lowest cell point
		double distanceMax = 0;
		double distanceMin = std::numeric_limits<double>::max();
		vtkIdType highestPointOfCell = -1;
		vtkIdType lowestPointOfCell = -1;


		for (vtkIdType i = 0; i < cellPointIds->GetNumberOfIds(); i++) {
			double distance = vtkMath::Distance2BetweenPoints(topPoint, data.getVtkData()->GetPoint(cellPointIds->GetId(i)));
			if (distance < distanceMin) {
				distanceMin = distance;
				highestPointOfCell = cellPointIds->GetId(i);
			}
			if (distance > distanceMax) {
				distanceMax = distance;
				lowestPointOfCell = cellPointIds->GetId(i);
			}
		}


		double highestPointOfCellCoord[3] = { 0, 0, 0 };
		highestPointOfCellCoord[0] = data.getVtkData()->GetPoint(highestPointOfCell)[0];
		highestPointOfCellCoord[1] = data.getVtkData()->GetPoint(highestPointOfCell)[1];
		highestPointOfCellCoord[2] = data.getVtkData()->GetPoint(highestPointOfCell)[2];

		double lowestPointOfCellCoord[3] = { 0, 0, 0 };
		lowestPointOfCellCoord[0] = data.getVtkData()->GetPoint(lowestPointOfCell)[0];
		lowestPointOfCellCoord[1] = data.getVtkData()->GetPoint(lowestPointOfCell)[1];
		lowestPointOfCellCoord[2] = data.getVtkData()->GetPoint(lowestPointOfCell)[2];

		set<vtkIdType> tempViewedCells;
		do {
			vtkIdType size = currentNeighborCellIds->GetNumberOfIds();
			for (vtkIdType i = 0; i < size; i++) {
				vtkIdType currentId = currentNeighborCellIds->GetId(i);

				double currentPoint[3] = { 0, 0, 0 };
				currentPoint[0] = centrePoints->GetPoint(currentId)[0];
				currentPoint[1] = centrePoints->GetPoint(currentId)[1];
				currentPoint[2] = centrePoints->GetPoint(currentId)[2];

				double UpperPrePlane = vtkPlane::Evaluate(normalVectorPlanePre, highestPointOfCellCoord, currentPoint);
				double UpperPostPlane = vtkPlane::Evaluate(normalVectorPlanePost, highestPointOfCellCoord, currentPoint);
				double LowerPrePlane = vtkPlane::Evaluate(normalVectorPlanePre, lowestPointOfCellCoord, currentPoint);
				double LowerPostPlane = vtkPlane::Evaluate(normalVectorPlanePost, lowestPointOfCellCoord, currentPoint);
				double leftPlane = vtkPlane::Evaluate(normalVectorLeftPlane, preParentPoint, currentPoint);
				double rightPlane = vtkPlane::Evaluate(normalVectorRightPlane, postParentPoint, currentPoint);


				if (!((leftPlane < 0) && (rightPlane >= 0))) {
					currentNeighborCellIds->DeleteId(currentId);
				}

				if (angleBetweenPlane < 180) {
					if (!((LowerPrePlane < 0) && (LowerPostPlane < 0))) {
						currentNeighborCellIds->DeleteId(currentId);
					}

					if (!((UpperPrePlane >= 0) || (UpperPostPlane >= 0))) {
						currentNeighborCellIds->DeleteId(currentId);
					}
				}
				else {
					if (!((LowerPrePlane < 0) || (LowerPostPlane < 0))) {
						currentNeighborCellIds->DeleteId(currentId);
					}

					if (!((UpperPrePlane >= 0) && (UpperPostPlane >= 0))) {
						currentNeighborCellIds->DeleteId(currentId);
					}
				}
				if ((viewedCells.find(currentId) != viewedCells.end()) ||
					(tempViewedCells.find(currentId) != tempViewedCells.end()) ||
					(data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(currentId, 0) != 0)) {
					currentNeighborCellIds->DeleteId(currentId);
				}
				tempViewedCells.insert(currentId);
			}


			for (vtkIdType i = 0; i < currentNeighborCellIds->GetNumberOfIds(); i++) {
				vtkIdType currentId = currentNeighborCellIds->GetId(i);

				double tempPoint[3] = { 0, 0, 0 };
				tempPoint[0] = tempArray->GetComponent(parentID, 0) + tempArray->GetComponent(currentId, 0);
				tempPoint[1] = tempArray->GetComponent(parentID, 1) + tempArray->GetComponent(currentId, 1);
				tempPoint[2] = tempArray->GetComponent(parentID, 2) + tempArray->GetComponent(currentId, 2);
				tempArray->SetTuple(currentId, tempPoint);
				changedCells.insert(currentId);
			}

			vtkSmartPointer<vtkIdList> newNeighborCellIds = vtkSmartPointer<vtkIdList>::New();

			for (vtkIdType i = 0; i < currentNeighborCellIds->GetNumberOfIds(); i++) {
				vtkIdType currentCellIds = currentNeighborCellIds->GetId(i);

				vtkSmartPointer<vtkIdList> cellPointIds = vtkSmartPointer<vtkIdList>::New();
				data.getVtkData()->GetCellPoints(currentCellIds, cellPointIds);

				for (vtkIdType i = 0; i < cellPointIds->GetNumberOfIds(); i++) {
					vtkSmartPointer<vtkIdList> parentCellPointIdList = vtkSmartPointer<vtkIdList>::New();
					parentCellPointIdList->InsertNextId(cellPointIds->GetId(i));

					vtkSmartPointer<vtkIdList> tempNeighborCellIds = vtkSmartPointer<vtkIdList>::New();
					data.getVtkData()->GetCellNeighbors(currentCellIds, parentCellPointIdList, tempNeighborCellIds);

					for (vtkIdType j = 0; j < tempNeighborCellIds->GetNumberOfIds(); j++) {
						if ((tempViewedCells.find(tempNeighborCellIds->GetId(j)) == tempViewedCells.end()) &&
							(viewedCells.find(tempNeighborCellIds->GetId(j)) == viewedCells.end()))
							newNeighborCellIds->InsertUniqueId(tempNeighborCellIds->GetId(j));
					} // all new neighbour cells
				}
			}
			currentNeighborCellIds->Reset();
			currentNeighborCellIds = newNeighborCellIds;
		} while (currentNeighborCellIds->GetNumberOfIds() > 0);
	}

	vtkSmartPointer<vtkIdList> returnIds = vtkSmartPointer<vtkIdList>::New();

	for (set<vtkIdType>::iterator iter = changedCells.begin(); iter != changedCells.end(); iter++) {
		returnIds->InsertNextId(*iter);
		if (!data.getDoubleLayer()) {
			if (data.getVtkData()->GetCellData()->GetArray("StatusEndo")->GetComponent(*iter, 0) == 0) {
				double *tempPoint = tempArray->GetTuple(*iter);
				vtkMath::Normalize(tempPoint);

				if (data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetComponent(*iter, 0) +
					data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetComponent(*iter, 1) +
					data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetComponent(*iter, 2) == 0) {
					data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->SetTuple(*iter, tempPoint);
				}
				data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->SetTuple(*iter, tempPoint);
				data.getVtkData()->GetCellData()->GetArray("basePathEndo")->SetComponent(*iter, 0, 1);
				data.getVtkData()->GetCellData()->GetArray("StatusEndo")->SetComponent(*iter, 0, 1);
			}
		}
		if (data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(*iter, 0) == 0) {
			double *tempPoint = tempArray->GetTuple(*iter);
			vtkMath::Normalize(tempPoint);

			if (data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetComponent(*iter, 0) +
				data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetComponent(*iter, 1) +
				data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetComponent(*iter, 2) == 0) {
				data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->SetTuple(*iter, tempPoint);
			}
			data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->SetTuple(*iter, tempPoint);
			data.getVtkData()->GetCellData()->GetArray("basePath")->SetComponent(*iter, 0, 1);
			data.getVtkData()->GetCellData()->GetArray("Status")->SetComponent(*iter, 0, 1);
		}
	}

	vtkSmartPointer<vtkIdList> borderPlane = Methods::add2Paths(returnIds, path);
	return borderPlane;
} // Methods::boundaryLayerAppend

/*! Grow/ interpolate a path to a point in the mesh and in one tissue class.
 \param data Pointer to the orginal mesh.
 \param border Calcuated path or border in the mesh.
 \param point The Point to that the path would be interpolated as Vector (x,y,z)
 \param inMaterial Tissue class in that the fiber a orientation could be interpolated.
 \retrun A vtkIdList with the cell-ids in that the fiber was interpolated.
 */
vtkSmartPointer<vtkIdList> Methods::growPathtoPoint(DataFormat &data, vtkSmartPointer<vtkIdList> border, std::vector<double> point, Material::Mat inMaterial) {
	return growPathtoPoint(data, border, point, inMaterial, inMaterial);
}

/*! Grow/ interpolate a path to a point in the mesh and in two tissue class.
 \param data Pointer to the orginal mesh.
 \param border Calcuated path or border in the mesh.
 \param point The Point to that the path would be interpolated as Vector (x,y,z)
 \param inMaterial1 First tissue class in that the fiber a orientation could be interpolated.
 \param inMaterial2 Second tissue class in that the fiber a orientation could be interpolated.
 \retrun A vtkIdList with the cell-ids in that the fiber was interpolated.
 */
vtkSmartPointer<vtkIdList> Methods::growPathtoPoint(DataFormat &data, vtkSmartPointer<vtkIdList> border, std::vector<double> point, Material::Mat inMaterial1, Material::Mat inMaterial2) {
	vtkSmartPointer<vtkPoints> centrePoints = data.getCentrePoints()->GetPoints();
	vtkSmartPointer<vtkUnstructuredGrid> centrePointsUGrid = data.getCentrePoints();

	double topPoint[3] = { point.at(0), point.at(1), point.at(2) };
	vtkIdType topPointId = centrePointsUGrid->FindPoint(point.at(0), point.at(1), point.at(2));

	vtkSmartPointer<vtkDoubleArray> tempArray = vtkSmartPointer<vtkDoubleArray>::New();

	tempArray->SetNumberOfComponents(3);
	tempArray->SetNumberOfTuples(centrePoints->GetNumberOfPoints());
	tempArray->FillComponent(0, 0);
	tempArray->FillComponent(1, 0);
	tempArray->FillComponent(2, 0);

	double maxDistancetoPoint = getMaxDistance(data, border, point);

	for (vtkIdType i = 0; i < border->GetNumberOfIds(); i++) {
		vtkIdType currentId = border->GetId(i);
		double tempPoint[3];
		tempPoint[0] = data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetComponent(currentId, 0);
		tempPoint[1] = data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetComponent(currentId, 1);
		tempPoint[2] = data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetComponent(currentId, 2);
		tempArray->InsertTuple(currentId, tempPoint);
	}

	vtkSmartPointer<vtkDoubleArray> cells = vtkSmartPointer<vtkDoubleArray>::New();
	cells->SetNumberOfComponents(2);
	cells->SetComponentName(0, "CellID");
	cells->SetComponentName(1, "Distance to the top point");

	set<vtkIdType> viewedCellList;

	bool end = false;
	set<vtkIdType> cellIdList;
	cellIdList.insert(topPointId);


	while (!end) {
		set<vtkIdType> newcellIdList;
		for (set<vtkIdType>::iterator iter = cellIdList.begin(); iter != cellIdList.end(); iter++) {
			vtkIdType currentCellIds = *iter;
			vtkSmartPointer<vtkIdList> cellPointIds = vtkSmartPointer<vtkIdList>::New();

			data.getVtkData()->GetCellPoints(currentCellIds, cellPointIds);

			for (vtkIdType i = 0; i < cellPointIds->GetNumberOfIds(); i++) {
				vtkSmartPointer<vtkIdList> CellPointIdList = vtkSmartPointer<vtkIdList>::New();
				CellPointIdList->InsertNextId(cellPointIds->GetId(i));

				vtkSmartPointer<vtkIdList> tempNeighborCellIds = vtkSmartPointer<vtkIdList>::New();
				data.getVtkData()->GetCellNeighbors(currentCellIds, CellPointIdList, tempNeighborCellIds);

				for (vtkIdType j = 0; j < tempNeighborCellIds->GetNumberOfIds(); j++) {
					if ((viewedCellList.find(tempNeighborCellIds->GetId(j)) == viewedCellList.end()) &&
						((data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(tempNeighborCellIds->GetId(j), 0) == 0) ||
						(!data.getDoubleLayer() &&
							(data.getVtkData()->GetCellData()->GetArray("StatusEndo")->GetComponent(tempNeighborCellIds->GetId(j), 0) == 0))) &&
							(sqrt(abs(vtkMath::Distance2BetweenPoints(topPoint, data.getCentrePoints()->GetPoint(tempNeighborCellIds->GetId(j)))))
								<= maxDistancetoPoint)
						) {
						double cellsArray[2] = { 0, 0 };
						cellsArray[0] = tempNeighborCellIds->GetId(j);
						cellsArray[1] = sqrt(abs(vtkMath::Distance2BetweenPoints(topPoint, centrePoints->GetPoint(tempNeighborCellIds->GetId(j)))));
						cells->InsertNextTuple(cellsArray);
						viewedCellList.insert(tempNeighborCellIds->GetId(j));
						newcellIdList.insert(tempNeighborCellIds->GetId(j));
					} // all new neighbour cells
				}
			}
		}
		if (newcellIdList.size() == 0) {
			end = true;
		}
		cellIdList.clear();
		cellIdList = newcellIdList;
	}

	vtkSortDataArray::SortArrayByComponent(cells, 1);

	vtkSmartPointer<vtkIdList> returnIds = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkIdList> diffVectorNull = vtkSmartPointer<vtkIdList>::New();

	for (vtkIdType j = cells->GetNumberOfTuples() - 1; j >= 0; j--) {
		vtkIdType currentCellIds = cells->vtkDataArray::GetComponent(j, 0);
		vtkSmartPointer<vtkIdList> cellPointIds = vtkSmartPointer<vtkIdList>::New();
		vtkSmartPointer<vtkIdList> neighbours = vtkSmartPointer<vtkIdList>::New();

		data.getVtkData()->GetCellPoints(currentCellIds, cellPointIds);

		for (vtkIdType i = 0; i < cellPointIds->GetNumberOfIds(); i++) {
			vtkSmartPointer<vtkIdList> CellPointIdList = vtkSmartPointer<vtkIdList>::New();
			CellPointIdList->InsertNextId(cellPointIds->GetId(i));

			vtkSmartPointer<vtkIdList> tempNeighborCellIds = vtkSmartPointer<vtkIdList>::New();
			data.getVtkData()->GetCellNeighbors(currentCellIds, CellPointIdList, tempNeighborCellIds);

			for (vtkIdType j = 0; j < tempNeighborCellIds->GetNumberOfIds(); j++) {
				neighbours->InsertUniqueId(tempNeighborCellIds->GetId(j));
			} // all new neighbour cells
		}

		double tempPoint[3] = { 0, 0, 0 };
		tempPoint[0] = tempArray->GetComponent(currentCellIds, 0);
		tempPoint[1] = tempArray->GetComponent(currentCellIds, 1);
		tempPoint[2] = tempArray->GetComponent(currentCellIds, 2);

		for (vtkIdType i = 0; i < neighbours->GetNumberOfIds(); i++) {
			vtkIdType neighbourID = neighbours->GetId(i);

			if (tempArray->GetComponent(neighbourID, 0) +
				tempArray->GetComponent(neighbourID, 1) + tempArray->GetComponent(neighbourID, 2) != 0) {
				tempPoint[0] += tempArray->GetComponent(neighbourID, 0);
				tempPoint[1] += tempArray->GetComponent(neighbourID, 1);
				tempPoint[2] += tempArray->GetComponent(neighbourID, 2);
			}
		}
		vtkMath::Normalize(tempPoint);
		tempArray->InsertTuple(currentCellIds, tempPoint);

		if (!data.getDoubleLayer()) {
			int material = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(currentCellIds, 0);

			if ((data.getVtkData()->GetCellData()->GetArray("StatusEndo")->GetComponent(currentCellIds, 0) == 0) &&
				((material == inMaterial1) || (material == inMaterial2))) {
				if ((tempArray->GetTuple(currentCellIds)[0] != 0) ||
					(tempArray->GetTuple(currentCellIds)[1] != 0) ||
					(tempArray->GetTuple(currentCellIds)[2] != 0)) {
					double *tempPoint = tempArray->GetTuple(currentCellIds);
					data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->InsertTuple(currentCellIds, tempPoint);
					data.getVtkData()->GetCellData()->GetArray("StatusEndo")->SetComponent(currentCellIds, 0, 3);
					returnIds->InsertNextId(currentCellIds);
				}
				else {
					diffVectorNull->InsertNextId(currentCellIds);
				}
			}
		}

		int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(currentCellIds, 0);

		if ((data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(currentCellIds, 0) == 0) &&
			((material == inMaterial1) || (material == inMaterial2))) {
			if ((tempArray->GetTuple(currentCellIds)[0] != 0) ||
				(tempArray->GetTuple(currentCellIds)[1] != 0) ||
				(tempArray->GetTuple(currentCellIds)[2] != 0)) {
				double *tempPoint = tempArray->GetTuple(currentCellIds);
				data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->InsertTuple(currentCellIds, tempPoint);
				data.getVtkData()->GetCellData()->GetArray("Status")->SetComponent(currentCellIds, 0, 3);
				returnIds->InsertNextId(currentCellIds);
			}
			else {
				diffVectorNull->InsertNextId(currentCellIds);
			}
		}
	}


	for (vtkIdType i = 0; i < diffVectorNull->GetNumberOfIds(); i++) {
		vtkIdType currentCellId = diffVectorNull->GetId(i);
		double tempPoint[3] = { 0, 0, 0 };
		tempPoint[0] = tempArray->GetComponent(currentCellId, 0);
		tempPoint[1] = tempArray->GetComponent(currentCellId, 1);
		tempPoint[2] = tempArray->GetComponent(currentCellId, 2);

		vtkSmartPointer<vtkIdList> neighbourCells = getCellNeighbours(data, currentCellId);

		for (vtkIdType i = 0; i < neighbourCells->GetNumberOfIds(); i++) {
			vtkIdType neighbourID = neighbourCells->GetId(i);

			if (tempArray->GetComponent(neighbourID, 0) +
				tempArray->GetComponent(neighbourID, 1) + tempArray->GetComponent(neighbourID, 2) != 0) {
				tempPoint[0] += tempArray->GetComponent(neighbourID, 0);
				tempPoint[1] += tempArray->GetComponent(neighbourID, 1);
				tempPoint[2] += tempArray->GetComponent(neighbourID, 2);
			}
		}

		vtkMath::Normalize(tempPoint);
		tempArray->InsertTuple(currentCellId, tempPoint);
		if (!data.getDoubleLayer()) {
			if ((tempArray->GetTuple(currentCellId)[0] != 0) && (tempArray->GetTuple(currentCellId)[1] != 0) &&
				(tempArray->GetTuple(currentCellId)[2] != 0)) {
				data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->InsertTuple(currentCellId, tempPoint);
				data.getVtkData()->GetCellData()->GetArray("StatusEndo")->SetComponent(currentCellId, 0, 3);
				returnIds->InsertNextId(currentCellId);
			}
		}

		if ((tempArray->GetTuple(currentCellId)[0] != 0) && (tempArray->GetTuple(currentCellId)[1] != 0) &&
			(tempArray->GetTuple(currentCellId)[2] != 0)) {
			data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->InsertTuple(currentCellId, tempPoint);
			data.getVtkData()->GetCellData()->GetArray("Status")->SetComponent(currentCellId, 0, 3);
			returnIds->InsertNextId(currentCellId);
		}
	}


	vtkSmartPointer<vtkIdList> returncells = Methods::add2Paths(returnIds, border);
	return returncells;
} // Methods::growPathtoPoint

/*! Grow the crista terminal path in the mesh and in three tissue class.
 \param data Pointer to the orginal mesh.
 \param path Founded path in the mesh.
 \param inMaterial1 First tissue class in that the fiber a orientation could be interpolated.
 \param inMaterial2 Second tissue class in that the fiber orientation could be interpolated.
 \param inMaterial3 Third tissue class in that the fiber orientation could be interpolated.
 \retrun A vtkIdList with the cell-ids of the grown crista terminal path.
 */
vtkSmartPointer<vtkIdList> Methods::growCTPath(DataFormat &data, vtkSmartPointer<vtkIdList> path, double minRadius, Material::Mat inMaterial1, Material::Mat inMaterial2, Material::Mat inMaterial3) {
	double radius = minRadius + 1.32;

	vtkSmartPointer<vtkPoints> centrePoints = data.getCentrePoints()->GetPoints();
	vtkSmartPointer<vtkUnstructuredGrid> centrePointsUGrid = data.getCentrePoints();
	vtkSmartPointer<vtkPointLocator> centrePointLocator = vtkSmartPointer<vtkPointLocator>::New();

	centrePointLocator->SetDataSet(centrePointsUGrid);
	centrePointLocator->BuildLocator();

	vtkSmartPointer<vtkDoubleArray> tempDiffVec = vtkDoubleArray::SafeDownCast(data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector"));

	vtkSmartPointer<vtkDoubleArray> tempArray = vtkSmartPointer<vtkDoubleArray>::New();
	tempArray->SetNumberOfComponents(3);
	tempArray->SetNumberOfTuples(centrePoints->GetNumberOfPoints());
	tempArray->FillComponent(0, 0);
	tempArray->FillComponent(1, 0);
	tempArray->FillComponent(2, 0);

	set<vtkIdType> changedId;
	set<vtkIdType> pathSearch;

	for (vtkIdType i = 0; i < path->GetNumberOfIds(); i++) {
		pathSearch.insert(path->GetId(i));
	}

	double pathLength = Methods::pathLength(data, path);
	double distance = 0;

	vtkIdType parentId = path->GetId(0);
	vtkSmartPointer<vtkIdList> neighbors = vtkSmartPointer<vtkIdList>::New();

	centrePointLocator->FindPointsWithinRadius((radius), centrePoints->GetPoint(parentId), neighbors);

	for (vtkIdType j = 0; j < neighbors->GetNumberOfIds(); j++) {
		vtkIdType currentId = neighbors->GetId(j);
		if (pathSearch.find(currentId) == pathSearch.end()) {
			tempArray->InsertComponent(currentId, 0, (tempDiffVec->vtkDataArray::GetComponent(parentId, 0) + tempArray->GetComponent(currentId, 0)));

			tempArray->InsertComponent(currentId, 1, (tempDiffVec->vtkDataArray::GetComponent(parentId, 1) + tempArray->GetComponent(currentId, 1)));

			tempArray->InsertComponent(currentId, 2, (tempDiffVec->vtkDataArray::GetComponent(parentId, 2) + tempArray->GetComponent(currentId, 2)));
			changedId.insert(currentId);
		}
	} // for start point


	for (vtkIdType i = 1; i < path->GetNumberOfIds(); i++) {
		vtkIdType lastId = path->GetId(i - 1);
		vtkIdType parentId = path->GetId(i);

		double x[3] = { 0, 0, 0 };
		x[0] = data.getCentrePoints()->GetPoint(lastId)[0];
		x[1] = data.getCentrePoints()->GetPoint(lastId)[1];
		x[2] = data.getCentrePoints()->GetPoint(lastId)[2];

		double y[3] = { 0, 0, 0 };
		y[0] = data.getCentrePoints()->GetPoint(parentId)[0];
		y[1] = data.getCentrePoints()->GetPoint(parentId)[1];
		y[2] = data.getCentrePoints()->GetPoint(parentId)[2];

		distance += sqrt(abs(vtkMath::Distance2BetweenPoints(x, y)));

		double width = radius - 1.32*(distance / pathLength);

		vtkSmartPointer<vtkIdList> neighbors = vtkSmartPointer<vtkIdList>::New();
		centrePointLocator->FindPointsWithinRadius(width, centrePoints->GetPoint(parentId), neighbors);

		for (vtkIdType j = 0; j < neighbors->GetNumberOfIds(); j++) {
			vtkIdType currentId = neighbors->GetId(j);
			if (pathSearch.find(currentId) == pathSearch.end()) {
				tempArray->InsertComponent(currentId, 0, (tempDiffVec->vtkDataArray::GetComponent(parentId, 0) + tempArray->GetComponent(currentId, 0)));

				tempArray->InsertComponent(currentId, 1, (tempDiffVec->vtkDataArray::GetComponent(parentId, 1) + tempArray->GetComponent(currentId, 1)));

				tempArray->InsertComponent(currentId, 2, (tempDiffVec->vtkDataArray::GetComponent(parentId, 2) + tempArray->GetComponent(currentId, 2)));
				changedId.insert(currentId);
			}
		}
	}

	vtkSmartPointer<vtkIdList> returnIds = vtkSmartPointer<vtkIdList>::New();

	for (set<vtkIdType>::iterator iter = changedId.begin(); iter != changedId.end(); iter++) {
		if (!data.getDoubleLayer()) {
			int material1 = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(*iter, 0);
			if ((data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(*iter, 0) == 0) &&
				((material1 == inMaterial1) || (material1 == inMaterial2) || (material1 == inMaterial3))) {
				double *tempPoint = tempArray->GetTuple(*iter);
				vtkMath::Normalize(tempPoint);
				data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->InsertTuple(*iter, tempPoint);
				data.getVtkData()->GetCellData()->GetArray("Status")->SetComponent(*iter, 0, 3);
				returnIds->InsertNextId(*iter);
			}
			if ((data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(*iter, 0) != 0) &&
				((material1 == inMaterial1) || (material1 == inMaterial2) || (material1 == inMaterial3))) {
				double currentVector[3];
				currentVector[0] = tempArray->GetTuple(*iter)[0];
				currentVector[1] = tempArray->GetTuple(*iter)[1];
				currentVector[2] = tempArray->GetTuple(*iter)[2];
				vtkMath::Normalize(currentVector);

				double oldDifferenceVector[3];
				oldDifferenceVector[0] = data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetTuple(*iter)[0];
				oldDifferenceVector[1] = data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetTuple(*iter)[1];
				oldDifferenceVector[2] = data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetTuple(*iter)[2];

				double newVector[3];
				vtkMath::Add(currentVector, oldDifferenceVector, newVector);
				vtkMath::Normalize(newVector);

				data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->InsertTuple(*iter, newVector);
				data.getVtkData()->GetCellData()->GetArray("Status")->SetComponent(*iter, 0, 3);
				returnIds->InsertNextId(*iter);
			}

			int material2 = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(*iter, 0);
			if ((data.getVtkData()->GetCellData()->GetArray("StatusEndo")->GetComponent(*iter, 0) == 0) &&
				((material2 == inMaterial1) || (material2 == inMaterial2) || (material2 == inMaterial3))) {
				double *tempPoint = tempArray->GetTuple(*iter);
				vtkMath::Normalize(tempPoint);
				data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->InsertTuple(*iter, tempPoint);
				data.getVtkData()->GetCellData()->GetArray("StatusEndo")->SetComponent(*iter, 0, 3);
				returnIds->InsertNextId(*iter);
			}
			if ((data.getVtkData()->GetCellData()->GetArray("StatusEndo")->GetComponent(*iter, 0) != 0) &&
				((material1 == inMaterial1) || (material1 == inMaterial2) || (material1 == inMaterial3))) {
				double currentVector[3];
				currentVector[0] = tempArray->GetTuple(*iter)[0];
				currentVector[1] = tempArray->GetTuple(*iter)[1];
				currentVector[2] = tempArray->GetTuple(*iter)[2];
				vtkMath::Normalize(currentVector);

				double oldDifferenceVector[3];
				oldDifferenceVector[0] = data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->GetTuple(*iter)[0];
				oldDifferenceVector[1] = data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->GetTuple(*iter)[1];
				oldDifferenceVector[2] = data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->GetTuple(*iter)[2];

				double newVector[3];
				vtkMath::Add(currentVector, oldDifferenceVector, newVector);
				vtkMath::Normalize(newVector);

				data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->InsertTuple(*iter, newVector);
				data.getVtkData()->GetCellData()->GetArray("StatusEndo")->SetComponent(*iter, 0, 3);
				returnIds->InsertNextId(*iter);
			}
		}
		else {
			int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(*iter, 0);
			if ((data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(*iter, 0) == 0) &&
				((material == inMaterial1) || (material == inMaterial2) || (material == inMaterial3))) {
				double *tempPoint = tempArray->GetTuple(*iter);
				vtkMath::Normalize(tempPoint);
				data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->InsertTuple(*iter, tempPoint);
				data.getVtkData()->GetCellData()->GetArray("Status")->SetComponent(*iter, 0, 3);
				returnIds->InsertNextId(*iter);
			}
			if ((data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(*iter, 0) != 0) &&
				((material == inMaterial1) || (material == inMaterial2) || (material == inMaterial3))) {
				double currentVector[3];
				currentVector[0] = tempArray->GetTuple(*iter)[0];
				currentVector[1] = tempArray->GetTuple(*iter)[1];
				currentVector[2] = tempArray->GetTuple(*iter)[2];
				vtkMath::Normalize(currentVector);

				double oldDifferenceVector[3];
				oldDifferenceVector[0] = data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetTuple(*iter)[0];
				oldDifferenceVector[1] = data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetTuple(*iter)[1];
				oldDifferenceVector[2] = data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetTuple(*iter)[2];

				double newVector[3];
				vtkMath::Add(currentVector, oldDifferenceVector, newVector);
				vtkMath::Normalize(newVector);

				data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->InsertTuple(*iter, newVector);
				data.getVtkData()->GetCellData()->GetArray("Status")->SetComponent(*iter, 0, 3);
				returnIds->InsertNextId(*iter);
			}
		}
	}

	vtkSmartPointer<vtkIdList> returnList = Methods::add2Paths(returnIds, path);
	return returnList;
} // Methods::growCTPath

/*! Grow transmural the crista terminal path in the mesh and in three tissue class.
 \param data Pointer to the orginal mesh.
 \param path Founded path in the mesh.
 \param inMaterial1 First tissue class in that the fiber a orientation could be interpolated.
 \param inMaterial2 Second tissue class in that the fiber orientation could be interpolated.
 \param inMaterial3 Third tissue class in that the fiber orientation could be interpolated.
 \param atrium In which atrium the path would be grown.
 \retrun A vtkIdList with the cell-ids of the grown crista terminal path.
 */
vtkSmartPointer<vtkIdList> Methods::growCTPathTransmural(DataFormat &data, vtkSmartPointer<vtkIdList> path, Material::Mat inMaterial1, Material::Mat inMaterial2, Material::Mat inMaterial3, string atrium) {
	vtkSmartPointer<vtkDoubleArray> tempDiffVec = vtkDoubleArray::SafeDownCast(data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector"));

	DataFormat surface;

	surface.setVtkData(data.getSurfacewithNormals());
	surface.setCentrePoints(data.getCentrePointsSurface());

	vtkSmartPointer<vtkDoubleArray> tempArray = vtkSmartPointer<vtkDoubleArray>::New();
	tempArray->SetNumberOfComponents(3);
	tempArray->SetNumberOfTuples(data.getCentrePoints()->GetNumberOfPoints());
	tempArray->FillComponent(0, 0);
	tempArray->FillComponent(1, 0);
	tempArray->FillComponent(2, 0);

	set<vtkIdType> changedId;
	vtkSmartPointer<vtkIdList> returnIds = vtkSmartPointer<vtkIdList>::New();

	double pathLength = Methods::pathLength(data, path);
	double distance = 0;

	for (vtkIdType i = 0; i < path->GetNumberOfIds(); i++) {
		vtkIdType parentId = path->GetId(i);

		if (i != 0) {
			vtkIdType lastId = path->GetId(i - 1);
			vtkIdType parentId = path->GetId(i);

			double x[3] = { 0, 0, 0 };
			x[0] = data.getCentrePoints()->GetPoint(lastId)[0];
			x[1] = data.getCentrePoints()->GetPoint(lastId)[1];
			x[2] = data.getCentrePoints()->GetPoint(lastId)[2];

			double y[3] = { 0, 0, 0 };
			y[0] = data.getCentrePoints()->GetPoint(parentId)[0];
			y[1] = data.getCentrePoints()->GetPoint(parentId)[1];
			y[2] = data.getCentrePoints()->GetPoint(parentId)[2];

			distance += sqrt(abs(vtkMath::Distance2BetweenPoints(x, y)));
		}

		vector<double> pathPoint = { data.getCentrePoints()->GetPoint(parentId)[0], data.getCentrePoints()->GetPoint(parentId)[1], data.getCentrePoints()->GetPoint(parentId)[2] };
		vector<double> pillPoint1;
		vector<double> pillPoint2;


		if (atrium.compare("right") == 0) {
			pillPoint1 = { Methods::findClosedPointinMaterialInRightEpi(surface, pathPoint) };
			pillPoint2 = { Methods::findClosedPointinMaterialInRightEndo(surface, pathPoint) };
		}
		else if (atrium.compare("left") == 0) {
			pillPoint1 = { Methods::findClosedPointinMaterialInLeftEpi(surface, pathPoint) };
			pillPoint2 = { Methods::findClosedPointinMaterialInLeftEndo(surface, pathPoint) };
		}


		double width = 3.96 - 1.32*(distance / pathLength);

		DataFormat pillData;
		pillData = getPill(pillPoint1, pillPoint2, 0.2, width);
		vtkSmartPointer<vtkPolyData> pill = vtkPolyData::SafeDownCast(pillData.getVtkData());

		vtkSmartPointer<vtkSelectEnclosedPoints> selectEnclosedPoints = vtkSmartPointer<vtkSelectEnclosedPoints>::New();
		selectEnclosedPoints->Initialize(pill);
		selectEnclosedPoints->SetTolerance(0);

		vtkSmartPointer<vtkIdList> cells = getCellsinRadius(data, pillPoint1, ConfigFiberorientation::maxWitdthAtrialWall);

		for (vtkIdType j = 0; j < cells->GetNumberOfIds(); j++) {
			vtkIdType currentId = cells->GetId(j);

			double *point = data.getCentrePoints()->GetPoint(currentId);
			int inside = selectEnclosedPoints->IsInsideSurface(point);

			if (inside == 1) {
				int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(currentId, 0);
				if ((material == inMaterial1) || (material == inMaterial2) || (material == inMaterial3)) {
					if ((tempDiffVec->GetComponent(parentId, 0) != 0) ||
						(tempDiffVec->GetComponent(parentId, 1) != 0) || (tempDiffVec->GetComponent(parentId, 2) != 0)) {
						tempArray->InsertComponent(currentId, 0, (tempDiffVec->GetComponent(parentId, 0) + tempArray->GetComponent(currentId, 0)));
						tempArray->InsertComponent(currentId, 1, (tempDiffVec->GetComponent(parentId, 1) + tempArray->GetComponent(currentId, 1)));
						tempArray->InsertComponent(currentId, 2, (tempDiffVec->GetComponent(parentId, 2) + tempArray->GetComponent(currentId, 2)));
					}
					changedId.insert(currentId);
				}
			}
		}
	}

	for (set<vtkIdType>::iterator iter = changedId.begin(); iter != changedId.end(); iter++) {
		int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(*iter, 0);
		if ((material == inMaterial1) || (material == inMaterial2) || (material == inMaterial3)) {
			double *tempPoint = tempArray->GetTuple(*iter);
			vtkMath::Normalize(tempPoint);
			data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->InsertTuple(*iter, tempPoint);
			data.getVtkData()->GetCellData()->GetArray("Status")->SetComponent(*iter, 0, 3);
			returnIds->InsertUniqueId(*iter);
		}
	}

	vtkSmartPointer<vtkIdList> returnList = Methods::add2Paths(returnIds, path);
	return returnList;
} // Methods::growCTPathTransmural

/*! Interpolate the region between two pectinate paths in the mesh and in one tissue class. The seed point for the
 interpolation region would be automatic calcuated.
 \param data Pointer to the orginal mesh.
 \param path1 Founded first pectinat path in the mesh.
 \param path2 Founded second pectinat path in the mesh.
 \param planePoint Help point for define the plane.
 \param inMaterial Tissue class in that the fiber a orientation could be interpolated.
 \retrun A vtkIdList with the cell-ids of the interpolate region.
 */
vtkSmartPointer<vtkIdList> Methods::growBetweenPecti(DataFormat &data, vtkSmartPointer<vtkIdList> path1, vtkSmartPointer<vtkIdList> path2, vector<double> planePoint, Material::Mat inMaterial) {
	vector<double> seedPoint = getSeedPoint(data, path1, path2, inMaterial, Material::Pektinatmuskeln, Material::Tricuspid_Valve_Ring, Material::Crista_Terminalis_2, 50, "right");

	return growBetweenPecti(data, path1, path2, planePoint, inMaterial, seedPoint);
}

/*! Interpolate the region between two pectinate paths in the mesh and in one tissue class.
 \param data Pointer to the orginal mesh.
 \param path1 Founded first pectinat path in the mesh.
 \param path2 Founded second pectinat path in the mesh.
 \param planePoint Help point for define the plane.
 \param inMaterial Tissue class in that the fiber a orientation could be interpolated.
 \param seedPoint Seed point for the interpolation region.
 \retrun A vtkIdList with the cell-ids of the interpolate region.
 */
vtkSmartPointer<vtkIdList> Methods::growBetweenPecti(DataFormat &data, vtkSmartPointer<vtkIdList> path1, vtkSmartPointer<vtkIdList> path2, vector<double> planePoint, Material::Mat inMaterial, vector<double> seedPoint) {
	vtkSmartPointer<vtkPoints> centrePoints = data.getCentrePoints()->GetPoints();
	vtkSmartPointer<vtkUnstructuredGrid> centrePointsUGrid = data.getCentrePoints();


	double startPointPlane1Array[3];

	centrePoints->GetPoint(path1->GetId(0), startPointPlane1Array);
	double targetPointPlane1Array[3];
	centrePoints->GetPoint(path1->GetId(path1->GetNumberOfIds() - 1), targetPointPlane1Array);
	double planePointArray[3] = { planePoint.at(0), planePoint.at(1), planePoint.at(2) };
	double normalVectorPlane1[3];
	double vectorStartTargetpointPlane1[3];
	double vectorPlaneStartpointPlane1[3];

	vtkMath::Subtract(targetPointPlane1Array, startPointPlane1Array, vectorStartTargetpointPlane1);
	vtkMath::Subtract(planePointArray, startPointPlane1Array, vectorPlaneStartpointPlane1);
	vtkMath::Cross(vectorStartTargetpointPlane1, vectorPlaneStartpointPlane1, normalVectorPlane1);
	vtkMath::Normalize(normalVectorPlane1);

	double startPointPlane2Array[3];
	centrePoints->GetPoint(path2->GetId(0), startPointPlane2Array);
	double targetPointPlane2Array[3];
	centrePoints->GetPoint(path2->GetId(path2->GetNumberOfIds() - 1), targetPointPlane2Array);

	double normalVectorPlane2[3];
	double vectorStartTargetpoint2Plane[3];
	double vectorPlaneStartpoint2Plane[3];

	vtkMath::Subtract(targetPointPlane2Array, startPointPlane2Array, vectorStartTargetpoint2Plane);
	vtkMath::Subtract(planePointArray, startPointPlane2Array, vectorPlaneStartpoint2Plane);
	vtkMath::Cross(vectorStartTargetpoint2Plane, vectorPlaneStartpoint2Plane, normalVectorPlane2);
	vtkMath::Normalize(normalVectorPlane2);

	double normalVectorPlaneStartPoints[3];
	double vectorStartStartpoint[3];
	double vectorStartCentrepoint[3];

	vtkMath::Subtract(startPointPlane2Array, startPointPlane1Array, vectorStartStartpoint);
	vtkMath::Subtract(planePointArray, startPointPlane1Array, vectorStartCentrepoint);
	vtkMath::Cross(vectorStartCentrepoint, vectorStartStartpoint, normalVectorPlaneStartPoints);
	vtkMath::Normalize(normalVectorPlaneStartPoints);

	double normalVectorPlaneTargetPoints[3];
	double vectorTargetTargetpoint[3];
	double vectorTargetCentrepoint[3];

	vtkMath::Subtract(targetPointPlane2Array, targetPointPlane1Array, vectorTargetTargetpoint);
	vtkMath::Subtract(planePointArray, targetPointPlane1Array, vectorTargetCentrepoint);
	vtkMath::Cross(vectorTargetCentrepoint, vectorTargetTargetpoint, normalVectorPlaneTargetPoints);
	vtkMath::Normalize(normalVectorPlaneTargetPoints);

	vtkSmartPointer<vtkUnstructuredGrid> path1SearchList = vtkSmartPointer<vtkUnstructuredGrid>::New();
	vtkSmartPointer<vtkPoints> pointsPath1 = vtkSmartPointer<vtkPoints>::New();


	for (vtkIdType i = 0; i < path1->GetNumberOfIds(); i++) {
		double point[3] = { centrePoints->GetPoint(path1->GetId(i))[0], centrePoints->GetPoint(path1->GetId(i))[1], centrePoints->GetPoint(path1->GetId(i))
		[2] };
		pointsPath1->InsertNextPoint(point);
	}
	path1SearchList->SetPoints(pointsPath1);

	vtkSmartPointer<vtkUnstructuredGrid> path2SearchList = vtkSmartPointer<vtkUnstructuredGrid>::New();
	vtkSmartPointer<vtkPoints> pointsPath2 = vtkSmartPointer<vtkPoints>::New();

	for (vtkIdType i = 0; i < path2->GetNumberOfIds(); i++) {
		double point[3] = { centrePoints->GetPoint(path2->GetId(i))[0], centrePoints->GetPoint(path2->GetId(i))[1], centrePoints->GetPoint(path2->GetId(i))
		[2] };
		pointsPath2->InsertNextPoint(point);
	}
	path2SearchList->SetPoints(pointsPath2);
	set<vtkIdType> viewedCells;


	vtkIdType seedPointId = centrePointsUGrid->FindPoint(seedPoint.at(0), seedPoint.at(1), seedPoint.at(2));

	vtkSmartPointer<vtkIdList> cellIdList = vtkSmartPointer<vtkIdList>::New();
	cellIdList->InsertNextId(seedPointId);

	set<vtkIdType> returnList;

	do {
		vtkSmartPointer<vtkIdList> newcellIdList = vtkSmartPointer<vtkIdList>::New();

		for (vtkIdType k = 0; k < cellIdList->GetNumberOfIds(); k++) {
			vtkIdType currentCellIds = cellIdList->GetId(k);
			vtkSmartPointer<vtkIdList> cellPointIds = vtkSmartPointer<vtkIdList>::New();

			data.getVtkData()->GetCellPoints(currentCellIds, cellPointIds);

			for (vtkIdType i = 0; i < cellPointIds->GetNumberOfIds(); i++) {
				vtkSmartPointer<vtkIdList> CellPointIdList = vtkSmartPointer<vtkIdList>::New();
				CellPointIdList->InsertNextId(cellPointIds->GetId(i));

				vtkSmartPointer<vtkIdList> tempNeighborCellIds = vtkSmartPointer<vtkIdList>::New();
				data.getVtkData()->GetCellNeighbors(currentCellIds, CellPointIdList, tempNeighborCellIds);

				for (vtkIdType j = 0; j < tempNeighborCellIds->GetNumberOfIds(); j++) {
					if (viewedCells.find(tempNeighborCellIds->GetId(j)) == viewedCells.end()) {
						newcellIdList->InsertNextId(tempNeighborCellIds->GetId(j));
					}
					viewedCells.insert(tempNeighborCellIds->GetId(j));
				}
			}
		}

		vtkIdType size = newcellIdList->GetNumberOfIds();

		for (vtkIdType i = 0; i < size; i++) {
			vtkIdType currentId = newcellIdList->GetId(i);

			double currentPoint[3] = { 0, 0, 0 };
			currentPoint[0] = centrePoints->GetPoint(currentId)[0];
			currentPoint[1] = centrePoints->GetPoint(currentId)[1];
			currentPoint[2] = centrePoints->GetPoint(currentId)[2];

			double plane1 = vtkPlane::Evaluate(normalVectorPlane1, startPointPlane1Array, currentPoint);
			double plane2 = vtkPlane::Evaluate(normalVectorPlane2, startPointPlane2Array, currentPoint);
			double planeStart = vtkPlane::Evaluate(normalVectorPlaneStartPoints, startPointPlane1Array, currentPoint);
			double planetarget = vtkPlane::Evaluate(normalVectorPlaneTargetPoints, targetPointPlane1Array, currentPoint);


			if (!((plane1 >= 0) && (plane2 < 0) && (planeStart >= 0) && (planetarget <= 0))) {
				newcellIdList->DeleteId(currentId);
			}
		}

		for (vtkIdType j = 0; j < newcellIdList->GetNumberOfIds(); j++) {
			vtkIdType currentId = newcellIdList->GetId(j);
			returnList.insert(currentId);
			if (!data.getDoubleLayer()) {
				int material1 = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(currentId, 0);
				if ((data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(currentId, 0) == 0) &&
					(material1 == inMaterial)) {
					double currentPoint[3] = { 0, 0, 0 };
					currentPoint[0] = centrePoints->GetPoint(currentId)[0];
					currentPoint[1] = centrePoints->GetPoint(currentId)[1];
					currentPoint[2] = centrePoints->GetPoint(currentId)[2];

					vtkIdType closestPosPath1 = path1SearchList->FindPoint(currentPoint);
					vtkIdType closestIdPath1 = path1->GetId(closestPosPath1);

					double closestPointPath1[3] = { 0, 0, 0 };
					closestPointPath1[0] = centrePointsUGrid->GetPoint(closestIdPath1)[0];
					closestPointPath1[1] = centrePointsUGrid->GetPoint(closestIdPath1)[1];
					closestPointPath1[2] = centrePointsUGrid->GetPoint(closestIdPath1)[2];

					double distancePath1 = sqrt(abs(vtkMath::Distance2BetweenPoints(currentPoint, closestPointPath1)));

					vtkIdType closestPosPath2 = path2SearchList->FindPoint(currentPoint);
					vtkIdType closestIdPath2 = path2->GetId(closestPosPath2);

					double closestPointPath2[3] = { 0, 0, 0 };
					closestPointPath2[0] = centrePointsUGrid->GetPoint(closestIdPath2)[0];
					closestPointPath2[1] = centrePointsUGrid->GetPoint(closestIdPath2)[1];
					closestPointPath2[2] = centrePointsUGrid->GetPoint(closestIdPath2)[2];

					double distancePath2 = sqrt(abs(vtkMath::Distance2BetweenPoints(currentPoint, closestPointPath2)));

					vector<double> vectorPath1(3);
					vectorPath1.at(0) = data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath1)[0];
					vectorPath1.at(1) = data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath1)[1];
					vectorPath1.at(2) = data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath1)[2];


					vector<double> vectorPath2(3);
					vectorPath2.at(0) = data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath2)[0];
					vectorPath2.at(1) = data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath2)[1];
					vectorPath2.at(2) = data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath2)[2];

					vector<double> newOrientationVector = AveragingOrientation::averageOrientation(vectorPath1, vectorPath2, distancePath1, distancePath2);

					double newDiffVector[3] = { newOrientationVector.at(0), newOrientationVector.at(1), newOrientationVector.at(2) };
					vtkMath::Normalize(newDiffVector);

					data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->InsertTuple(currentId, newDiffVector);
					data.getVtkData()->GetCellData()->GetArray("Status")->SetComponent(currentId, 0, 3);
				}
				int material2 = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(currentId, 0);
				if ((data.getVtkData()->GetCellData()->GetArray("StatusEndo")->GetComponent(currentId, 0) == 0) &&
					(material2 == inMaterial)) {
					double currentPoint[3] = { 0, 0, 0 };
					currentPoint[0] = centrePoints->GetPoint(currentId)[0];
					currentPoint[1] = centrePoints->GetPoint(currentId)[1];
					currentPoint[2] = centrePoints->GetPoint(currentId)[2];

					vtkIdType closestPosPath1 = path1SearchList->FindPoint(currentPoint);
					vtkIdType closestIdPath1 = path1->GetId(closestPosPath1);

					double closestPointPath1[3] = { 0, 0, 0 };
					closestPointPath1[0] = centrePointsUGrid->GetPoint(closestIdPath1)[0];
					closestPointPath1[1] = centrePointsUGrid->GetPoint(closestIdPath1)[1];
					closestPointPath1[2] = centrePointsUGrid->GetPoint(closestIdPath1)[2];

					double distancePath1 = sqrt(abs(vtkMath::Distance2BetweenPoints(currentPoint, closestPointPath1)));

					vtkIdType closestPosPath2 = path2SearchList->FindPoint(currentPoint);
					vtkIdType closestIdPath2 = path2->GetId(closestPosPath2);

					double closestPointPath2[3] = { 0, 0, 0 };
					closestPointPath2[0] = centrePointsUGrid->GetPoint(closestIdPath2)[0];
					closestPointPath2[1] = centrePointsUGrid->GetPoint(closestIdPath2)[1];
					closestPointPath2[2] = centrePointsUGrid->GetPoint(closestIdPath2)[2];

					double distancePath2 = sqrt(abs(vtkMath::Distance2BetweenPoints(currentPoint, closestPointPath2)));

					vector<double> vectorPath1(3);
					vectorPath1.at(0) = data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath1)[0];
					vectorPath1.at(1) = data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath1)[1];
					vectorPath1.at(2) = data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath1)[2];


					vector<double> vectorPath2(3);
					vectorPath2.at(0) = data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath2)[0];
					vectorPath2.at(1) = data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath2)[1];
					vectorPath2.at(2) = data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath2)[2];

					vector<double> newOrientationVector = AveragingOrientation::averageOrientation(vectorPath1, vectorPath2, distancePath1, distancePath2);

					double newDiffVector[3] = { newOrientationVector.at(0), newOrientationVector.at(1), newOrientationVector.at(2) };
					vtkMath::Normalize(newDiffVector);

					data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->InsertTuple(currentId, newDiffVector);
					data.getVtkData()->GetCellData()->GetArray("StatusEndo")->SetComponent(currentId, 0, 3);
				}
			}
			else {
				int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(currentId, 0);
				if ((data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(currentId, 0) == 0) && (material == inMaterial)) {
					double currentPoint[3] = { 0, 0, 0 };
					currentPoint[0] = centrePoints->GetPoint(currentId)[0];
					currentPoint[1] = centrePoints->GetPoint(currentId)[1];
					currentPoint[2] = centrePoints->GetPoint(currentId)[2];

					vtkIdType closestPosPath1 = path1SearchList->FindPoint(currentPoint);
					vtkIdType closestIdPath1 = path1->GetId(closestPosPath1);

					double closestPointPath1[3] = { 0, 0, 0 };
					closestPointPath1[0] = centrePointsUGrid->GetPoint(closestIdPath1)[0];
					closestPointPath1[1] = centrePointsUGrid->GetPoint(closestIdPath1)[1];
					closestPointPath1[2] = centrePointsUGrid->GetPoint(closestIdPath1)[2];

					double distancePath1 = sqrt(abs(vtkMath::Distance2BetweenPoints(currentPoint, closestPointPath1)));

					vtkIdType closestPosPath2 = path2SearchList->FindPoint(currentPoint);
					vtkIdType closestIdPath2 = path2->GetId(closestPosPath2);

					double closestPointPath2[3] = { 0, 0, 0 };
					closestPointPath2[0] = centrePointsUGrid->GetPoint(closestIdPath2)[0];
					closestPointPath2[1] = centrePointsUGrid->GetPoint(closestIdPath2)[1];
					closestPointPath2[2] = centrePointsUGrid->GetPoint(closestIdPath2)[2];

					double distancePath2 = sqrt(abs(vtkMath::Distance2BetweenPoints(currentPoint, closestPointPath2)));

					vector<double> vectorPath1(3);
					vectorPath1.at(0) = data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath1)[0];
					vectorPath1.at(1) = data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath1)[1];
					vectorPath1.at(2) = data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath1)[2];


					vector<double> vectorPath2(3);
					vectorPath2.at(0) = data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath2)[0];
					vectorPath2.at(1) = data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath2)[1];
					vectorPath2.at(2) = data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath2)[2];

					vector<double> newOrientationVector = AveragingOrientation::averageOrientation(vectorPath1, vectorPath2, distancePath1, distancePath2);

					double newDiffVector[3] = { newOrientationVector.at(0), newOrientationVector.at(1), newOrientationVector.at(2) };
					vtkMath::Normalize(newDiffVector);

					data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->InsertTuple(currentId, newDiffVector);
					data.getVtkData()->GetCellData()->GetArray("Status")->SetComponent(currentId, 0, 3);
				}
			}
		}


		cellIdList->Reset();
		cellIdList = newcellIdList;
	} while (cellIdList->GetNumberOfIds() != 0);

	vtkSmartPointer<vtkIdList> returnIdList = vtkSmartPointer<vtkIdList>::New();

	for (set<vtkIdType>::iterator iter = returnList.begin(); iter != returnList.end(); iter++) {
		returnIdList->InsertNextId(*iter);
	}

	return returnIdList;
} // Methods::growBetweenPecti

/*!Close the fiber orientation in a region.
 \param data Pointer to the orginal mesh.
 \param seedPoint Seed point for the interpolation region.
 \retrun A vtkIdList with the cell-ids of the interpolate region.
 */
vtkSmartPointer<vtkIdList> Methods::regionGrow(DataFormat &data, std::vector<double> seedPoint) {
	double seedPointArray[3] = { seedPoint.at(0), seedPoint.at(1), seedPoint.at(2) };

	vtkSmartPointer<vtkPoints> centrePoints = data.getCentrePoints()->GetPoints();
	vtkSmartPointer<vtkUnstructuredGrid> centrePointsUGrid = data.getCentrePoints();

	vtkIdType seedPointID = centrePointsUGrid->FindPoint(seedPointArray);

	vtkSmartPointer<vtkDoubleArray> cells = vtkSmartPointer<vtkDoubleArray>::New();

	cells->SetNumberOfComponents(2);
	cells->SetComponentName(0, "CellID");
	cells->SetComponentName(1, "distance to centre point");

	set<vtkIdType> viewedCellList;

	bool end = false;

	vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
	set<vtkIdType> cellIdList;

	cellIdList.insert(seedPointID);

	if (data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(seedPointID, 0) == 0) {
		double cellsArray[2] = { 0, 0 };
		cellsArray[0] = seedPointID;
		cells->InsertNextTuple(cellsArray);
		cellIds->InsertNextId(seedPointID);
	}

	viewedCellList.insert(seedPointID);
	unordered_map<int, vtkSmartPointer<vtkIdList>> containerMaterialCellIDs;

	while (!end) {
		set<vtkIdType> newcellIdList;

		for (set<vtkIdType>::iterator iter = cellIdList.begin(); iter != cellIdList.end(); iter++) {
			vtkIdType currentCellIds = *iter;
			vtkSmartPointer<vtkIdList> cellPointIds = vtkSmartPointer<vtkIdList>::New();

			data.getVtkData()->GetCellPoints(currentCellIds, cellPointIds);

			for (vtkIdType i = 0; i < cellPointIds->GetNumberOfIds(); i++) {
				vtkSmartPointer<vtkIdList> CellPointIdList = vtkSmartPointer<vtkIdList>::New();
				CellPointIdList->InsertNextId(cellPointIds->GetId(i));

				vtkSmartPointer<vtkIdList> tempNeighborCellIds = vtkSmartPointer<vtkIdList>::New();
				data.getVtkData()->GetCellNeighbors(currentCellIds, CellPointIdList, tempNeighborCellIds);

				for (vtkIdType j = 0; j < tempNeighborCellIds->GetNumberOfIds(); j++) {
					if (viewedCellList.find(tempNeighborCellIds->GetId(j)) == viewedCellList.end()) {
						int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(tempNeighborCellIds->GetId(j), 0);

						if (Material::isInLeft(material) || Material::isInRight(material)) {
							if (data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(tempNeighborCellIds->GetId(j), 0) == 0) {
								double cellsArray[2] = { 0, 0 };
								cellsArray[0] = tempNeighborCellIds->GetId(j);
								cells->InsertNextTuple(cellsArray);
								newcellIdList.insert(tempNeighborCellIds->GetId(j));
								cellIds->InsertNextId(tempNeighborCellIds->GetId(j));
							}
							else if ((material != Material::Vorhof_links) &&
								(material != Material::Vorhof_links_Endo) &&
								(material != Material::Vorhof_rechts) &&
								(data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetComponent(tempNeighborCellIds->GetId(j), 0) +
									data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetComponent(tempNeighborCellIds->GetId(j), 1) +
									data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetComponent(tempNeighborCellIds->GetId(j), 2) != 0)
								) {
								if (containerMaterialCellIDs.find(material) == containerMaterialCellIDs.end()) {
									vtkSmartPointer<vtkIdList> containerCellIds = vtkSmartPointer<vtkIdList>::New();
									containerCellIds->InsertNextId(tempNeighborCellIds->GetId(j));
									pair<int, vtkSmartPointer<vtkIdList>> insert(material, containerCellIds);
									containerMaterialCellIDs.insert(insert);
								}
								else {
									vtkSmartPointer<vtkIdList> containerCellIds = containerMaterialCellIDs.at(material);
									containerCellIds->InsertNextId(tempNeighborCellIds->GetId(j));
									containerMaterialCellIDs.at(material) = containerCellIds;
								}
							}
						}
					}

					viewedCellList.insert(tempNeighborCellIds->GetId(j));
				}
			}
		}
		if (newcellIdList.size() == 0) {
			end = true;
		}
		cellIdList.clear();
		cellIdList = newcellIdList;
	}

	vector<double> centroid = getCentroid(data, cellIds);

	double centroidArray[3] = { centroid.at(0), centroid.at(1), centroid.at(2) };

	for (vtkIdType i = 0; i < cells->GetNumberOfTuples(); i++) {
		double distance = sqrt(abs(vtkMath::Distance2BetweenPoints(centroidArray, centrePoints->GetPoint(cells->vtkDataArray::GetComponent(i, 0)))));
		cells->vtkDataArray::SetComponent(i, 1, distance);
	}


	vtkSortDataArray::SortArrayByComponent(cells, 1);

	unordered_map<int, vtkSmartPointer<vtkUnstructuredGrid>> containerMaterialSearchPath;
	unordered_map<int, vector<double>> mainOrientation;

	if (containerMaterialCellIDs.size() > 0) {
		if (containerMaterialCellIDs.size() < 2) {
			unordered_map<int, vtkSmartPointer<vtkIdList>>::const_iterator iter = containerMaterialCellIDs.begin();
			vtkSmartPointer<vtkUnstructuredGrid> SearchList = vtkSmartPointer<vtkUnstructuredGrid>::New();
			vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

			for (vtkIdType i = 0; i < iter->second->GetNumberOfIds(); i++) {
				points->InsertNextPoint(centrePoints->GetPoint(iter->second->GetId(i)));
			}
			SearchList->SetPoints(points);

			for (vtkIdType j = cells->GetNumberOfTuples() - 1; j >= 0; j--) {
				vtkIdType currentId = cells->vtkDataArray::GetComponent(j, 0);

				double currentPoint[3] = { 0, 0, 0 };
				currentPoint[0] = centrePoints->GetPoint(currentId)[0];
				currentPoint[1] = centrePoints->GetPoint(currentId)[1];
				currentPoint[2] = centrePoints->GetPoint(currentId)[2];

				vtkIdType closestPosPath = SearchList->FindPoint(currentPoint);
				vtkIdType closestIdPoint = iter->second->GetId(closestPosPath);

				double *newDiffVector = data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetTuple(closestIdPoint);
				vtkMath::Normalize(newDiffVector);

				data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->InsertTuple(currentId, newDiffVector);
				data.getVtkData()->GetCellData()->GetArray("Status")->SetComponent(currentId, 0, 3);
				data.getVtkData()->GetCellData()->GetArray("Material")->SetComponent(currentId, 0, iter->first);
			}
		}
		else {
			for (unordered_map<int, vtkSmartPointer<vtkIdList>>::const_iterator iter = containerMaterialCellIDs.begin();
				iter != containerMaterialCellIDs.end(); iter++) {
				vtkSmartPointer<vtkUnstructuredGrid> SearchList = vtkSmartPointer<vtkUnstructuredGrid>::New();
				vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

				double sumMainOrientationArray[3] = { 0, 0, 0 };

				for (vtkIdType i = 0; i < iter->second->GetNumberOfIds(); i++) {
					points->InsertNextPoint(centrePoints->GetPoint(iter->second->GetId(i)));
					double *orientationArray = data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetTuple(iter->second->GetId(i));
					vtkMath::Add(sumMainOrientationArray, orientationArray, sumMainOrientationArray);
				}
				SearchList->SetPoints(points);

				vtkMath::Normalize(sumMainOrientationArray);
				vector<double> sumMainOrientation = { sumMainOrientationArray[0], sumMainOrientationArray[1], sumMainOrientationArray[2] };


				pair<int, vector<double>> insertMainOrientation(iter->first, sumMainOrientation);
				mainOrientation.insert(insertMainOrientation);

				pair<int, vtkSmartPointer<vtkUnstructuredGrid>> insert(iter->first, SearchList);
				containerMaterialSearchPath.insert(insert);
			}

			for (vtkIdType j = cells->GetNumberOfTuples() - 1; j >= 0; j--) {
				vtkIdType currentId = cells->vtkDataArray::GetComponent(j, 0);

				double currentPoint[3] = { 0, 0, 0 };
				currentPoint[0] = centrePoints->GetPoint(currentId)[0];
				currentPoint[1] = centrePoints->GetPoint(currentId)[1];
				currentPoint[2] = centrePoints->GetPoint(currentId)[2];

				list<tuple<int, double, vector<double>, vector<double>>> differenceVectorList;
				double sumOfDistance = 0;
				int newMaterial = -1;

				unordered_map<int, vector<double>>::const_iterator iterMainOrientation = mainOrientation.begin();
				for (unordered_map<int, vtkSmartPointer<vtkUnstructuredGrid>>::const_iterator iter = containerMaterialSearchPath.begin();
					iter != containerMaterialSearchPath.end(); iter++) {
					vtkIdType closestPosPath = iter->second->FindPoint(currentPoint);
					vtkIdType closestIdPoint = containerMaterialCellIDs.find(iter->first)->second->GetId(closestPosPath);

					double closestPoint[3] = { 0, 0, 0 };
					closestPoint[0] = centrePoints->GetPoint(closestIdPoint)[0];
					closestPoint[1] = centrePoints->GetPoint(closestIdPoint)[1];
					closestPoint[2] = centrePoints->GetPoint(closestIdPoint)[2];

					double distance = sqrt(abs(vtkMath::Distance2BetweenPoints(currentPoint, closestPoint)));
					double *materialDiffArray = data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetTuple(closestIdPoint);

					vector<double> materialDiffVector = { materialDiffArray[0], materialDiffArray[1], materialDiffArray[2] };

					tuple<int, double, vector<double>, vector<double>> tupleDiffVector(iter->first, distance, materialDiffVector, iterMainOrientation->second);
					iterMainOrientation++;

					differenceVectorList.push_back(tupleDiffVector);
				}

				double minDistance = std::numeric_limits<double>::max();

				for (list<tuple<int, double, vector<double>, vector<double>>>::iterator iter = differenceVectorList.begin();
					iter != differenceVectorList.end(); iter++) {
					if (minDistance > get<1>(*iter)) {
						minDistance = get<1>(*iter);
						newMaterial = get<0>(*iter);
					}

					if (get<1>(*iter) > 0) {
						get<1>(*iter) = 1 / get<1>(*iter);
					}
					else {
						get<1>(*iter) = 1;
					}
					sumOfDistance += get<1>(*iter);
				}


				for (list<tuple<int, double, vector<double>, vector<double>>>::iterator iter = differenceVectorList.begin();
					iter != differenceVectorList.end(); iter++) {
					get<1>(*iter) = get<1>(*iter) / sumOfDistance;
				}

				vector<double> newOrientationVector = AveragingOrientation::AveragingNVectors(differenceVectorList);
				double newDiffVector[3] = { newOrientationVector.at(0), newOrientationVector.at(1), newOrientationVector.at(2) };

				data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->InsertTuple(currentId, newDiffVector);
				data.getVtkData()->GetCellData()->GetArray("Status")->SetComponent(currentId, 0, 3);
				data.getVtkData()->GetCellData()->GetArray("Material")->SetComponent(currentId, 0, newMaterial);
			}
		}
	}
	return cellIds;
} // Methods::regionGrow

/*!Close the fiber orientation at the free bridges.
 \param data Pointer to the orginal mesh.
 \param seedPoint Seed point for the interpolation region.
 \param freeBridgeMaterial Tissue class of the free bridge material.
 \retrun A vtkIdList with the cell-ids of the interpolate region.
 */
vtkSmartPointer<vtkIdList> Methods::regionGrowFreeBridge(DataFormat &data, std::vector<double> seedPoint, int freeBridgeMaterial) {
	double seedPointArray[3] = { seedPoint.at(0), seedPoint.at(1), seedPoint.at(2) };

	vtkSmartPointer<vtkPoints> centrePoints = data.getCentrePoints()->GetPoints();
	vtkSmartPointer<vtkUnstructuredGrid> centrePointsUGrid = data.getCentrePoints();

	vtkIdType seedPointID = centrePointsUGrid->FindPoint(seedPointArray);

	vtkSmartPointer<vtkDoubleArray> cells = vtkSmartPointer<vtkDoubleArray>::New();

	cells->SetNumberOfComponents(2);
	cells->SetComponentName(0, "CellID");
	cells->SetComponentName(1, "distance to centre point");

	set<vtkIdType> viewedCellList;

	bool end = false;

	vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
	set<vtkIdType> cellIdList;

	cellIdList.insert(seedPointID);

	if (data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(seedPointID, 0) == 0) {
		double cellsArray[2] = { 0, 0 };
		cellsArray[0] = seedPointID;
		cells->InsertNextTuple(cellsArray);
		cellIds->InsertNextId(seedPointID);
	}

	viewedCellList.insert(seedPointID);
	unordered_map<int, vtkSmartPointer<vtkIdList>> containerMaterialCellIDs;

	while (!end) {
		set<vtkIdType> newcellIdList;

		for (set<vtkIdType>::iterator iter = cellIdList.begin(); iter != cellIdList.end(); iter++) {
			vtkIdType currentCellIds = *iter;
			vtkSmartPointer<vtkIdList> cellPointIds = vtkSmartPointer<vtkIdList>::New();

			data.getVtkData()->GetCellPoints(currentCellIds, cellPointIds);

			for (vtkIdType i = 0; i < cellPointIds->GetNumberOfIds(); i++) {
				vtkSmartPointer<vtkIdList> CellPointIdList = vtkSmartPointer<vtkIdList>::New();
				CellPointIdList->InsertNextId(cellPointIds->GetId(i));

				vtkSmartPointer<vtkIdList> tempNeighborCellIds = vtkSmartPointer<vtkIdList>::New();
				data.getVtkData()->GetCellNeighbors(currentCellIds, CellPointIdList, tempNeighborCellIds);

				for (vtkIdType j = 0; j < tempNeighborCellIds->GetNumberOfIds(); j++) {
					if (viewedCellList.find(tempNeighborCellIds->GetId(j)) == viewedCellList.end()) {
						int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(tempNeighborCellIds->GetId(j), 0);

						if (Material::isInLeft(material) || Material::isInRight(material) || (material == freeBridgeMaterial)) {
							if (data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(tempNeighborCellIds->GetId(j), 0) == 0) {
								double cellsArray[2] = { 0, 0 };
								cellsArray[0] = tempNeighborCellIds->GetId(j);
								cells->InsertNextTuple(cellsArray);
								newcellIdList.insert(tempNeighborCellIds->GetId(j));
								cellIds->InsertNextId(tempNeighborCellIds->GetId(j));
							}
							else if ((material != Material::Vorhof_links) &&
								(material != Material::Vorhof_links_Endo) &&
								(material != Material::Vorhof_rechts) &&
								(data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetComponent(tempNeighborCellIds->GetId(j), 0) +
									data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetComponent(tempNeighborCellIds->GetId(j), 1) +
									data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetComponent(tempNeighborCellIds->GetId(j), 2) != 0)
								) {
								if (containerMaterialCellIDs.find(material) == containerMaterialCellIDs.end()) {
									vtkSmartPointer<vtkIdList> containerCellIds = vtkSmartPointer<vtkIdList>::New();
									containerCellIds->InsertNextId(tempNeighborCellIds->GetId(j));
									pair<int, vtkSmartPointer<vtkIdList>> insert(material, containerCellIds);
									containerMaterialCellIDs.insert(insert);
								}
								else {
									vtkSmartPointer<vtkIdList> containerCellIds = containerMaterialCellIDs.at(material);
									containerCellIds->InsertNextId(tempNeighborCellIds->GetId(j));
									containerMaterialCellIDs.at(material) = containerCellIds;
								}
							}
						}
					}

					viewedCellList.insert(tempNeighborCellIds->GetId(j));
				}
			}
		}
		if (newcellIdList.size() == 0) {
			end = true;
		}
		cellIdList.clear();
		cellIdList = newcellIdList;
	}

	vector<double> centroid = getCentroid(data, cellIds);

	double centroidArray[3] = { centroid.at(0), centroid.at(1), centroid.at(2) };

	for (vtkIdType i = 0; i < cells->GetNumberOfTuples(); i++) {
		double distance = sqrt(abs(vtkMath::Distance2BetweenPoints(centroidArray, centrePoints->GetPoint(cells->vtkDataArray::GetComponent(i, 0)))));
		cells->vtkDataArray::SetComponent(i, 1, distance);
	}


	vtkSortDataArray::SortArrayByComponent(cells, 1);

	unordered_map<int, vtkSmartPointer<vtkUnstructuredGrid>> containerMaterialSearchPath;
	unordered_map<int, vector<double>> mainOrientation;

	if (containerMaterialCellIDs.size() > 0) {
		if (containerMaterialCellIDs.size() < 2) {
			unordered_map<int, vtkSmartPointer<vtkIdList>>::const_iterator iter = containerMaterialCellIDs.begin();
			vtkSmartPointer<vtkUnstructuredGrid> SearchList = vtkSmartPointer<vtkUnstructuredGrid>::New();
			vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

			for (vtkIdType i = 0; i < iter->second->GetNumberOfIds(); i++) {
				points->InsertNextPoint(centrePoints->GetPoint(iter->second->GetId(i)));
			}
			SearchList->SetPoints(points);

			for (vtkIdType j = cells->GetNumberOfTuples() - 1; j >= 0; j--) {
				vtkIdType currentId = cells->vtkDataArray::GetComponent(j, 0);

				double currentPoint[3] = { 0, 0, 0 };
				currentPoint[0] = centrePoints->GetPoint(currentId)[0];
				currentPoint[1] = centrePoints->GetPoint(currentId)[1];
				currentPoint[2] = centrePoints->GetPoint(currentId)[2];

				vtkIdType closestPosPath = SearchList->FindPoint(currentPoint);
				vtkIdType closestIdPoint = iter->second->GetId(closestPosPath);

				double *newDiffVector = data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetTuple(closestIdPoint);
				vtkMath::Normalize(newDiffVector);

				data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->InsertTuple(currentId, newDiffVector);
				data.getVtkData()->GetCellData()->GetArray("Status")->SetComponent(currentId, 0, 3);
				data.getVtkData()->GetCellData()->GetArray("Material")->SetComponent(currentId, 0, iter->first);
			}
		}
		else {
			for (unordered_map<int, vtkSmartPointer<vtkIdList>>::const_iterator iter = containerMaterialCellIDs.begin();
				iter != containerMaterialCellIDs.end(); iter++) {
				vtkSmartPointer<vtkUnstructuredGrid> SearchList = vtkSmartPointer<vtkUnstructuredGrid>::New();
				vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

				double sumMainOrientationArray[3] = { 0, 0, 0 };

				for (vtkIdType i = 0; i < iter->second->GetNumberOfIds(); i++) {
					points->InsertNextPoint(centrePoints->GetPoint(iter->second->GetId(i)));
					double *orientationArray = data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetTuple(iter->second->GetId(i));
					vtkMath::Add(sumMainOrientationArray, orientationArray, sumMainOrientationArray);
				}
				SearchList->SetPoints(points);

				vtkMath::Normalize(sumMainOrientationArray);
				vector<double> sumMainOrientation = { sumMainOrientationArray[0], sumMainOrientationArray[1], sumMainOrientationArray[2] };


				pair<int, vector<double>> insertMainOrientation(iter->first, sumMainOrientation);
				mainOrientation.insert(insertMainOrientation);

				pair<int, vtkSmartPointer<vtkUnstructuredGrid>> insert(iter->first, SearchList);
				containerMaterialSearchPath.insert(insert);
			}

			for (vtkIdType j = cells->GetNumberOfTuples() - 1; j >= 0; j--) {
				vtkIdType currentId = cells->vtkDataArray::GetComponent(j, 0);

				double currentPoint[3] = { 0, 0, 0 };
				currentPoint[0] = centrePoints->GetPoint(currentId)[0];
				currentPoint[1] = centrePoints->GetPoint(currentId)[1];
				currentPoint[2] = centrePoints->GetPoint(currentId)[2];

				list<tuple<int, double, vector<double>, vector<double>>> differenceVectorList;
				double sumOfDistance = 0;
				int newMaterial = -1;

				unordered_map<int, vector<double>>::const_iterator iterMainOrientation = mainOrientation.begin();
				for (unordered_map<int, vtkSmartPointer<vtkUnstructuredGrid>>::const_iterator iter = containerMaterialSearchPath.begin();
					iter != containerMaterialSearchPath.end(); iter++) {
					vtkIdType closestPosPath = iter->second->FindPoint(currentPoint);
					vtkIdType closestIdPoint = containerMaterialCellIDs.find(iter->first)->second->GetId(closestPosPath);

					double closestPoint[3] = { 0, 0, 0 };
					closestPoint[0] = centrePoints->GetPoint(closestIdPoint)[0];
					closestPoint[1] = centrePoints->GetPoint(closestIdPoint)[1];
					closestPoint[2] = centrePoints->GetPoint(closestIdPoint)[2];

					double distance = sqrt(abs(vtkMath::Distance2BetweenPoints(currentPoint, closestPoint)));
					double *materialDiffArray = data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetTuple(closestIdPoint);

					vector<double> materialDiffVector = { materialDiffArray[0], materialDiffArray[1], materialDiffArray[2] };

					tuple<int, double, vector<double>, vector<double>> tupleDiffVector(iter->first, distance, materialDiffVector, iterMainOrientation->second);
					iterMainOrientation++;

					differenceVectorList.push_back(tupleDiffVector);
				}

				double minDistance = std::numeric_limits<double>::max();

				for (list<tuple<int, double, vector<double>, vector<double>>>::iterator iter = differenceVectorList.begin();
					iter != differenceVectorList.end(); iter++) {
					if (minDistance > get<1>(*iter)) {
						minDistance = get<1>(*iter);
						newMaterial = get<0>(*iter);
					}

					if (get<1>(*iter) > 0) {
						get<1>(*iter) = 1 / get<1>(*iter);
					}
					else {
						get<1>(*iter) = 1;
					}
					sumOfDistance += get<1>(*iter);
				}


				for (list<tuple<int, double, vector<double>, vector<double>>>::iterator iter = differenceVectorList.begin();
					iter != differenceVectorList.end(); iter++) {
					get<1>(*iter) = get<1>(*iter) / sumOfDistance;
				}

				vector<double> newOrientationVector = AveragingOrientation::AveragingNVectors(differenceVectorList);
				double newDiffVector[3] = { newOrientationVector.at(0), newOrientationVector.at(1), newOrientationVector.at(2) };

				data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->InsertTuple(currentId, newDiffVector);
				data.getVtkData()->GetCellData()->GetArray("Status")->SetComponent(currentId, 0, 3);
				data.getVtkData()->GetCellData()->GetArray("Material")->SetComponent(currentId, 0, newMaterial);
			}
		}
	}
	return cellIds;
} // Methods::regionGrowFreeBridge

/*! Grow the only the tissue class in a region without a fiber orientation.
 \param data Pointer to the orginal mesh.
 \param seedPoint Seed point for the interpolation region.
 */

void Methods::regionGrowMaterial(DataFormat &data, std::vector<double> seedPoint) {
	double seedPointArray[3] = { seedPoint.at(0), seedPoint.at(1), seedPoint.at(2) };

	vtkSmartPointer<vtkPoints> centrePoints = data.getCentrePoints()->GetPoints();
	vtkSmartPointer<vtkUnstructuredGrid> centrePointsUGrid = data.getCentrePoints();

	vtkSmartPointer<vtkUnstructuredGrid> searchPath = vtkSmartPointer<vtkUnstructuredGrid>::New();
	vtkSmartPointer<vtkPoints> pointsSearchPath = vtkSmartPointer<vtkPoints>::New();

	vtkIdType seedPointID = centrePointsUGrid->FindPoint(seedPointArray);

	vtkSmartPointer<vtkDoubleArray> cells = vtkSmartPointer<vtkDoubleArray>::New();

	cells->SetNumberOfComponents(2);
	cells->SetComponentName(0, "CellID");
	cells->SetComponentName(1, "distance to centre point");

	set<vtkIdType> viewedCellList;

	bool end = false;

	vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
	set<vtkIdType> cellIdList;

	cellIdList.insert(seedPointID);

	while (!end) {
		set<vtkIdType> newcellIdList;

		for (set<vtkIdType>::iterator iter = cellIdList.begin(); iter != cellIdList.end(); iter++) {
			vtkIdType currentCellIds = *iter;
			vtkSmartPointer<vtkIdList> cellPointIds = vtkSmartPointer<vtkIdList>::New();

			data.getVtkData()->GetCellPoints(currentCellIds, cellPointIds);

			for (vtkIdType i = 0; i < cellPointIds->GetNumberOfIds(); i++) {
				vtkSmartPointer<vtkIdList> CellPointIdList = vtkSmartPointer<vtkIdList>::New();
				CellPointIdList->InsertNextId(cellPointIds->GetId(i));

				vtkSmartPointer<vtkIdList> tempNeighborCellIds = vtkSmartPointer<vtkIdList>::New();
				data.getVtkData()->GetCellNeighbors(currentCellIds, CellPointIdList, tempNeighborCellIds);

				for (vtkIdType j = 0; j < tempNeighborCellIds->GetNumberOfIds(); j++) {
					if (viewedCellList.find(tempNeighborCellIds->GetId(j)) == viewedCellList.end()) {
						int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(tempNeighborCellIds->GetId(j), 0);

						if (Material::undef == material) {
							double cellsArray[2] = { 0, 0 };
							cellsArray[0] = tempNeighborCellIds->GetId(j);
							cells->InsertNextTuple(cellsArray);
							newcellIdList.insert(tempNeighborCellIds->GetId(j));
							cellIds->InsertNextId(tempNeighborCellIds->GetId(j));
						}
						else {
							pointsSearchPath->InsertNextPoint(data.getCentrePoints()->GetPoint(tempNeighborCellIds->GetId(j)));
						}
						viewedCellList.insert(tempNeighborCellIds->GetId(j));
					}
				}
			}
		}
		if (newcellIdList.size() == 0) {
			end = true;
		}
		cellIdList.clear();
		cellIdList = newcellIdList;
	}

	searchPath->SetPoints(pointsSearchPath);

	vector<double> centroid = getCentroid(data, cellIds);

	double centroidArray[3] = { centroid.at(0), centroid.at(1), centroid.at(2) };

	for (vtkIdType i = 0; i < cells->GetNumberOfTuples(); i++) {
		double distance = sqrt(abs(vtkMath::Distance2BetweenPoints(centroidArray, centrePoints->GetPoint(cells->vtkDataArray::GetComponent(i, 0)))));
		cells->vtkDataArray::SetComponent(i, 1, distance);
	}

	vtkSortDataArray::SortArrayByComponent(cells, 1);

	for (vtkIdType j = cells->GetNumberOfTuples() - 1; j >= 0; j--) {
		vtkIdType currentId = cells->vtkDataArray::GetComponent(j, 0);
		vtkIdType searchPathId = searchPath->FindPoint(data.getCentrePoints()->GetPoint(currentId));
		vtkIdType dataPointId = data.getCentrePoints()->FindPoint(searchPath->GetPoint(searchPathId));

		double material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(dataPointId, 0);
		data.getVtkData()->GetCellData()->GetArray("Material")->SetComponent(currentId, 0, material);
	}
} // Methods::regionGrowMaterial

/*! Grow the only the tissue class in a region without a fiber orientation.
 \param data Pointer to the orginal mesh.
 \param seedPoint Seed point for the interpolation region.
 */

void Methods::regionGrowMaterial(DataFormat &data, std::vector<double> seedPoint, Material::Mat fromMaterial1, Material::Mat fromMaterial2, Material::Mat inMaterial) {
	double seedPointArray[3] = { seedPoint.at(0), seedPoint.at(1), seedPoint.at(2) };

	vtkSmartPointer<vtkPoints> centrePoints = data.getCentrePoints()->GetPoints();
	vtkSmartPointer<vtkUnstructuredGrid> centrePointsUGrid = data.getCentrePoints();

	vtkSmartPointer<vtkUnstructuredGrid> searchPath = vtkSmartPointer<vtkUnstructuredGrid>::New();
	vtkSmartPointer<vtkPoints> pointsSearchPath = vtkSmartPointer<vtkPoints>::New();

	vtkIdType seedPointID = centrePointsUGrid->FindPoint(seedPointArray);

	vtkSmartPointer<vtkDoubleArray> cells = vtkSmartPointer<vtkDoubleArray>::New();

	cells->SetNumberOfComponents(2);
	cells->SetComponentName(0, "CellID");
	cells->SetComponentName(1, "distance to centre point");

	set<vtkIdType> viewedCellList;

	bool end = false;

	vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
	set<vtkIdType> cellIdList;

	cellIdList.insert(seedPointID);

	while (!end) {
		set<vtkIdType> newcellIdList;

		for (set<vtkIdType>::iterator iter = cellIdList.begin(); iter != cellIdList.end(); iter++) {
			vtkIdType currentCellIds = *iter;
			vtkSmartPointer<vtkIdList> cellPointIds = vtkSmartPointer<vtkIdList>::New();

			data.getVtkData()->GetCellPoints(currentCellIds, cellPointIds);

			for (vtkIdType i = 0; i < cellPointIds->GetNumberOfIds(); i++) {
				vtkSmartPointer<vtkIdList> CellPointIdList = vtkSmartPointer<vtkIdList>::New();
				CellPointIdList->InsertNextId(cellPointIds->GetId(i));

				vtkSmartPointer<vtkIdList> tempNeighborCellIds = vtkSmartPointer<vtkIdList>::New();
				data.getVtkData()->GetCellNeighbors(currentCellIds, CellPointIdList, tempNeighborCellIds);

				for (vtkIdType j = 0; j < tempNeighborCellIds->GetNumberOfIds(); j++) {
					if (viewedCellList.find(tempNeighborCellIds->GetId(j)) == viewedCellList.end()) {
						int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(tempNeighborCellIds->GetId(j), 0);

						if (inMaterial == material) {
							double cellsArray[2] = { 0, 0 };
							cellsArray[0] = tempNeighborCellIds->GetId(j);
							cells->InsertNextTuple(cellsArray);
							newcellIdList.insert(tempNeighborCellIds->GetId(j));
							cellIds->InsertNextId(tempNeighborCellIds->GetId(j));
						}
						else if ((fromMaterial1 == material) || (fromMaterial2 == material)) {
							pointsSearchPath->InsertNextPoint(data.getCentrePoints()->GetPoint(tempNeighborCellIds->GetId(j)));
						}
						viewedCellList.insert(tempNeighborCellIds->GetId(j));
					}
				}
			}
		}
		if (newcellIdList.size() == 0) {
			end = true;
		}
		cellIdList.clear();
		cellIdList = newcellIdList;
	}

	searchPath->SetPoints(pointsSearchPath);

	vector<double> centroid = getCentroid(data, cellIds);

	double centroidArray[3] = { centroid.at(0), centroid.at(1), centroid.at(2) };

	for (vtkIdType i = 0; i < cells->GetNumberOfTuples(); i++) {
		double distance = sqrt(abs(vtkMath::Distance2BetweenPoints(centroidArray, centrePoints->GetPoint(cells->vtkDataArray::GetComponent(i, 0)))));
		cells->vtkDataArray::SetComponent(i, 1, distance);
	}

	vtkSortDataArray::SortArrayByComponent(cells, 1);

	for (vtkIdType j = cells->GetNumberOfTuples() - 1; j >= 0; j--) {
		vtkIdType currentId = cells->vtkDataArray::GetComponent(j, 0);


		vtkIdType searchPathId = searchPath->FindPoint(data.getCentrePoints()->GetPoint(currentId));
		vtkIdType dataPointId = data.getCentrePoints()->FindPoint(searchPath->GetPoint(searchPathId));

		double material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(dataPointId, 0);
		data.getVtkData()->GetCellData()->GetArray("Material")->SetComponent(currentId, 0, material);
		data.getVtkData()->GetCellData()->GetArray("Status")->SetComponent(currentId, 0, 5);
	}
} // Methods::regionGrowMaterial

/*! Define a bool equation. If the first element of the tupele is euqal than the first element of the second tupel ->
 than true otherwise false. */
bool equal_differenceVectorList(tuple<int, double, vector<double>, vector<double>> frist, tuple<int, double, vector<double>, vector<double>> second)
{
	return get<0>(frist) == get<0>(second);
}

/*! Close the fiber orientation in a region at a surface mesh.
 \param data Pointer to the orginal mesh.
 \param seedPoint Seed point for the interpolation region.
 \param distanceEndoEpi Fictive distance between the endocard and epicard in the left atrial.
 \retrun A vtkIdList with the cell-ids of the interpolate region.
 */
vtkSmartPointer<vtkIdList> Methods::regionGrowOneSurface(DataFormat &data, std::vector<double> seedPoint, double distanceEndoEpi) {
	double seedPointArray[3] = { seedPoint.at(0), seedPoint.at(1), seedPoint.at(2) };

	vtkSmartPointer<vtkPoints> centrePoints = data.getCentrePoints()->GetPoints();

	vtkIdType seedPointID = data.getCentrePoints()->FindPoint(seedPointArray);

	vtkSmartPointer<vtkDoubleArray> cells = vtkSmartPointer<vtkDoubleArray>::New();

	cells->SetNumberOfComponents(4);
	cells->SetComponentName(0, "CellID");
	cells->SetComponentName(1, "distance to centre point");
	cells->SetComponentName(2, "Epi");
	cells->SetComponentName(3, "Endo");

	set<vtkIdType> viewedCellList;

	bool end = false;

	vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
	set<vtkIdType> cellIdList;

	cellIdList.insert(seedPointID);

	if ((data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(seedPointID, 0) == 0) ||
		(data.getVtkData()->GetCellData()->GetArray("StatusEndo")->GetComponent(seedPointID, 0) == 0)) {
		double cellsArray[4] = { 0, 0, 0, 0 };
		cellsArray[0] = seedPointID;
		if (data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(seedPointID, 0) == 0) {
			cellsArray[2] = 1;
		}
		if (data.getVtkData()->GetCellData()->GetArray("StatusEndo")->GetComponent(seedPointID, 0) == 0) {
			cellsArray[3] = 1;
		}
		cells->InsertNextTuple(cellsArray);
		cellIds->InsertNextId(seedPointID);
	}

	viewedCellList.insert(seedPointID);
	unordered_map<int, vtkSmartPointer<vtkIdList>> containerMaterialCellIDs;
	unordered_map<int, vtkSmartPointer<vtkIdList>> containerMaterialCellIDsEndo;


	while (!end) {
		set<vtkIdType> newcellIdList;

		for (set<vtkIdType>::iterator iter = cellIdList.begin(); iter != cellIdList.end(); iter++) {
			vtkIdType currentCellIds = *iter;
			vtkSmartPointer<vtkIdList> cellPointIds = vtkSmartPointer<vtkIdList>::New();

			data.getVtkData()->GetCellPoints(currentCellIds, cellPointIds);

			for (vtkIdType i = 0; i < cellPointIds->GetNumberOfIds(); i++) {
				vtkSmartPointer<vtkIdList> CellPointIdList = vtkSmartPointer<vtkIdList>::New();
				CellPointIdList->InsertNextId(cellPointIds->GetId(i));

				vtkSmartPointer<vtkIdList> tempNeighborCellIds = vtkSmartPointer<vtkIdList>::New();
				data.getVtkData()->GetCellNeighbors(currentCellIds, CellPointIdList, tempNeighborCellIds);

				for (vtkIdType j = 0; j < tempNeighborCellIds->GetNumberOfIds(); j++) {
					if (viewedCellList.find(tempNeighborCellIds->GetId(j)) == viewedCellList.end()) {
						int material1 = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(tempNeighborCellIds->GetId(j), 0);
						int material2 = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(tempNeighborCellIds->GetId(j), 0);

						if (Material::isInLeft(material1) || Material::isInRight(material1) || Material::isInLeft(material2) ||
							Material::isInRight(material2)) {
							if ((data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(tempNeighborCellIds->GetId(j), 0) == 0) ||
								(data.getVtkData()->GetCellData()->GetArray("StatusEndo")->GetComponent(tempNeighborCellIds->GetId(j), 0) == 0)) {
								double cellsArray[4] = { 0, 0, 0, 0 };
								cellsArray[0] = tempNeighborCellIds->GetId(j);
								if (data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(tempNeighborCellIds->GetId(j), 0) == 0) {
									cellsArray[2] = 1;
								}
								if (data.getVtkData()->GetCellData()->GetArray("StatusEndo")->GetComponent(tempNeighborCellIds->GetId(j), 0) == 0) {
									cellsArray[3] = 1;
								}
								cells->InsertNextTuple(cellsArray);
								newcellIdList.insert(tempNeighborCellIds->GetId(j));
								cellIds->InsertNextId(tempNeighborCellIds->GetId(j));
							}


							if ((material1 != Material::Vorhof_links) &&
								(material1 != Material::Vorhof_links_Endo) &&
								(material1 != Material::Vorhof_rechts) &&
								(data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(tempNeighborCellIds->GetId(j), 0) != 0) &&
								(data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetComponent(tempNeighborCellIds->
									GetId(j), 0) +
									data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetComponent(tempNeighborCellIds->
										GetId(j), 1) +
									data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetComponent(tempNeighborCellIds->
										GetId(j), 2) != 0)
								) {
								if (containerMaterialCellIDs.find(material1) == containerMaterialCellIDs.end()) {
									vtkSmartPointer<vtkIdList> containerCellIds = vtkSmartPointer<vtkIdList>::New();
									containerCellIds->InsertNextId(tempNeighborCellIds->GetId(j));
									pair<int, vtkSmartPointer<vtkIdList>> insert(material1, containerCellIds);
									containerMaterialCellIDs.insert(insert);
								}
								else {
									vtkSmartPointer<vtkIdList> containerCellIds = containerMaterialCellIDs.at(material1);
									containerCellIds->InsertUniqueId(tempNeighborCellIds->GetId(j));
									containerMaterialCellIDs.at(material1) = containerCellIds;
								}
							}
							if ((material2 != Material::Vorhof_links) &&
								(material2 != Material::Vorhof_links_Endo) &&
								(material2 != Material::Vorhof_rechts) &&
								(data.getVtkData()->GetCellData()->GetArray("StatusEndo")->GetComponent(tempNeighborCellIds->GetId(j), 0) != 0) &&
								(data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->GetComponent(tempNeighborCellIds
									->GetId(j), 0) +
									data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->GetComponent(tempNeighborCellIds
										->GetId(j), 1) +
									data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->GetComponent(tempNeighborCellIds
										->GetId(j), 2) != 0)
								) {
								if (containerMaterialCellIDsEndo.find(material2) == containerMaterialCellIDsEndo.end()) {
									vtkSmartPointer<vtkIdList> containerCellIdsEndo = vtkSmartPointer<vtkIdList>::New();
									containerCellIdsEndo->InsertNextId(tempNeighborCellIds->GetId(j));
									pair<int, vtkSmartPointer<vtkIdList>> insert(material2, containerCellIdsEndo);
									containerMaterialCellIDsEndo.insert(insert);
								}
								else {
									vtkSmartPointer<vtkIdList> containerCellIdsEndo = containerMaterialCellIDsEndo.at(material2);
									containerCellIdsEndo->InsertUniqueId(tempNeighborCellIds->GetId(j));
									containerMaterialCellIDsEndo.at(material2) = containerCellIdsEndo;
								}
							}
						}
					}

					viewedCellList.insert(tempNeighborCellIds->GetId(j));
				}
			}
		}
		if (newcellIdList.size() == 0) {
			end = true;
		}
		cellIdList.clear();
		cellIdList = newcellIdList;
	}

	vector<double> centroid = getCentroid(data, cellIds);

	double centroidArray[3] = { centroid.at(0), centroid.at(1), centroid.at(2) };

	for (vtkIdType i = 0; i < cells->GetNumberOfTuples(); i++) {
		double distance = sqrt(abs(vtkMath::Distance2BetweenPoints(centroidArray, centrePoints->GetPoint(cells->vtkDataArray::GetComponent(i, 0)))));
		cells->vtkDataArray::SetComponent(i, 1, distance);
	}

	vtkSortDataArray::SortArrayByComponent(cells, 1);

	unordered_map<int, vtkSmartPointer<vtkUnstructuredGrid>> containerMaterialSearchPath;
	unordered_map<int, vtkSmartPointer<vtkUnstructuredGrid>> containerMaterialSearchPathEndo;

	unordered_map<int, vector<double>> mainOrientation;
	unordered_map<int, vector<double>> mainOrientationEndo;

	if ((containerMaterialCellIDs.size() > 0) || (containerMaterialCellIDsEndo.size() > 0)) {
		for (unordered_map<int, vtkSmartPointer<vtkIdList>>::const_iterator iter = containerMaterialCellIDs.begin();
			iter != containerMaterialCellIDs.end(); iter++) {
			vtkSmartPointer<vtkUnstructuredGrid> SearchList = vtkSmartPointer<vtkUnstructuredGrid>::New();
			vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

			double sumMainOrientationArray[3] = { 0, 0, 0 };

			for (vtkIdType i = 0; i < iter->second->GetNumberOfIds(); i++) {
				points->InsertNextPoint(centrePoints->GetPoint(iter->second->GetId(i)));
				double *orientationArray = data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetTuple(iter->second->GetId(i));
				vtkMath::Add(sumMainOrientationArray, orientationArray, sumMainOrientationArray);
			}
			SearchList->SetPoints(points);

			vtkMath::Normalize(sumMainOrientationArray);
			vector<double> sumMainOrientation = { sumMainOrientationArray[0], sumMainOrientationArray[1], sumMainOrientationArray[2] };

			pair<int, vector<double>> insertMainOrientation(iter->first, sumMainOrientation);
			mainOrientation.insert(insertMainOrientation);

			pair<int, vtkSmartPointer<vtkUnstructuredGrid>> insert(iter->first, SearchList);
			containerMaterialSearchPath.insert(insert);
		}


		for (unordered_map<int, vtkSmartPointer<vtkIdList>>::const_iterator iter = containerMaterialCellIDsEndo.begin();
			iter != containerMaterialCellIDsEndo.end(); iter++) {
			vtkSmartPointer<vtkUnstructuredGrid> SearchList = vtkSmartPointer<vtkUnstructuredGrid>::New();
			vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

			double sumMainOrientationArray[3] = { 0, 0, 0 };

			for (vtkIdType i = 0; i < iter->second->GetNumberOfIds(); i++) {
				points->InsertNextPoint(centrePoints->GetPoint(iter->second->GetId(i)));
				double *orientationArray = data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->GetTuple(iter->second->GetId(i));
				vtkMath::Add(sumMainOrientationArray, orientationArray, sumMainOrientationArray);
			}
			SearchList->SetPoints(points);

			vtkMath::Normalize(sumMainOrientationArray);
			vector<double> sumMainOrientation = { sumMainOrientationArray[0], sumMainOrientationArray[1], sumMainOrientationArray[2] };

			pair<int, vector<double>> insertMainOrientation(iter->first, sumMainOrientation);
			mainOrientationEndo.insert(insertMainOrientation);

			pair<int, vtkSmartPointer<vtkUnstructuredGrid>> insert(iter->first, SearchList);
			containerMaterialSearchPathEndo.insert(insert);
		}


		for (vtkIdType j = cells->GetNumberOfTuples() - 1; j >= 0; j--) {
			if (cells->vtkDataArray::GetComponent(j, 2) == 1) { // Epi
				vtkIdType currentId = cells->vtkDataArray::GetComponent(j, 0);

				double currentPoint[3] = { 0, 0, 0 };
				currentPoint[0] = centrePoints->GetPoint(currentId)[0];
				currentPoint[1] = centrePoints->GetPoint(currentId)[1];
				currentPoint[2] = centrePoints->GetPoint(currentId)[2];

				list<tuple<int, double, vector<double>, vector<double>>> differenceVectorList;

				unordered_map<int, vector<double>>::const_iterator iterMainOriatnation = mainOrientation.begin();
				for (unordered_map<int, vtkSmartPointer<vtkUnstructuredGrid>>::const_iterator iter = containerMaterialSearchPath.begin();
					iter != containerMaterialSearchPath.end(); iter++) {
					vtkIdType closestPosPath = iter->second->FindPoint(currentPoint);
					vtkIdType closestIdPoint = containerMaterialCellIDs.find(iter->first)->second->GetId(closestPosPath);

					double closestPoint[3] = { 0, 0, 0 };
					closestPoint[0] = centrePoints->GetPoint(closestIdPoint)[0];
					closestPoint[1] = centrePoints->GetPoint(closestIdPoint)[1];
					closestPoint[2] = centrePoints->GetPoint(closestIdPoint)[2];

					double distance = sqrt(abs(vtkMath::Distance2BetweenPoints(currentPoint, closestPoint)));
					double *materialDiffArray = data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetTuple(closestIdPoint);

					vector<double> materialDiffVector = { materialDiffArray[0], materialDiffArray[1], materialDiffArray[2] };

					tuple<int, double, vector<double>, vector<double>> tupleDiffVector(iter->first, distance, materialDiffVector, iterMainOriatnation->second);
					iterMainOriatnation++;

					differenceVectorList.push_back(tupleDiffVector);
				}

				unordered_map<int, vector<double>>::const_iterator iterMainOriatnationEndo = mainOrientationEndo.begin();
				for (unordered_map<int, vtkSmartPointer<vtkUnstructuredGrid>>::const_iterator iter = containerMaterialSearchPathEndo.begin();
					iter != containerMaterialSearchPathEndo.end(); iter++) {
					vtkIdType closestPosPath = iter->second->FindPoint(currentPoint);
					vtkIdType closestIdPoint = containerMaterialCellIDsEndo.find(iter->first)->second->GetId(closestPosPath);


					double distance = calcDistancetoPointinEndo(data, currentId, closestIdPoint, distanceEndoEpi);
					double *materialDiffArray = data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->GetTuple(closestIdPoint);

					vector<double> materialDiffVector = { materialDiffArray[0], materialDiffArray[1], materialDiffArray[2] };
					tuple<int, double, vector<double>, vector<double>> tupleDiffVector(iter->first, distance, materialDiffVector, iterMainOriatnationEndo->second);
					iterMainOriatnationEndo++;
					differenceVectorList.push_back(tupleDiffVector);
				}

				differenceVectorList.unique(equal_differenceVectorList);

				int newMaterial = -1;
				double minDistance = std::numeric_limits<double>::max();
				double sumOfDistance = 0;

				for (list<tuple<int, double, vector<double>, vector<double>>>::iterator iter = differenceVectorList.begin();
					iter != differenceVectorList.end(); iter++) {
					if (minDistance > get<1>(*iter)) {
						minDistance = get<1>(*iter);
						newMaterial = get<0>(*iter);
					}
					if (get<1>(*iter) > 0) {
						get<1>(*iter) = 1 / get<1>(*iter);
					}
					else {
						get<1>(*iter) = 1;
					}
					sumOfDistance += get<1>(*iter);
				}


				for (list<tuple<int, double, vector<double>, vector<double>>>::iterator iter = differenceVectorList.begin();
					iter != differenceVectorList.end(); iter++) {
					get<1>(*iter) = get<1>(*iter) / sumOfDistance;
				}

				vector<double> newOrientationVector = AveragingOrientation::AveragingNVectors(differenceVectorList);
				double newDiffVector[3] = { newOrientationVector.at(0), newOrientationVector.at(1), newOrientationVector.at(2) };

				data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->InsertTuple(currentId, newDiffVector);
				data.getVtkData()->GetCellData()->GetArray("Status")->SetComponent(currentId, 0, 3);
				data.getVtkData()->GetCellData()->GetArray("Material")->SetComponent(currentId, 0, newMaterial);
			}

			if (cells->vtkDataArray::GetComponent(j, 3) == 1) { // Endo
				vtkIdType currentId = cells->vtkDataArray::GetComponent(j, 0);

				double currentPoint[3] = { 0, 0, 0 };
				currentPoint[0] = centrePoints->GetPoint(currentId)[0];
				currentPoint[1] = centrePoints->GetPoint(currentId)[1];
				currentPoint[2] = centrePoints->GetPoint(currentId)[2];

				list<tuple<int, double, vector<double>, vector<double>>> differenceVectorList;

				unordered_map<int, vector<double>>::const_iterator iterMainOriatnation = mainOrientation.begin();
				for (unordered_map<int, vtkSmartPointer<vtkUnstructuredGrid>>::const_iterator iter = containerMaterialSearchPath.begin();
					iter != containerMaterialSearchPath.end(); iter++) {
					vtkIdType closestPosPath = iter->second->FindPoint(currentPoint);
					vtkIdType closestIdPoint = containerMaterialCellIDs.find(iter->first)->second->GetId(closestPosPath);

					double distance = calcDistancetoPointinEpi(data, currentId, closestIdPoint, distanceEndoEpi);
					double *materialDiffArray = data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetTuple(closestIdPoint);

					vector<double> materialDiffVector = { materialDiffArray[0], materialDiffArray[1], materialDiffArray[2] };
					tuple<int, double, vector<double>, vector<double>> tupleDiffVector(iter->first, distance, materialDiffVector, iterMainOriatnation->second);
					differenceVectorList.push_back(tupleDiffVector);
				}

				unordered_map<int, vector<double>>::const_iterator iterMainOriatnationEndo = mainOrientationEndo.begin();
				for (unordered_map<int, vtkSmartPointer<vtkUnstructuredGrid>>::const_iterator iter = containerMaterialSearchPathEndo.begin();
					iter != containerMaterialSearchPathEndo.end(); iter++) {
					vtkIdType closestPosPath = iter->second->FindPoint(currentPoint);
					vtkIdType closestIdPoint = containerMaterialCellIDsEndo.find(iter->first)->second->GetId(closestPosPath);

					double closestPoint[3] = { 0, 0, 0 };
					closestPoint[0] = centrePoints->GetPoint(closestIdPoint)[0];
					closestPoint[1] = centrePoints->GetPoint(closestIdPoint)[1];
					closestPoint[2] = centrePoints->GetPoint(closestIdPoint)[2];

					double distance = sqrt(abs(vtkMath::Distance2BetweenPoints(currentPoint, closestPoint)));
					double *materialDiffArray = data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->GetTuple(closestIdPoint);

					vector<double> materialDiffVector = { materialDiffArray[0], materialDiffArray[1], materialDiffArray[2] };
					tuple<int, double, vector<double>, vector<double>> tupleDiffVector(iter->first, distance, materialDiffVector, iterMainOriatnationEndo->second);
					iterMainOriatnationEndo++;

					differenceVectorList.push_back(tupleDiffVector);
				}

				int newMaterial = -1;
				double minDistance = std::numeric_limits<double>::max();
				double sumOfDistance = 0;

				differenceVectorList.unique(equal_differenceVectorList);

				for (list<tuple<int, double, vector<double>, vector<double>>>::iterator iter = differenceVectorList.begin();
					iter != differenceVectorList.end(); iter++) {
					if (minDistance > get<1>(*iter)) {
						minDistance = get<1>(*iter);
						newMaterial = get<0>(*iter);
					}

					if (get<1>(*iter) > 0) {
						get<1>(*iter) = 1 / get<1>(*iter);
					}
					else {
						get<1>(*iter) = 1;
					}
					sumOfDistance += get<1>(*iter);
				}


				for (list<tuple<int, double, vector<double>, vector<double>>>::iterator iter = differenceVectorList.begin();
					iter != differenceVectorList.end(); iter++) {
					get<1>(*iter) = get<1>(*iter) / sumOfDistance;
				}

				vector<double> newOrientationVector = AveragingOrientation::AveragingNVectors(differenceVectorList);
				double newDiffVector[3] = { newOrientationVector.at(0), newOrientationVector.at(1), newOrientationVector.at(2) };

				data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->InsertTuple(currentId, newDiffVector);
				data.getVtkData()->GetCellData()->GetArray("StatusEndo")->SetComponent(currentId, 0, 3);
				data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->SetComponent(currentId, 0, newMaterial);
			}
		}
	}
	else {
		for (vtkIdType i = 0; i < cells->GetNumberOfTuples(); i++) {
			vtkIdType currentId = cells->vtkDataArray::GetComponent(i, 0);
			vtkSmartPointer<vtkIdList> cellNeighbours = getCellNeighbours(data, currentId);

			if (cells->vtkDataArray::GetComponent(i, 2) == 1) { // Epi
				int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(cellNeighbours->GetId(0), 0);
				data.getVtkData()->GetCellData()->GetArray("Material")->SetComponent(currentId, 0, material);
			}


			if (cells->vtkDataArray::GetComponent(i, 3) == 1) { // Endo
				int material = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(cellNeighbours->GetId(0), 0);
				data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->SetComponent(currentId, 0, material);
			}
		}
	}

	return cellIds;
} // Methods::regionGrowOneSurface

/*! Close the fiber orientation in a the bridge at a surface mesh.
 \param data Pointer to the orginal mesh.
 \param seedPoint Seed point for the interpolation region.
 \param distanceEndoEpi Fictive distance between the endocard and epicard in the left atrial.
 \param freeBridgeMaterial Tissue class of the free bridge material.
 \retrun A vtkIdList with the cell-ids of the interpolate region.
 */
vtkSmartPointer<vtkIdList> Methods::regionGrowOneSurfaceFreeBridge(DataFormat &data, std::vector<double> seedPoint, double distanceEndoEpi, int freeBridgeMaterial) {
	double seedPointArray[3] = { seedPoint.at(0), seedPoint.at(1), seedPoint.at(2) };

	vtkSmartPointer<vtkPoints> centrePoints = data.getCentrePoints()->GetPoints();

	vtkIdType seedPointID = data.getCentrePoints()->FindPoint(seedPointArray);

	vtkSmartPointer<vtkDoubleArray> cells = vtkSmartPointer<vtkDoubleArray>::New();

	cells->SetNumberOfComponents(4);
	cells->SetComponentName(0, "CellID");
	cells->SetComponentName(1, "distance to centre point");
	cells->SetComponentName(2, "Epi");
	cells->SetComponentName(3, "Endo");

	set<vtkIdType> viewedCellList;

	bool end = false;

	vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
	set<vtkIdType> cellIdList;

	cellIdList.insert(seedPointID);

	if ((data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(seedPointID, 0) == 0) ||
		(data.getVtkData()->GetCellData()->GetArray("StatusEndo")->GetComponent(seedPointID, 0) == 0)) {
		double cellsArray[4] = { 0, 0, 0, 0 };
		cellsArray[0] = seedPointID;
		if (data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(seedPointID, 0) == 0) {
			cellsArray[2] = 1;
		}
		if (data.getVtkData()->GetCellData()->GetArray("StatusEndo")->GetComponent(seedPointID, 0) == 0) {
			cellsArray[3] = 1;
		}
		cells->InsertNextTuple(cellsArray);
		cellIds->InsertNextId(seedPointID);
	}

	viewedCellList.insert(seedPointID);
	unordered_map<int, vtkSmartPointer<vtkIdList>> containerMaterialCellIDs;
	unordered_map<int, vtkSmartPointer<vtkIdList>> containerMaterialCellIDsEndo;


	while (!end) {
		set<vtkIdType> newcellIdList;

		for (set<vtkIdType>::iterator iter = cellIdList.begin(); iter != cellIdList.end(); iter++) {
			vtkIdType currentCellIds = *iter;
			vtkSmartPointer<vtkIdList> cellPointIds = vtkSmartPointer<vtkIdList>::New();

			data.getVtkData()->GetCellPoints(currentCellIds, cellPointIds);

			for (vtkIdType i = 0; i < cellPointIds->GetNumberOfIds(); i++) {
				vtkSmartPointer<vtkIdList> CellPointIdList = vtkSmartPointer<vtkIdList>::New();
				CellPointIdList->InsertNextId(cellPointIds->GetId(i));

				vtkSmartPointer<vtkIdList> tempNeighborCellIds = vtkSmartPointer<vtkIdList>::New();
				data.getVtkData()->GetCellNeighbors(currentCellIds, CellPointIdList, tempNeighborCellIds);

				for (vtkIdType j = 0; j < tempNeighborCellIds->GetNumberOfIds(); j++) {
					if (viewedCellList.find(tempNeighborCellIds->GetId(j)) == viewedCellList.end()) {
						int material1 = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(tempNeighborCellIds->GetId(j), 0);
						int material2 = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(tempNeighborCellIds->GetId(j), 0);

						if (Material::isInLeft(material1) || Material::isInRight(material1) || Material::isInLeft(material2) ||
							Material::isInRight(material2) || (material1 == freeBridgeMaterial) ||
							(material2 == freeBridgeMaterial)) {
							if ((data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(tempNeighborCellIds->GetId(j), 0) == 0) ||
								(data.getVtkData()->GetCellData()->GetArray("StatusEndo")->GetComponent(tempNeighborCellIds->GetId(j), 0) == 0)) {
								double cellsArray[4] = { 0, 0, 0, 0 };
								cellsArray[0] = tempNeighborCellIds->GetId(j);
								if (data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(tempNeighborCellIds->GetId(j), 0) == 0) {
									cellsArray[2] = 1;
								}
								if (data.getVtkData()->GetCellData()->GetArray("StatusEndo")->GetComponent(tempNeighborCellIds->GetId(j), 0) == 0) {
									cellsArray[3] = 1;
								}
								cells->InsertNextTuple(cellsArray);
								newcellIdList.insert(tempNeighborCellIds->GetId(j));
								cellIds->InsertNextId(tempNeighborCellIds->GetId(j));
							}


							if ((material1 != Material::Vorhof_links) &&
								(material1 != Material::Vorhof_links_Endo) &&
								(material1 != Material::Vorhof_rechts) &&
								(data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(tempNeighborCellIds->GetId(j), 0) != 0) &&
								(data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetComponent(tempNeighborCellIds->
									GetId(j), 0) +
									data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetComponent(tempNeighborCellIds->
										GetId(j), 1) +
									data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetComponent(tempNeighborCellIds->
										GetId(j), 2) != 0)
								) {
								if (containerMaterialCellIDs.find(material1) == containerMaterialCellIDs.end()) {
									vtkSmartPointer<vtkIdList> containerCellIds = vtkSmartPointer<vtkIdList>::New();
									containerCellIds->InsertNextId(tempNeighborCellIds->GetId(j));
									pair<int, vtkSmartPointer<vtkIdList>> insert(material1, containerCellIds);
									containerMaterialCellIDs.insert(insert);
								}
								else {
									vtkSmartPointer<vtkIdList> containerCellIds = containerMaterialCellIDs.at(material1);
									containerCellIds->InsertUniqueId(tempNeighborCellIds->GetId(j));
									containerMaterialCellIDs.at(material1) = containerCellIds;
								}
							}
							if ((material2 != Material::Vorhof_links) &&
								(material2 != Material::Vorhof_links_Endo) &&
								(material2 != Material::Vorhof_rechts) &&
								(data.getVtkData()->GetCellData()->GetArray("StatusEndo")->GetComponent(tempNeighborCellIds->GetId(j), 0) != 0) &&
								(data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->GetComponent(tempNeighborCellIds
									->GetId(j), 0) +
									data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->GetComponent(tempNeighborCellIds
										->GetId(j), 1) +
									data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->GetComponent(tempNeighborCellIds
										->GetId(j), 2) != 0)
								) {
								if (containerMaterialCellIDsEndo.find(material2) == containerMaterialCellIDsEndo.end()) {
									vtkSmartPointer<vtkIdList> containerCellIdsEndo = vtkSmartPointer<vtkIdList>::New();
									containerCellIdsEndo->InsertNextId(tempNeighborCellIds->GetId(j));
									pair<int, vtkSmartPointer<vtkIdList>> insert(material2, containerCellIdsEndo);
									containerMaterialCellIDsEndo.insert(insert);
								}
								else {
									vtkSmartPointer<vtkIdList> containerCellIdsEndo = containerMaterialCellIDsEndo.at(material2);
									containerCellIdsEndo->InsertUniqueId(tempNeighborCellIds->GetId(j));
									containerMaterialCellIDsEndo.at(material2) = containerCellIdsEndo;
								}
							}
						}
					}

					viewedCellList.insert(tempNeighborCellIds->GetId(j));
				}
			}
		}
		if (newcellIdList.size() == 0) {
			end = true;
		}
		cellIdList.clear();
		cellIdList = newcellIdList;
	}

	vector<double> centroid = getCentroid(data, cellIds);

	double centroidArray[3] = { centroid.at(0), centroid.at(1), centroid.at(2) };

	for (vtkIdType i = 0; i < cells->GetNumberOfTuples(); i++) {
		double distance = sqrt(abs(vtkMath::Distance2BetweenPoints(centroidArray, centrePoints->GetPoint(cells->vtkDataArray::GetComponent(i, 0)))));
		cells->vtkDataArray::SetComponent(i, 1, distance);
	}

	vtkSortDataArray::SortArrayByComponent(cells, 1);

	unordered_map<int, vtkSmartPointer<vtkUnstructuredGrid>> containerMaterialSearchPath;
	unordered_map<int, vtkSmartPointer<vtkUnstructuredGrid>> containerMaterialSearchPathEndo;

	unordered_map<int, vector<double>> mainOrientation;
	unordered_map<int, vector<double>> mainOrientationEndo;

	if ((containerMaterialCellIDs.size() > 0) || (containerMaterialCellIDsEndo.size() > 0)) {
		for (unordered_map<int, vtkSmartPointer<vtkIdList>>::const_iterator iter = containerMaterialCellIDs.begin();
			iter != containerMaterialCellIDs.end(); iter++) {
			vtkSmartPointer<vtkUnstructuredGrid> SearchList = vtkSmartPointer<vtkUnstructuredGrid>::New();
			vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

			double sumMainOrientationArray[3] = { 0, 0, 0 };

			for (vtkIdType i = 0; i < iter->second->GetNumberOfIds(); i++) {
				points->InsertNextPoint(centrePoints->GetPoint(iter->second->GetId(i)));
				double *orientationArray = data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetTuple(iter->second->GetId(i));
				vtkMath::Add(sumMainOrientationArray, orientationArray, sumMainOrientationArray);
			}
			SearchList->SetPoints(points);

			vtkMath::Normalize(sumMainOrientationArray);
			vector<double> sumMainOrientation = { sumMainOrientationArray[0], sumMainOrientationArray[1], sumMainOrientationArray[2] };

			pair<int, vector<double>> insertMainOrientation(iter->first, sumMainOrientation);
			mainOrientation.insert(insertMainOrientation);

			pair<int, vtkSmartPointer<vtkUnstructuredGrid>> insert(iter->first, SearchList);
			containerMaterialSearchPath.insert(insert);
		}


		for (unordered_map<int, vtkSmartPointer<vtkIdList>>::const_iterator iter = containerMaterialCellIDsEndo.begin();
			iter != containerMaterialCellIDsEndo.end(); iter++) {
			vtkSmartPointer<vtkUnstructuredGrid> SearchList = vtkSmartPointer<vtkUnstructuredGrid>::New();
			vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

			double sumMainOrientationArray[3] = { 0, 0, 0 };

			for (vtkIdType i = 0; i < iter->second->GetNumberOfIds(); i++) {
				points->InsertNextPoint(centrePoints->GetPoint(iter->second->GetId(i)));
				double *orientationArray = data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->GetTuple(iter->second->GetId(i));
				vtkMath::Add(sumMainOrientationArray, orientationArray, sumMainOrientationArray);
			}
			SearchList->SetPoints(points);

			vtkMath::Normalize(sumMainOrientationArray);
			vector<double> sumMainOrientation = { sumMainOrientationArray[0], sumMainOrientationArray[1], sumMainOrientationArray[2] };

			pair<int, vector<double>> insertMainOrientation(iter->first, sumMainOrientation);
			mainOrientationEndo.insert(insertMainOrientation);

			pair<int, vtkSmartPointer<vtkUnstructuredGrid>> insert(iter->first, SearchList);
			containerMaterialSearchPathEndo.insert(insert);
		}


		for (vtkIdType j = cells->GetNumberOfTuples() - 1; j >= 0; j--) {
			if (cells->vtkDataArray::GetComponent(j, 2) == 1) { // Epi
				vtkIdType currentId = cells->vtkDataArray::GetComponent(j, 0);

				double currentPoint[3] = { 0, 0, 0 };
				currentPoint[0] = centrePoints->GetPoint(currentId)[0];
				currentPoint[1] = centrePoints->GetPoint(currentId)[1];
				currentPoint[2] = centrePoints->GetPoint(currentId)[2];

				list<tuple<int, double, vector<double>, vector<double>>> differenceVectorList;

				unordered_map<int, vector<double>>::const_iterator iterMainOriatnation = mainOrientation.begin();
				for (unordered_map<int, vtkSmartPointer<vtkUnstructuredGrid>>::const_iterator iter = containerMaterialSearchPath.begin();
					iter != containerMaterialSearchPath.end(); iter++) {
					vtkIdType closestPosPath = iter->second->FindPoint(currentPoint);
					vtkIdType closestIdPoint = containerMaterialCellIDs.find(iter->first)->second->GetId(closestPosPath);

					double closestPoint[3] = { 0, 0, 0 };
					closestPoint[0] = centrePoints->GetPoint(closestIdPoint)[0];
					closestPoint[1] = centrePoints->GetPoint(closestIdPoint)[1];
					closestPoint[2] = centrePoints->GetPoint(closestIdPoint)[2];

					double distance = sqrt(abs(vtkMath::Distance2BetweenPoints(currentPoint, closestPoint)));
					double *materialDiffArray = data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetTuple(closestIdPoint);

					vector<double> materialDiffVector = { materialDiffArray[0], materialDiffArray[1], materialDiffArray[2] };

					tuple<int, double, vector<double>, vector<double>> tupleDiffVector(iter->first, distance, materialDiffVector, iterMainOriatnation->second);
					iterMainOriatnation++;

					differenceVectorList.push_back(tupleDiffVector);
				}

				unordered_map<int, vector<double>>::const_iterator iterMainOriatnationEndo = mainOrientationEndo.begin();
				for (unordered_map<int, vtkSmartPointer<vtkUnstructuredGrid>>::const_iterator iter = containerMaterialSearchPathEndo.begin();
					iter != containerMaterialSearchPathEndo.end(); iter++) {
					vtkIdType closestPosPath = iter->second->FindPoint(currentPoint);
					vtkIdType closestIdPoint = containerMaterialCellIDsEndo.find(iter->first)->second->GetId(closestPosPath);


					double distance = calcDistancetoPointinEndo(data, currentId, closestIdPoint, distanceEndoEpi);
					double *materialDiffArray = data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->GetTuple(closestIdPoint);

					vector<double> materialDiffVector = { materialDiffArray[0], materialDiffArray[1], materialDiffArray[2] };
					tuple<int, double, vector<double>, vector<double>> tupleDiffVector(iter->first, distance, materialDiffVector, iterMainOriatnationEndo->second);
					iterMainOriatnationEndo++;
					differenceVectorList.push_back(tupleDiffVector);
				}

				differenceVectorList.unique(equal_differenceVectorList);

				int newMaterial = -1;
				double minDistance = std::numeric_limits<double>::max();
				double sumOfDistance = 0;

				for (list<tuple<int, double, vector<double>, vector<double>>>::iterator iter = differenceVectorList.begin();
					iter != differenceVectorList.end(); iter++) {
					if (minDistance > get<1>(*iter)) {
						minDistance = get<1>(*iter);
						newMaterial = get<0>(*iter);
					}
					if (get<1>(*iter) > 0) {
						get<1>(*iter) = 1 / get<1>(*iter);
					}
					else {
						get<1>(*iter) = 1;
					}
					sumOfDistance += get<1>(*iter);
				}

				for (list<tuple<int, double, vector<double>, vector<double>>>::iterator iter = differenceVectorList.begin();
					iter != differenceVectorList.end(); iter++) {
					get<1>(*iter) = get<1>(*iter) / sumOfDistance;
				}

				vector<double> newOrientationVector = AveragingOrientation::AveragingNVectors(differenceVectorList);
				double newDiffVector[3] = { newOrientationVector.at(0), newOrientationVector.at(1), newOrientationVector.at(2) };

				data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->InsertTuple(currentId, newDiffVector);
				data.getVtkData()->GetCellData()->GetArray("Status")->SetComponent(currentId, 0, 3);
				data.getVtkData()->GetCellData()->GetArray("Material")->SetComponent(currentId, 0, newMaterial);
			}

			if (cells->vtkDataArray::GetComponent(j, 3) == 1) { // Endo
				vtkIdType currentId = cells->vtkDataArray::GetComponent(j, 0);

				double currentPoint[3] = { 0, 0, 0 };
				currentPoint[0] = centrePoints->GetPoint(currentId)[0];
				currentPoint[1] = centrePoints->GetPoint(currentId)[1];
				currentPoint[2] = centrePoints->GetPoint(currentId)[2];

				list<tuple<int, double, vector<double>, vector<double>>> differenceVectorList;

				unordered_map<int, vector<double>>::const_iterator iterMainOriatnation = mainOrientation.begin();
				for (unordered_map<int, vtkSmartPointer<vtkUnstructuredGrid>>::const_iterator iter = containerMaterialSearchPath.begin();
					iter != containerMaterialSearchPath.end(); iter++) {
					vtkIdType closestPosPath = iter->second->FindPoint(currentPoint);
					vtkIdType closestIdPoint = containerMaterialCellIDs.find(iter->first)->second->GetId(closestPosPath);

					double distance = calcDistancetoPointinEpi(data, currentId, closestIdPoint, distanceEndoEpi);
					double *materialDiffArray = data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetTuple(closestIdPoint);

					vector<double> materialDiffVector = { materialDiffArray[0], materialDiffArray[1], materialDiffArray[2] };
					tuple<int, double, vector<double>, vector<double>> tupleDiffVector(iter->first, distance, materialDiffVector, iterMainOriatnation->second);
					differenceVectorList.push_back(tupleDiffVector);
				}

				unordered_map<int, vector<double>>::const_iterator iterMainOriatnationEndo = mainOrientationEndo.begin();
				for (unordered_map<int, vtkSmartPointer<vtkUnstructuredGrid>>::const_iterator iter = containerMaterialSearchPathEndo.begin();
					iter != containerMaterialSearchPathEndo.end(); iter++) {
					vtkIdType closestPosPath = iter->second->FindPoint(currentPoint);
					vtkIdType closestIdPoint = containerMaterialCellIDsEndo.find(iter->first)->second->GetId(closestPosPath);

					double closestPoint[3] = { 0, 0, 0 };
					closestPoint[0] = centrePoints->GetPoint(closestIdPoint)[0];
					closestPoint[1] = centrePoints->GetPoint(closestIdPoint)[1];
					closestPoint[2] = centrePoints->GetPoint(closestIdPoint)[2];

					double distance = sqrt(abs(vtkMath::Distance2BetweenPoints(currentPoint, closestPoint)));
					double *materialDiffArray = data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->GetTuple(closestIdPoint);

					vector<double> materialDiffVector = { materialDiffArray[0], materialDiffArray[1], materialDiffArray[2] };
					tuple<int, double, vector<double>, vector<double>> tupleDiffVector(iter->first, distance, materialDiffVector, iterMainOriatnationEndo->second);
					iterMainOriatnationEndo++;

					differenceVectorList.push_back(tupleDiffVector);
				}

				int newMaterial = -1;
				double minDistance = std::numeric_limits<double>::max();
				double sumOfDistance = 0;

				differenceVectorList.unique(equal_differenceVectorList);

				for (list<tuple<int, double, vector<double>, vector<double>>>::iterator iter = differenceVectorList.begin();
					iter != differenceVectorList.end(); iter++) {
					if (minDistance > get<1>(*iter)) {
						minDistance = get<1>(*iter);
						newMaterial = get<0>(*iter);
					}

					if (get<1>(*iter) > 0) {
						get<1>(*iter) = 1 / get<1>(*iter);
					}
					else {
						get<1>(*iter) = 1;
					}
					sumOfDistance += get<1>(*iter);
				}


				for (list<tuple<int, double, vector<double>, vector<double>>>::iterator iter = differenceVectorList.begin();
					iter != differenceVectorList.end(); iter++) {
					get<1>(*iter) = get<1>(*iter) / sumOfDistance;
				}

				vector<double> newOrientationVector = AveragingOrientation::AveragingNVectors(differenceVectorList);
				double newDiffVector[3] = { newOrientationVector.at(0), newOrientationVector.at(1), newOrientationVector.at(2) };

				data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->InsertTuple(currentId, newDiffVector);
				data.getVtkData()->GetCellData()->GetArray("StatusEndo")->SetComponent(currentId, 0, 3);
				data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->SetComponent(currentId, 0, newMaterial);
			}
		}
	}
	else {
		for (vtkIdType i = 0; i < cells->GetNumberOfTuples(); i++) {
			vtkIdType currentId = cells->vtkDataArray::GetComponent(i, 0);
			vtkSmartPointer<vtkIdList> cellNeighbours = getCellNeighbours(data, currentId);

			if (cells->vtkDataArray::GetComponent(i, 2) == 1) { // Epi
				int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(cellNeighbours->GetId(0), 0);
				data.getVtkData()->GetCellData()->GetArray("Material")->SetComponent(currentId, 0, material);
			}


			if (cells->vtkDataArray::GetComponent(i, 3) == 1) { // Endo
				int material = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(cellNeighbours->GetId(0), 0);
				data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->SetComponent(currentId, 0, material);
			}
		}
	}

	return cellIds;
} // Methods::regionGrowOneSurfaceFreeBridge

/*! Interpolated the fiber orientation in the region between two paths in a tissue class.
 \param data Pointer to the orginal mesh.
 \param path1 First path with fiber orientation.
 \param path2 Second path with fiber orientation.
 \param inMaterial Tissue class in the fieber orientation could be interpolated.
 \param seedPoint Seed point for the interpolation region.
 \retrun A vtkIdList with the cell-ids of the interpolate region.
 */
vtkSmartPointer<vtkIdList> Methods::regionInterpolateBetween2PathwithSeed(DataFormat &data, vtkSmartPointer<vtkIdList> path1, vtkSmartPointer<vtkIdList> path2, Material::Mat inMaterial, vector<double> seedPoint) {
	return regionInterpolateBetween2PathwithSeed(data, path1, path2, inMaterial, inMaterial, seedPoint);
}

/*! Interpolated the fiber orientation in the region between two paths in two tissue class.
 \param data Pointer to the orginal mesh.
 \param path1 First path with fiber orientation.
 \param path2 Second path with fiber orientation.
 \param inMaterial1 Tissue class in the fieber orientation could be interpolated.
 \param inMaterial2 Tissue class in the fieber orientation could be interpolated.
 \param seedPoint Seed point for the interpolation region.
 \retrun A vtkIdList with the cell-ids of the interpolate region.
 */
vtkSmartPointer<vtkIdList> Methods::regionInterpolateBetween2PathwithSeed(DataFormat &data, vtkSmartPointer<vtkIdList> path1, vtkSmartPointer<vtkIdList> path2, Material::Mat inMaterial1, Material::Mat inMaterial2, vector<double> seedPoint) {
	double seedPointArray[3] = { seedPoint.at(0), seedPoint.at(1), seedPoint.at(2) };

	vtkSmartPointer<vtkPoints> centrePoints = data.getCentrePoints()->GetPoints();
	vtkSmartPointer<vtkUnstructuredGrid> centrePointsUGrid = data.getCentrePoints();

	vtkSmartPointer<vtkUnstructuredGrid> path1SearchList = vtkSmartPointer<vtkUnstructuredGrid>::New();
	vtkSmartPointer<vtkPoints> pointsPath1 = vtkSmartPointer<vtkPoints>::New();

	for (vtkIdType i = 0; i < path1->GetNumberOfIds(); i++) {
		double point[3] = { centrePoints->GetPoint(path1->GetId(i))[0], centrePoints->GetPoint(path1->GetId(i))[1], centrePoints->GetPoint(path1->GetId(i))
		[2] };
		pointsPath1->InsertNextPoint(point);
	}
	path1SearchList->SetPoints(pointsPath1);

	vtkSmartPointer<vtkUnstructuredGrid> path2SearchList = vtkSmartPointer<vtkUnstructuredGrid>::New();
	vtkSmartPointer<vtkPoints> pointsPath2 = vtkSmartPointer<vtkPoints>::New();

	for (vtkIdType i = 0; i < path2->GetNumberOfIds(); i++) {
		double point[3] = { centrePoints->GetPoint(path2->GetId(i))[0], centrePoints->GetPoint(path2->GetId(i))[1], centrePoints->GetPoint(path2->GetId(i))
		[2] };
		pointsPath2->InsertNextPoint(point);
	}
	path2SearchList->SetPoints(pointsPath2);

	vtkIdType seedPointID = centrePointsUGrid->FindPoint(seedPointArray);

	vtkSmartPointer<vtkIdList> cellIdList = vtkSmartPointer<vtkIdList>::New();
	cellIdList->InsertNextId(seedPointID);

	set<vtkIdType> returnList;
	set<vtkIdType> viewedCells;

	do {
		vtkSmartPointer<vtkIdList> newcellIdList = vtkSmartPointer<vtkIdList>::New();

		for (vtkIdType k = 0; k < cellIdList->GetNumberOfIds(); k++) {
			vtkIdType currentCellIds = cellIdList->GetId(k);
			vtkSmartPointer<vtkIdList> cellPointIds = vtkSmartPointer<vtkIdList>::New();

			data.getVtkData()->GetCellPoints(currentCellIds, cellPointIds);

			for (vtkIdType i = 0; i < cellPointIds->GetNumberOfIds(); i++) {
				vtkSmartPointer<vtkIdList> CellPointIdList = vtkSmartPointer<vtkIdList>::New();
				CellPointIdList->InsertNextId(cellPointIds->GetId(i));

				vtkSmartPointer<vtkIdList> neighbors = vtkSmartPointer<vtkIdList>::New();
				data.getVtkData()->GetCellNeighbors(currentCellIds, CellPointIdList, neighbors);

				for (vtkIdType j = 0; j < neighbors->GetNumberOfIds(); j++) {
					vtkIdType currentId = neighbors->GetId(j);
					if (!data.getDoubleLayer()) {
						int material = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(currentId, 0);
						if ((viewedCells.find(currentId) == viewedCells.end()) &&
							((material == inMaterial1) || (material == inMaterial2))) {
							newcellIdList->InsertNextId(currentId);
						}
					}

					int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(currentId, 0);
					if ((viewedCells.find(currentId) == viewedCells.end()) &&
						((material == inMaterial1) || (material == inMaterial2))) {
						newcellIdList->InsertNextId(currentId);
					}
					viewedCells.insert(currentId);
				}
			}
		}

		for (vtkIdType j = 0; j < newcellIdList->GetNumberOfIds(); j++) {
			vtkIdType currentId = newcellIdList->GetId(j);

			returnList.insert(currentId);


			int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(currentId, 0);
			if (((material == inMaterial1) || (material == inMaterial2)) &&
				(data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(currentId, 0) == 0)) {
				double currentPoint[3] = { 0, 0, 0 };
				currentPoint[0] = centrePoints->GetPoint(currentId)[0];
				currentPoint[1] = centrePoints->GetPoint(currentId)[1];
				currentPoint[2] = centrePoints->GetPoint(currentId)[2];

				vtkIdType closestPosPath1 = path1SearchList->FindPoint(currentPoint);
				vtkIdType closestIdPath1 = path1->GetId(closestPosPath1);

				double closestPointPath1[3] = { 0, 0, 0 };
				closestPointPath1[0] = centrePointsUGrid->GetPoint(closestIdPath1)[0];
				closestPointPath1[1] = centrePointsUGrid->GetPoint(closestIdPath1)[1];
				closestPointPath1[2] = centrePointsUGrid->GetPoint(closestIdPath1)[2];

				double distancePath1 = sqrt(abs(vtkMath::Distance2BetweenPoints(currentPoint, closestPointPath1)));

				vtkIdType closestPosPath2 = path2SearchList->FindPoint(currentPoint);
				vtkIdType closestIdPath2 = path2->GetId(closestPosPath2);

				double closestPointPath2[3] = { 0, 0, 0 };
				closestPointPath2[0] = centrePointsUGrid->GetPoint(closestIdPath2)[0];
				closestPointPath2[1] = centrePointsUGrid->GetPoint(closestIdPath2)[1];
				closestPointPath2[2] = centrePointsUGrid->GetPoint(closestIdPath2)[2];

				double distancePath2 = sqrt(abs(vtkMath::Distance2BetweenPoints(currentPoint, closestPointPath2)));

				vector<double> vectorPath1(3);
				vectorPath1.at(0) = data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath1)[0];
				vectorPath1.at(1) = data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath1)[1];
				vectorPath1.at(2) = data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath1)[2];

				vector<double> vectorPath2(3);
				vectorPath2.at(0) = data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath2)[0];
				vectorPath2.at(1) = data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath2)[1];
				vectorPath2.at(2) = data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath2)[2];

				vector<double> newOrientationVector = AveragingOrientation::averageOrientation(vectorPath1, vectorPath2, distancePath1, distancePath2);

				double newDiffVector[3] = { newOrientationVector.at(0), newOrientationVector.at(1), newOrientationVector.at(2) };
				vtkMath::Normalize(newDiffVector);

				data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->InsertTuple(currentId, newDiffVector);
				data.getVtkData()->GetCellData()->GetArray("Status")->SetComponent(currentId, 0, 3);

				if (!data.getDoubleLayer() && ((material == inMaterial1) || (material == inMaterial2)) &&
					(data.getVtkData()->GetCellData()->GetArray("StatusEndo")->GetComponent(currentId, 0) == 0)) {
					data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->InsertTuple(currentId, newDiffVector);
					data.getVtkData()->GetCellData()->GetArray("StatusEndo")->SetComponent(currentId, 0, 3);
				}
			}

			if (!data.getDoubleLayer()) {
				int material = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(currentId, 0);
				if (((material == inMaterial1) || (material == inMaterial2)) &&
					(data.getVtkData()->GetCellData()->GetArray("StatusEndo")->GetComponent(currentId, 0) == 0)) {
					double currentPoint[3] = { 0, 0, 0 };
					currentPoint[0] = centrePoints->GetPoint(currentId)[0];
					currentPoint[1] = centrePoints->GetPoint(currentId)[1];
					currentPoint[2] = centrePoints->GetPoint(currentId)[2];

					vtkIdType closestPosPath1 = path1SearchList->FindPoint(currentPoint);
					vtkIdType closestIdPath1 = path1->GetId(closestPosPath1);

					double closestPointPath1[3] = { 0, 0, 0 };
					closestPointPath1[0] = centrePointsUGrid->GetPoint(closestIdPath1)[0];
					closestPointPath1[1] = centrePointsUGrid->GetPoint(closestIdPath1)[1];
					closestPointPath1[2] = centrePointsUGrid->GetPoint(closestIdPath1)[2];

					double distancePath1 = sqrt(abs(vtkMath::Distance2BetweenPoints(currentPoint, closestPointPath1)));

					vtkIdType closestPosPath2 = path2SearchList->FindPoint(currentPoint);
					vtkIdType closestIdPath2 = path2->GetId(closestPosPath2);

					double closestPointPath2[3] = { 0, 0, 0 };
					closestPointPath2[0] = centrePointsUGrid->GetPoint(closestIdPath2)[0];
					closestPointPath2[1] = centrePointsUGrid->GetPoint(closestIdPath2)[1];
					closestPointPath2[2] = centrePointsUGrid->GetPoint(closestIdPath2)[2];

					double distancePath2 = sqrt(abs(vtkMath::Distance2BetweenPoints(currentPoint, closestPointPath2)));

					vector<double> vectorPath1(3);
					vectorPath1.at(0) = data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath1)[0];
					vectorPath1.at(1) = data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath1)[1];
					vectorPath1.at(2) = data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath1)[2];

					vector<double> vectorPath2(3);
					vectorPath2.at(0) = data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath2)[0];
					vectorPath2.at(1) = data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath2)[1];
					vectorPath2.at(2) = data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath2)[2];

					vector<double> newOrientationVector = AveragingOrientation::averageOrientation(vectorPath1, vectorPath2, distancePath1, distancePath2);

					double newDiffVector[3] = { newOrientationVector.at(0), newOrientationVector.at(1), newOrientationVector.at(2) };
					vtkMath::Normalize(newDiffVector);

					data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->InsertTuple(currentId, newDiffVector);
					data.getVtkData()->GetCellData()->GetArray("StatusEndo")->SetComponent(currentId, 0, 3);
				}
			}
		}

		cellIdList->Reset();
		cellIdList = newcellIdList;
	} while (cellIdList->GetNumberOfIds() != 0);

	vtkSmartPointer<vtkIdList> returnIdList = vtkSmartPointer<vtkIdList>::New();

	for (set<vtkIdType>::iterator iter = returnList.begin(); iter != returnList.end(); iter++) {
		returnIdList->InsertNextId(*iter);
	}

	return returnIdList;
} // Methods::regionInterpolateBetween2PathwithSeed

/*! Interpolated the fiber orientation in the region between two paths.
 \param data Pointer to the orginal mesh.
 \param path1 First path with fiber orientation.
 \param path2 Second path with fiber orientation.
 \param inMaterial Tissue class in the fieber orientation could be interpolated.
 \param atrium Atium in that the fiber would be interpolated.
 \retrun A vtkIdList with the cell-ids of the interpolate region.
 */
vtkSmartPointer<vtkIdList> Methods::regionInterpolateBetween2Path(DataFormat &data, vtkSmartPointer<vtkIdList> path1, vtkSmartPointer<vtkIdList> path2, Material::Mat inMaterial, string atrium) {
	return regionInterpolateBetween2Path(data, path1, path2, inMaterial, inMaterial, atrium);
}

/*! Interpolated the fiber orientation in the region between two paths.
 \param data Pointer to the orginal mesh.
 \param path1 First path with fiber orientation.
 \param path2 Second path with fiber orientation.
 \param inMaterial1 First tissue class in the fieber orientation could be interpolated.
 \param inMaterial2 Second tissue class in the fieber orientation could be interpolated.
 \param atrium Atium in that the fiber would be interpolated.
 \retrun A vtkIdList with the cell-ids of the interpolate region.
 */
vtkSmartPointer<vtkIdList> Methods::regionInterpolateBetween2Path(DataFormat &data, vtkSmartPointer<vtkIdList> path1, vtkSmartPointer<vtkIdList> path2, Material::Mat inMaterial1, Material::Mat inMaterial2, string atrium) {
	vector<double> seedPoint = getSeedPoint(data, path1, path2, inMaterial1, inMaterial2, atrium);

	return regionInterpolateBetween2PathwithSeed(data, path1, path2, inMaterial1, inMaterial2, seedPoint);
}

/*! Interpolated the fiber orientation in the region between two paths. The fiber would be interpolate in the direction
 of the first direction point of plane 1 but not in the driection of direction point 2 of plane 2. .
 \param data Pointer to the orginal mesh.
 \param path1 First path with fiber orientation.
 \param path2 Second path with fiber orientation.
 \param inMaterial Tissue class in the fieber orientation could be interpolated.
 \param seedPoint Seed point of the region in that the fiber would be interpolated.
 \param point1Plane1 Point 1 of the first plane as Vector (x,y,z).
 \param point2Plane1 Point 2 of the first plane as Vector (x,y,z).
 \param point3Plane1 Point 3 of the first plane as Vector (x,y,z).
 \param directionPoint1 Direction point of the first plane as Vector (x,y,z) in that the fiber would be interpolated.
 \param point1Plane2 Point 1 of the second plane as Vector (x,y,z).
 \param point2Plane2 Point 2 of the second plane as Vector (x,y,z).
 \param point3Plane2 Point 3 of the second plane as Vector (x,y,z).
 \param directionPoint2 Direction point of the second plane as Vector (x,y,z) in that the fiber would be not
 interpolated.
 \retrun A vtkIdList with the cell-ids of the interpolate region.
 */
vtkSmartPointer<vtkIdList> Methods::regionGrowBetween2PathwithDifferentOrientationInDirectionAndNotInDirection(DataFormat &data, vtkSmartPointer<vtkIdList> path1, vtkSmartPointer<vtkIdList> path2, Material::Mat inMaterial, std::vector<double> seedPoint, std::vector<double> point1Plane1, std::vector<double> point2Plane1, std::vector<double> point3Plane1, std::vector<double> directionPoint1, std::vector<double> point1Plane2, std::vector<double> point2Plane2, std::vector<double> point3Plane2, std::vector<double> directionPoint2) {
	double seedPointArray[3] = { seedPoint.at(0), seedPoint.at(1), seedPoint.at(2) };
	double point1Plane1Array[3] = { point1Plane1.at(0), point1Plane1.at(1), point1Plane1.at(2) };
	double point2Plane1Array[3] = { point2Plane1.at(0), point2Plane1.at(1), point2Plane1.at(2) };
	double point3Plane1Array[3] = { point3Plane1.at(0), point3Plane1.at(1), point3Plane1.at(2) };
	double directionPoint1Array[3] = { directionPoint1.at(0), directionPoint1.at(1), directionPoint1.at(2) };

	double normalPlane1[3];
	double vector1Point1Point2[3];
	double vector1point1Point3[3];

	vtkMath::Subtract(point2Plane1Array, point1Plane1Array, vector1Point1Point2);
	vtkMath::Subtract(point3Plane1Array, point1Plane1Array, vector1point1Point3);
	vtkMath::Cross(vector1Point1Point2, vector1point1Point3, normalPlane1);
	vtkMath::Normalize(normalPlane1);

	double evalPlaneEvalPoint = vtkPlane::Evaluate(normalPlane1, point1Plane1Array, directionPoint1Array);


	double point1Plane2Array[3] = { point1Plane2.at(0), point1Plane2.at(1), point1Plane2.at(2) };
	double point2Plane2Array[3] = { point2Plane2.at(0), point2Plane2.at(1), point2Plane2.at(2) };
	double point3Plane2Array[3] = { point3Plane2.at(0), point3Plane2.at(1), point3Plane2.at(2) };
	double directionPoint2Array[3] = { directionPoint2.at(0), directionPoint2.at(1), directionPoint2.at(2) };

	double normalPlane2[3];
	double vector2Point1Point2[3];
	double vector2point1Point3[3];

	vtkMath::Subtract(point2Plane2Array, point1Plane2Array, vector2Point1Point2);
	vtkMath::Subtract(point3Plane2Array, point1Plane2Array, vector2point1Point3);
	vtkMath::Cross(vector2Point1Point2, vector2point1Point3, normalPlane2);
	vtkMath::Normalize(normalPlane2);

	double evalPlane2EvalPoint = vtkPlane::Evaluate(normalPlane2, point1Plane2Array, directionPoint2Array);

	vtkSmartPointer<vtkPoints> centrePoints = data.getCentrePoints()->GetPoints();
	vtkSmartPointer<vtkUnstructuredGrid> centrePointsUGrid = data.getCentrePoints();

	vtkSmartPointer<vtkUnstructuredGrid> path1SearchList = vtkSmartPointer<vtkUnstructuredGrid>::New();
	vtkSmartPointer<vtkPoints> pointsPath1 = vtkSmartPointer<vtkPoints>::New();

	for (vtkIdType i = 0; i < path1->GetNumberOfIds(); i++) {
		double point[3] = { centrePoints->GetPoint(path1->GetId(i))[0], centrePoints->GetPoint(path1->GetId(i))[1], centrePoints->GetPoint(path1->GetId(i))
		[2] };
		pointsPath1->InsertNextPoint(point);
	}
	path1SearchList->SetPoints(pointsPath1);

	vtkSmartPointer<vtkUnstructuredGrid> path2SearchList = vtkSmartPointer<vtkUnstructuredGrid>::New();
	vtkSmartPointer<vtkPoints> pointsPath2 = vtkSmartPointer<vtkPoints>::New();

	for (vtkIdType i = 0; i < path2->GetNumberOfIds(); i++) {
		double point[3] = { centrePoints->GetPoint(path2->GetId(i))[0], centrePoints->GetPoint(path2->GetId(i))[1], centrePoints->GetPoint(path2->GetId(i))
		[2] };
		pointsPath2->InsertNextPoint(point);
	}
	path2SearchList->SetPoints(pointsPath2);


	// find a starting Point in the allowed material
	vtkIdType seedPointID = centrePointsUGrid->FindPoint(seedPointArray);

	vtkSmartPointer<vtkIdList> cellIdList = vtkSmartPointer<vtkIdList>::New();

	set<vtkIdType> returnList;
	set<vtkIdType> viewedCells;
	cellIdList->InsertNextId(seedPointID);

	do {
		vtkSmartPointer<vtkIdList> newcellIdList = vtkSmartPointer<vtkIdList>::New();

		for (vtkIdType k = 0; k < cellIdList->GetNumberOfIds(); k++) {
			vtkIdType currentCellIds = cellIdList->GetId(k);
			vtkSmartPointer<vtkIdList> cellPointIds = vtkSmartPointer<vtkIdList>::New();

			data.getVtkData()->GetCellPoints(currentCellIds, cellPointIds);

			for (vtkIdType i = 0; i < cellPointIds->GetNumberOfIds(); i++) {
				vtkSmartPointer<vtkIdList> CellPointIdList = vtkSmartPointer<vtkIdList>::New();
				CellPointIdList->InsertNextId(cellPointIds->GetId(i));

				vtkSmartPointer<vtkIdList> neighbors = vtkSmartPointer<vtkIdList>::New();
				data.getVtkData()->GetCellNeighbors(currentCellIds, CellPointIdList, neighbors);

				for (vtkIdType j = 0; j < neighbors->GetNumberOfIds(); j++) {
					vtkIdType currentId = neighbors->GetId(j);

					int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(currentId, 0);
					double evalPlane1 = vtkPlane::Evaluate(normalPlane1, point1Plane1Array, centrePoints->GetPoint(currentId));
					double evalPlane2 = vtkPlane::Evaluate(normalPlane2, point1Plane2Array, centrePoints->GetPoint(currentId));

					if ((viewedCells.find(currentId) == viewedCells.end()) && (inMaterial == material) &&
						(((evalPlane1 <= 0) && (evalPlaneEvalPoint <= 0)) ||
						((evalPlane1 >= 0) && (evalPlaneEvalPoint >= 0))) &&
						!(((evalPlane2 >= 0) && (evalPlane2EvalPoint >= 0)) ||
						((evalPlane2 <= 0) && (evalPlane2EvalPoint <= 0)))) {
						newcellIdList->InsertNextId(currentId);
					}
					viewedCells.insert(currentId);
				}
			}
		}

		for (vtkIdType j = 0; j < newcellIdList->GetNumberOfIds(); j++) {
			vtkIdType currentId = newcellIdList->GetId(j);
			returnList.insert(currentId);
			if (!data.getDoubleLayer()) {
				if (data.getVtkData()->GetCellData()->GetArray("StatusEndo")->GetComponent(currentId, 0) == 0) {
					double currentPoint[3] = { 0, 0, 0 };
					currentPoint[0] = centrePoints->GetPoint(currentId)[0];
					currentPoint[1] = centrePoints->GetPoint(currentId)[1];
					currentPoint[2] = centrePoints->GetPoint(currentId)[2];

					vtkIdType closestPosPath1 = path1SearchList->FindPoint(currentPoint);
					vtkIdType closestIdPath1 = path1->GetId(closestPosPath1);

					double closestPointPath1[3] = { 0, 0, 0 };
					closestPointPath1[0] = centrePointsUGrid->GetPoint(closestIdPath1)[0];
					closestPointPath1[1] = centrePointsUGrid->GetPoint(closestIdPath1)[1];
					closestPointPath1[2] = centrePointsUGrid->GetPoint(closestIdPath1)[2];

					double distancePath1 = sqrt(abs(vtkMath::Distance2BetweenPoints(currentPoint, closestPointPath1)));

					vtkIdType closestPosPath2 = path2SearchList->FindPoint(currentPoint);
					vtkIdType closestIdPath2 = path2->GetId(closestPosPath2);

					double closestPointPath2[3] = { 0, 0, 0 };
					closestPointPath2[0] = centrePointsUGrid->GetPoint(closestIdPath2)[0];
					closestPointPath2[1] = centrePointsUGrid->GetPoint(closestIdPath2)[1];
					closestPointPath2[2] = centrePointsUGrid->GetPoint(closestIdPath2)[2];

					double distancePath2 = sqrt(abs(vtkMath::Distance2BetweenPoints(currentPoint, closestPointPath2)));

					double weightPath2 = distancePath1 / (distancePath1 + distancePath2);
					double weightPath1 = 1 - weightPath2;

					double vectorPath1[3] = {
					data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath1)[0], data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath1)[1], data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath1)[2]
					};

					double vectorPath2[3] = {
					data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath2)[0], data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath2)[1], data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath2)[2]
					};

					vtkMath::MultiplyScalar(vectorPath1, weightPath1);
					vtkMath::MultiplyScalar(vectorPath2, -1 * weightPath2);
					double newDiffVector[3] = { 0, 0, 0 };

					vtkMath::Add(vectorPath1, vectorPath2, newDiffVector);
					vtkMath::Normalize(newDiffVector);

					data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->InsertTuple(currentId, newDiffVector);
					data.getVtkData()->GetCellData()->GetArray("StatusEndo")->SetComponent(currentId, 0, 3);
				}
			}
			if (data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(currentId, 0) == 0) {
				double currentPoint[3] = { 0, 0, 0 };
				currentPoint[0] = centrePoints->GetPoint(currentId)[0];
				currentPoint[1] = centrePoints->GetPoint(currentId)[1];
				currentPoint[2] = centrePoints->GetPoint(currentId)[2];

				vtkIdType closestPosPath1 = path1SearchList->FindPoint(currentPoint);
				vtkIdType closestIdPath1 = path1->GetId(closestPosPath1);

				double closestPointPath1[3] = { 0, 0, 0 };
				closestPointPath1[0] = centrePointsUGrid->GetPoint(closestIdPath1)[0];
				closestPointPath1[1] = centrePointsUGrid->GetPoint(closestIdPath1)[1];
				closestPointPath1[2] = centrePointsUGrid->GetPoint(closestIdPath1)[2];

				double distancePath1 = sqrt(abs(vtkMath::Distance2BetweenPoints(currentPoint, closestPointPath1)));

				vtkIdType closestPosPath2 = path2SearchList->FindPoint(currentPoint);
				vtkIdType closestIdPath2 = path2->GetId(closestPosPath2);

				double closestPointPath2[3] = { 0, 0, 0 };
				closestPointPath2[0] = centrePointsUGrid->GetPoint(closestIdPath2)[0];
				closestPointPath2[1] = centrePointsUGrid->GetPoint(closestIdPath2)[1];
				closestPointPath2[2] = centrePointsUGrid->GetPoint(closestIdPath2)[2];

				double distancePath2 = sqrt(abs(vtkMath::Distance2BetweenPoints(currentPoint, closestPointPath2)));


				double weightPath2 = distancePath1 / (distancePath1 + distancePath2);
				double weightPath1 = 1 - weightPath2;

				double vectorPath1[3] = {
				data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath1)[0], data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath1)[1], data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath1)[2]
				};

				double vectorPath2[3] = {
				data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath2)[0], data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath2)[1], data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath2)[2]
				};

				vtkMath::MultiplyScalar(vectorPath1, weightPath1);
				vtkMath::MultiplyScalar(vectorPath2, -1 * weightPath2);
				double newDiffVector[3] = { 0, 0, 0 };

				vtkMath::Add(vectorPath1, vectorPath2, newDiffVector);
				vtkMath::Normalize(newDiffVector);


				data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->InsertTuple(currentId, newDiffVector);
				data.getVtkData()->GetCellData()->GetArray("Status")->SetComponent(currentId, 0, 3);
			}
		}


		cellIdList->Reset();
		cellIdList = newcellIdList;
	} while (cellIdList->GetNumberOfIds() != 0);

	vtkSmartPointer<vtkIdList> returnIdList = vtkSmartPointer<vtkIdList>::New();

	for (set<vtkIdType>::iterator iter = returnList.begin(); iter != returnList.end(); iter++) {
		returnIdList->InsertNextId(*iter);
	}

	return returnIdList;
} // Methods::regionGrowBetween2PathwithDifferentOrientationInDirectionAndNotInDirection

/*! Interpolated the fiber orientation in the region between two paths. The fiber would be interpolate in the direction
 of the first direction point of plane 1.
 \param data Pointer to the orginal mesh.
 \param path1 First path with fiber orientation.
 \param path2 Second path with fiber orientation.
 \param inMaterial Tissue class in the fieber orientation could be interpolated.
 \param seedPoint Seed point of the region in that the fiber would be interpolated.
 \param point1Plane Point 1 of the first plane as Vector (x,y,z).
 \param point2Plane Point 2 of the first plane as Vector (x,y,z).
 \param point3Plane Point 3 of the first plane as Vector (x,y,z).
 \param directionPoint Direction point of the first plane as Vector (x,y,z) in that the fiber would be interpolated.
 \retrun A vtkIdList with the cell-ids of the interpolate region.
 */
vtkSmartPointer<vtkIdList> Methods::regionInterpolateBetween2PathsInDirection(DataFormat &data, vtkSmartPointer<vtkIdList> path1, vtkSmartPointer<vtkIdList> path2, Material::Mat inMaterial, std::vector<double> seedPoint, std::vector<double> point1Plane, std::vector<double> point2Plane, std::vector<double> point3Plane, std::vector<double> directionPoint)
{
	double seedPointArray[3] = { seedPoint.at(0), seedPoint.at(1), seedPoint.at(2) };
	double point1PlaneArray[3] = { point1Plane.at(0), point1Plane.at(1), point1Plane.at(2) };
	double point2PlaneArray[3] = { point2Plane.at(0), point2Plane.at(1), point2Plane.at(2) };
	double point3PlaneArray[3] = { point3Plane.at(0), point3Plane.at(1), point3Plane.at(2) };
	double directionPointArray[3] = { directionPoint.at(0), directionPoint.at(1), directionPoint.at(2) };

	double normalPlane[3];
	double vectorPoint1Point2[3];
	double vectorpoint1Point3[3];

	vtkMath::Subtract(point2PlaneArray, point1PlaneArray, vectorPoint1Point2);
	vtkMath::Subtract(point3PlaneArray, point1PlaneArray, vectorpoint1Point3);
	vtkMath::Cross(vectorPoint1Point2, vectorpoint1Point3, normalPlane);
	vtkMath::Normalize(normalPlane);

	double evalPlaneEvalPoint = vtkPlane::Evaluate(normalPlane, point1PlaneArray, directionPointArray);

	vtkSmartPointer<vtkPoints> centrePoints = data.getCentrePoints()->GetPoints();
	vtkSmartPointer<vtkUnstructuredGrid> centrePointsUGrid = data.getCentrePoints();

	vtkSmartPointer<vtkUnstructuredGrid> path1SearchList = vtkSmartPointer<vtkUnstructuredGrid>::New();
	vtkSmartPointer<vtkPoints> pointsPath1 = vtkSmartPointer<vtkPoints>::New();

	for (vtkIdType i = 0; i < path1->GetNumberOfIds(); i++) {
		double point[3] = { centrePoints->GetPoint(path1->GetId(i))[0], centrePoints->GetPoint(path1->GetId(i))[1], centrePoints->GetPoint(path1->GetId(i))
		[2] };
		pointsPath1->InsertNextPoint(point);
	}
	path1SearchList->SetPoints(pointsPath1);

	vtkSmartPointer<vtkUnstructuredGrid> path2SearchList = vtkSmartPointer<vtkUnstructuredGrid>::New();
	vtkSmartPointer<vtkPoints> pointsPath2 = vtkSmartPointer<vtkPoints>::New();

	for (vtkIdType i = 0; i < path2->GetNumberOfIds(); i++) {
		double point[3] = { centrePoints->GetPoint(path2->GetId(i))[0], centrePoints->GetPoint(path2->GetId(i))[1], centrePoints->GetPoint(path2->GetId(i))
		[2] };
		pointsPath2->InsertNextPoint(point);
	}
	path2SearchList->SetPoints(pointsPath2);

	// find a starting Point in the allowed material
	vtkIdType seedPointID = centrePointsUGrid->FindPoint(seedPointArray);

	vtkSmartPointer<vtkIdList> cellIdList = vtkSmartPointer<vtkIdList>::New();


	set<vtkIdType> returnList;
	set<vtkIdType> viewedCells;
	cellIdList->InsertNextId(seedPointID);

	do {
		vtkSmartPointer<vtkIdList> newcellIdList = vtkSmartPointer<vtkIdList>::New();

		for (vtkIdType k = 0; k < cellIdList->GetNumberOfIds(); k++) {
			vtkIdType currentCellIds = cellIdList->GetId(k);
			vtkSmartPointer<vtkIdList> cellPointIds = vtkSmartPointer<vtkIdList>::New();

			data.getVtkData()->GetCellPoints(currentCellIds, cellPointIds);

			for (vtkIdType i = 0; i < cellPointIds->GetNumberOfIds(); i++) {
				vtkSmartPointer<vtkIdList> CellPointIdList = vtkSmartPointer<vtkIdList>::New();
				CellPointIdList->InsertNextId(cellPointIds->GetId(i));

				vtkSmartPointer<vtkIdList> neighbors = vtkSmartPointer<vtkIdList>::New();
				data.getVtkData()->GetCellNeighbors(currentCellIds, CellPointIdList, neighbors);

				for (vtkIdType j = 0; j < neighbors->GetNumberOfIds(); j++) {
					vtkIdType currentId = neighbors->GetId(j);

					int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(currentId, 0);
					double evalPlane = vtkPlane::Evaluate(normalPlane, point1PlaneArray, centrePoints->GetPoint(currentId));
					if ((viewedCells.find(currentId) == viewedCells.end()) && (inMaterial == material) &&
						(((evalPlane <= 0) && (evalPlaneEvalPoint <= 0)) || ((evalPlane >= 0) && (evalPlaneEvalPoint >= 0)))) {
						newcellIdList->InsertNextId(currentId);
					}
					viewedCells.insert(currentId);
				}
			}
		}

		for (vtkIdType j = 0; j < newcellIdList->GetNumberOfIds(); j++) {
			vtkIdType currentId = newcellIdList->GetId(j);
			returnList.insert(currentId);
			if (!data.getDoubleLayer()) {
				if (data.getVtkData()->GetCellData()->GetArray("StatusEndo")->GetComponent(currentId, 0) == 0) {
					double currentPoint[3] = { 0, 0, 0 };
					currentPoint[0] = centrePoints->GetPoint(currentId)[0];
					currentPoint[1] = centrePoints->GetPoint(currentId)[1];
					currentPoint[2] = centrePoints->GetPoint(currentId)[2];

					vtkIdType closestPosPath1 = path1SearchList->FindPoint(currentPoint);
					vtkIdType closestIdPath1 = path1->GetId(closestPosPath1);

					double closestPointPath1[3] = { 0, 0, 0 };
					closestPointPath1[0] = centrePointsUGrid->GetPoint(closestIdPath1)[0];
					closestPointPath1[1] = centrePointsUGrid->GetPoint(closestIdPath1)[1];
					closestPointPath1[2] = centrePointsUGrid->GetPoint(closestIdPath1)[2];

					double distancePath1 = sqrt(abs(vtkMath::Distance2BetweenPoints(currentPoint, closestPointPath1)));

					vtkIdType closestPosPath2 = path2SearchList->FindPoint(currentPoint);
					vtkIdType closestIdPath2 = path2->GetId(closestPosPath2);

					double closestPointPath2[3] = { 0, 0, 0 };
					closestPointPath2[0] = centrePointsUGrid->GetPoint(closestIdPath2)[0];
					closestPointPath2[1] = centrePointsUGrid->GetPoint(closestIdPath2)[1];
					closestPointPath2[2] = centrePointsUGrid->GetPoint(closestIdPath2)[2];

					double distancePath2 = sqrt(abs(vtkMath::Distance2BetweenPoints(currentPoint, closestPointPath2)));

					vector<double> vectorPath1(3);
					vectorPath1.at(0) = data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath1)[0];
					vectorPath1.at(1) = data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath1)[1];
					vectorPath1.at(2) = data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath1)[2];


					vector<double> vectorPath2(3);
					vectorPath2.at(0) = data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath2)[0];
					vectorPath2.at(1) = data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath2)[1];
					vectorPath2.at(2) = data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath2)[2];

					vector<double> newOrientationVector = AveragingOrientation::averageOrientation(vectorPath1, vectorPath2, distancePath1, distancePath2);

					double newDiffVector[3] = { newOrientationVector.at(0), newOrientationVector.at(1), newOrientationVector.at(2) };
					vtkMath::Normalize(newDiffVector);


					data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->InsertTuple(currentId, newDiffVector);
					data.getVtkData()->GetCellData()->GetArray("StatusEndo")->SetComponent(currentId, 0, 3);
				}
			}
			if (data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(currentId, 0) == 0) {
				double currentPoint[3] = { 0, 0, 0 };
				currentPoint[0] = centrePoints->GetPoint(currentId)[0];
				currentPoint[1] = centrePoints->GetPoint(currentId)[1];
				currentPoint[2] = centrePoints->GetPoint(currentId)[2];

				vtkIdType closestPosPath1 = path1SearchList->FindPoint(currentPoint);
				vtkIdType closestIdPath1 = path1->GetId(closestPosPath1);

				double closestPointPath1[3] = { 0, 0, 0 };
				closestPointPath1[0] = centrePointsUGrid->GetPoint(closestIdPath1)[0];
				closestPointPath1[1] = centrePointsUGrid->GetPoint(closestIdPath1)[1];
				closestPointPath1[2] = centrePointsUGrid->GetPoint(closestIdPath1)[2];

				double distancePath1 = sqrt(abs(vtkMath::Distance2BetweenPoints(currentPoint, closestPointPath1)));

				vtkIdType closestPosPath2 = path2SearchList->FindPoint(currentPoint);
				vtkIdType closestIdPath2 = path2->GetId(closestPosPath2);

				double closestPointPath2[3] = { 0, 0, 0 };
				closestPointPath2[0] = centrePointsUGrid->GetPoint(closestIdPath2)[0];
				closestPointPath2[1] = centrePointsUGrid->GetPoint(closestIdPath2)[1];
				closestPointPath2[2] = centrePointsUGrid->GetPoint(closestIdPath2)[2];

				double distancePath2 = sqrt(abs(vtkMath::Distance2BetweenPoints(currentPoint, closestPointPath2)));

				vector<double> vectorPath1(3);
				vectorPath1.at(0) = data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath1)[0];
				vectorPath1.at(1) = data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath1)[1];
				vectorPath1.at(2) = data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath1)[2];


				vector<double> vectorPath2(3);
				vectorPath2.at(0) = data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath2)[0];
				vectorPath2.at(1) = data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath2)[1];
				vectorPath2.at(2) = data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath2)[2];

				vector<double> newOrientationVector = AveragingOrientation::averageOrientation(vectorPath1, vectorPath2, distancePath1, distancePath2);

				double newDiffVector[3] = { newOrientationVector.at(0), newOrientationVector.at(1), newOrientationVector.at(2) };
				vtkMath::Normalize(newDiffVector);


				data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->InsertTuple(currentId, newDiffVector);
				data.getVtkData()->GetCellData()->GetArray("Status")->SetComponent(currentId, 0, 3);
			}
		}


		cellIdList->Reset();
		cellIdList = newcellIdList;
	} while (cellIdList->GetNumberOfIds() != 0);

	vtkSmartPointer<vtkIdList> returnIdList = vtkSmartPointer<vtkIdList>::New();

	for (set<vtkIdType>::iterator iter = returnList.begin(); iter != returnList.end(); iter++) {
		returnIdList->InsertNextId(*iter);
	}

	return returnIdList;
} // Methods::regionInterpolateBetween2PathsInDirection

/*! Grow the fiber orientation of a path in the region.
 \param data Pointer to the orginal mesh.
 \param path Path with fiber orientation.
 \param region Region in that the path would grown.
 \param inMaterial1 First tissue class in the fieber orientation could be grown.
 \param inMaterial2 Second tissue class in the fieber orientation could be grown.
 \retrun A vtkIdList with the cell-ids.
 */
vtkSmartPointer<vtkIdList> Methods::growPathInRegion(DataFormat &data, vtkSmartPointer<vtkIdList> path, vtkSmartPointer<vtkIdList> region, Material::Mat inMaterial1, Material::Mat inMaterial2) {
	vtkSmartPointer<vtkPoints> centrePoints = data.getCentrePoints()->GetPoints();

	vtkSmartPointer<vtkUnstructuredGrid> pathSearchList = vtkSmartPointer<vtkUnstructuredGrid>::New();
	vtkSmartPointer<vtkPoints> pointsPath = vtkSmartPointer<vtkPoints>::New();

	for (vtkIdType i = 0; i < path->GetNumberOfIds(); i++) {
		double point[3] = { centrePoints->GetPoint(path->GetId(i))[0], centrePoints->GetPoint(path->GetId(i))[1], centrePoints->GetPoint(path->GetId(i))
		[2] };
		pointsPath->InsertNextPoint(point);
	}
	pathSearchList->SetPoints(pointsPath);

	set<vtkIdType> cells;

	for (vtkIdType j = 0; j < region->GetNumberOfIds(); j++) {
		if (!data.getDoubleLayer()) {
			int material = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(region->GetId(j), 0);
			if ((data.getVtkData()->GetCellData()->GetArray("StatusEndo")->GetComponent(region->GetId(j), 0) == 0) &&
				((inMaterial1 == material) || (inMaterial2 == material))) {
				cells.insert(region->GetId(j));
			}
		}

		int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(region->GetId(j), 0);
		if ((data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(region->GetId(j), 0) == 0) &&
			((inMaterial1 == material) || (inMaterial2 == material))) {
			cells.insert(region->GetId(j));
		}
	}

	vtkSmartPointer<vtkIdList> returnIdList = vtkSmartPointer<vtkIdList>::New();

	for (set<vtkIdType>::iterator iter = cells.begin(); iter != cells.end(); iter++) {
		vtkIdType currentId = *iter;

		double currentPoint[3] = { 0, 0, 0 };
		currentPoint[0] = centrePoints->GetPoint(currentId)[0];
		currentPoint[1] = centrePoints->GetPoint(currentId)[1];
		currentPoint[2] = centrePoints->GetPoint(currentId)[2];

		vtkIdType closestPosPath = pathSearchList->FindPoint(currentPoint);
		vtkIdType closestIdPath = path->GetId(closestPosPath);

		double *closestPointPath1 = data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath);

		if (!data.getDoubleLayer()) {
			data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->InsertTuple(currentId, closestPointPath1);
			data.getVtkData()->GetCellData()->GetArray("StatusEndo")->SetComponent(currentId, 0, 3);
		}

		data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->InsertTuple(currentId, closestPointPath1);
		data.getVtkData()->GetCellData()->GetArray("Status")->SetComponent(currentId, 0, 3);

		returnIdList->InsertNextId(*iter);
	}


	return returnIdList;
} // Methods::growPathInRegion

/*! Grow the fiber orientation of a path in the region.
 \param data Pointer to the orginal mesh.
 \param path Path with fiber orientation.
 \param inMaterial1 First tissue class in the fieber orientation could be grown.
 \param inMaterial12Second tissue class in the fieber orientation could be grown.
 \param seedPoint Seed point of the region in that the fiber would be grown.
 \retrun A vtkIdList with the cell-ids.
 */
vtkSmartPointer<vtkIdList> Methods::growPathInRegion(DataFormat &data, vtkSmartPointer<vtkIdList> path, Material::Mat inMaterial1, Material::Mat inMaterial2, vector<double> seedPoint) {
	double seedPointArray[3] = { seedPoint.at(0), seedPoint.at(1), seedPoint.at(2) };

	vtkSmartPointer<vtkPoints> centrePoints = data.getCentrePoints()->GetPoints();
	vtkSmartPointer<vtkUnstructuredGrid> centrePointsUGrid = data.getCentrePoints();


	vtkSmartPointer<vtkUnstructuredGrid> pathSearchList = vtkSmartPointer<vtkUnstructuredGrid>::New();
	vtkSmartPointer<vtkPoints> pointsPath = vtkSmartPointer<vtkPoints>::New();

	for (vtkIdType i = 0; i < path->GetNumberOfIds(); i++) {
		double point[3] = { centrePoints->GetPoint(path->GetId(i))[0], centrePoints->GetPoint(path->GetId(i))[1], centrePoints->GetPoint(path->GetId(i))
		[2] };
		pointsPath->InsertNextPoint(point);
	}
	pathSearchList->SetPoints(pointsPath);


	// find a starting Point in the allowed material
	set<vtkIdType> viewedCellList;

	vtkIdType seedPointID = centrePointsUGrid->FindPoint(seedPointArray);
	set<vtkIdType> cellIdList;

	cellIdList.insert(seedPointID);

	bool end = false;

	while (!end) {
		set<vtkIdType> newcellIdList;
		for (set<vtkIdType>::iterator iter = cellIdList.begin(); iter != cellIdList.end(); iter++) {
			vtkIdType currentCellIds = *iter;
			vtkSmartPointer<vtkIdList> cellPointIds = vtkSmartPointer<vtkIdList>::New();

			data.getVtkData()->GetCellPoints(currentCellIds, cellPointIds);

			for (vtkIdType i = 0; i < cellPointIds->GetNumberOfIds(); i++) {
				vtkSmartPointer<vtkIdList> CellPointIdList = vtkSmartPointer<vtkIdList>::New();
				CellPointIdList->InsertNextId(cellPointIds->GetId(i));

				vtkSmartPointer<vtkIdList> tempNeighborCellIds = vtkSmartPointer<vtkIdList>::New();
				data.getVtkData()->GetCellNeighbors(currentCellIds, CellPointIdList, tempNeighborCellIds);

				for (vtkIdType j = 0; j < tempNeighborCellIds->GetNumberOfIds(); j++) {
					if (!data.getDoubleLayer()) {
						int material = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(tempNeighborCellIds->GetId(j), 0);

						if ((viewedCellList.find(tempNeighborCellIds->GetId(j)) == viewedCellList.end()) &&
							(data.getVtkData()->GetCellData()->GetArray("StatusEndo")->GetComponent(tempNeighborCellIds->GetId(j), 0) == 0) &&
							((inMaterial1 == material) || (inMaterial2 == material))
							) {
							viewedCellList.insert(tempNeighborCellIds->GetId(j));
							newcellIdList.insert(tempNeighborCellIds->GetId(j));
						}
					}
					int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(tempNeighborCellIds->GetId(j), 0);

					if ((viewedCellList.find(tempNeighborCellIds->GetId(j)) == viewedCellList.end()) &&
						(data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(tempNeighborCellIds->GetId(j), 0) == 0) &&
						((inMaterial1 == material) || (inMaterial2 == material))
						) {
						viewedCellList.insert(tempNeighborCellIds->GetId(j));
						newcellIdList.insert(tempNeighborCellIds->GetId(j));
					} // alle nachbarzellen neuen
				}
			}
		}
		if (newcellIdList.size() == 0) {
			end = true;
		}
		cellIdList.clear();
		cellIdList = newcellIdList;
	}

	vtkSmartPointer<vtkIdList> returnIdList = vtkSmartPointer<vtkIdList>::New();

	for (set<vtkIdType>::iterator iter = viewedCellList.begin(); iter != viewedCellList.end(); iter++) {
		vtkIdType currentId = *iter;

		double currentPoint[3] = { 0, 0, 0 };
		currentPoint[0] = centrePoints->GetPoint(currentId)[0];
		currentPoint[1] = centrePoints->GetPoint(currentId)[1];
		currentPoint[2] = centrePoints->GetPoint(currentId)[2];

		vtkIdType closestPosPath = pathSearchList->FindPoint(currentPoint);
		vtkIdType closestIdPath = path->GetId(closestPosPath);

		if (!data.getDoubleLayer()) {
			double *closestPointPath1 = data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath);
			data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->InsertTuple(currentId, closestPointPath1);
			data.getVtkData()->GetCellData()->GetArray("StatusEndo")->SetComponent(currentId, 0, 3);
		}

		double *closestPointPath1 = data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath);
		data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->InsertTuple(currentId, closestPointPath1);
		data.getVtkData()->GetCellData()->GetArray("Status")->SetComponent(currentId, 0, 3);

		returnIdList->InsertNextId(*iter);
	}


	return returnIdList;
} // Methods::growPathInRegion

/*! Grow the fiber orientation of a path in the region.
 \param data Pointer to the orginal mesh.
 \param path Path with fiber orientation.
 \param inMaterial1 First tissue class in the fieber orientation could be grown.
 \param inMaterial12Second tissue class in the fieber orientation could be grown.
 \param seedPoint Seed point of the region in that the fiber would be grown.
 \retrun A vtkIdList with the cell-ids.
 */
vtkSmartPointer<vtkIdList> Methods::growPathInRegion(DataFormat &data, vtkSmartPointer<vtkIdList> path, double width, Material::Mat inMaterial1, Material::Mat inMaterial2, vector<double> seedPoint) {
	double seedPointArray[3] = { seedPoint.at(0), seedPoint.at(1), seedPoint.at(2) };

	vtkSmartPointer<vtkPoints> centrePoints = data.getCentrePoints()->GetPoints();
	vtkSmartPointer<vtkUnstructuredGrid> centrePointsUGrid = data.getCentrePoints();


	vtkSmartPointer<vtkUnstructuredGrid> pathSearchList = vtkSmartPointer<vtkUnstructuredGrid>::New();
	vtkSmartPointer<vtkPoints> pointsPath = vtkSmartPointer<vtkPoints>::New();

	for (vtkIdType i = 0; i < path->GetNumberOfIds(); i++) {
		double point[3] = { centrePoints->GetPoint(path->GetId(i))[0], centrePoints->GetPoint(path->GetId(i))[1], centrePoints->GetPoint(path->GetId(i))
		[2] };
		pointsPath->InsertNextPoint(point);
	}
	pathSearchList->SetPoints(pointsPath);

	vtkSmartPointer<vtkPointLocator> centrePointLocator = vtkSmartPointer<vtkPointLocator>::New();
	centrePointLocator->SetDataSet(centrePointsUGrid);
	centrePointLocator->BuildLocator();

	vtkSmartPointer<vtkDoubleArray> tempArray = vtkSmartPointer<vtkDoubleArray>::New();
	tempArray->SetNumberOfComponents(3);
	tempArray->SetNumberOfTuples(centrePoints->GetNumberOfPoints());
	tempArray->FillComponent(0, 0);
	tempArray->FillComponent(1, 0);
	tempArray->FillComponent(2, 0);

	set<vtkIdType> region;

	for (vtkIdType j = 0; j < path->GetNumberOfIds(); j++) {
		vtkIdType parentId = path->GetId(j);

		vtkSmartPointer<vtkIdList> neighbors = vtkSmartPointer<vtkIdList>::New();
		centrePointLocator->FindPointsWithinRadius(width, centrePoints->GetPoint(parentId), neighbors);

		for (vtkIdType k = 0; k < neighbors->GetNumberOfIds(); k++) {
			if (region.find(neighbors->GetId(k)) == region.end()) {
				region.insert(neighbors->GetId(k));
			}
		}
	}


	// find a starting Point in the allowed material
	set<vtkIdType> viewedCellList;

	vtkIdType seedPointID = centrePointsUGrid->FindPoint(seedPointArray);
	set<vtkIdType> cellIdList;

	cellIdList.insert(seedPointID);

	bool end = false;

	while (!end) {
		set<vtkIdType> newcellIdList;
		for (set<vtkIdType>::iterator iter = cellIdList.begin(); iter != cellIdList.end(); iter++) {
			vtkIdType currentCellIds = *iter;
			vtkSmartPointer<vtkIdList> cellPointIds = vtkSmartPointer<vtkIdList>::New();

			data.getVtkData()->GetCellPoints(currentCellIds, cellPointIds);

			for (vtkIdType i = 0; i < cellPointIds->GetNumberOfIds(); i++) {
				vtkSmartPointer<vtkIdList> CellPointIdList = vtkSmartPointer<vtkIdList>::New();
				CellPointIdList->InsertNextId(cellPointIds->GetId(i));

				vtkSmartPointer<vtkIdList> tempNeighborCellIds = vtkSmartPointer<vtkIdList>::New();
				data.getVtkData()->GetCellNeighbors(currentCellIds, CellPointIdList, tempNeighborCellIds);

				for (vtkIdType j = 0; j < tempNeighborCellIds->GetNumberOfIds(); j++) {
					if (!data.getDoubleLayer()) {
						int material = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(tempNeighborCellIds->GetId(j), 0);

						if ((viewedCellList.find(tempNeighborCellIds->GetId(j)) == viewedCellList.end()) &&
							(data.getVtkData()->GetCellData()->GetArray("StatusEndo")->GetComponent(tempNeighborCellIds->GetId(j), 0) == 0) &&
							((inMaterial1 == material) || (inMaterial2 == material))
							) {
							viewedCellList.insert(tempNeighborCellIds->GetId(j));
							newcellIdList.insert(tempNeighborCellIds->GetId(j));
						}
					}
					int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(tempNeighborCellIds->GetId(j), 0);

					if ((viewedCellList.find(tempNeighborCellIds->GetId(j)) == viewedCellList.end()) &&
						(data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(tempNeighborCellIds->GetId(j), 0) == 0) &&
						((inMaterial1 == material) || (inMaterial2 == material))
						) {
						viewedCellList.insert(tempNeighborCellIds->GetId(j));
						newcellIdList.insert(tempNeighborCellIds->GetId(j));
					} // alle nachbarzellen neuen
				}
			}
		}
		if (newcellIdList.size() == 0) {
			end = true;
		}
		cellIdList.clear();
		cellIdList = newcellIdList;
	}

	vtkSmartPointer<vtkIdList> returnIdList = vtkSmartPointer<vtkIdList>::New();

	for (set<vtkIdType>::iterator iter = viewedCellList.begin(); iter != viewedCellList.end(); iter++) {
		vtkIdType currentId = *iter;

		if (region.find(currentId) != region.end()) {
			double currentPoint[3] = { 0, 0, 0 };
			currentPoint[0] = centrePoints->GetPoint(currentId)[0];
			currentPoint[1] = centrePoints->GetPoint(currentId)[1];
			currentPoint[2] = centrePoints->GetPoint(currentId)[2];

			vtkIdType closestPosPath = pathSearchList->FindPoint(currentPoint);
			vtkIdType closestIdPath = path->GetId(closestPosPath);

			if (!data.getDoubleLayer()) {
				double *closestPointPath1 = data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath);
				data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->InsertTuple(currentId, closestPointPath1);
				data.getVtkData()->GetCellData()->GetArray("StatusEndo")->SetComponent(currentId, 0, 3);
			}

			double *closestPointPath1 = data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(closestIdPath);
			data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->InsertTuple(currentId, closestPointPath1);
			data.getVtkData()->GetCellData()->GetArray("Status")->SetComponent(currentId, 0, 3);

			returnIdList->InsertNextId(*iter);
		}
	}
	returnIdList = add2Paths(returnIdList, path);

	return returnIdList;
} // Methods::growPathInRegion

/*! Grow the fiber orientation from the vein base to the top of the vein on a surface mesh. The fiber grow only in
 define region.
 \param data Pointer to the orginal mesh.
 \param path Path with fiber orientation.
 \param inMaterial1 First tissue class in the fieber orientation could be grown.
 \param inMaterial2 Second tissue class in the fieber orientation could be grown.
 \param materialVein Tissue class of the vein oufire.
 \param seedPoint Seed point of the region in that the fiber would be grown.
 \param distancePoint Outer point of the vein .
 \param distanceEndoEpi Fictive distance between endo- and epicard in the left atrial.
 \param regionList Region in that the fiber could be grown.
 \retrun A vtkIdList with the cell-ids.
 */
vtkSmartPointer<vtkIdList> Methods::growVeinInRegionOneSurface(DataFormat &data, Material::Mat inMaterial1, Material::Mat inMaterial2, Material::Mat materialVein, vector<double> distancePoint, double distanceEndoEpi, vtkSmartPointer<vtkIdList> regionList) {
	double distancePointArray[3] = { distancePoint.at(0), distancePoint.at(1), distancePoint.at(2) };

	vtkSmartPointer<vtkPoints> centrePoints = data.getCentrePoints()->GetPoints();

	vtkIdType seedPointID = Methods::findClosedPointIdOnPathWithoutOrientationInMaterial(data, distancePoint, regionList, inMaterial1, inMaterial2);

	vtkSmartPointer<vtkDoubleArray> cells = vtkSmartPointer<vtkDoubleArray>::New();

	cells->SetNumberOfComponents(4);
	cells->SetComponentName(0, "CellID");
	cells->SetComponentName(1, "distance to centre point");
	cells->SetComponentName(2, "Epi");
	cells->SetComponentName(3, "Endo");

	set<vtkIdType> viewedCellList;

	bool end = false;

	vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
	set<vtkIdType> cellIdList;

	cellIdList.insert(seedPointID);

	if ((data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(seedPointID, 0) == 0) ||
		(data.getVtkData()->GetCellData()->GetArray("StatusEndo")->GetComponent(seedPointID, 0) == 0)) {
		double cellsArray[4] = { 0, 0, 0, 0 };
		cellsArray[0] = seedPointID;
		if (data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(seedPointID, 0) == 0) {
			cellsArray[2] = 1;
		}
		if (data.getVtkData()->GetCellData()->GetArray("StatusEndo")->GetComponent(seedPointID, 0) == 0) {
			cellsArray[3] = 1;
		}
		cells->InsertNextTuple(cellsArray);
		cellIds->InsertNextId(seedPointID);
	}

	set<vtkIdType> region;

	for (vtkIdType i = 0; i < regionList->GetNumberOfIds(); i++) {
		region.insert(regionList->GetId(i));
	}

	viewedCellList.insert(seedPointID);
	unordered_map<int, vtkSmartPointer<vtkIdList>> containerMaterialCellIDs;
	unordered_map<int, vtkSmartPointer<vtkIdList>> containerMaterialCellIDsEndo;


	while (!end) {
		set<vtkIdType> newcellIdList;

		for (set<vtkIdType>::iterator iter = cellIdList.begin(); iter != cellIdList.end(); iter++) {
			vtkIdType currentCellIds = *iter;
			vtkSmartPointer<vtkIdList> cellPointIds = vtkSmartPointer<vtkIdList>::New();

			data.getVtkData()->GetCellPoints(currentCellIds, cellPointIds);

			for (vtkIdType i = 0; i < cellPointIds->GetNumberOfIds(); i++) {
				vtkSmartPointer<vtkIdList> CellPointIdList = vtkSmartPointer<vtkIdList>::New();
				CellPointIdList->InsertNextId(cellPointIds->GetId(i));

				vtkSmartPointer<vtkIdList> tempNeighborCellIds = vtkSmartPointer<vtkIdList>::New();
				data.getVtkData()->GetCellNeighbors(currentCellIds, CellPointIdList, tempNeighborCellIds);

				for (vtkIdType j = 0; j < tempNeighborCellIds->GetNumberOfIds(); j++) {
					if ((viewedCellList.find(tempNeighborCellIds->GetId(j)) == viewedCellList.end()) &&
						(region.find(tempNeighborCellIds->GetId(j)) != region.end())) {
						int material1 = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(tempNeighborCellIds->GetId(j), 0);
						int material2 = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(tempNeighborCellIds->GetId(j), 0);

						if (Material::isInLeft(material1) || Material::isInRight(material1) || Material::isInLeft(material2) ||
							Material::isInRight(material2)) {
							if ((data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(tempNeighborCellIds->GetId(j), 0) == 0) ||
								(data.getVtkData()->GetCellData()->GetArray("StatusEndo")->GetComponent(tempNeighborCellIds->GetId(j), 0) == 0)) {
								double cellsArray[4] = { 0, 0, 0, 0 };
								cellsArray[0] = tempNeighborCellIds->GetId(j);
								if (data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(tempNeighborCellIds->GetId(j), 0) == 0) {
									cellsArray[2] = 1;
								}
								if (data.getVtkData()->GetCellData()->GetArray("StatusEndo")->GetComponent(tempNeighborCellIds->GetId(j), 0) == 0) {
									cellsArray[3] = 1;
								}
								cells->InsertNextTuple(cellsArray);
								newcellIdList.insert(tempNeighborCellIds->GetId(j));
								cellIds->InsertNextId(tempNeighborCellIds->GetId(j));
							}


							if ((material1 != Material::Vorhof_links) &&
								(material1 != Material::Vorhof_links_Endo) &&
								(material1 != Material::Vorhof_rechts) &&
								(data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(tempNeighborCellIds->GetId(j), 0) != 0) &&
								(data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetComponent(tempNeighborCellIds->
									GetId(j), 0) +
									data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetComponent(tempNeighborCellIds->
										GetId(j), 1) +
									data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetComponent(tempNeighborCellIds->
										GetId(j), 2) != 0)
								) {
								if (containerMaterialCellIDs.find(material1) == containerMaterialCellIDs.end()) {
									vtkSmartPointer<vtkIdList> containerCellIds = vtkSmartPointer<vtkIdList>::New();
									containerCellIds->InsertNextId(tempNeighborCellIds->GetId(j));
									pair<int, vtkSmartPointer<vtkIdList>> insert(material1, containerCellIds);
									containerMaterialCellIDs.insert(insert);
								}
								else {
									vtkSmartPointer<vtkIdList> containerCellIds = containerMaterialCellIDs.at(material1);
									containerCellIds->InsertUniqueId(tempNeighborCellIds->GetId(j));
									containerMaterialCellIDs.at(material1) = containerCellIds;
								}
							}
							if ((material2 != Material::Vorhof_links) &&
								(material2 != Material::Vorhof_links_Endo) &&
								(material2 != Material::Vorhof_rechts) &&
								(data.getVtkData()->GetCellData()->GetArray("StatusEndo")->GetComponent(tempNeighborCellIds->GetId(j), 0) != 0) &&
								(data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->GetComponent(tempNeighborCellIds
									->GetId(j), 0) +
									data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->GetComponent(tempNeighborCellIds
										->GetId(j), 1) +
									data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->GetComponent(tempNeighborCellIds
										->GetId(j), 2) != 0)
								) {
								if (containerMaterialCellIDsEndo.find(material2) == containerMaterialCellIDsEndo.end()) {
									vtkSmartPointer<vtkIdList> containerCellIdsEndo = vtkSmartPointer<vtkIdList>::New();
									containerCellIdsEndo->InsertNextId(tempNeighborCellIds->GetId(j));
									pair<int, vtkSmartPointer<vtkIdList>> insert(material2, containerCellIdsEndo);
									containerMaterialCellIDsEndo.insert(insert);
								}
								else {
									vtkSmartPointer<vtkIdList> containerCellIdsEndo = containerMaterialCellIDsEndo.at(material2);
									containerCellIdsEndo->InsertUniqueId(tempNeighborCellIds->GetId(j));
									containerMaterialCellIDsEndo.at(material2) = containerCellIdsEndo;
								}
							}
						}
					}

					viewedCellList.insert(tempNeighborCellIds->GetId(j));
				}
			}
		}
		if (newcellIdList.size() == 0) {
			end = true;
		}
		cellIdList.clear();
		cellIdList = newcellIdList;
	}

	for (vtkIdType i = 0; i < cells->GetNumberOfTuples(); i++) {
		double distance = sqrt(abs(vtkMath::Distance2BetweenPoints(distancePointArray, centrePoints->GetPoint(cells->vtkDataArray::GetComponent(i, 0)))));
		cells->vtkDataArray::SetComponent(i, 1, distance);
	}

	vtkSortDataArray::SortArrayByComponent(cells, 1);

	unordered_map<int, vtkSmartPointer<vtkUnstructuredGrid>> containerMaterialSearchPath;
	unordered_map<int, vtkSmartPointer<vtkUnstructuredGrid>> containerMaterialSearchPathEndo;


	if ((containerMaterialCellIDs.size() > 0) || (containerMaterialCellIDsEndo.size() > 0)) {
		for (unordered_map<int, vtkSmartPointer<vtkIdList>>::const_iterator iter = containerMaterialCellIDs.begin();
			iter != containerMaterialCellIDs.end(); iter++) {
			vtkSmartPointer<vtkUnstructuredGrid> SearchList = vtkSmartPointer<vtkUnstructuredGrid>::New();
			vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();


			for (vtkIdType i = 0; i < iter->second->GetNumberOfIds(); i++) {
				points->InsertNextPoint(centrePoints->GetPoint(iter->second->GetId(i)));
			}
			SearchList->SetPoints(points);

			pair<int, vtkSmartPointer<vtkUnstructuredGrid>> insert(iter->first, SearchList);
			containerMaterialSearchPath.insert(insert);
		}


		for (unordered_map<int, vtkSmartPointer<vtkIdList>>::const_iterator iter = containerMaterialCellIDsEndo.begin();
			iter != containerMaterialCellIDsEndo.end(); iter++) {
			vtkSmartPointer<vtkUnstructuredGrid> SearchList = vtkSmartPointer<vtkUnstructuredGrid>::New();
			vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();


			for (vtkIdType i = 0; i < iter->second->GetNumberOfIds(); i++) {
				points->InsertNextPoint(centrePoints->GetPoint(iter->second->GetId(i)));
			}
			SearchList->SetPoints(points);

			pair<int, vtkSmartPointer<vtkUnstructuredGrid>> insert(iter->first, SearchList);
			containerMaterialSearchPathEndo.insert(insert);
		}


		for (vtkIdType j = cells->GetNumberOfTuples() - 1; j >= 0; j--) {
			if (cells->vtkDataArray::GetComponent(j, 2) == 1) { // Epi
				vtkIdType currentId = cells->vtkDataArray::GetComponent(j, 0);

				double currentPoint[3] = { 0, 0, 0 };
				currentPoint[0] = centrePoints->GetPoint(currentId)[0];
				currentPoint[1] = centrePoints->GetPoint(currentId)[1];
				currentPoint[2] = centrePoints->GetPoint(currentId)[2];

				list<tuple<int, double, vector<double>, vector<double>>> differenceVectorList;

				for (unordered_map<int, vtkSmartPointer<vtkUnstructuredGrid>>::const_iterator iter = containerMaterialSearchPath.begin();
					iter != containerMaterialSearchPath.end(); iter++) {
					vtkIdType closestPosPath = iter->second->FindPoint(currentPoint);
					vtkIdType closestIdPoint = containerMaterialCellIDs.find(iter->first)->second->GetId(closestPosPath);

					double closestPoint[3] = { 0, 0, 0 };
					closestPoint[0] = centrePoints->GetPoint(closestIdPoint)[0];
					closestPoint[1] = centrePoints->GetPoint(closestIdPoint)[1];
					closestPoint[2] = centrePoints->GetPoint(closestIdPoint)[2];

					double distance = sqrt(abs(vtkMath::Distance2BetweenPoints(currentPoint, closestPoint)));
					double *materialDiffArray = data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetTuple(closestIdPoint);

					vector<double> materialDiffVector = { materialDiffArray[0], materialDiffArray[1], materialDiffArray[2] };
					tuple<int, double, vector<double>, vector<double>> tupleDiffVector(iter->first, distance, materialDiffVector, materialDiffVector);

					differenceVectorList.push_back(tupleDiffVector);
				}

				for (unordered_map<int, vtkSmartPointer<vtkUnstructuredGrid>>::const_iterator iter = containerMaterialSearchPathEndo.begin();
					iter != containerMaterialSearchPathEndo.end(); iter++) {
					vtkIdType closestPosPath = iter->second->FindPoint(currentPoint);
					vtkIdType closestIdPoint = containerMaterialCellIDsEndo.find(iter->first)->second->GetId(closestPosPath);


					double distance = calcDistancetoPointinEndo(data, currentId, closestIdPoint, distanceEndoEpi);
					double *materialDiffArray = data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->GetTuple(closestIdPoint);

					vector<double> materialDiffVector = { materialDiffArray[0], materialDiffArray[1], materialDiffArray[2] };
					tuple<int, double, vector<double>, vector<double>> tupleDiffVector(iter->first, distance, materialDiffVector, materialDiffVector);

					differenceVectorList.push_back(tupleDiffVector);
				}

				double minDistance = std::numeric_limits<double>::max();
				double sumOfDistance = 0;

				for (list<tuple<int, double, vector<double>, vector<double>>>::iterator iter = differenceVectorList.begin();
					iter != differenceVectorList.end(); iter++) {
					if (minDistance > get<1>(*iter)) {
						minDistance = get<1>(*iter);
					}
					if (get<1>(*iter) > 0) {
						get<1>(*iter) = 1 / get<1>(*iter);
					}
					else {
						get<1>(*iter) = 1;
					}
					sumOfDistance += get<1>(*iter);
				}


				for (list<tuple<int, double, vector<double>, vector<double>>>::iterator iter = differenceVectorList.begin();
					iter != differenceVectorList.end(); iter++) {
					get<1>(*iter) = get<1>(*iter) / sumOfDistance;
				}

				vector<double> newOrientationVector = AveragingOrientation::AveragingNVectors(differenceVectorList);
				double newDiffVector[3] = { newOrientationVector.at(0), newOrientationVector.at(1), newOrientationVector.at(2) };

				data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->InsertTuple(currentId, newDiffVector);
				data.getVtkData()->GetCellData()->GetArray("Status")->SetComponent(currentId, 0, 3);
				data.getVtkData()->GetCellData()->GetArray("Material")->SetComponent(currentId, 0, materialVein);
			}

			if (cells->vtkDataArray::GetComponent(j, 3) == 1) { // Endo
				vtkIdType currentId = cells->vtkDataArray::GetComponent(j, 0);

				double currentPoint[3] = { 0, 0, 0 };
				currentPoint[0] = centrePoints->GetPoint(currentId)[0];
				currentPoint[1] = centrePoints->GetPoint(currentId)[1];
				currentPoint[2] = centrePoints->GetPoint(currentId)[2];

				list<tuple<int, double, vector<double>, vector<double>>> differenceVectorList;
				for (unordered_map<int, vtkSmartPointer<vtkUnstructuredGrid>>::const_iterator iter = containerMaterialSearchPath.begin();
					iter != containerMaterialSearchPath.end(); iter++) {
					vtkIdType closestPosPath = iter->second->FindPoint(currentPoint);
					vtkIdType closestIdPoint = containerMaterialCellIDs.find(iter->first)->second->GetId(closestPosPath);

					double distance = calcDistancetoPointinEpi(data, currentId, closestIdPoint, distanceEndoEpi);
					double *materialDiffArray = data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetTuple(closestIdPoint);

					vector<double> materialDiffVector = { materialDiffArray[0], materialDiffArray[1], materialDiffArray[2] };
					tuple<int, double, vector<double>, vector<double>> tupleDiffVector(iter->first, distance, materialDiffVector, materialDiffVector);
					differenceVectorList.push_back(tupleDiffVector);
				}

				for (unordered_map<int, vtkSmartPointer<vtkUnstructuredGrid>>::const_iterator iter = containerMaterialSearchPathEndo.begin();
					iter != containerMaterialSearchPathEndo.end(); iter++) {
					vtkIdType closestPosPath = iter->second->FindPoint(currentPoint);
					vtkIdType closestIdPoint = containerMaterialCellIDsEndo.find(iter->first)->second->GetId(closestPosPath);

					double closestPoint[3] = { 0, 0, 0 };
					closestPoint[0] = centrePoints->GetPoint(closestIdPoint)[0];
					closestPoint[1] = centrePoints->GetPoint(closestIdPoint)[1];
					closestPoint[2] = centrePoints->GetPoint(closestIdPoint)[2];

					double distance = sqrt(abs(vtkMath::Distance2BetweenPoints(currentPoint, closestPoint)));
					double *materialDiffArray = data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->GetTuple(closestIdPoint);

					vector<double> materialDiffVector = { materialDiffArray[0], materialDiffArray[1], materialDiffArray[2] };
					tuple<int, double, vector<double>, vector<double>> tupleDiffVector(iter->first, distance, materialDiffVector, materialDiffVector);

					differenceVectorList.push_back(tupleDiffVector);
				}

				double minDistance = std::numeric_limits<double>::max();
				double sumOfDistance = 0;

				for (list<tuple<int, double, vector<double>, vector<double>>>::iterator iter = differenceVectorList.begin();
					iter != differenceVectorList.end(); iter++) {
					if (minDistance > get<1>(*iter)) {
						minDistance = get<1>(*iter);
					}

					if (get<1>(*iter) > 0) {
						get<1>(*iter) = 1 / get<1>(*iter);
					}
					else {
						get<1>(*iter) = 1;
					}
					sumOfDistance += get<1>(*iter);
				}


				for (list<tuple<int, double, vector<double>, vector<double>>>::iterator iter = differenceVectorList.begin();
					iter != differenceVectorList.end(); iter++) {
					get<1>(*iter) = get<1>(*iter) / sumOfDistance;
				}

				vector<double> newOrientationVector = AveragingOrientation::AveragingNVectors(differenceVectorList);
				double newDiffVector[3] = { newOrientationVector.at(0), newOrientationVector.at(1), newOrientationVector.at(2) };

				data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->InsertTuple(currentId, newDiffVector);
				data.getVtkData()->GetCellData()->GetArray("StatusEndo")->SetComponent(currentId, 0, 3);
				data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->SetComponent(currentId, 0, materialVein);
			}
		}
	}
	else {
		for (vtkIdType i = 0; i < cells->GetNumberOfTuples(); i++) {
			vtkIdType currentId = cells->vtkDataArray::GetComponent(i, 0);
			vtkSmartPointer<vtkIdList> cellNeighbours = getCellNeighbours(data, currentId);

			if (cells->vtkDataArray::GetComponent(i, 2) == 1) { // Epi
				int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(cellNeighbours->GetId(0), 0);
				data.getVtkData()->GetCellData()->GetArray("Material")->SetComponent(currentId, 0, material);
			}


			if (cells->vtkDataArray::GetComponent(i, 3) == 1) { // Endo
				int material = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(cellNeighbours->GetId(0), 0);
				data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->SetComponent(currentId, 0, material);
			}
		}
	}

	return cellIds;
} // Methods::growVeinInRegionOneSurface

/*! Grow the fiber orientation from the vein base to the top of the vein. The fiber grow only in define region.
 \param data Pointer to the orginal mesh.
 \param path Path with fiber orientation.
 \param inMaterial1 First tissue class in the fieber orientation could be grown.
 \param inMaterial2 Second tissue class in the fieber orientation could be grown.
 \param materialVein Tissue class of the vein oufire.
 \param distancePoint Outer point of the vein .
 \param regionList Region in that the fiber could be grown.
 \retrun A vtkIdList with the cell-ids.
 */
vtkSmartPointer<vtkIdList> Methods::growVeinInRegion(DataFormat &data, Material::Mat inMaterial1, Material::Mat inMaterial2, Material::Mat materialVein, vector<double> distancePoint, vtkSmartPointer<vtkIdList> regionList) {

	vtkIdType seedPointID = Methods::findClosedPointIdOnPathWithoutOrientationInMaterial(data, distancePoint, regionList, inMaterial1, inMaterial2);

	vtkSmartPointer<vtkDoubleArray> cells = vtkSmartPointer<vtkDoubleArray>::New();

	cells->SetNumberOfComponents(2);
	cells->SetComponentName(0, "CellID");
	cells->SetComponentName(1, "distance to centre point");


	set<vtkIdType> viewedCellList;

	bool end = false;

	set<vtkIdType> cellIdsWithoutOrienation;
	set<vtkIdType> cellIdList;

	cellIdList.insert(seedPointID);

	if (data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(seedPointID, 0) == 0) {
		double cellsArray[2] = { 0, 0 };
		cellsArray[0] = seedPointID;
		cells->InsertNextTuple(cellsArray);
		cellIdsWithoutOrienation.insert(seedPointID);
	}

	viewedCellList.insert(seedPointID);
	vtkSmartPointer<vtkIdList> cellIdsWithOrienation = vtkSmartPointer<vtkIdList>::New();
	set<vtkIdType> viewedCellIdsWithOrienation;

	set<vtkIdType> region;

	for (vtkIdType i = 0; i < regionList->GetNumberOfIds(); i++) {
		region.insert(regionList->GetId(i));
	}

	while (!end) {
		set<vtkIdType> newcellIdList;

		for (set<vtkIdType>::iterator iter = cellIdList.begin(); iter != cellIdList.end(); iter++) {
			vtkIdType currentCellIds = *iter;
			vtkSmartPointer<vtkIdList> cellPointIds = vtkSmartPointer<vtkIdList>::New();

			data.getVtkData()->GetCellPoints(currentCellIds, cellPointIds);

			for (vtkIdType i = 0; i < cellPointIds->GetNumberOfIds(); i++) {
				vtkSmartPointer<vtkIdList> CellPointIdList = vtkSmartPointer<vtkIdList>::New();
				CellPointIdList->InsertNextId(cellPointIds->GetId(i));

				vtkSmartPointer<vtkIdList> tempNeighborCellIds = vtkSmartPointer<vtkIdList>::New();
				data.getVtkData()->GetCellNeighbors(currentCellIds, CellPointIdList, tempNeighborCellIds);

				for (vtkIdType j = 0; j < tempNeighborCellIds->GetNumberOfIds(); j++) {
					if ((viewedCellList.find(tempNeighborCellIds->GetId(j)) == viewedCellList.end()) &&
						(region.find(tempNeighborCellIds->GetId(j)) != region.end())) {
						int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(tempNeighborCellIds->GetId(j), 0);

						if (Material::isInLeft(material) || Material::isInRight(material)) {
							if (data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(tempNeighborCellIds->GetId(j), 0) == 0) {
								double cellsArray[2] = { 0, 0 };
								cellsArray[0] = tempNeighborCellIds->GetId(j);
								cells->InsertNextTuple(cellsArray);
								newcellIdList.insert(tempNeighborCellIds->GetId(j));
								cellIdsWithoutOrienation.insert(tempNeighborCellIds->GetId(j));
							}
							else if ((material != Material::Vorhof_links) &&
								(material != Material::Vorhof_links_Endo) &&
								(material != Material::Vorhof_rechts) &&
								(data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetComponent(tempNeighborCellIds->GetId(j), 0) +
									data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetComponent(tempNeighborCellIds->GetId(j), 1) +
									data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetComponent(tempNeighborCellIds->GetId(j), 2) != 0)
								&& (viewedCellIdsWithOrienation.find(tempNeighborCellIds->GetId(j)) == viewedCellIdsWithOrienation.end())
								) {
								cellIdsWithOrienation->InsertNextId(tempNeighborCellIds->GetId(j));
								viewedCellIdsWithOrienation.insert(tempNeighborCellIds->GetId(j));
							}
						}
					}
					viewedCellList.insert(tempNeighborCellIds->GetId(j));
				}
			}
		}
		if (newcellIdList.size() == 0) {
			end = true;
		}
		cellIdList.clear();
		cellIdList = newcellIdList;
	}

	vtkSmartPointer<vtkIdList> returnIdList = vtkSmartPointer<vtkIdList>::New();

	if (cellIdsWithOrienation->GetNumberOfIds() > 0) {
		cellIdsWithOrienation = Methods::growPathThroughLayerInRegion(data, cellIdsWithOrienation, cellIdsWithoutOrienation, 4, inMaterial1, inMaterial2);

		returnIdList = growPathtoPoint(data, cellIdsWithOrienation, distancePoint, inMaterial1, inMaterial2);
		setMaterial(data, returnIdList, materialVein, inMaterial1, inMaterial2);
	}
	return returnIdList;
} // Methods::growVeinInRegion

/*! Grow the fiber orientation from a path with a fix width.
 \param data Pointer to the orginal mesh.
 \param path Path with fiber orientation.
 \param width The width which the path would be grown.
 \retrun A vtkIdList with the cell-ids.
 */
vtkSmartPointer<vtkIdList> Methods::growPath(DataFormat &data, vtkSmartPointer<vtkIdList> path, double width) {
	vtkSmartPointer<vtkPoints> centrePoints = data.getCentrePoints()->GetPoints();
	vtkSmartPointer<vtkUnstructuredGrid> centrePointsUGrid = data.getCentrePoints();
	vtkSmartPointer<vtkPointLocator> centrePointLocator = vtkSmartPointer<vtkPointLocator>::New();

	centrePointLocator->SetDataSet(centrePointsUGrid);
	centrePointLocator->BuildLocator();

	vtkSmartPointer<vtkDoubleArray> tempDiffVec = vtkDoubleArray::SafeDownCast(data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector"));

	vtkSmartPointer<vtkDoubleArray> tempArray = vtkSmartPointer<vtkDoubleArray>::New();
	tempArray->SetNumberOfComponents(3);
	tempArray->SetNumberOfTuples(centrePoints->GetNumberOfPoints());
	tempArray->FillComponent(0, 0);
	tempArray->FillComponent(1, 0);
	tempArray->FillComponent(2, 0);

	set<vtkIdType> changedId;
	set<vtkIdType> searchPath;

	for (vtkIdType i = 0; i < path->GetNumberOfIds(); i++) {
		searchPath.insert(path->GetId(i));
	}


	for (set<vtkIdType>::iterator iter = searchPath.begin(); iter != searchPath.end(); iter++) {
		vtkIdType parentId = *iter;
		vtkSmartPointer<vtkIdList> neighbors = vtkSmartPointer<vtkIdList>::New();

		centrePointLocator->FindPointsWithinRadius(width, centrePoints->GetPoint(parentId), neighbors);

		for (vtkIdType j = 0; j < neighbors->GetNumberOfIds(); j++) {
			vtkIdType currentId = neighbors->GetId(j);
			if (searchPath.find(currentId) == searchPath.end()) {
				tempArray->InsertComponent(currentId, 0, (tempDiffVec->vtkDataArray::GetComponent(parentId, 0) + tempArray->GetComponent(currentId, 0)));

				tempArray->InsertComponent(currentId, 1, (tempDiffVec->vtkDataArray::GetComponent(parentId, 1) + tempArray->GetComponent(currentId, 1)));

				tempArray->InsertComponent(currentId, 2, (tempDiffVec->vtkDataArray::GetComponent(parentId, 2) + tempArray->GetComponent(currentId, 2)));
				changedId.insert(currentId);
			}
		}
	}

	vtkSmartPointer<vtkIdList> returnIds = vtkSmartPointer<vtkIdList>::New();

	for (set<vtkIdType>::iterator iter = changedId.begin(); iter != changedId.end(); iter++) {
		int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(*iter, 0);
		if (Material::isInLeft(material) || Material::isInRight(material)) {
			double *tempPoint = tempArray->GetTuple(*iter);
			vtkMath::Normalize(tempPoint);
			data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->InsertTuple(*iter, tempPoint);
			data.getVtkData()->GetCellData()->GetArray("Status")->SetComponent(*iter, 0, 3);
			returnIds->InsertNextId(*iter);
		}
	}
	vtkSmartPointer<vtkIdList> returnList = Methods::add2Paths(returnIds, path);
	return returnList;
} // Methods::growPath

/*! Grow the fiber orientation from a path with a fix width in the left atrial.
 \param data Pointer to the orginal mesh.
 \param path Path with fiber orientation.
 \param width The width which the path would be grown.
 \retrun A vtkIdList with the cell-ids.
 */
vtkSmartPointer<vtkIdList> Methods::growPathInLeft(DataFormat &data, vtkSmartPointer<vtkIdList> path, double width) {
	vtkSmartPointer<vtkPoints> centrePoints = data.getCentrePoints()->GetPoints();
	vtkSmartPointer<vtkUnstructuredGrid> centrePointsUGrid = data.getCentrePoints();
	vtkSmartPointer<vtkPointLocator> centrePointLocator = vtkSmartPointer<vtkPointLocator>::New();

	centrePointLocator->SetDataSet(centrePointsUGrid);
	centrePointLocator->BuildLocator();

	vtkSmartPointer<vtkDoubleArray> tempDiffVec = vtkDoubleArray::SafeDownCast(data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector"));

	vtkSmartPointer<vtkDoubleArray> tempArray = vtkSmartPointer<vtkDoubleArray>::New();
	tempArray->SetNumberOfComponents(3);
	tempArray->SetNumberOfTuples(centrePoints->GetNumberOfPoints());
	tempArray->FillComponent(0, 0);
	tempArray->FillComponent(1, 0);
	tempArray->FillComponent(2, 0);

	set<vtkIdType> changedId;
	set<vtkIdType> searchPath;

	for (vtkIdType i = 0; i < path->GetNumberOfIds(); i++) {
		searchPath.insert(path->GetId(i));
	}


	for (set<vtkIdType>::iterator iter = searchPath.begin(); iter != searchPath.end(); iter++) {
		vtkIdType parentId = *iter;
		vtkSmartPointer<vtkIdList> neighbors = vtkSmartPointer<vtkIdList>::New();

		centrePointLocator->FindPointsWithinRadius(width, centrePoints->GetPoint(parentId), neighbors);

		for (vtkIdType j = 0; j < neighbors->GetNumberOfIds(); j++) {
			vtkIdType currentId = neighbors->GetId(j);
			if (searchPath.find(currentId) == searchPath.end()) {
				tempArray->InsertComponent(currentId, 0, (tempDiffVec->vtkDataArray::GetComponent(parentId, 0) + tempArray->GetComponent(currentId, 0)));

				tempArray->InsertComponent(currentId, 1, (tempDiffVec->vtkDataArray::GetComponent(parentId, 1) + tempArray->GetComponent(currentId, 1)));

				tempArray->InsertComponent(currentId, 2, (tempDiffVec->vtkDataArray::GetComponent(parentId, 2) + tempArray->GetComponent(currentId, 2)));
				changedId.insert(currentId);
			}
		}
	}

	vtkSmartPointer<vtkIdList> returnIds = vtkSmartPointer<vtkIdList>::New();

	for (set<vtkIdType>::iterator iter = changedId.begin(); iter != changedId.end(); iter++) {
		if (!data.getDoubleLayer()) {
			int material = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(*iter, 0);
			if (Material::isInLeft(material)) {
				double *tempPoint = tempArray->GetTuple(*iter);
				vtkMath::Normalize(tempPoint);
				data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->InsertTuple(*iter, tempPoint);
				data.getVtkData()->GetCellData()->GetArray("StatusEndo")->SetComponent(*iter, 0, 3);
				returnIds->InsertNextId(*iter);
			}
		}
		int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(*iter, 0);
		if (Material::isInLeft(material)) {
			double *tempPoint = tempArray->GetTuple(*iter);
			vtkMath::Normalize(tempPoint);
			data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->InsertTuple(*iter, tempPoint);
			data.getVtkData()->GetCellData()->GetArray("Status")->SetComponent(*iter, 0, 3);
			returnIds->InsertNextId(*iter);
		}
	}
	vtkSmartPointer<vtkIdList> returnList = Methods::add2Paths(returnIds, path);
	return returnList;
} // Methods::growPathInLeft

/*! Grow the fiber orientation from a path with a fix width in the right atrial.
 \param data Pointer to the orginal mesh.
 \param path Path with fiber orientation.
 \param width The width which the path would be grown.
 \retrun A vtkIdList with the cell-ids.
 */
vtkSmartPointer<vtkIdList> Methods::growPathInRight(DataFormat &data, vtkSmartPointer<vtkIdList> path, double width) {
	vtkSmartPointer<vtkPoints> centrePoints = data.getCentrePoints()->GetPoints();
	vtkSmartPointer<vtkUnstructuredGrid> centrePointsUGrid = data.getCentrePoints();
	vtkSmartPointer<vtkPointLocator> centrePointLocator = vtkSmartPointer<vtkPointLocator>::New();

	centrePointLocator->SetDataSet(centrePointsUGrid);
	centrePointLocator->BuildLocator();

	vtkSmartPointer<vtkDoubleArray> tempDiffVec = vtkDoubleArray::SafeDownCast(data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector"));

	vtkSmartPointer<vtkDoubleArray> tempArray = vtkSmartPointer<vtkDoubleArray>::New();
	tempArray->SetNumberOfComponents(3);
	tempArray->SetNumberOfTuples(centrePoints->GetNumberOfPoints());
	tempArray->FillComponent(0, 0);
	tempArray->FillComponent(1, 0);
	tempArray->FillComponent(2, 0);

	set<vtkIdType> changedId;
	set<vtkIdType> searchPath;

	for (vtkIdType i = 0; i < path->GetNumberOfIds(); i++) {
		searchPath.insert(path->GetId(i));
	}


	for (set<vtkIdType>::iterator iter = searchPath.begin(); iter != searchPath.end(); iter++) {
		vtkIdType parentId = *iter;
		vtkSmartPointer<vtkIdList> neighbors = vtkSmartPointer<vtkIdList>::New();

		centrePointLocator->FindPointsWithinRadius(width, centrePoints->GetPoint(parentId), neighbors);

		for (vtkIdType j = 0; j < neighbors->GetNumberOfIds(); j++) {
			vtkIdType currentId = neighbors->GetId(j);
			if (searchPath.find(currentId) == searchPath.end()) {
				tempArray->InsertComponent(currentId, 0, (tempDiffVec->vtkDataArray::GetComponent(parentId, 0) + tempArray->GetComponent(currentId, 0)));

				tempArray->InsertComponent(currentId, 1, (tempDiffVec->vtkDataArray::GetComponent(parentId, 1) + tempArray->GetComponent(currentId, 1)));

				tempArray->InsertComponent(currentId, 2, (tempDiffVec->vtkDataArray::GetComponent(parentId, 2) + tempArray->GetComponent(currentId, 2)));
				changedId.insert(currentId);
			}
		}
	}

	vtkSmartPointer<vtkIdList> returnIds = vtkSmartPointer<vtkIdList>::New();

	for (set<vtkIdType>::iterator iter = changedId.begin(); iter != changedId.end(); iter++) {
		if (!data.getDoubleLayer()) {
			int material = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(*iter, 0);
			if (Material::isInRight(material)) {
				double *tempPoint = tempArray->GetTuple(*iter);
				vtkMath::Normalize(tempPoint);
				data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->InsertTuple(*iter, tempPoint);
				data.getVtkData()->GetCellData()->GetArray("StatusEndo")->SetComponent(*iter, 0, 3);
				returnIds->InsertNextId(*iter);
			}
		}


		int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(*iter, 0);
		if (Material::isInRight(material)) {
			double *tempPoint = tempArray->GetTuple(*iter);
			vtkMath::Normalize(tempPoint);
			data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->InsertTuple(*iter, tempPoint);
			data.getVtkData()->GetCellData()->GetArray("Status")->SetComponent(*iter, 0, 3);
			returnIds->InsertNextId(*iter);
		}
	}
	vtkSmartPointer<vtkIdList> returnList = Methods::add2Paths(returnIds, path);
	return returnList;
} // Methods::growPathInRight

/*! Grow the fiber orientation from a path form one material with a fix width in anther material.
 \param data Pointer to the orginal mesh.
 \param path Path with fiber orientation.
 \param width The width which the path would be grown.
 \param outMaterial Tissue class from that the fieber orientation would be grown.
 \param inMaterial Tissue class in the fieber orientation could be grown.
 \retrun A vtkIdList with the cell-ids.
 */
vtkSmartPointer<vtkIdList> Methods::growPath(DataFormat &data, vtkSmartPointer<vtkIdList> path, double width, Material::Mat outMaterial, Material::Mat inMaterial) {
	return growPath(data, path, width, outMaterial, inMaterial, inMaterial);
}

/*! Grow the fiber orientation from a path form one material with a fix width in two other material.
 \param data Pointer to the orginal mesh.
 \param path Path with fiber orientation.
 \param width The width which the path would be grown.
 \param outMaterial Tissue class from that the fieber orientation would be grown.
 \param inMaterial1 First tissue class in the fieber orientation could be grown.
 \param inMaterial2 Second tissue class in the fieber orientation could be grown.
 \retrun A vtkIdList with the cell-ids.
 */
vtkSmartPointer<vtkIdList> Methods::growPath(DataFormat &data, vtkSmartPointer<vtkIdList> path, double width, Material::Mat outMaterial, Material::Mat inMaterial1, Material::Mat inMaterial2) {
	vtkSmartPointer<vtkPoints> centrePoints = data.getCentrePoints()->GetPoints();
	vtkSmartPointer<vtkUnstructuredGrid> centrePointsUGrid = data.getCentrePoints();
	vtkSmartPointer<vtkPointLocator> centrePointLocator = vtkSmartPointer<vtkPointLocator>::New();

	centrePointLocator->SetDataSet(centrePointsUGrid);
	centrePointLocator->BuildLocator();

	vtkSmartPointer<vtkIdList> growedPathEndo = vtkSmartPointer<vtkIdList>::New();

	if (!data.getDoubleLayer()) {
		if ((Material::isInLeftEndo(inMaterial1) || Material::isInLeftEndo(inMaterial2)) &&
			!Material::isInLeftEndo(outMaterial)) {
			double widthEndo = sqrt(pow(2.2, 2) + pow(width, 2));
			if ((Material::isInLeftEndo(inMaterial1) && Material::isInLeftEndo(inMaterial2))) {
				growedPathEndo = growPathEndoInSurface(data, path, widthEndo, inMaterial1, inMaterial2);
				inMaterial1 = Material::undef;
				inMaterial2 = Material::undef;
			}
			else if (Material::isInLeftEndo(inMaterial1) && !Material::isInLeftEndo(inMaterial2)) {
				growedPathEndo = growPathEndoInSurface(data, path, widthEndo, inMaterial1, inMaterial1);
				inMaterial1 = Material::undef;
			}
			else if (!Material::isInLeftEndo(inMaterial1) && Material::isInLeftEndo(inMaterial2)) {
				growedPathEndo = growPathEndoInSurface(data, path, widthEndo, inMaterial2, inMaterial2);
				inMaterial2 = Material::undef;
			}
		}
	}

	vtkSmartPointer<vtkDoubleArray> tempDiffVec = vtkDoubleArray::SafeDownCast(data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector"));

	vtkSmartPointer<vtkDoubleArray> tempArray = vtkSmartPointer<vtkDoubleArray>::New();
	tempArray->SetNumberOfComponents(3);
	tempArray->SetNumberOfTuples(centrePoints->GetNumberOfPoints());
	tempArray->FillComponent(0, 0);
	tempArray->FillComponent(1, 0);
	tempArray->FillComponent(2, 0);

	set<vtkIdType> changedId;

	for (vtkIdType j = 0; j < path->GetNumberOfIds(); j++) {
		vtkIdType parentId = path->GetId(j);

		vtkSmartPointer<vtkIdList> neighbors = vtkSmartPointer<vtkIdList>::New();
		centrePointLocator->FindPointsWithinRadius(width, centrePoints->GetPoint(parentId), neighbors);

		for (vtkIdType k = 0; k < neighbors->GetNumberOfIds(); k++) {
			vtkIdType currentId = neighbors->GetId(k);

			if ((tempDiffVec->GetComponent(parentId, 0) != 0) ||
				(tempDiffVec->GetComponent(parentId, 1) != 0) || (tempDiffVec->GetComponent(parentId, 2) != 0)) {
				tempArray->InsertComponent(currentId, 0, (tempDiffVec->GetComponent(parentId, 0) + tempArray->GetComponent(currentId, 0)));
				tempArray->InsertComponent(currentId, 1, (tempDiffVec->GetComponent(parentId, 1) + tempArray->GetComponent(currentId, 1)));
				tempArray->InsertComponent(currentId, 2, (tempDiffVec->GetComponent(parentId, 2) + tempArray->GetComponent(currentId, 2)));
				changedId.insert(currentId);
			}
		}
	}
	vtkSmartPointer<vtkIdList> returnIds = vtkSmartPointer<vtkIdList>::New();

	for (set<vtkIdType>::iterator iter = changedId.begin(); iter != changedId.end(); iter++) {
		if (!data.getDoubleLayer()) {
			int material = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(*iter, 0);
			if ((data.getVtkData()->GetCellData()->GetArray("StatusEndo")->GetComponent(*iter, 0) == 0) &&
				((material == inMaterial1) || (material == inMaterial2))) {
				double *tempPoint = tempArray->GetTuple(*iter);
				vtkMath::Normalize(tempPoint);
				data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->InsertTuple(*iter, tempPoint);
				data.getVtkData()->GetCellData()->GetArray("StatusEndo")->SetComponent(*iter, 0, 3);
				returnIds->InsertUniqueId(*iter);
			}
		}

		int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(*iter, 0);
		if ((data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(*iter, 0) == 0) &&
			((material == inMaterial1) || (material == inMaterial2))) {
			double *tempPoint = tempArray->GetTuple(*iter);
			vtkMath::Normalize(tempPoint);
			data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->InsertTuple(*iter, tempPoint);
			data.getVtkData()->GetCellData()->GetArray("Status")->SetComponent(*iter, 0, 3);
			returnIds->InsertUniqueId(*iter);
		}
	}

	vtkSmartPointer<vtkIdList> returnList = Methods::add2Paths(returnIds, path);
	if (!data.getDoubleLayer()) {
		returnList = Methods::add2Paths(returnList, growedPathEndo);
	}

	return returnList;
} // Methods::growPath

/*! Grow the fiber orientation from a path form one material with a fix width in two other material.
 \param data Pointer to the orginal mesh.
 \param path Path with fiber orientation.
 \param width The width which the path would be grown.
 \param outMaterial Tissue class from that the fieber orientation would be grown.
 \param inMaterial1 First tissue class in the fieber orientation could be grown.
 \param inMaterial2 Second tissue class in the fieber orientation could be grown.
 \retrun A vtkIdList with the cell-ids.
 */
vtkSmartPointer<vtkIdList> Methods::growPathThroughLayerInRegion(DataFormat &data, vtkSmartPointer<vtkIdList> path, set<vtkIdType> region, double width, Material::Mat inMaterial1, Material::Mat inMaterial2) {
	vtkSmartPointer<vtkPoints> centrePoints = data.getCentrePoints()->GetPoints();
	vtkSmartPointer<vtkUnstructuredGrid> centrePointsUGrid = data.getCentrePoints();
	vtkSmartPointer<vtkPointLocator> centrePointLocator = vtkSmartPointer<vtkPointLocator>::New();

	centrePointLocator->SetDataSet(centrePointsUGrid);
	centrePointLocator->BuildLocator();

	vtkSmartPointer<vtkDoubleArray> tempDiffVec = vtkDoubleArray::SafeDownCast(data.getVtkData()->GetCellData()->GetArray("DifferenceVector"));

	vtkSmartPointer<vtkDoubleArray> tempArray = vtkSmartPointer<vtkDoubleArray>::New();
	tempArray->SetNumberOfComponents(3);
	tempArray->SetNumberOfTuples(centrePoints->GetNumberOfPoints());
	tempArray->FillComponent(0, 0);
	tempArray->FillComponent(1, 0);
	tempArray->FillComponent(2, 0);

	set<vtkIdType> changedId;

	for (vtkIdType j = 0; j < path->GetNumberOfIds(); j++) {
		vtkIdType parentId = path->GetId(j);

		vtkSmartPointer<vtkIdList> neighbors = vtkSmartPointer<vtkIdList>::New();
		centrePointLocator->FindPointsWithinRadius(width, centrePoints->GetPoint(parentId), neighbors);

		for (vtkIdType k = 0; k < neighbors->GetNumberOfIds(); k++) {
			vtkIdType currentId = neighbors->GetId(k);


			if (((tempDiffVec->GetComponent(parentId, 0) != 0) ||
				(tempDiffVec->GetComponent(parentId, 1) != 0) || (tempDiffVec->GetComponent(parentId, 2) != 0)) &&
				(region.find(currentId) != region.end())) {
				tempArray->InsertComponent(currentId, 0, (tempDiffVec->GetComponent(parentId, 0) + tempArray->GetComponent(currentId, 0)));
				tempArray->InsertComponent(currentId, 1, (tempDiffVec->GetComponent(parentId, 1) + tempArray->GetComponent(currentId, 1)));
				tempArray->InsertComponent(currentId, 2, (tempDiffVec->GetComponent(parentId, 2) + tempArray->GetComponent(currentId, 2)));
				changedId.insert(currentId);
			}
		}
	}
	vtkSmartPointer<vtkIdList> returnIds = vtkSmartPointer<vtkIdList>::New();

	for (set<vtkIdType>::iterator iter = changedId.begin(); iter != changedId.end(); iter++) {
		int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(*iter, 0);
		if ((data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(*iter, 0) == 0) &&
			((material == inMaterial1) || (material == inMaterial2))) {
			double *tempPoint = tempArray->GetTuple(*iter);
			vtkMath::Normalize(tempPoint);
			data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->InsertTuple(*iter, tempPoint);
			data.getVtkData()->GetCellData()->GetArray("Status")->SetComponent(*iter, 0, 3);
			returnIds->InsertUniqueId(*iter);
		}
	}

	vtkSmartPointer<vtkIdList> returnList = Methods::add2Paths(returnIds, path);

	return returnList;
} // Methods::growPathThroughLayerInRegion

/*! Grow the fiber orientation from a path with a fix width in two other material in the endocard in a surface mesh.
 \param data Pointer to the orginal mesh.
 \param path Path with fiber orientation.
 \param width The width which the path would be grown.
 \param inMaterial1 First tissue class in the fieber orientation could be grown.
 \param inMaterial2 Second tissue class in the fieber orientation could be grown.
 \retrun A vtkIdList with the cell-ids.
 */
vtkSmartPointer<vtkIdList> Methods::growPathEndoInSurface(DataFormat &data, vtkSmartPointer<vtkIdList> path, double width, Material::Mat inMaterial1, Material::Mat inMaterial2) {
	vtkSmartPointer<vtkPoints> centrePoints = data.getCentrePoints()->GetPoints();
	vtkSmartPointer<vtkUnstructuredGrid> centrePointsUGrid = data.getCentrePoints();
	vtkSmartPointer<vtkPointLocator> centrePointLocator = vtkSmartPointer<vtkPointLocator>::New();

	centrePointLocator->SetDataSet(centrePointsUGrid);
	centrePointLocator->BuildLocator();

	vtkSmartPointer<vtkDoubleArray> tempDiffVec = vtkDoubleArray::SafeDownCast(data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector"));

	vtkSmartPointer<vtkDoubleArray> tempArray = vtkSmartPointer<vtkDoubleArray>::New();
	tempArray->SetNumberOfComponents(3);
	tempArray->SetNumberOfTuples(centrePoints->GetNumberOfPoints());
	tempArray->FillComponent(0, 0);
	tempArray->FillComponent(1, 0);
	tempArray->FillComponent(2, 0);

	set<vtkIdType> changedId;
	set<vtkIdType> searchPath;

	for (vtkIdType i = 0; i < path->GetNumberOfIds(); i++) {
		searchPath.insert(path->GetId(i));
	}


	for (vtkIdType j = 0; j < path->GetNumberOfIds(); j++) {
		vtkIdType parentId = path->GetId(j);

		double tempVec[3] = { 0, 0, 0 };

		if ((data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(parentId, 0) != 1) ||
			(!data.getDoubleLayer() &&
			(data.getVtkData()->GetCellData()->GetArray("StatusEndo")->GetComponent(parentId, 0) != 1))) {
			if (j < path->GetNumberOfIds() - 1) {
				vtkIdType nextPointId = path->GetId(j + 1);
				double currentPoint[3] = { centrePoints->GetPoint(parentId)[0], centrePoints->GetPoint(parentId)[1], centrePoints->GetPoint(parentId)[2] };
				double nextPoint[3] = { centrePoints->GetPoint(nextPointId)[0], centrePoints->GetPoint(nextPointId)[1], centrePoints->GetPoint(nextPointId)[2] };

				vtkMath::Subtract(nextPoint, currentPoint, tempVec);
				vtkMath::Normalize(tempVec);
			}
			else if (j == path->GetNumberOfIds() - 1) {
				vtkIdType nextPointId = path->GetId(j - 1);
				double currentPoint[3] = { centrePoints->GetPoint(parentId)[0], centrePoints->GetPoint(parentId)[1], centrePoints->GetPoint(parentId)[2] };
				double lastPoint[3] = { centrePoints->GetPoint(nextPointId)[0], centrePoints->GetPoint(nextPointId)[1], centrePoints->GetPoint(nextPointId)[2] };

				vtkMath::Subtract(currentPoint, lastPoint, tempVec);
				vtkMath::Normalize(tempVec);
			}
		}

		vtkSmartPointer<vtkIdList> neighbors = vtkSmartPointer<vtkIdList>::New();
		centrePointLocator->FindPointsWithinRadius(width, centrePoints->GetPoint(parentId), neighbors);

		for (vtkIdType k = 0; k < neighbors->GetNumberOfIds(); k++) {
			vtkIdType currentId = neighbors->GetId(k);

			if (!data.getDoubleLayer()) {
				if (data.getVtkData()->GetCellData()->GetArray("StatusEndo")->GetComponent(parentId, 0) == 1) {
					if ((data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(currentId, 0) == 0) ||
						(data.getVtkData()->GetCellData()->GetArray("StatusEndo")->GetComponent(currentId, 0) == 0)) {
						tempArray->InsertComponent(currentId, 0, (tempDiffVec->vtkDataArray::GetComponent(parentId, 0) +
							tempArray->GetComponent(currentId, 0)));
						tempArray->InsertComponent(currentId, 1, (tempDiffVec->vtkDataArray::GetComponent(parentId, 1) +
							tempArray->GetComponent(currentId, 1)));
						tempArray->InsertComponent(currentId, 2, (tempDiffVec->vtkDataArray::GetComponent(parentId, 2) +
							tempArray->GetComponent(currentId, 2)));
						changedId.insert(currentId);
					}
				}
				else if ((data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(currentId, 0) == 0) ||
					(data.getVtkData()->GetCellData()->GetArray("StatusEndo")->GetComponent(currentId, 0) == 0)) {
					tempArray->InsertComponent(currentId, 0, tempVec[0] + tempArray->GetComponent(currentId, 0));
					tempArray->InsertComponent(currentId, 1, tempVec[1] + tempArray->GetComponent(currentId, 1));
					tempArray->InsertComponent(currentId, 2, tempVec[2] + tempArray->GetComponent(currentId, 2));
					changedId.insert(currentId);
				}
			}
			if (data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(parentId, 0) == 1) {
				if ((data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(currentId, 0) == 0) ||
					(!data.getDoubleLayer() &&
					(data.getVtkData()->GetCellData()->GetArray("StatusEndo")->GetComponent(currentId, 0) == 0))) {
					tempArray->InsertComponent(currentId, 0, (tempDiffVec->vtkDataArray::GetComponent(parentId, 0) +
						tempArray->GetComponent(currentId, 0)));
					tempArray->InsertComponent(currentId, 1, (tempDiffVec->vtkDataArray::GetComponent(parentId, 1) +
						tempArray->GetComponent(currentId, 1)));
					tempArray->InsertComponent(currentId, 2, (tempDiffVec->vtkDataArray::GetComponent(parentId, 2) +
						tempArray->GetComponent(currentId, 2)));
					changedId.insert(currentId);
				}
			}
			else if ((data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(currentId, 0) == 0) ||
				(!data.getDoubleLayer() &&
				(data.getVtkData()->GetCellData()->GetArray("StatusEndo")->GetComponent(currentId, 0) == 0))) {
				tempArray->InsertComponent(currentId, 0, tempVec[0] + tempArray->GetComponent(currentId, 0));
				tempArray->InsertComponent(currentId, 1, tempVec[1] + tempArray->GetComponent(currentId, 1));
				tempArray->InsertComponent(currentId, 2, tempVec[2] + tempArray->GetComponent(currentId, 2));
				changedId.insert(currentId);
			}
		}
	}
	vtkSmartPointer<vtkIdList> returnIds = vtkSmartPointer<vtkIdList>::New();

	for (set<vtkIdType>::iterator iter = changedId.begin(); iter != changedId.end(); iter++) {
		if (!data.getDoubleLayer()) {
			int material = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(*iter, 0);
			if ((data.getVtkData()->GetCellData()->GetArray("StatusEndo")->GetComponent(*iter, 0) == 0) &&
				((material == inMaterial1) || (material == inMaterial2))) {
				double *tempPoint = tempArray->GetTuple(*iter);
				vtkMath::Normalize(tempPoint);
				data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->InsertTuple(*iter, tempPoint);
				data.getVtkData()->GetCellData()->GetArray("StatusEndo")->SetComponent(*iter, 0, 3);
				returnIds->InsertUniqueId(*iter);
			}
		}

		int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(*iter, 0);
		if ((data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(*iter, 0) == 0) &&
			((material == inMaterial1) || (material == inMaterial2))) {
			double *tempPoint = tempArray->GetTuple(*iter);
			vtkMath::Normalize(tempPoint);
			data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->InsertTuple(*iter, tempPoint);
			data.getVtkData()->GetCellData()->GetArray("Status")->SetComponent(*iter, 0, 3);
			returnIds->InsertUniqueId(*iter);
		}
	}

	vtkSmartPointer<vtkIdList> returnList = Methods::add2Paths(returnIds, path);
	return returnList;
} // Methods::growPathEndoInSurface

/*! Get the cell id in a radius arround a path in a material.
 \param data Pointer to the orginal mesh.
 \param path Path with fiber orientation.
 \param width The radius in which the cell would be return.
 \param inMaterial Tissue class in that the cells must have.
 \retrun A vtkIdList with the cell-ids.
 */
vtkSmartPointer<vtkIdList> Methods::getCellsArroundPathwithStatusZero(DataFormat &data, vtkSmartPointer<vtkIdList> path, double width, Material::Mat inMaterial) {
	vtkSmartPointer<vtkPoints> centrePoints = data.getCentrePoints()->GetPoints();
	vtkSmartPointer<vtkUnstructuredGrid> centrePointsUGrid = data.getCentrePoints();
	vtkSmartPointer<vtkPointLocator> centrePointLocator = vtkSmartPointer<vtkPointLocator>::New();

	centrePointLocator->SetDataSet(centrePointsUGrid);
	centrePointLocator->BuildLocator();

	vtkSmartPointer<vtkIdList> returnIds = vtkSmartPointer<vtkIdList>::New();

	for (vtkIdType i = 0; i < path->GetNumberOfIds(); i++) {
		vtkIdType parentId = path->GetId(i);
		vtkSmartPointer<vtkIdList> neighbors = vtkSmartPointer<vtkIdList>::New();

		centrePointLocator->FindPointsWithinRadius(width, centrePoints->GetPoint(parentId), neighbors);

		for (vtkIdType j = 0; j < neighbors->GetNumberOfIds(); j++) {
			vtkIdType currentId = neighbors->GetId(j);

			if (!data.getDoubleLayer()) {
				if ((data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(currentId, 0) == inMaterial) &&
					(data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(currentId, 0) == 0)) {
					returnIds->InsertUniqueId(currentId);
				}
			}
			if ((data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(currentId, 0) == inMaterial) &&
				(data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(currentId, 0) == 0)) {
				returnIds->InsertUniqueId(currentId);
			}
		}
	}
	return returnIds;
} // Methods::getCellsArroundPathwithStatusZero

vtkSmartPointer<vtkIdList> Methods::getCellsArroundPath(DataFormat &data, vtkSmartPointer<vtkIdList> path, double width, Material::Mat inMaterial) {
	vtkSmartPointer<vtkPoints> centrePoints = data.getCentrePoints()->GetPoints();
	vtkSmartPointer<vtkUnstructuredGrid> centrePointsUGrid = data.getCentrePoints();
	vtkSmartPointer<vtkPointLocator> centrePointLocator = vtkSmartPointer<vtkPointLocator>::New();

	centrePointLocator->SetDataSet(centrePointsUGrid);
	centrePointLocator->BuildLocator();

	vtkSmartPointer<vtkIdList> returnIds = vtkSmartPointer<vtkIdList>::New();

	for (vtkIdType i = 0; i < path->GetNumberOfIds(); i++) {
		vtkIdType parentId = path->GetId(i);
		vtkSmartPointer<vtkIdList> neighbors = vtkSmartPointer<vtkIdList>::New();

		centrePointLocator->FindPointsWithinRadius(width, centrePoints->GetPoint(parentId), neighbors);

		for (vtkIdType j = 0; j < neighbors->GetNumberOfIds(); j++) {
			vtkIdType currentId = neighbors->GetId(j);

			if (!data.getDoubleLayer()) {
				if (data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(currentId, 0) == inMaterial) {
					returnIds->InsertUniqueId(currentId);
				}
			}
			if (data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(currentId, 0) == inMaterial) {
				returnIds->InsertUniqueId(currentId);
			}
		}
	}
	return returnIds;
} // Methods::getCellsArroundPath

// Methods::getCellsArroundPath

/*! Get the cell id in a radius arround a path in a atrium.
 \param data Pointer to the orginal mesh.
 \param path Path with fiber orientation.
 \param width The radius in which the cell would be return.
 \param atrium Atrium in that the cells lie.
 \retrun A vtkIdList with the cell-ids.
 */
vtkSmartPointer<vtkIdList> Methods::getCellsArroundPath(DataFormat &data, vtkSmartPointer<vtkIdList> path, double width, string atrium) {
	vtkSmartPointer<vtkPoints> centrePoints = data.getCentrePoints()->GetPoints();
	vtkSmartPointer<vtkUnstructuredGrid> centrePointsUGrid = data.getCentrePoints();
	vtkSmartPointer<vtkPointLocator> centrePointLocator = vtkSmartPointer<vtkPointLocator>::New();

	centrePointLocator->SetDataSet(centrePointsUGrid);
	centrePointLocator->BuildLocator();

	vtkSmartPointer<vtkIdList> returnIds = vtkSmartPointer<vtkIdList>::New();

	for (vtkIdType i = 0; i < path->GetNumberOfIds(); i++) {
		vtkIdType parentId = path->GetId(i);
		vtkSmartPointer<vtkIdList> neighbors = vtkSmartPointer<vtkIdList>::New();

		centrePointLocator->FindPointsWithinRadius(width, centrePoints->GetPoint(parentId), neighbors);

		for (vtkIdType j = 0; j < neighbors->GetNumberOfIds(); j++) {
			vtkIdType currentId = neighbors->GetId(j);
			int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(currentId, 0);
			if (atrium.compare("right") == 0) {
				if (Material::isInRight(material)) {
					returnIds->InsertUniqueId(currentId);
				}
			}
			else if (atrium.compare("left") == 0) {
				if (Material::isInLeft(material)) {
					returnIds->InsertUniqueId(currentId);
				}
			}
			else if (Material::isInLeft(material) || Material::isInRight(material)) {
				returnIds->InsertUniqueId(currentId);
			}
		}
	}
	return returnIds;
} // Methods::getCellsArroundPath


/*! Calculate the fiber angle from the mesh as normal theta and phi.
 \param data Pointer to the orginal mesh.
 */
void Methods::setThetaPhi(DataFormat &data) {
	for (vtkIdType i = 0; i < data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetNumberOfTuples(); i++) {
		int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(i, 0);
		if (Material::isInLeft(material) || Material::isInRight(material)) {
			vector<double> angles = angle(data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetComponent(i, 0), data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetComponent(i, 1), data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetComponent(i, 2));

			double angleThetaPhiArray[2] = { 0, 0 };
			angleThetaPhiArray[0] = angles.at(0);
			angleThetaPhiArray[1] = angles.at(1);

			data.getVtkData()->GetCellData()->GetArray("ThetaPhi")->InsertTuple(i, angleThetaPhiArray);
			data.getVtkData()->GetCellData()->GetArray("Theta")->InsertComponent(i, 0, angles.at(0));
			data.getVtkData()->GetCellData()->GetArray("Phi")->InsertComponent(i, 0, angles.at(1));
			data.getVtkData()->GetCellData()->GetArray("Status")->SetComponent(i, 0, 4);
		}
		if (!data.getDoubleLayer()) {
			int materialEndo = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(i, 0);
			if (Material::isInLeft(materialEndo) || Material::isInRight(materialEndo)) {
				vector<double> anglesEndo = angle(data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->GetComponent(i, 0), data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->GetComponent(i, 1), data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->GetComponent(i, 2));

				double angleThetaPhiArrayEndo[2] = { 0, 0 };
				angleThetaPhiArrayEndo[0] = anglesEndo.at(0);
				angleThetaPhiArrayEndo[1] = anglesEndo.at(1);

				data.getVtkData()->GetCellData()->GetArray("ThetaPhiEndo")->InsertTuple(i, angleThetaPhiArrayEndo);
				data.getVtkData()->GetCellData()->GetArray("ThetaEndo")->InsertComponent(i, 0, anglesEndo.at(0));
				data.getVtkData()->GetCellData()->GetArray("PhiEndo")->InsertComponent(i, 0, anglesEndo.at(1));
				data.getVtkData()->GetCellData()->GetArray("StatusEndo")->SetComponent(i, 0, 4);
			}
		}
	}
} // Methods::setThetaPhi


/*! Set the tissue form the cells.
 \param data Pointer to the orginal mesh.
 \param path Cell which material would be set.
 \param setMaterial Tissue class as that the cell would be marked.
 */
void Methods::setMaterial(DataFormat &data, vtkSmartPointer<vtkIdList> path, Material::Mat setMaterial) {
	for (vtkIdType i = 0; i < path->GetNumberOfIds(); i++) {
		if (!data.getDoubleLayer()) {
			data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->SetComponent(path->GetId(i), 0, setMaterial);
			data.getVtkData()->GetCellData()->GetArray("StatusEndo")->SetComponent(path->GetId(i), 0, 5);
		}

		data.getVtkData()->GetCellData()->GetArray("Material")->SetComponent(path->GetId(i), 0, setMaterial);
		data.getVtkData()->GetCellData()->GetArray("Status")->SetComponent(path->GetId(i), 0, 5);
	}
}

/*! Set the tissue form the cells in a dataArray.
 \param data Pointer to the orginal mesh.
 \param path Cell which material would be set.
 \param cellArray Cellarray Name which would be modify.
 \param setMaterial Tissue class as that the cell would be marked.
 */
void Methods::setMaterial(DataFormat &data, vtkSmartPointer<vtkIdList> path, const char *cellArray, int setMaterial) {
	for (vtkIdType i = 0; i < path->GetNumberOfIds(); i++) {
		data.getVtkData()->GetCellData()->GetArray(cellArray)->SetComponent(path->GetId(i), 0, setMaterial);
	}
}

/*! Set the tissue form the cells in the left atrial.
 \param data Pointer to the orginal mesh.
 \param path Cell which material would be set .
 \param setMaterial Tissue class as that the cell would be marked.
 */
void Methods::setMaterialinLeft(DataFormat &data, vtkSmartPointer<vtkIdList> path, Material::Mat setMaterial) {
	for (vtkIdType i = 0; i < path->GetNumberOfIds(); i++) {
		if (!data.getDoubleLayer()) {
			int material = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(path->GetId(i), 0);

			if (Material::isInLeft(material)) {
				data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->SetComponent(path->GetId(i), 0, setMaterial);
				data.getVtkData()->GetCellData()->GetArray("StatusEndo")->SetComponent(path->GetId(i), 0, 5);
			}
		}

		int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(path->GetId(i), 0);
		if (Material::isInLeft(material)) {
			data.getVtkData()->GetCellData()->GetArray("Material")->SetComponent(path->GetId(i), 0, setMaterial);
			data.getVtkData()->GetCellData()->GetArray("Status")->SetComponent(path->GetId(i), 0, 5);
		}
	}
}

/*! Set the tissue form the cells in the right atrial.
 \param data Pointer to the orginal mesh.
 \param path Cell which material would be set.
 \param setMaterial Tissue class as that the cell would be marked.
 */
void Methods::setMaterialinRight(DataFormat &data, vtkSmartPointer<vtkIdList> path, Material::Mat setMaterial) {
	for (vtkIdType i = 0; i < path->GetNumberOfIds(); i++) {
		if (!data.getDoubleLayer()) {
			int material = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(path->GetId(i), 0);

			if (Material::isInRight(material)) {
				data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->SetComponent(path->GetId(i), 0, setMaterial);
				data.getVtkData()->GetCellData()->GetArray("StatusEndo")->SetComponent(path->GetId(i), 0, 5);
			}
		}

		int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(path->GetId(i), 0);

		if (Material::isInRight(material)) {
			data.getVtkData()->GetCellData()->GetArray("Material")->SetComponent(path->GetId(i), 0, setMaterial);
			data.getVtkData()->GetCellData()->GetArray("Status")->SetComponent(path->GetId(i), 0, 5);
		}
	}
}

/*! Set the tissue form the cells only in one material.
 \param data Pointer to the orginal mesh.
 \param path Cell which material would be set.
 \param setMaterial Tissue class as that the cell would be marked.
 \param inMaterial Tissue class in that the cell would be marked.
 */
void Methods::setMaterial(DataFormat &data, vtkSmartPointer<vtkIdList> path, Material::Mat setMaterial, Material::Mat inMaterial) {
	Methods::setMaterial(data, path, setMaterial, inMaterial, inMaterial, inMaterial);
}

/*! Set the tissue form the cells only in two material.
 \param data Pointer to the orginal mesh.
 \param path Cell which material would be set.
 \param setMaterial Tissue class as that the cell would be marked.
 \param inMaterial1 Frirst tissue class in that the cell would be marked.
 \param inMaterial2 Second tissue class in that the cell would be marked.

 */
void Methods::setMaterial(DataFormat &data, vtkSmartPointer<vtkIdList> path, Material::Mat setMaterial, Material::Mat inMaterial1, Material::Mat inMaterial2) {
	Methods::setMaterial(data, path, setMaterial, inMaterial1, inMaterial2, inMaterial1);
}

/*! Set the tissue form the cells only in three material.
 \param data Pointer to the orginal mesh.
 \param path Cell which material would be set.
 \param setMaterial Tissue class as that the cell would be marked.
 \param inMaterial1 First tissue class in that the cell would be marked.
 \param inMaterial2 Second tissue class in that the cell would be marked.
 \param inMaterial3 Third tissue class in that the cell would be marked.

 */
void Methods::setMaterial(DataFormat &data, vtkSmartPointer<vtkIdList> path, Material::Mat setMaterial, Material::Mat inMaterial1, Material::Mat inMaterial2, Material::Mat inMaterial3) {
	for (vtkIdType i = 0; i < path->GetNumberOfIds(); i++) {
		if (!data.getDoubleLayer()) {
			int material = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(path->GetId(i), 0);
			if ((material == inMaterial1) || (material == inMaterial2) || (material == inMaterial3)) {
				data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->SetComponent(path->GetId(i), 0, setMaterial);
				data.getVtkData()->GetCellData()->GetArray("StatusEndo")->SetComponent(path->GetId(i), 0, 5);
			}
		}

		int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(path->GetId(i), 0);
		if ((material == inMaterial1) || (material == inMaterial2) || (material == inMaterial3)) {
			data.getVtkData()->GetCellData()->GetArray("Material")->SetComponent(path->GetId(i), 0, setMaterial);
			data.getVtkData()->GetCellData()->GetArray("Status")->SetComponent(path->GetId(i), 0, 5);
		}
	}
}

/*! Set the tissue form the cells only in one material.
 \param data Pointer to the orginal mesh.
 \param path Cell which material would be set.
 \param setMaterial Tissue class as that the cell would be marked.
 \param inMaterial Tissue class in that the cell would be marked.
 */
void Methods::setMaterial(DataFormat &data, vtkSmartPointer<vtkIdList> path, int setMaterial, Material::Mat inMaterial) {
	for (vtkIdType i = 0; i < path->GetNumberOfIds(); i++) {
		if (!data.getDoubleLayer()) {
			int material = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(path->GetId(i), 0);
			if (material == inMaterial) {
				data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->SetComponent(path->GetId(i), 0, setMaterial);
				data.getVtkData()->GetCellData()->GetArray("StatusEndo")->SetComponent(path->GetId(i), 0, 5);
			}
		}

		int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(path->GetId(i), 0);
		if (material == inMaterial) {
			data.getVtkData()->GetCellData()->GetArray("Material")->SetComponent(path->GetId(i), 0, setMaterial);
			data.getVtkData()->GetCellData()->GetArray("Status")->SetComponent(path->GetId(i), 0, 5);
		}
	}
}

/*! Set the tissue form the cells as the nerreast different material.
 \param data Pointer to the orginal mesh.
 \param cells Cell which material would be set.
 */

 /*! Set the tissue form the cells only in one material and one direction of a plane.
  \param data Pointer to the orginal mesh.
  \param path Cell which material would be set.
  \param setMaterial Tissue class as that the cell would be marked.
  \param inMaterial Tissue class in that the cell would be marked.
  \param point1Plane1 Point 1 of the first plane as Vector (x,y,z).
  \param point2Plane1 Point 2 of the first plane as Vector (x,y,z).
  \param point3Plane1 Point 3 of the first plane as Vector (x,y,z).
  \param directionPoint1 Direction point of the first plane as Vector (x,y,z) in that the cells do not be marked.
  \param point1Plane2 Point 1 of the second plane as Vector (x,y,z).
  \param point2Plane2 Point 2 of the second plane as Vector (x,y,z).
  \param point3Plane2 Point 3 of the second plane as Vector (x,y,z).
  \param directionPoint2 Direction point of the second plane as Vector (x,y,z) in that the cells do not be marked.
  \param seedPoint Seed point of the region that would e marked.
  */
void Methods::setMaterialinInMaterialNotinDirectionandDirection(DataFormat &data, Material::Mat setMaterial, Material::Mat inMaterial, std::vector<double> point1Plane1, std::vector<double> point2Plane1, std::vector<double> point3Plane1, std::vector<double> directionPointPlane1, std::vector<double> point1Plane2, std::vector<double> point2Plane2, std::vector<double> point3Plane2, std::vector<double> directionPointPlane2, vector<double> seedPoint) {
	double point1Plane1Array[3] = { point1Plane1.at(0), point1Plane1.at(1), point1Plane1.at(2) };
	double point2Plane1Array[3] = { point2Plane1.at(0), point2Plane1.at(1), point2Plane1.at(2) };
	double point3Plane1Array[3] = { point3Plane1.at(0), point3Plane1.at(1), point3Plane1.at(2) };
	double directionPointPlane1Array[3] = { directionPointPlane1.at(0), directionPointPlane1.at(1), directionPointPlane1.at(2) };

	double normalPlane1[3];
	double vector1Point1Point2[3];
	double vector1point1Point3[3];

	vtkMath::Subtract(point2Plane1Array, point1Plane1Array, vector1Point1Point2);
	vtkMath::Subtract(point3Plane1Array, point1Plane1Array, vector1point1Point3);
	vtkMath::Cross(vector1Point1Point2, vector1point1Point3, normalPlane1);
	vtkMath::Normalize(normalPlane1);

	double evalPlane1EqualPoint = vtkPlane::Evaluate(normalPlane1, point1Plane1Array, directionPointPlane1Array);


	double point1Plane2Array[3] = { point1Plane2.at(0), point1Plane2.at(1), point1Plane2.at(2) };
	double point2Plane2Array[3] = { point2Plane2.at(0), point2Plane2.at(1), point2Plane2.at(2) };
	double point3Plane2Array[3] = { point3Plane2.at(0), point3Plane2.at(1), point3Plane2.at(2) };
	double directionPointPlane2Array[3] = { directionPointPlane2.at(0), directionPointPlane2.at(1), directionPointPlane2.at(2) };

	double normalPlane2[3];
	double vector2Point1Point2[3];
	double vector2point1Point3[3];

	vtkMath::Subtract(point2Plane2Array, point1Plane2Array, vector2Point1Point2);
	vtkMath::Subtract(point3Plane2Array, point1Plane2Array, vector2point1Point3);
	vtkMath::Cross(vector2Point1Point2, vector2point1Point3, normalPlane2);
	vtkMath::Normalize(normalPlane2);


	double evalPlane2EqualPoint = vtkPlane::Evaluate(normalPlane2, point1Plane2Array, directionPointPlane2Array);

	double seedPointArray[3] = { seedPoint.at(0), seedPoint.at(1), seedPoint.at(2) };

	vtkSmartPointer<vtkPoints> centrePoints = data.getCentrePoints()->GetPoints();
	vtkSmartPointer<vtkUnstructuredGrid> centrePointsUGrid = data.getCentrePoints();

	vtkIdType seedPointID = centrePointsUGrid->FindPoint(seedPointArray);

	set<vtkIdType> viewedCellList;

	bool end = false;

	set<vtkIdType> cellIdList;
	cellIdList.insert(seedPointID);

	while (!end) {
		set<vtkIdType> newcellIdList;
		for (set<vtkIdType>::iterator iter = cellIdList.begin(); iter != cellIdList.end(); iter++) {
			vtkIdType currentCellIds = *iter;
			vtkSmartPointer<vtkIdList> cellPointIds = vtkSmartPointer<vtkIdList>::New();

			data.getVtkData()->GetCellPoints(currentCellIds, cellPointIds);

			for (vtkIdType i = 0; i < cellPointIds->GetNumberOfIds(); i++) {
				vtkSmartPointer<vtkIdList> CellPointIdList = vtkSmartPointer<vtkIdList>::New();
				CellPointIdList->InsertNextId(cellPointIds->GetId(i));

				vtkSmartPointer<vtkIdList> tempNeighborCellIds = vtkSmartPointer<vtkIdList>::New();
				data.getVtkData()->GetCellNeighbors(currentCellIds, CellPointIdList, tempNeighborCellIds);

				for (vtkIdType j = 0; j < tempNeighborCellIds->GetNumberOfIds(); j++) {
					if (!data.getDoubleLayer()) {
						if ((viewedCellList.find(tempNeighborCellIds->GetId(j)) == viewedCellList.end()) &&
							(data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(tempNeighborCellIds->GetId(j), 0) == inMaterial)) {
							double evalPlane1 = vtkPlane::Evaluate(normalPlane1, point1Plane1Array, centrePoints->GetPoint(tempNeighborCellIds->GetId(j)));
							double evalPlane2 = vtkPlane::Evaluate(normalPlane2, point1Plane2Array, centrePoints->GetPoint(tempNeighborCellIds->GetId(j)));

							if ((((evalPlane1 <= 0) && (evalPlane1EqualPoint >= 0)) ||
								((evalPlane1 >= 0) && (evalPlane1EqualPoint <= 0))) &&
								(((evalPlane2 <= 0) && (evalPlane2EqualPoint >= 0)) ||
								((evalPlane2 >= 0) && (evalPlane2EqualPoint <= 0)))) {
								data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->SetComponent(tempNeighborCellIds->GetId(j), 0, setMaterial);
								data.getVtkData()->GetCellData()->GetArray("StatusEndo")->SetComponent(tempNeighborCellIds->GetId(j), 0, 5);

								newcellIdList.insert(tempNeighborCellIds->GetId(j));
							}
						}
					}

					if ((viewedCellList.find(tempNeighborCellIds->GetId(j)) == viewedCellList.end()) &&
						(data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(tempNeighborCellIds->GetId(j), 0) == inMaterial)) {
						double evalPlane1 = vtkPlane::Evaluate(normalPlane1, point1Plane1Array, centrePoints->GetPoint(tempNeighborCellIds->GetId(j)));
						double evalPlane2 = vtkPlane::Evaluate(normalPlane2, point1Plane2Array, centrePoints->GetPoint(tempNeighborCellIds->GetId(j)));


						if ((((evalPlane1 <= 0) && (evalPlane1EqualPoint >= 0)) ||
							((evalPlane1 >= 0) && (evalPlane1EqualPoint <= 0))) &&
							(((evalPlane2 <= 0) && (evalPlane2EqualPoint >= 0)) ||
							((evalPlane2 >= 0) && (evalPlane2EqualPoint <= 0)))) {
							data.getVtkData()->GetCellData()->GetArray("Material")->SetComponent(tempNeighborCellIds->GetId(j), 0, setMaterial);
							data.getVtkData()->GetCellData()->GetArray("Status")->SetComponent(tempNeighborCellIds->GetId(j), 0, 5);

							newcellIdList.insert(tempNeighborCellIds->GetId(j));
						}
					} // alle nachbarzellen neuen
					viewedCellList.insert(tempNeighborCellIds->GetId(j));
				}
			}
		}
		if (newcellIdList.size() == 0) {
			end = true;
		}
		cellIdList.clear();
		cellIdList = newcellIdList;
	}
} // Methods::setMaterialinInMaterialNotinDirectionandDirection

/*! Set the tissue form the cells only in one material and one direction of a plane.
 \param data Pointer to the orginal mesh.
 \param path Cell which material would be set.
 \param setMaterial Tissue class as that the cell would be marked.
 \param inMaterial Tissue class in that the cell would be marked.
 \param point1Plane Point 1 of the first plane as Vector (x,y,z).
 \param point2Plane Point 2 of the first plane as Vector (x,y,z).
 \param point3Plane Point 3 of the first plane as Vector (x,y,z).
 \param directionPoint Direction point of the first plane as Vector (x,y,z) in that the cells do not be marked.
 */
void Methods::setMaterialNotinDirection(DataFormat &data, vtkSmartPointer<vtkIdList> path, Material::Mat setMaterial, Material::Mat inMaterial, std::vector<double> point1Plane, std::vector<double> point2Plane, std::vector<double> point3Plane, std::vector<double> directionPoint) {
	vtkSmartPointer<vtkPoints> centrePoints = data.getCentrePoints()->GetPoints();


	double point1PlaneArray[3] = { point1Plane.at(0), point1Plane.at(1), point1Plane.at(2) };
	double point2PlaneArray[3] = { point2Plane.at(0), point2Plane.at(1), point2Plane.at(2) };
	double point3PlaneArray[3] = { point3Plane.at(0), point3Plane.at(1), point3Plane.at(2) };
	double directionPointArray[3] = { directionPoint.at(0), directionPoint.at(1), directionPoint.at(2) };

	double normalPlane[3];
	double vectorPoint1Point2[3];
	double vectorpoint1Point3[3];

	vtkMath::Subtract(point2PlaneArray, point1PlaneArray, vectorPoint1Point2);
	vtkMath::Subtract(point3PlaneArray, point1PlaneArray, vectorpoint1Point3);
	vtkMath::Cross(vectorPoint1Point2, vectorpoint1Point3, normalPlane);
	vtkMath::Normalize(normalPlane);

	double evalPlaneEqualPoint = vtkPlane::Evaluate(normalPlane, point1PlaneArray, directionPointArray);


	for (vtkIdType i = 0; i < path->GetNumberOfIds(); i++) {
		if (!data.getDoubleLayer()) {
			int material = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(path->GetId(i), 0);
			double evalPlane = vtkPlane::Evaluate(normalPlane, point1PlaneArray, centrePoints->GetPoint(path->GetId(i)));
			if ((inMaterial == material) &&
				(((evalPlane <= 0) && (evalPlaneEqualPoint > 0)) || ((evalPlane >= 0) && (evalPlaneEqualPoint < 0)))) {
				data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->SetComponent(path->GetId(i), 0, setMaterial);
			}
		}
		int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(path->GetId(i), 0);
		double evalPlane = vtkPlane::Evaluate(normalPlane, point1PlaneArray, centrePoints->GetPoint(path->GetId(i)));
		if ((inMaterial == material) &&
			(((evalPlane <= 0) && (evalPlaneEqualPoint > 0)) || ((evalPlane >= 0) && (evalPlaneEqualPoint < 0)))) {
			data.getVtkData()->GetCellData()->GetArray("Material")->SetComponent(path->GetId(i), 0, setMaterial);
		}
	}
} // Methods::setMaterialNotinDirection

/*! Calculate the length of a path
 \param data Pointer to the orginal mesh.
 \param path Path that length would be calculated.
 \return The length of a path.
 */
double Methods::pathLength(DataFormat &data, vtkSmartPointer<vtkIdList> path) {
	double distance = 0;
	vtkIdType lastId = path->GetId(0);

	for (vtkIdType i = 1; i < path->GetNumberOfIds(); i++) {
		vtkIdType currentId = path->GetId(i);
		double x[3] = { 0, 0, 0 };
		x[0] = data.getCentrePoints()->GetPoint(lastId)[0];
		x[1] = data.getCentrePoints()->GetPoint(lastId)[1];
		x[2] = data.getCentrePoints()->GetPoint(lastId)[2];

		double y[3] = { 0, 0, 0 };
		y[0] = data.getCentrePoints()->GetPoint(currentId)[0];
		y[1] = data.getCentrePoints()->GetPoint(currentId)[1];
		y[2] = data.getCentrePoints()->GetPoint(currentId)[2];

		distance += sqrt(abs(vtkMath::Distance2BetweenPoints(x, y)));
		lastId = currentId;
	}
	return distance;
}

/*! Return the cell-id at percent position on the path.
 \param data Pointer to the orginal mesh.
 \param path The hole path.
 \param percent Percent at that the cell-id would be return.
 \return The cell-id at the position on the path.
 */
vtkIdType Methods::pointIdAtPercentOfPath(DataFormat &data, vtkSmartPointer<vtkIdList> path, double percent) {
	double maxdistance = pathLength(data, path);

	double distance = maxdistance * percent / 100;


	double tempdistance = 0;

	vtkIdType lastId = path->GetId(0);
	vtkIdType currentId = -1;
	vtkIdType i = 1;

	do {
		if (i == path->GetNumberOfIds()) {
			return path->GetId(path->GetNumberOfIds() - 1);
		}

		currentId = path->GetId(i);
		double x[3] = { 0, 0, 0 };
		x[0] = data.getCentrePoints()->GetPoint(lastId)[0];
		x[1] = data.getCentrePoints()->GetPoint(lastId)[1];
		x[2] = data.getCentrePoints()->GetPoint(lastId)[2];

		double y[3] = { 0, 0, 0 };
		y[0] = data.getCentrePoints()->GetPoint(currentId)[0];
		y[1] = data.getCentrePoints()->GetPoint(currentId)[1];
		y[2] = data.getCentrePoints()->GetPoint(currentId)[2];

		tempdistance += sqrt(abs(vtkMath::Distance2BetweenPoints(x, y)));
		i++;


		lastId = currentId;
	} while (tempdistance < distance);

	return currentId;
} // Methods::pointIdAtPercentOfPath

/*! Return the center point of the cell at percent position on the path.
 \param data Pointer to the orginal mesh.
 \param path The hole path.
 \param percent Percent at that the center point of the cell would be return.
 \return The center point of the cell at the position on the path.
 */
vector<double> Methods::pointAtPercentOfPath(DataFormat &data, vtkSmartPointer<vtkIdList> path, double percent) {
	vtkIdType pointId = Methods::pointIdAtPercentOfPath(data, path, percent);

	// cout << "PointID: " << pointId << endl;
	vector<double> point(3);

	point.at(0) = data.getCentrePoints()->GetPoint(pointId)[0];
	point.at(1) = data.getCentrePoints()->GetPoint(pointId)[1];
	point.at(2) = data.getCentrePoints()->GetPoint(pointId)[2];
	return point;
}

/*! Return the cell-id at mm position on the path.
 \param data Pointer to the orginal mesh.
 \param path The hole path.
 \param distance Distance in mm at that the cell-id would be return.
 \return The cell-id at the position on the path.
 */
vtkIdType Methods::pointIdAtMMOfPath(DataFormat &data, vtkSmartPointer<vtkIdList> path, double distance) {
	vtkIdType lastId = path->GetId(0);
	vtkIdType currentId = -1;
	double tempdistance = 0;
	vtkIdType i = 1;

	do {
		currentId = path->GetId(i);

		double x[3] = { 0, 0, 0 };
		x[0] = data.getCentrePoints()->GetPoint(lastId)[0];
		x[1] = data.getCentrePoints()->GetPoint(lastId)[1];
		x[2] = data.getCentrePoints()->GetPoint(lastId)[2];

		double y[3] = { 0, 0, 0 };
		y[0] = data.getCentrePoints()->GetPoint(currentId)[0];
		y[1] = data.getCentrePoints()->GetPoint(currentId)[1];
		y[2] = data.getCentrePoints()->GetPoint(currentId)[2];

		tempdistance += sqrt(abs(vtkMath::Distance2BetweenPoints(x, y)));

		i++;
		lastId = currentId;
		if (i >= path->GetNumberOfIds()) {
			cout << "Distance > Pathlength; it will be return the last ID of the Path" << endl;
			return path->GetId(path->GetNumberOfIds() - 1);
		}
	} while (tempdistance < distance);

	return currentId;
} // Methods::pointIdAtMMOfPath

/*! Return the center point of the cell at mm position on the path.
 \param data Pointer to the orginal mesh.
 \param path The hole path.
 \param distance Distance in mm at that the center point of the cell would be return.
 \return The center point of the cell at the position on the path.
 */
vector<double> Methods::pointAtMMOfPath(DataFormat &data, vtkSmartPointer<vtkIdList> path, double distance) {
	vtkIdType pointId = Methods::pointIdAtMMOfPath(data, path, distance);

	vector<double> point(3);

	point.at(0) = data.getCentrePoints()->GetPoint(pointId)[0];
	point.at(1) = data.getCentrePoints()->GetPoint(pointId)[1];
	point.at(2) = data.getCentrePoints()->GetPoint(pointId)[2];
	return point;
}

/*! Calc x- percent of th path.
 \param data Pointer to the orginal mesh.
 \param path The hole path.
 \param Percent Percent of the path that would be return.
 \return Return x- percent of th path form the beginning.
 */
vtkSmartPointer<vtkIdList> Methods::getPercentOfPath(DataFormat &data, vtkSmartPointer<vtkIdList> path, double percent) {
	vtkSmartPointer<vtkIdList> newPath = vtkSmartPointer<vtkIdList>::New();

	vtkIdType targetPointId = Methods::pointIdAtPercentOfPath(data, path, percent);

	for (vtkIdType i = 0; i < path->GetNumberOfIds(); i++) {
		vtkIdType currentId = path->GetId(i);
		newPath->InsertNextId(currentId);
		if (currentId == targetPointId)
			return newPath;
	}
	newPath->Initialize();
	newPath->InsertNextId(-1);
	return newPath;
}

/*! Calc x- mm of th path.
 \param data Pointer to the orginal mesh.
 \param path The hole path.
 \param distance Distance in mm of the path that would be return.
 \return Return x- mm of th path form the beginning.
 */
vtkSmartPointer<vtkIdList> Methods::getMMOfPath(DataFormat &data, vtkSmartPointer<vtkIdList> path, double distance) {
	vtkSmartPointer<vtkIdList> newPath = vtkSmartPointer<vtkIdList>::New();

	vtkIdType targetPointId = Methods::pointIdAtMMOfPath(data, path, distance);


	for (vtkIdType i = 0; i < path->GetNumberOfIds(); i++) {
		vtkIdType currentId = path->GetId(i);
		newPath->InsertNextId(currentId);
		if (currentId == targetPointId)
			return newPath;
	}
	newPath->Initialize();
	newPath->InsertNextId(-1);
	return newPath;
}

/*! Invert the path and return x- mm of th path.
 \param data Pointer to the orginal mesh.
 \param path The hole path.
 \param distance Distance in mm of the path that would be return.
 \return Return x- mm of th path form the invert path.
 */
vtkSmartPointer<vtkIdList> Methods::getMMInverseOfPath(DataFormat &data, vtkSmartPointer<vtkIdList> path, double distance) {
	vtkSmartPointer<vtkIdList> newPath = vtkSmartPointer<vtkIdList>::New();
	double lengthPath = pathLength(data, path);

	vtkIdType targetPointId = Methods::pointIdAtMMOfPath(data, path, (lengthPath - distance));

	for (vtkIdType i = path->GetNumberOfIds() - 1; i >= 0; i--) {
		vtkIdType currentId = path->GetId(i);
		newPath->InsertNextId(currentId);
		if (currentId == targetPointId)
			return inversePath(newPath);
	}
	newPath->Initialize();
	newPath->InsertNextId(-1);
	return newPath;
}

/*! Get a Point in a define material on a path. The point lie after a min distance in the material.
 \param data Pointer to the orginal mesh.
 \param path The hole path.
 \param material Tissue class in that the point must lie.
 \param minDistanceInMaterial Min distance that the point lies in the material.
 \return Return a point (X,y,z) in a define material on a path or NAN if no point was found.
 */
vector<double> Methods::getPointAlongPathOfPointInMaterial(DataFormat &data, vtkSmartPointer<vtkIdList> path, int inMaterial1, double minDistanceInMaterial) {
	return getPointAlongPathOfPointInMaterial(data, path, inMaterial1, inMaterial1, minDistanceInMaterial);
}

/*! Get a Point in a define material on a path. The point lie after a min distance in the material.
 \param data Pointer to the orginal mesh.
 \param path The hole path.
 \param material1 First tissue class in that the point must lie.
 \param material2 Second tissue class in that the point must lie.
 \param minDistanceInMaterial Min distance that the point lies in the material.
 \return Return a point (X,y,z) in a define material on a path or NAN if no point was found.
 */
vector<double> Methods::getPointAlongPathOfPointInMaterial(DataFormat &data, vtkSmartPointer<vtkIdList> path, int inMaterial1, int inMaterial2, double minDistanceInMaterial) {
	bool end = false;

	vtkIdType i = 0;

	vector<double> point(3);

	point.at(0) = NAN;
	point.at(1) = NAN;
	point.at(2) = NAN;
	double lastPoint[3] = { data.getCentrePoints()->GetPoint(path->GetId(0))[0], data.getCentrePoints()->GetPoint(path->GetId(0))[1], data.getCentrePoints()->GetPoint(path->GetId(0))[2] };

	double distance = 0;
	double distanceEndo = 0;

	while (!end) {
		double currentPoint[3] = { data.getCentrePoints()->GetPoint(path->GetId(i))[0], data.getCentrePoints()->GetPoint(path->GetId(i))
		[1], data.getCentrePoints()->GetPoint(path->GetId(i))[2] };

		if (!data.getDoubleLayer()) {
			int currentMaterial1 = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(path->GetId(i), 0);
			int currentMaterial2 = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(path->GetId(i), 0);

			if ((currentMaterial1 == inMaterial1) || (currentMaterial1 == inMaterial2) ||
				(currentMaterial2 == inMaterial1) || (currentMaterial2 == inMaterial2)) {
				distanceEndo += sqrt(abs(vtkMath::Distance2BetweenPoints(lastPoint, currentPoint)));

				if (distanceEndo >= minDistanceInMaterial) {
					end = true;
					point.at(0) = data.getCentrePoints()->GetPoint(path->GetId(i))[0];
					point.at(1) = data.getCentrePoints()->GetPoint(path->GetId(i))[1];
					point.at(2) = data.getCentrePoints()->GetPoint(path->GetId(i))[2];
				}
			}

			lastPoint[0] = currentPoint[0];
			lastPoint[1] = currentPoint[1];
			lastPoint[2] = currentPoint[2];

			i++;

			if (i == path->GetNumberOfIds())
				end = true;
		}
		else {
			int currentMaterial = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(path->GetId(i), 0);
			if ((currentMaterial == inMaterial1) || (currentMaterial == inMaterial2)) {
				distance += sqrt(abs(vtkMath::Distance2BetweenPoints(lastPoint, currentPoint)));
				if (distance >= minDistanceInMaterial) {
					end = true;
					point.at(0) = data.getCentrePoints()->GetPoint(path->GetId(i))[0];
					point.at(1) = data.getCentrePoints()->GetPoint(path->GetId(i))[1];
					point.at(2) = data.getCentrePoints()->GetPoint(path->GetId(i))[2];
				}
			}

			lastPoint[0] = currentPoint[0];
			lastPoint[1] = currentPoint[1];
			lastPoint[2] = currentPoint[2];

			i++;
			if (i == path->GetNumberOfIds())
				end = true;
		}
	}

	if (::isnan(point.at(0)) && (distance < minDistanceInMaterial) && (distanceEndo < minDistanceInMaterial)) {
		if (!data.getDoubleLayer()) {
			int currentMaterial1 = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(path->GetId(path->GetNumberOfIds() - 1), 0);
			int currentMaterial2 = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(path->GetId(path->GetNumberOfIds() - 1), 0);

			if ((currentMaterial1 == inMaterial1) || (currentMaterial1 == inMaterial2) ||
				(currentMaterial2 == inMaterial1) || (currentMaterial2 == inMaterial2)) {
				point.at(0) = data.getCentrePoints()->GetPoint(path->GetId(path->GetNumberOfIds() - 1))[0];
				point.at(1) = data.getCentrePoints()->GetPoint(path->GetId(path->GetNumberOfIds() - 1))[1];
				point.at(2) = data.getCentrePoints()->GetPoint(path->GetId(path->GetNumberOfIds() - 1))[2];
			}
		}
		else {
			int currentMaterial = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(path->GetId(path->GetNumberOfIds() - 1), 0);
			if ((currentMaterial == inMaterial1) || (currentMaterial == inMaterial2)) {
				point.at(0) = data.getCentrePoints()->GetPoint(path->GetId(path->GetNumberOfIds() - 1))[0];
				point.at(1) = data.getCentrePoints()->GetPoint(path->GetId(path->GetNumberOfIds() - 1))[1];
				point.at(2) = data.getCentrePoints()->GetPoint(path->GetId(path->GetNumberOfIds() - 1))[2];
			}
		}
	}

	return point;
} // Methods::getPointAlongPathOfPointInMaterial


vector<double> Methods::getPointAlongPathOfPointInMaterial(DataFormat &data, vtkSmartPointer<vtkIdList> path, int inMaterial1, int inMaterial2, double offsetMin, double offsetMax) {
	bool end = false;

	vtkIdType i = 0;

	vector<double> point(3);

	point.at(0) = NAN;
	point.at(1) = NAN;
	point.at(2) = NAN;

	double lastPoint[3] = { data.getCentrePoints()->GetPoint(path->GetId(0))[0], data.getCentrePoints()->GetPoint(path->GetId(0))[1], data.getCentrePoints()->GetPoint(path->GetId(0))[2] };

	double distance = 0;
	while (!end) {
		double currentPoint[3] = { data.getCentrePoints()->GetPoint(path->GetId(i))[0], data.getCentrePoints()->GetPoint(path->GetId(i))
		[1], data.getCentrePoints()->GetPoint(path->GetId(i))[2] };
		distance += sqrt(abs(vtkMath::Distance2BetweenPoints(lastPoint, currentPoint)));
		if ((distance >= offsetMin) && (distance <= offsetMax)) {
			if (!data.getDoubleLayer()) {
				int currentMaterial1 = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(path->GetId(i), 0);
				int currentMaterial2 = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(path->GetId(i), 0);

				if ((currentMaterial1 == inMaterial1) || (currentMaterial1 == inMaterial2) ||
					(currentMaterial2 == inMaterial1) || (currentMaterial2 == inMaterial2)) {
					end = true;
					point.at(0) = data.getCentrePoints()->GetPoint(path->GetId(i))[0];
					point.at(1) = data.getCentrePoints()->GetPoint(path->GetId(i))[1];
					point.at(2) = data.getCentrePoints()->GetPoint(path->GetId(i))[2];
				}
			}
			else {
				int currentMaterial = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(path->GetId(i), 0);

				if ((currentMaterial == inMaterial1) || (currentMaterial == inMaterial2)) {
					end = true;
					point.at(0) = data.getCentrePoints()->GetPoint(path->GetId(i))[0];
					point.at(1) = data.getCentrePoints()->GetPoint(path->GetId(i))[1];
					point.at(2) = data.getCentrePoints()->GetPoint(path->GetId(i))[2];
				}
			}
		}
		else if (distance > offsetMax) {
			cout << "Punkt konnte nicht im Offest aus dem Matierial verschoben werden" << endl;
			end = true;
		}
		lastPoint[0] = currentPoint[0];
		lastPoint[1] = currentPoint[1];
		lastPoint[2] = currentPoint[2];

		i++;
	}
	return point;
} // Methods::getPointAlongPathOfPointInMaterial

/*! Sub a path form another path.
 \param path1 The first path.
 \param path2 The second path.
 \return path1 - path2.
 */
vtkSmartPointer<vtkIdList> Methods::subPathFromPath(vtkSmartPointer<vtkIdList> path1, vtkSmartPointer<vtkIdList> path2) {
	set<vtkIdType> path2Set;

	for (vtkIdType i = 0; i < path2->GetNumberOfIds(); i++) {
		path2Set.insert(path2->GetId(i));
	}


	vtkSmartPointer<vtkIdList> tempPath = vtkSmartPointer<vtkIdList>::New();
	for (vtkIdType i = 0; i < path1->GetNumberOfIds(); i++) {
		if (path2Set.find(path1->GetId(i)) == path2Set.end())
			tempPath->InsertNextId(path1->GetId(i));
	}

	return tempPath;
}

/*! Add a path to another path.
 \param path1 The first path.
 \param path2 The second path.
 \return path1 + path2.
 */
vtkSmartPointer<vtkIdList> Methods::add2Paths(vtkSmartPointer<vtkIdList> path1, vtkSmartPointer<vtkIdList> path2) {
	vtkSmartPointer<vtkIdList> tempPath = vtkSmartPointer<vtkIdList>::New();

	tempPath->DeepCopy(path1);

	for (vtkIdType i = 0; i < path2->GetNumberOfIds(); i++) {
		tempPath->InsertUniqueId(path2->GetId(i));
	}
	return tempPath;
}

/*! Add three paths.
 \param path1 The first path.
 \param path2 The second path.
 \param path3 The third path.
 \return path1 + path2 + path3.
 */
vtkSmartPointer<vtkIdList> Methods::add3Paths(vtkSmartPointer<vtkIdList> path1, vtkSmartPointer<vtkIdList> path2, vtkSmartPointer<vtkIdList> path3) {
	vtkSmartPointer<vtkIdList> tempPath = vtkSmartPointer<vtkIdList>::New();

	tempPath->DeepCopy(path1);

	for (vtkIdType i = 0; i < path2->GetNumberOfIds(); i++) {
		tempPath->InsertUniqueId(path2->GetId(i));
	}
	for (vtkIdType i = 0; i < path3->GetNumberOfIds(); i++) {
		tempPath->InsertUniqueId(path3->GetId(i));
	}
	return tempPath;
}

/*! Add four paths.
 \param path1 The first path.
 \param path2 The second path.
 \param path3 The third path.
 \param path4 The fourth path.
 \return path1 + path2 + path3 + path4.
 */
vtkSmartPointer<vtkIdList> Methods::add4Paths(vtkSmartPointer<vtkIdList> path1, vtkSmartPointer<vtkIdList> path2, vtkSmartPointer<vtkIdList> path3, vtkSmartPointer<vtkIdList> path4) {
	vtkSmartPointer<vtkIdList> tempPath = vtkSmartPointer<vtkIdList>::New();

	tempPath->DeepCopy(path1);

	for (vtkIdType i = 0; i < path2->GetNumberOfIds(); i++) {
		tempPath->InsertUniqueId(path2->GetId(i));
	}
	for (vtkIdType i = 0; i < path3->GetNumberOfIds(); i++) {
		tempPath->InsertUniqueId(path3->GetId(i));
	}
	for (vtkIdType i = 0; i < path4->GetNumberOfIds(); i++) {
		tempPath->InsertUniqueId(path4->GetId(i));
	}
	return tempPath;
}

/*! Add five paths.
 \param path1 The first path.
 \param path2 The second path.
 \param path3 The third path.
 \param path4 The fourth path.
 \param path5 The fifth path.
 \return path1 + path2 + path3 + path4 + path5.
 */
vtkSmartPointer<vtkIdList> Methods::add5Paths(vtkSmartPointer<vtkIdList> path1, vtkSmartPointer<vtkIdList> path2, vtkSmartPointer<vtkIdList> path3, vtkSmartPointer<vtkIdList> path4, vtkSmartPointer<vtkIdList> path5) {
	vtkSmartPointer<vtkIdList> tempPath = vtkSmartPointer<vtkIdList>::New();

	tempPath->DeepCopy(path1);

	for (vtkIdType i = 0; i < path2->GetNumberOfIds(); i++) {
		tempPath->InsertUniqueId(path2->GetId(i));
	}
	for (vtkIdType i = 0; i < path3->GetNumberOfIds(); i++) {
		tempPath->InsertUniqueId(path3->GetId(i));
	}
	for (vtkIdType i = 0; i < path4->GetNumberOfIds(); i++) {
		tempPath->InsertUniqueId(path4->GetId(i));
	}
	for (vtkIdType j = 0; j < path5->GetNumberOfIds(); j++) {
		tempPath->InsertUniqueId(path5->GetId(j));
	}
	return tempPath;
}

/*! Invert a path.
 \param path The path that would be invert.
 \return The invert path.
 */
vtkSmartPointer<vtkIdList> Methods::inversePath(vtkSmartPointer<vtkIdList> path) {
	vtkSmartPointer<vtkIdList> tempPath = vtkSmartPointer<vtkIdList>::New();

	for (vtkIdType i = path->GetNumberOfIds() - 1; i >= 0; i--) {
		tempPath->InsertNextId(path->GetId(i));
	}
	return tempPath;
}

/*! Calc the angle (theta, phi) of a vector.
 \param x X componet of the vector.
 \param y Y componet of the vector.
 \param z Z componet of the vector.
 \return The angle of the vector as a vector (theta, phi).
 */
vector<double> Methods::angle(double x, double y, double z) {
	vector<double> angles(2);

	double theta = acos(z);

	angles.at(0) = theta;

	double phi = (double)atan2(y, x);
	angles.at(1) = phi;

	return angles;
}

/*! Round a double number.
 \param figure The number the would be round.
 \param digits Digtis of the round number.
 \return The round number.
 */
double Methods::round(double figure, int digits) {
	double order = pow(10.0, digits);

	return (int)(figure*order + (figure > 0 ? 0.5 : -0.5)) / order;
}

/*! Test linear if an id is in the array.
 \param array The aray in which the id would be search.
 \param search ID which be searched.
 \return The position of the id in the array or -1.
 */
vtkIdType Methods::isIDInsideArrayLinear(vtkSmartPointer<vtkDataArray> array, vtkIdType search) {
	for (int i = 0; i < array->GetNumberOfTuples(); i++) {
		if (search == array->GetComponent(i, 0))
			return i;
	}
	return -1;
}

/*! Calc the angle between two vectors.
 \param vector1 First vector
 \param vector2 Second vecotr.
 \return The angle between the two vectors.
 */
double Methods::angleBetween2Vectors(vector<double> vector1, vector<double> vector2) { // it is possible to replace
 // this function in vtk6.2 with
 // the vtkMath function
 // AngleBetweenVectors
	double v1[3] = { vector1.at(0), vector1.at(1), vector1.at(2) };
	double v2[3] = { vector2.at(0), vector2.at(1), vector2.at(2) };

	double dotvec = vtkMath::Dot(v1, v2);

	dotvec = round(dotvec, 2);

	double angle = NAN;
	double crossvec[3] = { 0, 0, 0 };
	vtkMath::Cross(v1, v2, crossvec);

	angle = acos(dotvec) * 180 / M_PI;

	if (round(crossvec[2], 2) < 0)
		angle += 180;

	return angle;
}

/*! Find a seedpoint in the region between to paths at 50% in a define material.
 \param data Pointer to the orginal mesh.
 \param path1 The first path.
 \param path2 The scond path.
 \param inMaterial Tissue class in that the point must lie.
 \param atrium Atrium in which the seed point would be searched.
 \return Return a seed point as Vector (x,y,z)
 */
vector<double> Methods::getSeedPoint(DataFormat &data, vtkSmartPointer<vtkIdList> path1, vtkSmartPointer<vtkIdList> path2, Material::Mat inMaterial, string atrium) {
	return getSeedPoint(data, path1, path2, inMaterial, inMaterial, inMaterial, inMaterial, 50, atrium);
}

/*! Find a seedpoint in the region between to paths at 50% in a define material.
 \param data Pointer to the orginal mesh.
 \param path1 The first path.
 \param path2 The scond path.
 \param inMaterial1 First tissue class in that the point could be lie.
 \param inMaterial2 Second tissue class in that the point could be lie.
 \param atrium Atrium in which the seed point would be searched.
 \return Return a seed point as Vector (x,y,z)
 */
vector<double> Methods::getSeedPoint(DataFormat &data, vtkSmartPointer<vtkIdList> path1, vtkSmartPointer<vtkIdList> path2, Material::Mat inMaterial1, Material::Mat inMaterial2, string atrium) {
	return getSeedPoint(data, path1, path2, inMaterial1, inMaterial2, inMaterial1, inMaterial1, 50, atrium);
}

/*! Find a seedpoint in the region between to paths at 50% in a define material.
 \param data Pointer to the orginal mesh.
 \param path1 The first path.
 \param path2 The scond path.
 \param inMaterial1 First tissue class in that the point could be lie.
 \param inMaterial2 Second tissue class in that the point could be lie.
 \param inMaterial3 Third tissue class in that the point could be lie.
 \param atrium Atrium in which the seed point would be searched.
 \return Return a seed point as Vector (x,y,z)
 */
vector<double> Methods::getSeedPoint(DataFormat &data, vtkSmartPointer<vtkIdList> path1, vtkSmartPointer<vtkIdList> path2, Material::Mat inMaterial1, Material::Mat inMaterial2, Material::Mat inMaterial3, string atrium) {
	return getSeedPoint(data, path1, path2, inMaterial1, inMaterial2, inMaterial3, inMaterial1, 50, atrium);
}

/*! Find a seedpoint in the region between to paths at 50% in a define material.
 \param data Pointer to the orginal mesh.
 \param path1 The first path.
 \param path2 The scond path.
 \param inMaterial1 First tissue class in that the point could be lie.
 \param inMaterial2 Second tissue class in that the point could be lie.
 \param inMaterial3 Third tissue class in that the point could be lie.
 \param inMaterial4 Fourth tissue class in that the point could be lie.
 \param atrium Atrium in which the seed point would be searched.
 \return Return a seed point as Vector (x,y,z)
 */
vector<double> Methods::getSeedPoint(DataFormat &data, vtkSmartPointer<vtkIdList> path1, vtkSmartPointer<vtkIdList> path2, Material::Mat inMaterial1, Material::Mat inMaterial2, Material::Mat inMaterial3, Material::Mat inMaterial4, string atrium) {
	return getSeedPoint(data, path1, path2, inMaterial1, inMaterial2, inMaterial3, inMaterial4, 50, atrium);
}

/*! Find a seedpoint in the region between to paths at X-percent in a define material.
 \param data Pointer to the orginal mesh.
 \param path1 The first path.
 \param path2 The scond path.
 \param inMaterial Tissue class in that the point must lie.
 \param percent Between this position on the path the seed point would be searched.
 \param atrium Atrium in which the seed point would be searched.
 \return Return a seed point as Vector (x,y,z)
 */
vector<double> Methods::getSeedPoint(DataFormat &data, vtkSmartPointer<vtkIdList> path1, vtkSmartPointer<vtkIdList> path2, Material::Mat inMaterial, double percent, string atrium) {
	return getSeedPoint(data, path1, path2, inMaterial, inMaterial, inMaterial, inMaterial, percent, atrium);
}

/*! Find a seedpoint in the region between to paths at X-percent in a define material.
 \param data Pointer to the orginal mesh.
 \param path1 The first path.
 \param path2 The scond path.
 \param inMaterial1 First tissue class in that the point could be lie.
 \param inMaterial2 Second tissue class in that the point could be lie.
 \param percent Between this position on the path the seed point would be searched.
 \param atrium Atrium in which the seed point would be searched.
 \return Return a seed point as Vector (x,y,z)
 */
vector<double> Methods::getSeedPoint(DataFormat &data, vtkSmartPointer<vtkIdList> path1, vtkSmartPointer<vtkIdList> path2, Material::Mat inMaterial1, Material::Mat inMaterial2, double percent, string atrium) {
	return getSeedPoint(data, path1, path2, inMaterial1, inMaterial2, inMaterial1, inMaterial1, percent, atrium);
}

/*! Find a seedpoint in the region between to paths at X-percent in a define material.
 \param data Pointer to the orginal mesh.
 \param path1 The first path.
 \param path2 The scond path.
 \param inMaterial1 First tissue class in that the point could be lie.
 \param inMaterial2 Second tissue class in that the point could be lie.
 \param inMaterial3 Third tissue class in that the point could be lie.
 \param percent Between this position on the path the seed point would be searched.
 \param atrium Atrium in which the seed point would be searched.
 \return Return a seed point as Vector (x,y,z)
 */
vector<double> Methods::getSeedPoint(DataFormat &data, vtkSmartPointer<vtkIdList> path1, vtkSmartPointer<vtkIdList> path2, Material::Mat inMaterial1, Material::Mat inMaterial2, Material::Mat inMaterial3, double percent, string atrium) {
	return getSeedPoint(data, path1, path2, inMaterial1, inMaterial2, inMaterial3, inMaterial1, percent, atrium);
}

/*! Find a seedpoint in the region between to paths at X-percent in a define material.
 \param data Pointer to the orginal mesh.
 \param path1 The first path.
 \param path2 The scond path.
 \param inMaterial1 First tissue class in that the point could be lie.
 \param inMaterial2 Second tissue class in that the point could be lie.
 \param inMaterial3 Third tissue class in that the point could be lie.
 \param inMaterial4 Fourth tissue class in that the point could be lie.
 \param percent Between this position on the path the seed point would be searched.
 \param atrium Atrium in which the seed point would be searched.
 \return Return a seed point as Vector (x,y,z)
 */
vector<double> Methods::getSeedPoint(DataFormat &data, vtkSmartPointer<vtkIdList> path1, vtkSmartPointer<vtkIdList> path2, Material::Mat inMaterial1, Material::Mat inMaterial2, Material::Mat inMaterial3, Material::Mat inMaterial4, double percent, string atrium) {
	vector<double> vectorPointPath1 = pointAtPercentOfPath(data, path1, percent);
	vector<double> vectorPointPath2 = pointAtPercentOfPath(data, path2, percent);

	vector<double> seedPoint(3);
	vtkIdType seedPointId = -1;
	vtkSmartPointer<vtkIdList> temppath = Methods::pathSearch(data, vectorPointPath1, vectorPointPath2, inMaterial1, inMaterial2, inMaterial3, inMaterial4, atrium);

	if (temppath->GetNumberOfIds() > 0) {
		seedPointId = Methods::pointIdAtPercentOfPath(data, temppath, 50);
	}
	else {
		cout << "SeedPoint not found on Path" << endl;
		seedPointId = findClosedPointIdinMaterialwithoutOrientation(data, vectorPointPath1, inMaterial1);
	}

	seedPoint.at(0) = data.getCentrePoints()->GetPoint(seedPointId)[0];
	seedPoint.at(1) = data.getCentrePoints()->GetPoint(seedPointId)[1];
	seedPoint.at(2) = data.getCentrePoints()->GetPoint(seedPointId)[2];

	if (!data.getDoubleLayer()) {
		int material2 = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(seedPointId, 0);
		if ((material2 != inMaterial1) && (material2 != inMaterial2) && (material2 != inMaterial3) &&
			(material2 != inMaterial4)) {
			int material1 = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(seedPointId, 0);
			if ((material1 != inMaterial1) && (material1 != inMaterial2) && (material1 != inMaterial3) &&
				(material1 != inMaterial4)) {
				cout << "getseedPoint: material not correct" << endl;
			}
		}
	}
	else {
		int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(seedPointId, 0);
		if ((material != inMaterial1) && (material != inMaterial2) && (material != inMaterial3) &&
			(material != inMaterial4)) {
			cout << "getseedPoint: material not correct" << endl;
		}
	}
	return seedPoint;
} // Methods::getSeedPoint

/*! Change the status of cells in a material.
 \param data Pointer to the orginal mesh.
 \param path The cells which status would be change.
 \param status The status that would be set.
 \param inMaterial Tissue class in that the status change.
 */
void Methods::setStatus(DataFormat &data, vtkSmartPointer<vtkIdList> path, int status, Material::Mat inMaterial) {
	setStatus(data, path, status, inMaterial, inMaterial);
}

/*! Change the status of cells in a material.
 \param data Pointer to the orginal mesh.
 \param path The cells which status would be change.
 \param status The status that would be set.
 \param inMaterial1 First tissue class in that the status change.
 \param inMaterial2 Second tissue class in that the status change.
 */
void Methods::setStatus(DataFormat &data, vtkSmartPointer<vtkIdList> path, int status, Material::Mat inMaterial1, Material::Mat inMaterial2) {
	for (vtkIdType i = 0; i < path->GetNumberOfIds(); i++) {
		if (!data.getDoubleLayer()) {
			int material = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(path->GetId(i), 0);
			if ((material == inMaterial1) || (material == inMaterial2)) {
				data.getVtkData()->GetCellData()->GetArray("StatusEndo")->SetComponent(path->GetId(i), 0, status);
			}
		}

		int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(path->GetId(i), 0);
		if ((material == inMaterial1) || (material == inMaterial2)) {
			data.getVtkData()->GetCellData()->GetArray("Status")->SetComponent(path->GetId(i), 0, status);
		}
	}
}

/*! Calc the cenroid of cells.
 \param data Pointer to the orginal mesh.
 \param path The cell ids of that the cenroid would be calculated .
 \return The cenroid of the cells as vector (x,y,z)
 */
vector<double> Methods::getCentroid(DataFormat &data, vtkSmartPointer<vtkIdList> path) {
	vector<double> value = getCentroidandBounds(data, path);
	vector<double> returnValue = { value.at(0), value.at(1), value.at(2) };

	return returnValue;
}

/*! Calc the cenroid and the bounding box of cells.
 \param data Pointer to the orginal mesh.
 \param path The cell ids of that the cenroid and the bounding box would be calculated .
 \return The cenroid of the cells as vector (x,y,z, x min, x max, y min, y max, z min, z max)
 */
vector<double> Methods::getCentroidandBounds(DataFormat &data, vtkSmartPointer<vtkIdList> path) {
	vtkSmartPointer<vtkPoints> centrePoints = data.getCentrePoints()->GetPoints();

	double xmax = -std::numeric_limits<double>::max();
	double xmin = std::numeric_limits<double>::max();
	double ymax = -std::numeric_limits<double>::max();
	double ymin = std::numeric_limits<double>::max();
	double zmax = -std::numeric_limits<double>::max();
	double zmin = std::numeric_limits<double>::max();


	for (vtkIdType i = 0; i < path->GetNumberOfIds(); i++) {
		vtkIdType currentId = path->GetId(i);
		double *tempPoint = centrePoints->GetPoint(currentId);
		if (tempPoint[0] < xmin) {
			xmin = tempPoint[0];
		}
		if (tempPoint[0] > xmax) {
			xmax = tempPoint[0];
		}
		if (tempPoint[1] < ymin) {
			ymin = tempPoint[1];
		}
		if (tempPoint[1] > ymax) {
			ymax = tempPoint[1];
		}
		if (tempPoint[2] < zmin) {
			zmin = tempPoint[2];
		}
		if (tempPoint[2] > zmax) {
			zmax = tempPoint[2];
		}
	}
	vector<double> returnValue(9);
	returnValue.at(0) = (xmin + xmax) / 2;
	returnValue.at(1) = (ymin + ymax) / 2;
	returnValue.at(2) = (zmin + zmax) / 2;
	returnValue.at(3) = xmin;
	returnValue.at(4) = xmax;
	returnValue.at(5) = ymin;
	returnValue.at(6) = ymax;
	returnValue.at(7) = zmin;
	returnValue.at(8) = zmax;

	return returnValue;
} // Methods::getCentroidandBounds


/*! Disconneted region with material 1 form the material 2 and 3.
 \param data Pointer to the orginal mesh.
 \param inMaterial1 Tissue class that would be disconnected.
 \param inMaterial2 Tissue class that would be disconnected.
 \param inMaterial3 Tissue class that would be disconnected.
 */
void Methods::uncoupleMaterials(DataFormat &data, Material::Mat inMaterial1, Material::Mat inMaterial2, Material::Mat inMaterial3) {
	return uncoupleMaterials(data, inMaterial1, inMaterial1, inMaterial2, inMaterial3);
}

/*! Disconneted region with material 1 and 2 form the material 3 and 4.
 \param data Pointer to the orginal mesh.
 \param inMaterial1 Tissue class that would be disconnected.
 \param inMaterial2 Tissue class that would be disconnected.
 \param inMaterial3 Tissue class that would be disconnected.
 \param inMaterial4 Tissue class that would be disconnected.
 */
void Methods::uncoupleMaterials(DataFormat &data, Material::Mat inMaterial1, Material::Mat inMaterial2, Material::Mat inMaterial3, Material::Mat inMaterial4) {
	vtkSmartPointer<vtkUnstructuredGrid> dataUGrid;
	vtkSmartPointer<vtkPolyData> dataPoly;
	vtkSmartPointer<vtkPoints> points = data.getVtkData()->GetPoints();
	vtkSmartPointer<vtkIdList> removedCells = vtkSmartPointer<vtkIdList>::New();


	if (data.getInputType() == DataFormat::PossibleInputType::vtu) {
		dataUGrid = vtkUnstructuredGrid::SafeDownCast(data.getVtkData());
	}

	if (data.getInputType() == DataFormat::PossibleInputType::vtp) {
		dataPoly = vtkPolyData::SafeDownCast(data.getVtkData());
	}
	unordered_map<vtkIdType, vtkIdType> newPointIDs;

	for (vtkIdType i = 0; i < data.getVtkData()->GetNumberOfCells(); i++) {
		int materialParent = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(i, 0);

		if ((materialParent == inMaterial3) || (materialParent == inMaterial4)) {
			vtkSmartPointer<vtkIdList> cellPointIds = vtkSmartPointer<vtkIdList>::New();
			vtkSmartPointer<vtkIdList> currentNeighborCellIds = vtkSmartPointer<vtkIdList>::New();

			data.getVtkData()->GetCellPoints(i, cellPointIds);

			for (vtkIdType i = 0; i < cellPointIds->GetNumberOfIds(); i++) {
				vtkSmartPointer<vtkIdList> PointIdList = vtkSmartPointer<vtkIdList>::New();
				PointIdList->InsertNextId(cellPointIds->GetId(i));

				vtkSmartPointer<vtkIdList> tempNeighborCellIds = vtkSmartPointer<vtkIdList>::New();

				data.getVtkData()->GetCellNeighbors(i, PointIdList, tempNeighborCellIds);

				for (vtkIdType j = 0; j < tempNeighborCellIds->GetNumberOfIds(); j++) {
					currentNeighborCellIds->InsertUniqueId(tempNeighborCellIds->GetId(j));
				}
			}
			for (vtkIdType j = 0; j < currentNeighborCellIds->GetNumberOfIds(); j++) {
				vtkIdType currentCellID = currentNeighborCellIds->GetId(j);
				int materialCurrentCell = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(currentCellID, 0);
				if ((materialCurrentCell == inMaterial1) || (materialCurrentCell == inMaterial2)) {
					removedCells->InsertUniqueId(currentCellID);
					vtkSmartPointer<vtkIdList> currentCellPointIds = vtkSmartPointer<vtkIdList>::New();
					data.getVtkData()->GetCellPoints(currentCellID, currentCellPointIds);

					std::vector<vtkIdType> pointIds(data.getVtkData()->GetCell(currentCellID)->GetNumberOfPoints(), 0);
					for (int k = 0; k < data.getVtkData()->GetCell(currentCellID)->GetNumberOfPoints(); k++) {
						pointIds[k] = data.getVtkData()->GetCell(currentCellID)->GetPointId(k);
					}

					bool change = false;

					for (vtkIdType k = 0; k < currentCellPointIds->GetNumberOfIds(); k++) {
						for (vtkIdType l = 0; l < cellPointIds->GetNumberOfIds(); l++) {
							if (currentCellPointIds->GetId(k) == cellPointIds->GetId(l)) {
								change = true;

								unordered_map<vtkIdType, vtkIdType>::iterator iter;
								iter = newPointIDs.find(currentCellPointIds->GetId(k));
								if (iter != newPointIDs.end()) {
									pointIds[k] = newPointIDs.at(currentCellPointIds->GetId(k));
									if (data.getInputType() == DataFormat::PossibleInputType::vtu) {
										dataUGrid->GetCell(currentCellID)->GetPointIds()->SetId(k, newPointIDs.at(currentCellPointIds->GetId(k)));
									}

									if (data.getInputType() == DataFormat::PossibleInputType::vtp) {
										dataPoly->GetCell(currentCellID)->GetPointIds()->SetId(k, newPointIDs.at(currentCellPointIds->GetId(k)));
									}
								}
								else {
									double pointCoords[3] = {
									data.getVtkData()->GetPoint(currentCellPointIds->GetId(k))[0], data.getVtkData()->GetPoint(currentCellPointIds->GetId(k))[1], data.getVtkData()->GetPoint(currentCellPointIds->GetId(k))[2], };

									if (data.getInputType() == DataFormat::PossibleInputType::vtu) {
										vtkIdType newPointID = points->InsertNextPoint(pointCoords);
										pointIds[k] = newPointID;

										// cout << newPointID <<endl;
										std::pair<vtkIdType, vtkIdType> insertIDs(currentCellPointIds->GetId(k), newPointID);
										newPointIDs.insert(insertIDs);
									}

									if (data.getInputType() == DataFormat::PossibleInputType::vtp) {
										vtkIdType newPointID = points->InsertNextPoint(pointCoords);
										pointIds[k] = newPointID;
										std::pair<vtkIdType, vtkIdType> insertIDs(currentCellPointIds->GetId(k), newPointID);
										newPointIDs.insert(insertIDs);
									}
								}
							}
						}
					}

					if (change) {
						if (data.getInputType() == DataFormat::PossibleInputType::vtu) {
							dataUGrid->SetPoints(points);
							dataUGrid->ReplaceCell(currentCellID, (int)data.getVtkData()->GetCell(currentCellID)->GetNumberOfPoints(), &(pointIds[0]));
						}

						if (data.getInputType() == DataFormat::PossibleInputType::vtp) {
							dataPoly->SetPoints(points);
							dataPoly->ReplaceCell(currentCellID, (int)data.getVtkData()->GetCell(currentCellID)->GetNumberOfPoints(), &(pointIds[0]));
						}
					}
				}
			}
		}
	}
	if (data.getInputType() == DataFormat::PossibleInputType::vtu) {
		dataUGrid->BuildLinks();
		data.setVtkData(dataUGrid);
	}

	if (data.getInputType() == DataFormat::PossibleInputType::vtp) {
		dataPoly->BuildLinks();
		data.setVtkData(dataPoly);
	}
	vector<vector<double>> rmCells;
	for (vtkIdType i = 0; i < removedCells->GetNumberOfIds(); i++) {
		vector<double> currentCell = { data.getCentrePoints()->GetPoint(removedCells->GetId(i))[0], data.getCentrePoints()->GetPoint(removedCells->GetId(i))[1], data.getCentrePoints()->GetPoint(removedCells->GetId(i))[2] };
		rmCells.push_back(currentCell);
	}


	data.setRemoveCells(rmCells);
} // Methods::uncoupleMaterials

/*! Disconneted the rifht atrium from the left atrium.
 \param data Pointer to the orginal mesh.
 */
void Methods::uncoupleRightFromLeft(DataFormat &data) {
	vtkSmartPointer<vtkUnstructuredGrid> dataUGrid;
	vtkSmartPointer<vtkPolyData> dataPoly;
	vtkSmartPointer<vtkPoints> points = data.getVtkData()->GetPoints();
	vtkSmartPointer<vtkIdList> removedCells = vtkSmartPointer<vtkIdList>::New();


	if (data.getInputType() == DataFormat::PossibleInputType::vtu) {
		dataUGrid = vtkUnstructuredGrid::SafeDownCast(data.getVtkData());
	}

	if (data.getInputType() == DataFormat::PossibleInputType::vtp) {
		dataPoly = vtkPolyData::SafeDownCast(data.getVtkData());
	}
	unordered_map<vtkIdType, vtkIdType> newPointIDs;

	for (vtkIdType i = 0; i < data.getVtkData()->GetNumberOfCells(); i++) {
		int materialParent = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(i, 0);

		if (Material::isInLeftWithBlood(materialParent)) {
			vtkSmartPointer<vtkIdList> cellPointIds = vtkSmartPointer<vtkIdList>::New();
			vtkSmartPointer<vtkIdList> currentNeighborCellIds = vtkSmartPointer<vtkIdList>::New();

			data.getVtkData()->GetCellPoints(i, cellPointIds);

			for (vtkIdType i = 0; i < cellPointIds->GetNumberOfIds(); i++) {
				vtkSmartPointer<vtkIdList> PointIdList = vtkSmartPointer<vtkIdList>::New();
				PointIdList->InsertNextId(cellPointIds->GetId(i));

				vtkSmartPointer<vtkIdList> tempNeighborCellIds = vtkSmartPointer<vtkIdList>::New();

				data.getVtkData()->GetCellNeighbors(i, PointIdList, tempNeighborCellIds);

				for (vtkIdType j = 0; j < tempNeighborCellIds->GetNumberOfIds(); j++) {
					currentNeighborCellIds->InsertUniqueId(tempNeighborCellIds->GetId(j));
				}
			}
			for (vtkIdType j = 0; j < currentNeighborCellIds->GetNumberOfIds(); j++) {
				vtkIdType currentCellID = currentNeighborCellIds->GetId(j);
				int materialCurrentCell = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(currentCellID, 0);
				if (Material::isInRightWithBlood(materialCurrentCell)) {
					removedCells->InsertUniqueId(currentCellID);
					vtkSmartPointer<vtkIdList> currentCellPointIds = vtkSmartPointer<vtkIdList>::New();
					data.getVtkData()->GetCellPoints(currentCellID, currentCellPointIds);


					std::vector<vtkIdType> pointIds(data.getVtkData()->GetCell(currentCellID)->GetNumberOfPoints(), 0);
					for (int k = 0; k < data.getVtkData()->GetCell(currentCellID)->GetNumberOfPoints(); k++) {
						pointIds[k] = data.getVtkData()->GetCell(currentCellID)->GetPointId(k);
					}

					bool change = false;

					for (vtkIdType k = 0; k < currentCellPointIds->GetNumberOfIds(); k++) {
						for (vtkIdType l = 0; l < cellPointIds->GetNumberOfIds(); l++) {
							if (currentCellPointIds->GetId(k) == cellPointIds->GetId(l)) {
								change = true;

								unordered_map<vtkIdType, vtkIdType>::iterator iter;
								iter = newPointIDs.find(currentCellPointIds->GetId(k));
								if (iter != newPointIDs.end()) {
									pointIds[k] = newPointIDs.at(currentCellPointIds->GetId(k));
									if (data.getInputType() == DataFormat::PossibleInputType::vtu) {
										dataUGrid->GetCell(currentCellID)->GetPointIds()->SetId(k, newPointIDs.at(currentCellPointIds->GetId(k)));
									}

									if (data.getInputType() == DataFormat::PossibleInputType::vtp) {
										dataPoly->GetCell(currentCellID)->GetPointIds()->SetId(k, newPointIDs.at(currentCellPointIds->GetId(k)));
									}
								}
								else {
									double pointCoords[3] = {
									data.getVtkData()->GetPoint(currentCellPointIds->GetId(k))[0], data.getVtkData()->GetPoint(currentCellPointIds->GetId(k))[1], data.getVtkData()->GetPoint(currentCellPointIds->GetId(k))[2], };

									if (data.getInputType() == DataFormat::PossibleInputType::vtu) {
										vtkIdType newPointID = points->InsertNextPoint(pointCoords);
										pointIds[k] = newPointID;

										// cout << newPointID <<endl;
										std::pair<vtkIdType, vtkIdType> insertIDs(currentCellPointIds->GetId(k), newPointID);
										newPointIDs.insert(insertIDs);
									}

									if (data.getInputType() == DataFormat::PossibleInputType::vtp) {
										vtkIdType newPointID = points->InsertNextPoint(pointCoords);
										pointIds[k] = newPointID;
										std::pair<vtkIdType, vtkIdType> insertIDs(currentCellPointIds->GetId(k), newPointID);
										newPointIDs.insert(insertIDs);
									}
								}
							}
						}
					}

					if (change) {
						if (data.getInputType() == DataFormat::PossibleInputType::vtu) {
							dataUGrid->SetPoints(points);
							dataUGrid->ReplaceCell(currentCellID, (int)data.getVtkData()->GetCell(currentCellID)->GetNumberOfPoints(), &(pointIds[0]));
						}

						if (data.getInputType() == DataFormat::PossibleInputType::vtp) {
							dataPoly->SetPoints(points);
							dataPoly->ReplaceCell(currentCellID, (int)data.getVtkData()->GetCell(currentCellID)->GetNumberOfPoints(), &(pointIds[0]));
						}
					}
				}
			}
		}
	}
	if (data.getInputType() == DataFormat::PossibleInputType::vtu) {
		dataUGrid->BuildLinks();
		data.setVtkData(dataUGrid);
	}

	if (data.getInputType() == DataFormat::PossibleInputType::vtp) {
		dataPoly->BuildLinks();
		data.setVtkData(dataPoly);
	}

	vector<vector<double>> rmCells;
	for (vtkIdType i = 0; i < removedCells->GetNumberOfIds(); i++) {
		vector<double> currentCell = { data.getCentrePoints()->GetPoint(removedCells->GetId(i))[0], data.getCentrePoints()->GetPoint(removedCells->GetId(i))[1], data.getCentrePoints()->GetPoint(removedCells->GetId(i))[2] };
		rmCells.push_back(currentCell);
	}


	data.setRemoveCells(rmCells);
} // Methods::uncoupleRightFromLeft


/*! Get the max distance between a path and a point.
 \param data Pointer to the orginal mesh.
 \param path The Path
 \param point Point as vecotr (x,y,z)
 \return Max distance betwwen path and point.
 */
double Methods::getMaxDistance(DataFormat &data, vtkSmartPointer<vtkIdList> path, vector<double> point) {
	double pointArray[3] = { point.at(0), point.at(1), point.at(2) };
	double maxDistance = -1;
	vtkIdType maxID;

	for (vtkIdType i = 0; i < path->GetNumberOfIds(); i++) {
		double tempdistance = sqrt(abs(vtkMath::Distance2BetweenPoints(pointArray, data.getCentrePoints()->GetPoint(path->GetId(i)))));
		if (tempdistance > maxDistance) {
			maxDistance = tempdistance;
		}
	}
	return maxDistance;
}

/*! Get all cell neighbours form one cell that are connected to that.
 \param data Pointer to the orginal mesh.
 \param cellID The cell id
 \return List with the neighbour cells.
 */
vtkSmartPointer<vtkIdList> Methods::getAllCellNeighbours(DataFormat &data, vtkIdType cellID) {
	vtkSmartPointer<vtkIdList> neighbourCells = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkIdList> cellPointIds = vtkSmartPointer<vtkIdList>::New();

	bool end = false;
	set<vtkIdType> viewedCells;

	vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();

	cellIds->InsertNextId(cellID);

	while (!end) {
		vtkSmartPointer<vtkIdList> newCellIds = vtkSmartPointer<vtkIdList>::New();

		for (vtkIdType i = 0; i < cellIds->GetNumberOfIds(); i++) {
			vtkIdType cellId = cellIds->GetId(i);
			data.getVtkData()->GetCellPoints(cellId, cellPointIds);

			for (vtkIdType j = 0; j < cellPointIds->GetNumberOfIds(); j++) {
				vtkSmartPointer<vtkIdList> CellPointIdList = vtkSmartPointer<vtkIdList>::New();
				CellPointIdList->InsertNextId(cellPointIds->GetId(j));

				vtkSmartPointer<vtkIdList> tempNeighborCellIds = vtkSmartPointer<vtkIdList>::New();
				data.getVtkData()->GetCellNeighbors(cellId, CellPointIdList, tempNeighborCellIds);

				for (vtkIdType k = 0; k < tempNeighborCellIds->GetNumberOfIds(); k++) {
					if (viewedCells.find(tempNeighborCellIds->GetId(k)) == viewedCells.end()) {
						neighbourCells->InsertNextId(tempNeighborCellIds->GetId(k));
						newCellIds->InsertNextId(tempNeighborCellIds->GetId(k));
						viewedCells.insert(tempNeighborCellIds->GetId(k));
					}
				}
			}
		}
		cellIds = newCellIds;

		if (cellIds->GetNumberOfIds() == 0) {
			end = true;
		}
	}

	return neighbourCells;
} // Methods::getAllCellNeighbours

/*! Get the direct cell neighbours form one cell.
 \param data Pointer to the orginal mesh.
 \param cellID The cell id
 \return List with the neighbour cells.
 */
vtkSmartPointer<vtkIdList> Methods::getCellNeighbours(DataFormat &data, vtkIdType cellID) {
	vtkSmartPointer<vtkIdList> neighbourCells = vtkSmartPointer<vtkIdList>::New();

	vtkSmartPointer<vtkIdList> cellPointIds = vtkSmartPointer<vtkIdList>::New();

	data.getVtkData()->GetCellPoints(cellID, cellPointIds);

	for (vtkIdType i = 0; i < cellPointIds->GetNumberOfIds(); i++) {
		vtkSmartPointer<vtkIdList> CellPointIdList = vtkSmartPointer<vtkIdList>::New();
		CellPointIdList->InsertNextId(cellPointIds->GetId(i));

		vtkSmartPointer<vtkIdList> tempNeighborCellIds = vtkSmartPointer<vtkIdList>::New();
		data.getVtkData()->GetCellNeighbors(cellID, CellPointIdList, tempNeighborCellIds);

		for (vtkIdType j = 0; j < tempNeighborCellIds->GetNumberOfIds(); j++) {
			neighbourCells->InsertUniqueId(tempNeighborCellIds->GetId(j));
		}
	}

	return neighbourCells;
}

/*! Get all cells in a radius of a cell.
 \param data Pointer to the orginal mesh.
 \param cellID The cell id
 \param radius Radius in mm
 \return List with the neighbour cells.
 */
vtkSmartPointer<vtkIdList> Methods::getCellsinRadius(DataFormat &data, vtkIdType cellID, double radius) {
	vector<double> point = { data.getCentrePoints()->GetPoint(cellID)[0], data.getCentrePoints()->GetPoint(cellID)[1], data.getCentrePoints()->GetPoint(cellID)[2] };

	return getCellsinRadius(data, point, radius);
}

/*! Get all cells in a radius of a point.
 \param data Pointer to the orginal mesh.
 \param point Point as vector (x,y,z)
 \param radius Radius in mm
 \return List with the neighbour cells.
 */
vtkSmartPointer<vtkIdList> Methods::getCellsinRadius(DataFormat &data, vector<double> point, double radius) {
	vtkSmartPointer<vtkUnstructuredGrid> centrePointsUGrid = data.getCentrePoints();
	vtkSmartPointer<vtkPointLocator> centrePointLocator = vtkSmartPointer<vtkPointLocator>::New();

	centrePointLocator->SetDataSet(centrePointsUGrid);
	centrePointLocator->BuildLocator();

	double pointArray[3] = { point.at(0), point.at(1), point.at(2) };

	vtkSmartPointer<vtkIdList> neighbourCells = vtkSmartPointer<vtkIdList>::New();

	centrePointLocator->FindPointsWithinRadius(radius, pointArray, neighbourCells);

	return neighbourCells;
}

/*! Test if a point is in a material.
 \param data Pointer to the orginal mesh.
 \param point Point as vector (x,y,z)
 \param inMaterial Tissue class in tha the point must lie.
 \return true if it is or false.
 */
bool Methods::isPointinMaterial(DataFormat &data, vector<double> point, Material::Mat inMaterial) {
	double pointArray[3] = { point.at(0), point.at(1), point.at(2) };

	vtkIdType cellId = data.getCentrePoints()->FindPoint(pointArray);

	if (!data.getDoubleLayer()) {
		int material = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(cellId, 0);

		if (material == inMaterial) {
			cout << "The point: " << point.at(0) << ", " << point.at(1) << ", " << point.at(2) << " is in material: " <<
				inMaterial << endl;
			return true;
		}
	}
	int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(cellId, 0);

	if (material == inMaterial) {
		cout << "The point: " << point.at(0) << ", " << point.at(1) << ", " << point.at(2) << " is in material: " <<
			inMaterial << endl;
		return true;
	}
	cout << "The point: " << point.at(0) << ", " << point.at(1) << ", " << point.at(2) << " isn't in material: " <<
		inMaterial << endl;

	return false;
} // Methods::isPointinMaterial

/*!Find the closest point in a material without a fiber orientation.
 \param data Pointer to the orginal mesh.
 \param point Point as vector (x,y,z)
 \param inMaterial Tissue class in tha the point must lie.
 \return The point who would be found as vector (x,y,z).
 */
vector<double> Methods::findClosedPointinMaterialwithoutOrientation(DataFormat &data, vector<double> point, Material::Mat inMaterial) {
	return findClosedPointinMaterialwithoutOrientation(data, point, inMaterial, inMaterial);
}

/*!Find the closest point in a material without a fiber orientation.
 \param data Pointer to the orginal mesh.
 \param point Point as vector (x,y,z)
 \param inMaterial1 First tissue class in tha the point must lie.
 \param inMaterial2 Second tissue class in tha the point must lie.
 \return The point who would be found as vector (x,y,z).
 */
vector<double> Methods::findClosedPointinMaterialwithoutOrientation(DataFormat &data, vector<double> point, Material::Mat inMaterial1, Material::Mat inMaterial2) {
	double radius = 15;
	double pointArray[3] = { point.at(0), point.at(1), point.at(2) };

	double mindistance = std::numeric_limits<double>::max();

	vtkIdType cellId = data.getCentrePoints()->FindPoint(pointArray);

	if (!data.getDoubleLayer()) {
		int material = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(cellId, 0);

		if (((material == inMaterial1) || (material == inMaterial2)) &&
			(data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->GetComponent(cellId, 0) +
				data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->GetComponent(cellId, 1) +
				data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->GetComponent(cellId, 2) == 0)) {
			point = { data.getCentrePoints()->GetPoint(cellId)[0], data.getCentrePoints()->GetPoint(cellId)[1], data.getCentrePoints()->GetPoint(cellId)[2] };
			return point;
		}
	}
	int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(cellId, 0);

	if (((material == inMaterial1) || (material == inMaterial2)) &&
		(data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetComponent(cellId, 0) +
			data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetComponent(cellId, 1) +
			data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetComponent(cellId, 2) == 0)) {
		point = { data.getCentrePoints()->GetPoint(cellId)[0], data.getCentrePoints()->GetPoint(cellId)[1], data.getCentrePoints()->GetPoint(cellId)[2] };
		return point;
	}

	vtkSmartPointer<vtkIdList> neighbours = getCellsinRadius(data, cellId, radius);

	vtkIdType newCellId = -1;

	vector<double> newPoint(3);
	newPoint.at(0) = NAN;
	newPoint.at(1) = NAN;
	newPoint.at(2) = NAN;


	for (vtkIdType i = 0; i < neighbours->GetNumberOfIds(); i++) {
		double distance = sqrt(abs(vtkMath::Distance2BetweenPoints(pointArray, data.getCentrePoints()->GetPoint(neighbours->GetId(i)))));
		if (!data.getDoubleLayer()) {
			int material = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(neighbours->GetId(i), 0);

			if ((distance < mindistance) && ((material == inMaterial1) || (material == inMaterial2)) &&
				(data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->GetComponent(neighbours->GetId(i), 0) +
					data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->GetComponent(neighbours->GetId(i), 1) +
					data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->GetComponent(neighbours->GetId(i), 2) == 0)) {
				mindistance = distance;
				newCellId = neighbours->GetId(i);
			}
		}

		int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(neighbours->GetId(i), 0);

		if ((distance < mindistance) && ((material == inMaterial1) || (material == inMaterial2)) &&
			(data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetComponent(neighbours->GetId(i), 0) +
				data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetComponent(neighbours->GetId(i), 1) +
				data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetComponent(neighbours->GetId(i), 2) == 0)) {
			mindistance = distance;
			newCellId = neighbours->GetId(i);
		}
	}

	if (newCellId != -1) {
		double *newPointArray = data.getCentrePoints()->GetPoint(newCellId);
		newPoint.at(0) = newPointArray[0];
		newPoint.at(1) = newPointArray[1];
		newPoint.at(2) = newPointArray[2];
		cout << "The point: " << point.at(0) << ", " << point.at(1) << ", " << point.at(2) << " was moved in: " <<
			newPoint.at(0) << ", " << newPoint.at(1) << ", " << newPoint.at(2) << endl;
		return newPoint;
	}
	else {
		cout << "The point: " << point.at(0) << ", " << point.at(1) << ", " << point.at(2) << " can't be moved!!!" << endl;
		cout << "Point in Material " << inMaterial1 << "or " << inMaterial2 << " was not found!!!" << endl;
		return newPoint;
	}
} // Methods::findClosedPointinMaterialwithoutOrientation

/*!Find the closest cell in a material without a fiber orientation.
 \param data Pointer to the orginal mesh.
 \param point Point as vector (x,y,z)
 \param inMaterial Tissue class in tha the point must lie.
 \return The cell who would be found as cell id.
 */
vtkIdType Methods::findClosedPointIdinMaterialwithoutOrientation(DataFormat &data, vector<double> point, Material::Mat inMaterial) {
	return findClosedPointIdinMaterialwithoutOrientation(data, point, inMaterial, inMaterial);
}

/*!Find the closest cell in a material without a fiber orientation.
 \param data Pointer to the orginal mesh.
 \param point Point as vector (x,y,z)
 \param inMaterial1 First tissue class in tha the point must lie.
 \param inMaterial2 Second tissue class in tha the point must lie.
 \return The cell who would be found as cell id.
 */
vtkIdType Methods::findClosedPointIdinMaterialwithoutOrientation(DataFormat &data, vector<double> point, Material::Mat inMaterial1, Material::Mat inMaterial2) {
	double radius = 15;
	double pointArray[3] = { point.at(0), point.at(1), point.at(2) };

	double mindistance = std::numeric_limits<double>::max();

	vtkIdType cellId = data.getCentrePoints()->FindPoint(pointArray);

	if (!data.getDoubleLayer()) {
		int material = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(cellId, 0);

		if (((material == inMaterial1) || (material == inMaterial2)) &&
			(data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->GetComponent(cellId, 0) +
				data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->GetComponent(cellId, 1) +
				data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->GetComponent(cellId, 2) == 0)) {
			return cellId;
		}
	}
	int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(cellId, 0);

	if (((material == inMaterial1) || (material == inMaterial2)) &&
		(data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetComponent(cellId, 0) +
			data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetComponent(cellId, 1) +
			data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetComponent(cellId, 2) == 0)) {
		return cellId;
	}

	vtkSmartPointer<vtkIdList> neighbours = getCellsinRadius(data, cellId, radius);

	vtkIdType newCellId = -1;

	for (vtkIdType i = 0; i < neighbours->GetNumberOfIds(); i++) {
		double distance = sqrt(abs(vtkMath::Distance2BetweenPoints(pointArray, data.getCentrePoints()->GetPoint(neighbours->GetId(i)))));
		if (!data.getDoubleLayer()) {
			int material = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(neighbours->GetId(i), 0);

			if ((distance < mindistance) && ((material == inMaterial1) || (material == inMaterial2)) &&
				(data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->GetComponent(neighbours->GetId(i), 0) +
					data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->GetComponent(neighbours->GetId(i), 1) +
					data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->GetComponent(neighbours->GetId(i), 2) == 0)) {
				mindistance = distance;
				newCellId = neighbours->GetId(i);
			}
		}

		int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(neighbours->GetId(i), 0);

		if ((distance < mindistance) && ((material == inMaterial1) || (material == inMaterial2)) &&
			(data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetComponent(neighbours->GetId(i), 0) +
				data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetComponent(neighbours->GetId(i), 1) +
				data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetComponent(neighbours->GetId(i), 2) == 0)) {
			mindistance = distance;
			newCellId = neighbours->GetId(i);
		}
	}

	return newCellId;
} // Methods::findClosedPointIdinMaterialwithoutOrientation

/*!Find the closest point in the right epicard.
 \param data Pointer to the orginal mesh.
 \param point Point as vector (x,y,z)
 \return The point who would be found as vector (x,y,z).
 */
vector<double> Methods::findClosedPointinMaterialInRightEpi(DataFormat &data, vector<double> point) {
	vtkIdType newCellId = findClosedPointIdinMaterialInRightEpi(data, point);

	vector<double> newPoint(3);

	newPoint.at(0) = NAN;
	newPoint.at(1) = NAN;
	newPoint.at(2) = NAN;


	if (newCellId != -1) {
		double *newPointArray = data.getCentrePoints()->GetPoint(newCellId);
		newPoint.at(0) = newPointArray[0];
		newPoint.at(1) = newPointArray[1];
		newPoint.at(2) = newPointArray[2];
		cout << "The point: " << point.at(0) << ", " << point.at(1) << ", " << point.at(2) << " was moved in: " <<
			newPoint.at(0) << ", " << newPoint.at(1) << ", " << newPoint.at(2) << endl;
		return newPoint;
	}
	else {
		cout << "The point: " << point.at(0) << ", " << point.at(1) << ", " << point.at(2) << " can't be moved!!!" << endl;
		cout << "Point in right epi material was not found!!!" << endl;
		return newPoint;
	}
}

/*!Find the closest point in the left epicard.
 \param data Pointer to the orginal mesh.
 \param point Point as vector (x,y,z)
 \return The point who would be found as vector (x,y,z).
 */
vector<double> Methods::findClosedPointinMaterialInLeftEpi(DataFormat &data, vector<double> point) {
	vtkIdType newCellId = findClosedPointIdinMaterialInLeftEpi(data, point);

	vector<double> newPoint(3);

	newPoint.at(0) = NAN;
	newPoint.at(1) = NAN;
	newPoint.at(2) = NAN;


	if (newCellId != -1) {
		double *newPointArray = data.getCentrePoints()->GetPoint(newCellId);
		newPoint.at(0) = newPointArray[0];
		newPoint.at(1) = newPointArray[1];
		newPoint.at(2) = newPointArray[2];
		cout << "The point: " << point.at(0) << ", " << point.at(1) << ", " << point.at(2) << " was moved in: " <<
			newPoint.at(0) << ", " << newPoint.at(1) << ", " << newPoint.at(2) << endl;
		return newPoint;
	}
	else {
		cout << "The point: " << point.at(0) << ", " << point.at(1) << ", " << point.at(2) << " can't be moved!!!" << endl;
		cout << "Point in left material was not found!!!" << endl;
		return newPoint;
	}
}

/*!Find the closest point in the left endocard.
 \param data Pointer to the orginal mesh.
 \param point Point as vector (x,y,z)
 \return The point who would be found as vector (x,y,z).
 */
vector<double> Methods::findClosedPointinMaterialInLeftEndo(DataFormat &data, vector<double> point) {
	vtkIdType newCellId = findClosedPointIdinMaterialInLeftEndo(data, point);

	vector<double> newPoint(3);

	newPoint.at(0) = NAN;
	newPoint.at(1) = NAN;
	newPoint.at(2) = NAN;


	if (newCellId != -1) {
		double *newPointArray = data.getCentrePoints()->GetPoint(newCellId);
		newPoint.at(0) = newPointArray[0];
		newPoint.at(1) = newPointArray[1];
		newPoint.at(2) = newPointArray[2];
		cout << "The point: " << point.at(0) << ", " << point.at(1) << ", " << point.at(2) << " was moved in: " <<
			newPoint.at(0) << ", " << newPoint.at(1) << ", " << newPoint.at(2) << endl;
		return newPoint;
	}
	else {
		cout << "The point: " << point.at(0) << ", " << point.at(1) << ", " << point.at(2) << " can't be moved!!!" << endl;
		cout << "Point in left endo material was not found!!!" << endl;
		return newPoint;
	}
}

/*!Find the closest point in the right endocard.
 \param data Pointer to the orginal mesh.
 \param point Point as vector (x,y,z)
 \return The point who would be found as vector (x,y,z).
 */
vector<double> Methods::findClosedPointinMaterialInRightEndo(DataFormat &data, vector<double> point) {
	vtkIdType newCellId = findClosedPointIdinMaterialInRightEndo(data, point);

	vector<double> newPoint(3);

	newPoint.at(0) = NAN;
	newPoint.at(1) = NAN;
	newPoint.at(2) = NAN;


	if (newCellId != -1) {
		double *newPointArray = data.getCentrePoints()->GetPoint(newCellId);
		newPoint.at(0) = newPointArray[0];
		newPoint.at(1) = newPointArray[1];
		newPoint.at(2) = newPointArray[2];
		cout << "The point: " << point.at(0) << ", " << point.at(1) << ", " << point.at(2) << " was moved in: " <<
			newPoint.at(0) << ", " << newPoint.at(1) << ", " << newPoint.at(2) << endl;
		return newPoint;
	}
	else {
		cout << "The point: " << point.at(0) << ", " << point.at(1) << ", " << point.at(2) << " can't be moved!!!" << endl;
		cout << "Point in right endo material was not found!!!" << endl;

		return newPoint;
	}
}

/*!Find the closest point in a region of the right atrium.
 \param data Pointer to the orginal mesh.
 \param point Point as vector (x,y,z)
 \param region Region in that the point could lie.
 \return The point who would be found as vector (x,y,z).
 */
vector<double> Methods::findClosedPointinMaterialInRightInRegion(DataFormat &data, vector<double> point, vtkSmartPointer<vtkIdList> region) {
	vtkIdType newCellId = findClosedPointIdinMaterialInRightInRegion(data, point, region);

	vector<double> newPoint(3);

	newPoint.at(0) = NAN;
	newPoint.at(1) = NAN;
	newPoint.at(2) = NAN;


	if (newCellId != -1) {
		double *newPointArray = data.getCentrePoints()->GetPoint(newCellId);
		newPoint.at(0) = newPointArray[0];
		newPoint.at(1) = newPointArray[1];
		newPoint.at(2) = newPointArray[2];
		cout << "The point: " << point.at(0) << ", " << point.at(1) << ", " << point.at(2) << " was moved in: " <<
			newPoint.at(0) << ", " << newPoint.at(1) << ", " << newPoint.at(2) << endl;
		return newPoint;
	}
	else {
		cout << "The point: " << point.at(0) << ", " << point.at(1) << ", " << point.at(2) << " can't be moved!!!" << endl;
		cout << "Point in right material was not found!!!" << endl;
		return newPoint;
	}
}

/*!Find the closest point in a region of the left atrium.
 \param data Pointer to the orginal mesh.
 \param point Point as vector (x,y,z)
 \param region Region in that the point could lie.
 \return The point who would be found as vector (x,y,z).
 */
vector<double> Methods::findClosedPointinMaterialInLeftInRegion(DataFormat &data, vector<double> point, vtkSmartPointer<vtkIdList> region) {
	vtkIdType newCellId = findClosedPointIdinMaterialInLeftInRegion(data, point, region);

	vector<double> newPoint(3);

	newPoint.at(0) = NAN;
	newPoint.at(1) = NAN;
	newPoint.at(2) = NAN;

	if (newCellId != -1) {
		double *newPointArray = data.getCentrePoints()->GetPoint(newCellId);
		newPoint.at(0) = newPointArray[0];
		newPoint.at(1) = newPointArray[1];
		newPoint.at(2) = newPointArray[2];
		cout << "The point: " << point.at(0) << ", " << point.at(1) << ", " << point.at(2) << " was moved in: " <<
			newPoint.at(0) << ", " << newPoint.at(1) << ", " << newPoint.at(2) << endl;
		return newPoint;
	}
	else {
		cout << "The point: " << point.at(0) << ", " << point.at(1) << ", " << point.at(2) << " can't be moved!!!" << endl;
		cout << "Point in left material was not found!!!" << endl;
		return newPoint;
	}
}

/*!Find the closest cell in the right atrial.
 \param data Pointer to the orginal mesh.
 \param point Point as vector (x,y,z)
 \return The cell who would be found as cell id.
 */
vtkIdType Methods::findClosedPointIdinMaterialInRight(DataFormat &data, vector<double> point) {
	double radius = 15;
	double pointArray[3] = { point.at(0), point.at(1), point.at(2) };
	double mindistance = std::numeric_limits<double>::max();

	vtkIdType cellId = data.getCentrePoints()->FindPoint(pointArray);

	if (!data.getDoubleLayer()) {
		int material = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(cellId, 0);
		if (Material::isInRight(material)) {
			return cellId;
		}
	}
	int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(cellId, 0);
	if (Material::isInRight(material)) {
		return cellId;
	}

	vtkSmartPointer<vtkIdList> neighbours = getCellsinRadius(data, cellId, radius);

	vtkIdType newCellId = -1;


	for (vtkIdType i = 0; i < neighbours->GetNumberOfIds(); i++) {
		double distance = sqrt(abs(vtkMath::Distance2BetweenPoints(pointArray, data.getCentrePoints()->GetPoint(neighbours->GetId(i)))));
		if (!data.getDoubleLayer()) {
			int material = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(neighbours->GetId(i), 0);

			if ((distance < mindistance) && Material::isInRight(material)) {
				mindistance = distance;
				newCellId = neighbours->GetId(i);
			}
		}

		int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(neighbours->GetId(i), 0);

		if ((distance < mindistance) && Material::isInRight(material)) {
			mindistance = distance;
			newCellId = neighbours->GetId(i);
		}
	}

	if (newCellId != -1) {
		return newCellId;
	}
	else {
		cout << "Point in right material was not found!!!" << endl;
		return newCellId;
	}
} // Methods::findClosedPointIdinMaterialInRight

/*!Find the closest cell in the right epicard.
 \param data Pointer to the orginal mesh.
 \param point Point as vector (x,y,z)
 \return The cell who would be found as cell id.
 */
vtkIdType Methods::findClosedPointIdinMaterialInRightEpi(DataFormat &data, vector<double> point) {
	double radius = 15;
	double pointArray[3] = { point.at(0), point.at(1), point.at(2) };
	double mindistance = std::numeric_limits<double>::max();

	vtkIdType cellId = data.getCentrePoints()->FindPoint(pointArray);

	if (!data.getDoubleLayer()) {
		int material = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(cellId, 0);
		if (Material::isInRightWithoutEndo(material)) {
			return cellId;
		}
	}
	int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(cellId, 0);
	if (Material::isInRightWithoutEndo(material)) {
		return cellId;
	}

	vtkSmartPointer<vtkIdList> neighbours = getCellsinRadius(data, cellId, radius);

	vtkIdType newCellId = -1;


	for (vtkIdType i = 0; i < neighbours->GetNumberOfIds(); i++) {
		double distance = sqrt(abs(vtkMath::Distance2BetweenPoints(pointArray, data.getCentrePoints()->GetPoint(neighbours->GetId(i)))));
		if (!data.getDoubleLayer()) {
			int material = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(neighbours->GetId(i), 0);

			if ((distance < mindistance) && Material::isInRightWithoutEndo(material)) {
				mindistance = distance;
				newCellId = neighbours->GetId(i);
			}
		}

		int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(neighbours->GetId(i), 0);

		if ((distance < mindistance) && Material::isInRightWithoutEndo(material)) {
			mindistance = distance;
			newCellId = neighbours->GetId(i);
		}
	}

	if (newCellId != -1) {
		return newCellId;
	}
	else {
		cout << "Point in right material was not found!!!" << endl;
		return newCellId;
	}
} // Methods::findClosedPointIdinMaterialInRightEpi

/*!Find the closest cell in the epicard.
 \param data Pointer to the orginal mesh.
 \param point Point as vector (x,y,z)
 \return The cell who would be found as cell id.
 */
vtkIdType Methods::findClosedPointIdinMaterialInEpi(DataFormat &data, vector<double> point) {
	double radius = 15;
	double pointArray[3] = { point.at(0), point.at(1), point.at(2) };
	double mindistance = std::numeric_limits<double>::max();

	vtkIdType cellId = data.getCentrePoints()->FindPoint(pointArray);

	if (!data.getDoubleLayer()) {
		int material = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(cellId, 0);
		if (Material::isInRightWithoutEndo(material) || Material::isInLeftEpi(material)) {
			return cellId;
		}
	}
	int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(cellId, 0);
	if (Material::isInRightWithoutEndo(material) || Material::isInLeftEpi(material)) {
		return cellId;
	}

	vtkSmartPointer<vtkIdList> neighbours = getCellsinRadius(data, cellId, radius);

	vtkIdType newCellId = -1;


	for (vtkIdType i = 0; i < neighbours->GetNumberOfIds(); i++) {
		double distance = sqrt(abs(vtkMath::Distance2BetweenPoints(pointArray, data.getCentrePoints()->GetPoint(neighbours->GetId(i)))));
		if (!data.getDoubleLayer()) {
			int material = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(neighbours->GetId(i), 0);

			if ((distance < mindistance) && (Material::isInRightWithoutEndo(material) || Material::isInLeftEpi(material))) {
				mindistance = distance;
				newCellId = neighbours->GetId(i);
			}
		}

		int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(neighbours->GetId(i), 0);

		if ((distance < mindistance) && (Material::isInRightWithoutEndo(material) || Material::isInLeftEpi(material))) {
			mindistance = distance;
			newCellId = neighbours->GetId(i);
		}
	}

	if (newCellId != -1) {
		return newCellId;
	}
	else {
		return newCellId;
	}
} // Methods::findClosedPointIdinMaterialInEpi

/*!Find the closest cell in the endocard.
 \param data Pointer to the orginal mesh.
 \param point Point as vector (x,y,z)
 \return The cell who would be found as cell id.
 */
vtkIdType Methods::findClosedPointIdinMaterialInEndo(DataFormat &data, vector<double> point) {
	double radius = 15;
	double pointArray[3] = { point.at(0), point.at(1), point.at(2) };
	double mindistance = std::numeric_limits<double>::max();

	vtkIdType cellId = data.getCentrePoints()->FindPoint(pointArray);

	if (!data.getDoubleLayer()) {
		int material = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(cellId, 0);
		if ((material == Material::subendo_right_Atrial_Cardium) || Material::isInLeftEndo(material)) {
			return cellId;
		}
	}
	int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(cellId, 0);
	if ((material == Material::subendo_right_Atrial_Cardium) || Material::isInLeftEndo(material)) {
		return cellId;
	}

	vtkSmartPointer<vtkIdList> neighbours = getCellsinRadius(data, cellId, radius);

	vtkIdType newCellId = -1;


	for (vtkIdType i = 0; i < neighbours->GetNumberOfIds(); i++) {
		double distance = sqrt(abs(vtkMath::Distance2BetweenPoints(pointArray, data.getCentrePoints()->GetPoint(neighbours->GetId(i)))));
		if (!data.getDoubleLayer()) {
			int material = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(neighbours->GetId(i), 0);

			if ((distance < mindistance) &&
				((material == Material::subendo_right_Atrial_Cardium) || Material::isInLeftEndo(material))) {
				mindistance = distance;
				newCellId = neighbours->GetId(i);
			}
		}

		int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(neighbours->GetId(i), 0);

		if ((distance < mindistance) &&
			((material == Material::subendo_right_Atrial_Cardium) || Material::isInLeftEndo(material))) {
			mindistance = distance;
			newCellId = neighbours->GetId(i);
		}
	}

	if (newCellId != -1) {
		return newCellId;
	}
	else {
		return newCellId;
	}
} // Methods::findClosedPointIdinMaterialInEndo

/*!Find the closest cell in the left atrial.
 \param data Pointer to the orginal mesh.
 \param point Point as vector (x,y,z)
 \return The cell who would be found as cell id.
 */
vtkIdType Methods::findClosedPointIdinMaterialInLeft(DataFormat &data, vector<double> point) {
	double radius = 15;
	double pointArray[3] = { point.at(0), point.at(1), point.at(2) };
	double mindistance = std::numeric_limits<double>::max();


	vtkIdType cellId = data.getCentrePoints()->FindPoint(pointArray);

	if (!data.getDoubleLayer()) {
		int material = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(cellId, 0);
		if (Material::isInLeft(material)) {
			return cellId;
		}
	}
	int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(cellId, 0);
	if (Material::isInLeft(material)) {
		return cellId;
	}

	vtkSmartPointer<vtkIdList> neighbours = getCellsinRadius(data, cellId, radius);

	vtkIdType newCellId = -1;

	for (vtkIdType i = 0; i < neighbours->GetNumberOfIds(); i++) {
		double distance = sqrt(abs(vtkMath::Distance2BetweenPoints(pointArray, data.getCentrePoints()->GetPoint(neighbours->GetId(i)))));
		if (!data.getDoubleLayer()) {
			int material = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(neighbours->GetId(i), 0);
			if ((distance < mindistance) && Material::isInLeft(material)) {
				mindistance = distance;
				newCellId = neighbours->GetId(i);
			}
		}

		int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(neighbours->GetId(i), 0);
		if ((distance < mindistance) && Material::isInLeft(material)) {
			mindistance = distance;
			newCellId = neighbours->GetId(i);
		}
	}

	if (newCellId != -1) {
		return newCellId;
	}
	else {
		cout << "Point in left material was not found!!!" << endl;
		return newCellId;
	}
} // Methods::findClosedPointIdinMaterialInLeft

/*!Find the closest cell in the left epicard.
 \param data Pointer to the orginal mesh.
 \param point Point as vector (x,y,z)
 \return The cell who would be found as cell id.
 */
vtkIdType Methods::findClosedPointIdinMaterialInLeftEpi(DataFormat &data, vector<double> point) {
	double radius = 15;
	double pointArray[3] = { point.at(0), point.at(1), point.at(2) };
	double mindistance = std::numeric_limits<double>::max();


	vtkIdType cellId = data.getCentrePoints()->FindPoint(pointArray);

	int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(cellId, 0);

	if (Material::isInLeftEpi(material)) {
		return cellId;
	}

	vtkSmartPointer<vtkIdList> neighbours = getCellsinRadius(data, cellId, radius);

	vtkIdType newCellId = -1;

	for (vtkIdType i = 0; i < neighbours->GetNumberOfIds(); i++) {
		double distance = sqrt(abs(vtkMath::Distance2BetweenPoints(pointArray, data.getCentrePoints()->GetPoint(neighbours->GetId(i)))));

		int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(neighbours->GetId(i), 0);
		if ((distance < mindistance) && Material::isInLeftEpi(material)) {
			mindistance = distance;
			newCellId = neighbours->GetId(i);
		}
	}

	if (newCellId != -1) {
		return newCellId;
	}
	else {
		cout << "Point in left epi material was not found!!!" << endl;
		return newCellId;
	}
} // Methods::findClosedPointIdinMaterialInLeftEpi

/*!Find the closest cell in the left endocard.
 \param data Pointer to the orginal mesh.
 \param point Point as vector (x,y,z)
 \return The cell who would be found as cell id.
 */
vtkIdType Methods::findClosedPointIdinMaterialInLeftEndo(DataFormat &data, vector<double> point) {
	double radius = 15;
	double pointArray[3] = { point.at(0), point.at(1), point.at(2) };
	double mindistance = std::numeric_limits<double>::max();

	vtkIdType cellId = data.getCentrePoints()->FindPoint(pointArray);

	if (!data.getDoubleLayer()) {
		int material = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(cellId, 0);
		if (Material::isInLeftEndo(material)) {
			return cellId;
		}
	}

	int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(cellId, 0);
	if (Material::isInLeftEndo(material)) {
		return cellId;
	}

	vtkSmartPointer<vtkIdList> neighbours = getCellsinRadius(data, cellId, radius);

	vtkIdType newCellId = -1;

	for (vtkIdType i = 0; i < neighbours->GetNumberOfIds(); i++) {
		double distance = sqrt(abs(vtkMath::Distance2BetweenPoints(pointArray, data.getCentrePoints()->GetPoint(neighbours->GetId(i)))));

		if (!data.getDoubleLayer()) {
			int material = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(neighbours->GetId(i), 0);
			if ((distance < mindistance) && Material::isInLeftEndo(material)) {
				mindistance = distance;
				newCellId = neighbours->GetId(i);
			}
		}

		int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(neighbours->GetId(i), 0);
		if ((distance < mindistance) && Material::isInLeftEndo(material)) {
			mindistance = distance;
			newCellId = neighbours->GetId(i);
		}
	}

	if (newCellId != -1) {
		return newCellId;
	}
	else {
		cout << "Point in left endo material was not found!!!" << endl;
		return newCellId;
	}
} // Methods::findClosedPointIdinMaterialInLeftEndo

/*!Find the closest cell in the right endocard.
 \param data Pointer to the orginal mesh.
 \param point Point as vector (x,y,z)
 \return The cell who would be found as cell id.
 */
vtkIdType Methods::findClosedPointIdinMaterialInRightEndo(DataFormat &data, vector<double> point) {
	double radius = 15;
	double pointArray[3] = { point.at(0), point.at(1), point.at(2) };
	double mindistance = std::numeric_limits<double>::max();

	vtkIdType cellId = data.getCentrePoints()->FindPoint(pointArray);

	if (!data.getDoubleLayer()) {
		int material = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(cellId, 0);
		if (material == Material::subendo_right_Atrial_Cardium) {
			return cellId;
		}
	}

	int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(cellId, 0);
	if (material == Material::subendo_right_Atrial_Cardium) {
		return cellId;
	}

	vtkSmartPointer<vtkIdList> neighbours = getCellsinRadius(data, cellId, radius);

	vtkIdType newCellId = -1;

	for (vtkIdType i = 0; i < neighbours->GetNumberOfIds(); i++) {
		double distance = sqrt(abs(vtkMath::Distance2BetweenPoints(pointArray, data.getCentrePoints()->GetPoint(neighbours->GetId(i)))));

		if (!data.getDoubleLayer()) {
			int material = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(neighbours->GetId(i), 0);
			if ((distance < mindistance) && (material == Material::subendo_right_Atrial_Cardium)) {
				mindistance = distance;
				newCellId = neighbours->GetId(i);
			}
		}

		int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(neighbours->GetId(i), 0);
		if ((distance < mindistance) && (material == Material::subendo_right_Atrial_Cardium)) {
			mindistance = distance;
			newCellId = neighbours->GetId(i);
		}
	}

	if (newCellId != -1) {
		return newCellId;
	}
	else {
		cout << "Point in right endo material was not found!!!" << endl;
		return newCellId;
	}
} // Methods::findClosedPointIdinMaterialInRightEndo

/*!Find the closest cell in a region in the right atrial.
 \param data Pointer to the orginal mesh.
 \param point Point as vector (x,y,z)
 \param region Region in that the cell must lie.
 \return The cell who would be found as cell id.
 */
vtkIdType Methods::findClosedPointIdinMaterialInRightInRegion(DataFormat &data, vector<double> point, vtkSmartPointer<vtkIdList> region) {
	double pointArray[3] = { point.at(0), point.at(1), point.at(2) };
	double mindistance = std::numeric_limits<double>::max();

	vtkIdType newCellId = -1;

	for (vtkIdType i = 0; i < region->GetNumberOfIds(); i++) {
		double distance = sqrt(abs(vtkMath::Distance2BetweenPoints(pointArray, data.getCentrePoints()->GetPoint(region->GetId(i)))));
		if (!data.getDoubleLayer()) {
			int material = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(region->GetId(i), 0);

			if ((distance < mindistance) && Material::isInRightWithoutBridges(material)) {
				mindistance = distance;
				newCellId = region->GetId(i);
			}
		}

		int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(region->GetId(i), 0);

		if ((distance < mindistance) && Material::isInRightWithoutBridges(material)) {
			mindistance = distance;
			newCellId = region->GetId(i);
		}
	}

	if (newCellId != -1) {
		return newCellId;
	}
	else {
		cout << "Point in right material was not found!!!" << endl;
		return newCellId;
	}
} // Methods::findClosedPointIdinMaterialInRightInRegion

/*!Find the closest cell in a region in the left atrial.
 \param data Pointer to the orginal mesh.
 \param point Point as vector (x,y,z)
 \param region Region in that the cell must lie.
 \return The cell who would be found as cell id.
 */
vtkIdType Methods::findClosedPointIdinMaterialInLeftInRegion(DataFormat &data, vector<double> point, vtkSmartPointer<vtkIdList> region) {
	double pointArray[3] = { point.at(0), point.at(1), point.at(2) };
	double mindistance = std::numeric_limits<double>::max();

	vtkIdType newCellId = -1;

	for (vtkIdType i = 0; i < region->GetNumberOfIds(); i++) {
		double distance = sqrt(abs(vtkMath::Distance2BetweenPoints(pointArray, data.getCentrePoints()->GetPoint(region->GetId(i)))));
		if (!data.getDoubleLayer()) {
			int material = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(region->GetId(i), 0);
			if ((distance < mindistance) && Material::isInLeftWithoutBridge(material)) {
				mindistance = distance;
				newCellId = region->GetId(i);
			}
		}

		int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(region->GetId(i), 0);
		if ((distance < mindistance) && Material::isInLeftWithoutBridge(material)) {
			mindistance = distance;
			newCellId = region->GetId(i);
		}
	}

	if (newCellId != -1) {
		return newCellId;
	}
	else {
		cout << "Point in left material was not found!!!" << endl;
		return newCellId;
	}
} // Methods::findClosedPointIdinMaterialInLeftInRegion

/*!Find the closest point in a material in a radius of 15 mm.
 \param data Pointer to the orginal mesh.
 \param point Point as vector (x,y,z)
 \param inMaterial Tissue class in that the point must lie.
 \return The point who would be found as vector (x,y,z).
 */
vector<double> Methods::findClosedPointinMaterial(DataFormat &data, vector<double> point, Material::Mat inMaterial) {
	return findClosedPointinMaterialinRadius(data, point, inMaterial, inMaterial, 15);
}

/*!Find the closest point in a material in a radius of 15 mm.
 \param data Pointer to the orginal mesh.
 \param point Point as vector (x,y,z)
 \param inMaterial1 First tissue class in that the point must lie.
 \param inMaterial2 Second tissue class in that the point must lie.
 \return The point who would be found as vector (x,y,z).
 */
vector<double> Methods::findClosedPointinMaterial(DataFormat &data, vector<double> point, Material::Mat inMaterial1, Material::Mat inMaterial2) {
	return findClosedPointinMaterialinRadius(data, point, inMaterial1, inMaterial2, 15);
}

/*!Find the closest point in a material.
 \param data Pointer to the orginal mesh.
 \param point Point as vector (x,y,z)
 \param inMaterial Tissue class in that the point must lie.
 \param radius Radius in that the point would searched.
 \return The point who would be found as vector (x,y,z).
 */
vector<double> Methods::findClosedPointinMaterialinRadius(DataFormat &data, vector<double> point, Material::Mat inMaterial, double radius) {
	return findClosedPointinMaterialinRadius(data, point, inMaterial, inMaterial, radius);
}

/*!Find the closest point in a material.
 \param data Pointer to the orginal mesh.
 \param point Point as vector (x,y,z)
 \param inMaterial1 First tissue class in that the point must lie.
 \param inMaterial2 Second tissue class in that the point must lie.
 \param radius Radius in that the point would searched.
 \return The point who would be found as vector (x,y,z).
 */
vector<double> Methods::findClosedPointinMaterialinRadius(DataFormat &data, vector<double> point, Material::Mat inMaterial1, Material::Mat inMaterial2, double radius) {
	double pointArray[3] = { point.at(0), point.at(1), point.at(2) };
	double mindistance = std::numeric_limits<double>::max();

	vtkIdType cellId = data.getCentrePoints()->FindPoint(pointArray);

	if (!data.getDoubleLayer()) {
		int material = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(cellId, 0);
		if ((material == inMaterial1) || (material == inMaterial2)) {
			vector<double> closedPoint = { data.getCentrePoints()->GetPoint(cellId)[0], data.getCentrePoints()->GetPoint(cellId)[1], data.getCentrePoints()->GetPoint(cellId)[2] };
			return closedPoint;
		}
	}

	int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(cellId, 0);
	if ((material == inMaterial1) || (material == inMaterial2)) {
		vector<double> closedPoint = { data.getCentrePoints()->GetPoint(cellId)[0], data.getCentrePoints()->GetPoint(cellId)[1], data.getCentrePoints()->GetPoint(cellId)[2] };
		return closedPoint;
	}

	vtkSmartPointer<vtkIdList> neighbours = getCellsinRadius(data, cellId, radius);

	vtkIdType newCellId = -1;

	vector<double> newPoint(3);
	newPoint.at(0) = NAN;
	newPoint.at(1) = NAN;
	newPoint.at(2) = NAN;


	for (vtkIdType i = 0; i < neighbours->GetNumberOfIds(); i++) {
		double distance = sqrt(abs(vtkMath::Distance2BetweenPoints(pointArray, data.getCentrePoints()->GetPoint(neighbours->GetId(i)))));

		if (!data.getDoubleLayer()) {
			int material = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(neighbours->GetId(i), 0);
			if ((distance < mindistance) && ((material == inMaterial1) || (material == inMaterial2))) {
				mindistance = distance;
				newCellId = neighbours->GetId(i);
			}
		}

		int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(neighbours->GetId(i), 0);
		if ((distance < mindistance) && ((material == inMaterial1) || (material == inMaterial2))) {
			mindistance = distance;
			newCellId = neighbours->GetId(i);
		}
	}

	if (newCellId != -1) {
		double *newPointArray = data.getCentrePoints()->GetPoint(newCellId);
		newPoint.at(0) = newPointArray[0];
		newPoint.at(1) = newPointArray[1];
		newPoint.at(2) = newPointArray[2];
		return newPoint;
	}
	else {
		cout << "The point: " << point.at(0) << ", " << point.at(1) << ", " << point.at(2) << " can't be moved!!!" << endl;
		cout << "Point in Material " << inMaterial1 << " or " << inMaterial2 << " was not found!!!" << endl;
		return newPoint;
	}
} // Methods::findClosedPointinMaterialinRadius

/*!Find the closest point in direction of another point in a material.
 \param data Pointer to the orginal mesh.
 \param point Point as vector (x,y,z)
 \param toPoint Reference point
 \param minDistance The minimum distance to the old point.
 \param inMaterial Tissue class in that the point must lie.
 \return The point who would be found as vector (x,y,z).
 */
vector<double> Methods::findClosedPointinMaterialinDirectionPoint(DataFormat &data, vector<double> point, vector<double> toPoint, double minDistance, Material::Mat inMaterial) {
	return findClosedPointinMaterialinDirectionPoint(data, point, toPoint, minDistance, inMaterial, inMaterial);
}

/*!Find the closest point in direction of another point in a material.
 \param data Pointer to the orginal mesh.
 \param point Point as vector (x,y,z)
 \param toPoint Reference point
 \param minDistance The minimum distance to the old point.
 \param inMaterial1 First tissue class in that the point must lie.
 \param inMaterial2 Second tissue class in that the point must lie.
 \return The point who would be found as vector (x,y,z).
 */
vector<double> Methods::findClosedPointinMaterialinDirectionPoint(DataFormat &data, vector<double> point, vector<double> toPoint, double minDistance, Material::Mat inMaterial1, Material::Mat inMaterial2) {
	double radius = 15;
	double pointArray[3] = { point.at(0), point.at(1), point.at(2) };
	double toPointArray[3] = { toPoint.at(0), toPoint.at(1), toPoint.at(2) };

	double mininmalDistance = std::numeric_limits<double>::max();

	double toPointDistance = sqrt(abs(vtkMath::Distance2BetweenPoints(toPointArray, pointArray)));

	vtkIdType cellId = data.getCentrePoints()->FindPoint(pointArray);

	vtkSmartPointer<vtkIdList> neighbours = getCellsinRadius(data, cellId, radius);

	vtkIdType newCellId = -1;

	vector<double> newPoint(3);

	newPoint.at(0) = NAN;
	newPoint.at(1) = NAN;
	newPoint.at(2) = NAN;

	for (vtkIdType i = 0; i < neighbours->GetNumberOfIds(); i++) {
		double distance = sqrt(abs(vtkMath::Distance2BetweenPoints(pointArray, data.getCentrePoints()->GetPoint(neighbours->GetId(i)))));
		double distancetoPoint = sqrt(abs(vtkMath::Distance2BetweenPoints(toPointArray, data.getCentrePoints()->GetPoint(neighbours->GetId(i)))));
		if (!data.getDoubleLayer()) {
			int material = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(neighbours->GetId(i), 0);

			if ((distance < mininmalDistance) && (distance > minDistance) && (distancetoPoint < toPointDistance) &&
				((material == inMaterial1) || (material == inMaterial2)) &&
				(data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->GetComponent(neighbours->GetId(i), 0) +
					data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->GetComponent(neighbours->GetId(i), 1) +
					data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->GetComponent(neighbours->GetId(i), 2) == 0)) {
				mininmalDistance = distance;
				newCellId = neighbours->GetId(i);
			}
		}

		int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(neighbours->GetId(i), 0);
		if ((distancetoPoint < toPointDistance) && (distance < mininmalDistance) && (distance > minDistance) &&
			((material == inMaterial1) || (material == inMaterial2)) &&
			(data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetComponent(neighbours->GetId(i), 0) +
				data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetComponent(neighbours->GetId(i), 1) +
				data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetComponent(neighbours->GetId(i), 2) == 0)) {
			mininmalDistance = distance;
			newCellId = neighbours->GetId(i);
		}
	}

	if (newCellId != -1) {
		double *newPointArray = data.getCentrePoints()->GetPoint(newCellId);
		newPoint.at(0) = newPointArray[0];
		newPoint.at(1) = newPointArray[1];
		newPoint.at(2) = newPointArray[2];
		cout << "The point: " << point.at(0) << ", " << point.at(1) << ", " << point.at(2) << " was moved in: " <<
			newPoint.at(0) << ", " << newPoint.at(1) << ", " << newPoint.at(2) << endl;
		return newPoint;
	}
	else {
		cout << "The point: " << point.at(0) << ", " << point.at(1) << ", " << point.at(2) << " can't be moved!!!" << endl;
		cout << "Point in Material " << inMaterial1 << "or " << inMaterial2 << " was not found!!!" << endl;
		return newPoint;
	}
} // Methods::findClosedPointinMaterialinDirectionPoint

/*!Find the closest cell to another point in a define radius.
 \param data Pointer to the orginal mesh.
 \param point Point as vector (x,y,z)
 \param toPoint Reference point
 \param radius Define the searchregion arround the point.
 \return The cell id who would be found.
 */
vtkIdType Methods::findClosedCellIdinDirectionPoint(DataFormat &data, vector<double> point, vector<double> toPoint, double radius) {
	double pointArray[3] = { point.at(0), point.at(1), point.at(2) };
	double toPointArray[3] = { toPoint.at(0), toPoint.at(1), toPoint.at(2) };

	double mininmalDistance = std::numeric_limits<double>::max();

	double toPointDistance = sqrt(abs(vtkMath::Distance2BetweenPoints(toPointArray, pointArray)));

	vtkIdType cellId = data.getCentrePoints()->FindPoint(pointArray);

	vtkSmartPointer<vtkIdList> neighbours = getCellsinRadius(data, cellId, radius);

	vtkIdType newCellId = -1;

	for (vtkIdType i = 0; i < neighbours->GetNumberOfIds(); i++) {
		double distance = sqrt(abs(vtkMath::Distance2BetweenPoints(pointArray, data.getCentrePoints()->GetPoint(neighbours->GetId(i)))));
		double distancetoPoint = sqrt(abs(vtkMath::Distance2BetweenPoints(toPointArray, data.getCentrePoints()->GetPoint(neighbours->GetId(i)))));

		if ((distancetoPoint < toPointDistance) && (distance < mininmalDistance)) {
			mininmalDistance = distance;
			newCellId = neighbours->GetId(i);
		}
	}

	if (newCellId != -1) {
		return newCellId;
	}
	else {
		cout << "The point: " << point.at(0) << ", " << point.at(1) << ", " << point.at(2) << " can't be moved!!!" << endl;
		cout << "Point was not found!!!" << endl;
		return newCellId;
	}
} // Methods::findClosedCellIdinDirectionPoint

/*!Find the furtherst point to another point in a define radius and material.
 \param data Pointer to the orginal mesh.
 \param point Point as vector (x,y,z)
 \param toPoint Reference point
 \param inMaterial Tissue class in that the point must lie.
 \param radius Define the searchregion arround the point.
 \return The point who would be found as vector (x,y,z).
 */
vector<double> Methods::findFurthestPointinMaterialtoPoint(DataFormat &data, vector<double> point, vector<double> toPoint, Material::Mat inMaterial1, double radius) {
	double pointArray[3] = { point.at(0), point.at(1), point.at(2) };
	double toPointArray[3] = { toPoint.at(0), toPoint.at(1), toPoint.at(2) };

	double maxdistance = -std::numeric_limits<double>::max();

	vtkIdType cellId = data.getCentrePoints()->FindPoint(pointArray);

	vtkSmartPointer<vtkIdList> neighbours = getCellsinRadius(data, cellId, radius);

	vtkIdType newCellId = -1;

	vector<double> newPoint(3);

	newPoint.at(0) = NAN;
	newPoint.at(1) = NAN;
	newPoint.at(2) = NAN;

	for (vtkIdType i = 0; i < neighbours->GetNumberOfIds(); i++) {
		double distance = sqrt(abs(vtkMath::Distance2BetweenPoints(toPointArray, data.getCentrePoints()->GetPoint(neighbours->GetId(i)))));

		if (!data.getDoubleLayer()) {
			int material = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(neighbours->GetId(i), 0);
			if ((distance > maxdistance) && (material == inMaterial1) &&
				(data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->GetComponent(neighbours->GetId(i), 0) +
					data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->GetComponent(neighbours->GetId(i), 1) +
					data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->GetComponent(neighbours->GetId(i), 2) == 0)) {
				maxdistance = distance;
				newCellId = neighbours->GetId(i);
			}
		}

		int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(neighbours->GetId(i), 0);
		if ((distance > maxdistance) && (material == inMaterial1) &&
			(data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetComponent(neighbours->GetId(i), 0) +
				data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetComponent(neighbours->GetId(i), 1) +
				data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetComponent(neighbours->GetId(i), 2) == 0)) {
			maxdistance = distance;
			newCellId = neighbours->GetId(i);
		}
	}

	if (newCellId != -1) {
		double *newPointArray = data.getCentrePoints()->GetPoint(newCellId);
		newPoint.at(0) = newPointArray[0];
		newPoint.at(1) = newPointArray[1];
		newPoint.at(2) = newPointArray[2];
		cout << "The point: " << point.at(0) << ", " << point.at(1) << ", " << point.at(2) << " was moved in: " <<
			newPoint.at(0) << ", " << newPoint.at(1) << ", " << newPoint.at(2) << endl;
		return newPoint;
	}
	else {
		cout << "The point: " << point.at(0) << ", " << point.at(1) << ", " << point.at(2) << " can't be moved!!!" << endl;
		cout << "Point in Material " << inMaterial1 << " was not found!!!" << endl;
		return newPoint;
	}
} // Methods::findFurthestPointinMaterialtoPoint

/*!Find the closest point to another point in a define radius and in the left atrial.
 \param data Pointer to the orginal mesh.
 \param startPoint Start Point for the line as vector (x,y,z)
 \param targetPoint Target Point for the line as vector (x,y,z)
 \param radius Define the searchregion arround the line.
 \param resolution Resolution of the mesh.
 \return The point who would be found as vector (x,y,z).
 */
vector<double> Methods::findClosedPointinMaterialinRadiustoPointLeft(DataFormat &data, vector<double> startPoint, vector<double> targetPoint, double radius, double resolution) {

	double tempStartPointArray[3] = { startPoint.at(0), startPoint.at(1), startPoint.at(2) };
	vtkIdType startPointId = data.getCentrePoints()->FindPoint(tempStartPointArray);
	double startPointArray[3] = { data.getCentrePoints()->GetPoint(startPointId)[0], data.getCentrePoints()->GetPoint(startPointId)[1], data.getCentrePoints()->GetPoint(startPointId)[2] };

	double tempTargetPointArray[3] = { targetPoint.at(0), targetPoint.at(1), targetPoint.at(2) };
	vtkIdType targetPointId = data.getCentrePoints()->FindPoint(tempTargetPointArray);
	double targetPointArray[3] = { data.getCentrePoints()->GetPoint(targetPointId)[0], data.getCentrePoints()->GetPoint(targetPointId)[1], data.getCentrePoints()->GetPoint(targetPointId)[2] };

	vector<double> point = { targetPointArray[0], targetPointArray[1], targetPointArray[2] };

	int resolution2 = ceil(2 * radius*M_PI / 3 / resolution);

	vtkSmartPointer<vtkLineSource> lineSource = vtkSmartPointer<vtkLineSource>::New();

	lineSource->SetPoint1(startPointArray[0], startPointArray[1], startPointArray[2]);
	lineSource->SetPoint2(targetPointArray[0], targetPointArray[1], targetPointArray[2]);
	lineSource->Update();

	vtkSmartPointer<vtkTubeFilter> cylindertube = vtkSmartPointer<vtkTubeFilter>::New();
	cylindertube->SetInputConnection(lineSource->GetOutputPort());
	cylindertube->SetRadius(radius + resolution);
	cylindertube->CappingOn();
	cylindertube->SetNumberOfSides(resolution2);
	cylindertube->Update();

	vtkSmartPointer<vtkTriangleFilter> triangleFilter = vtkSmartPointer<vtkTriangleFilter>::New();
	triangleFilter->SetInputConnection(cylindertube->GetOutputPort());
	triangleFilter->Update();

	vtkSmartPointer<vtkPolyData> cylinder = triangleFilter->GetOutput();

	vtkSmartPointer<vtkCleanPolyData> cleanFilterCylinderBridge = vtkSmartPointer<vtkCleanPolyData>::New();
	cleanFilterCylinderBridge->SetInputData(cylinder);
	cleanFilterCylinderBridge->Update();
	cylinder = cleanFilterCylinderBridge->GetOutput();

	DataFormat cylinderData;
	cylinderData.setInputType(DataFormat::vtp);
	cylinderData.setVtkData(cylinder);
	cylinderData.setCentrePoints(Methods::calcCellCentroids(cylinder));

	vtkSmartPointer<vtkIdList> region = getCellsWhichWereInside(data, cylinderData);

	vector<double> closestPoint = findClosedPointinMaterialInLeftInRegion(data, point, region);

	return closestPoint;
} // Methods::findClosedPointinMaterialinRadiustoPointLeft

/*!Find the closest point to another point in a define radius and in the right atrial.
 \param data Pointer to the orginal mesh.
 \param startPoint Start Point for the line as vector (x,y,z)
 \param targetPoint Target Point for the line as vector (x,y,z)
 \param radius Define the searchregion arround the line.
 \param resolution Resolution of the mesh.
 \return The point who would be found as vector (x,y,z).
 */
vector<double> Methods::findClosedPointinMaterialinRadiustoPointRight(DataFormat &data, vector<double> startPoint, vector<double> targetPoint, double radius, double resolution) {

	double tempStartPointArray[3] = { startPoint.at(0), startPoint.at(1), startPoint.at(2) };
	vtkIdType startPointId = data.getCentrePoints()->FindPoint(tempStartPointArray);
	double startPointArray[3] = { data.getCentrePoints()->GetPoint(startPointId)[0], data.getCentrePoints()->GetPoint(startPointId)[1], data.getCentrePoints()->GetPoint(startPointId)[2] };

	double tempTargetPointArray[3] = { targetPoint.at(0), targetPoint.at(1), targetPoint.at(2) };
	vtkIdType targetPointId = data.getCentrePoints()->FindPoint(tempTargetPointArray);
	double targetPointArray[3] = { data.getCentrePoints()->GetPoint(targetPointId)[0], data.getCentrePoints()->GetPoint(targetPointId)[1], data.getCentrePoints()->GetPoint(targetPointId)[2] };

	vector<double> point = { targetPointArray[0], targetPointArray[1], targetPointArray[2] };

	int resolution2 = ceil(2 * radius*M_PI / 3 / resolution);

	vtkSmartPointer<vtkLineSource> lineSource = vtkSmartPointer<vtkLineSource>::New();

	lineSource->SetPoint1(startPointArray[0], startPointArray[1], startPointArray[2]);
	lineSource->SetPoint2(targetPointArray[0], targetPointArray[1], targetPointArray[2]);
	lineSource->Update();

	vtkSmartPointer<vtkTubeFilter> cylindertube = vtkSmartPointer<vtkTubeFilter>::New();
	cylindertube->SetInputConnection(lineSource->GetOutputPort());
	cylindertube->SetRadius(radius + resolution);
	cylindertube->CappingOn();
	cylindertube->SetNumberOfSides(resolution2);
	cylindertube->Update();

	vtkSmartPointer<vtkTriangleFilter> triangleFilter = vtkSmartPointer<vtkTriangleFilter>::New();
	triangleFilter->SetInputConnection(cylindertube->GetOutputPort());
	triangleFilter->Update();

	vtkSmartPointer<vtkPolyData> cylinder = triangleFilter->GetOutput();

	vtkSmartPointer<vtkCleanPolyData> cleanFilterCylinderBridge = vtkSmartPointer<vtkCleanPolyData>::New();
	cleanFilterCylinderBridge->SetInputData(cylinder);
	cleanFilterCylinderBridge->Update();
	cylinder = cleanFilterCylinderBridge->GetOutput();

	DataFormat cylinderData;
	cylinderData.setInputType(DataFormat::vtp);
	cylinderData.setVtkData(cylinder);
	cylinderData.setCentrePoints(Methods::calcCellCentroids(cylinder));

	vtkSmartPointer<vtkIdList> region = getCellsWhichWereInside(data, cylinderData);

	vector<double> closestPoint = findClosedPointinMaterialInRightInRegion(data, point, region);

	return closestPoint;
} // Methods::findClosedPointinMaterialinRadiustoPointRight

/*!Find the closest point on a path.
 \param data Pointer to the orginal mesh.
 \param point Reference point as vector (x,y,z)
 \param path Path on that the point must lie.
 \return The point who would be found as vector (x,y,z).
 */
vector<double> Methods::findClosedPointOnPath(DataFormat &data, vector<double> point, vtkSmartPointer<vtkIdList> path) {
	vtkSmartPointer<vtkUnstructuredGrid> pathSearchList = vtkSmartPointer<vtkUnstructuredGrid>::New();
	vtkSmartPointer<vtkPoints> pointsPath = vtkSmartPointer<vtkPoints>::New();

	for (vtkIdType i = 0; i < path->GetNumberOfIds(); i++) {
		double point[3] = { data.getCentrePoints()->GetPoint(path->GetId(i))[0], data.getCentrePoints()->GetPoint(path->GetId(i))[1], data.getCentrePoints()->GetPoint(path->GetId(i))[2] };
		pointsPath->InsertNextPoint(point);
	}
	pathSearchList->SetPoints(pointsPath);

	double pointArray[3] = { point.at(0), point.at(1), point.at(2) };
	vtkIdType pointId = pathSearchList->FindPoint(pointArray);

	double *newPointArray = data.getCentrePoints()->GetPoint(path->GetId(pointId));
	vector<double> newPoint = { newPointArray[0], newPointArray[1], newPointArray[2] };

	return newPoint;
}

/*!Find the closest cell on a path.
 \param data Pointer to the orginal mesh.
 \param point Reference point as vector (x,y,z)
 \param path Path on that the point must lie.
 \return The cell id who would be found.
 */
vtkIdType Methods::findClosedPointIdOnPath(DataFormat &data, vector<double> point, vtkSmartPointer<vtkIdList> path) {
	vtkSmartPointer<vtkUnstructuredGrid> pathSearchList = vtkSmartPointer<vtkUnstructuredGrid>::New();
	vtkSmartPointer<vtkPoints> pointsPath = vtkSmartPointer<vtkPoints>::New();

	for (vtkIdType i = 0; i < path->GetNumberOfIds(); i++) {
		double point[3] = { data.getCentrePoints()->GetPoint(path->GetId(i))[0], data.getCentrePoints()->GetPoint(path->GetId(i))[1], data.getCentrePoints()->GetPoint(path->GetId(i))[2] };
		pointsPath->InsertNextPoint(point);
	}
	pathSearchList->SetPoints(pointsPath);

	double pointArray[3] = { point.at(0), point.at(1), point.at(2) };
	vtkIdType pointId = pathSearchList->FindPoint(pointArray);

	return pointId;
}

/*!Find the closest cell on a path without an orientation in a material.
 \param data Pointer to the orginal mesh.
 \param point Reference point as vector (x,y,z)
 \param path Path on that the point must lie.
 \param inMaterial1 First tissue class in that the cell must lie.
 \param inMaterial2 Second tissue class in that the cell must lie.
 \return The cell id who would be found.
 */
vtkIdType Methods::findClosedPointIdOnPathWithoutOrientationInMaterial(DataFormat &data, vector<double> point, vtkSmartPointer<vtkIdList> path, Material::Mat inMaterial1, Material::Mat inMaterial2) {
	double pointArray[3] = { point.at(0), point.at(1), point.at(2) };
	double mindistance = std::numeric_limits<double>::max();
	vtkIdType newCellId = -1;

	for (vtkIdType i = 0; i < path->GetNumberOfIds(); i++) {
		double distance = vtkMath::Distance2BetweenPoints(pointArray, data.getCentrePoints()->GetPoint(path->GetId(i)));
		int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(path->GetId(i), 0);


		if ((distance < mindistance) && ((material == inMaterial1) || (material == inMaterial2)) &&
			(data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetComponent(path->GetId(i), 0) +
				data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetComponent(path->GetId(i), 1) +
				data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetComponent(path->GetId(i), 2) == 0)) {
			mindistance = distance;
			newCellId = path->GetId(i);
		}
	}

	return newCellId;
}

/*!Calc an out vein point.
 \param point1 Reference point as vector (x,y,z)
 \param point2 Reference point as vector (x,y,z)
 \param point3 Reference point as vector (x,y,z)
 \param centroid Centroid of the 3 points as vector (x,y,z)
 \param distanceToCentrepoint distance to the centroid.
 \return The point who would be found as vector (x,y,z).
 */
vector<double> Methods::getOuterVeinCentrepoint(vector<double> point1, vector<double> point2, vector<double> point3, vector<double> centroid, double distanceToCentrepoint) {
	double centroidArray[3] = { centroid.at(0), centroid.at(1), centroid.at(2) };
	double point1Array[3] = { point1.at(0), point1.at(1), point1.at(2) };
	double point2Array[3] = { point2.at(0), point2.at(1), point2.at(2) };
	double point3Array[3] = { point3.at(0), point3.at(1), point3.at(2) };
	double normalVectorPlane[3];
	double vectorPoint2Point1[3];
	double vectorPoint2Point3[3];

	vtkMath::Subtract(point1Array, point2Array, vectorPoint2Point1);
	vtkMath::Subtract(point3Array, point1Array, vectorPoint2Point3);
	vtkMath::Cross(vectorPoint2Point3, vectorPoint2Point1, normalVectorPlane);
	vtkMath::Normalize(normalVectorPlane);

	vtkMath::MultiplyScalar(normalVectorPlane, distanceToCentrepoint);
	double returnPointArray[3] = { 0, 0, 0 };
	vtkMath::Add(centroidArray, normalVectorPlane, returnPointArray);
	vector<double> returnPoint = { returnPointArray[0], returnPointArray[1], returnPointArray[2] };

	return returnPoint;
}

/*!Test if it is an outer vein.
 \param data Pointer to the orginal mesh.
 \param point1 Reference point as vector (x,y,z)
 \param point2 Reference point as vector (x,y,z)
 \param point3 Reference point as vector (x,y,z)
 \param centroid Centroid of the 3 points as vector (x,y,z)
 \param refPoint Refence Point to define outer site.
 \param radius Radius in that a cell would be searched.
 \param inMaterial1 First tissue class in that the cell must lie.
 \param inMaterial2 Second tissue class in that the cell must lie.
 \return True if there is an outer vein or not then false.
 */

bool Methods::isVein(DataFormat &data, vector<double> point1, vector<double> point2, vector<double> point3, vector<double> centroid, vector<double> refPoint, double radius, Material::Mat inMaterial1, Material::Mat inMaterial2) {
	double centroidArray[3] = { centroid.at(0), centroid.at(1), centroid.at(2) };
	double refPointArray[3] = { refPoint.at(0), refPoint.at(1), refPoint.at(2) };
	double point1Array[3] = { point1.at(0), point1.at(1), point1.at(2) };
	double point2Array[3] = { point2.at(0), point2.at(1), point2.at(2) };
	double point3Array[3] = { point3.at(0), point3.at(1), point3.at(2) };
	double normalVectorPlane[3];
	double vectorPoint2Point1[3];
	double vectorPoint3Point1[3];

	vtkMath::Subtract(point2Array, point1Array, vectorPoint2Point1);
	vtkMath::Subtract(point3Array, point1Array, vectorPoint3Point1);
	vtkMath::Cross(vectorPoint2Point1, vectorPoint3Point1, normalVectorPlane);
	vtkMath::Normalize(normalVectorPlane);

	double evalPlaneRefPoint = vtkPlane::Evaluate(normalVectorPlane, centroidArray, refPointArray);
	if (!signbit(evalPlaneRefPoint)) {
		vtkMath::MultiplyScalar(normalVectorPlane, -1);
	}

	double copyNormalvectorPlane[3] = { normalVectorPlane[0], normalVectorPlane[1], normalVectorPlane[2] };

	vtkMath::MultiplyScalar(copyNormalvectorPlane, radius / 2);
	double centrePoint[3];
	vtkMath::Add(centroidArray, copyNormalvectorPlane, centrePoint);

	vtkSmartPointer<vtkPointLocator> centrePointLocator = vtkSmartPointer<vtkPointLocator>::New();
	centrePointLocator->SetDataSet(data.getCentrePoints());
	centrePointLocator->BuildLocator();

	vtkSmartPointer<vtkIdList> cells = vtkSmartPointer<vtkIdList>::New();

	vtkSmartPointer<vtkIdList> possibleCells = vtkSmartPointer<vtkIdList>::New();
	centrePointLocator->FindPointsWithinRadius(radius, centrePoint, possibleCells);

	double originPlane1[3];
	vtkMath::Add(centrePoint, normalVectorPlane, originPlane1);

	double originPlane2[3];
	double negativeNormalvectorPlane[3] = { normalVectorPlane[0], normalVectorPlane[1], normalVectorPlane[2] };
	vtkMath::MultiplyScalar(negativeNormalvectorPlane, -1);
	vtkMath::Add(centrePoint, negativeNormalvectorPlane, originPlane2);

	for (vtkIdType i = 0; i < possibleCells->GetNumberOfIds(); i++) {
		double *currentPoint = data.getCentrePoints()->GetPoint(possibleCells->GetId(i));

		double evaluatePlane1 = vtkPlane::Evaluate(normalVectorPlane, originPlane1, currentPoint);
		double evaluatePlane2 = vtkPlane::Evaluate(normalVectorPlane, originPlane2, currentPoint);

		if ((evaluatePlane1 <= 0) && (evaluatePlane2 >= 0)) {
			cells->InsertNextId(possibleCells->GetId(i));
		}
	}

	for (vtkIdType j = 0; j < cells->GetNumberOfIds(); j++) {
		if (!data.getDoubleLayer()) {
			int material = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(cells->GetId(j), 0);
			if ((material == inMaterial1) || (material == inMaterial2)) {
				return true;
			}
		}

		int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(cells->GetId(j), 0);
		if ((material == inMaterial1) || (material == inMaterial2)) {
			return true;
		}
	}

	return false;
} // Methods::isVein

/*! Translate the points.
 \param points Point that must translate.
 \param vector Translation vector (x,y,z)
 \return Translated points.
 */
vtkSmartPointer<vtkPoints> Methods::movePoints(vtkSmartPointer<vtkPoints> points, vector<double> vector) {
	double vectorArray[3] = { vector.at(0), vector.at(1), vector.at(2) };

	for (vtkIdType i = 0; i < points->GetNumberOfPoints(); i++) {
		if (!::isnan(points->GetPoint(i)[0])) {
			double *currentPoint = points->GetPoint(i);
			vtkMath::Add(currentPoint, vectorArray, currentPoint);
			points->SetPoint(i, currentPoint);
		}
	}
	return points;
}

/*! Rotate the points.
 \param points Point that must translate.
 \param matrix 3x3 rotation matrix
 \return Rotated points.
 */
vtkSmartPointer<vtkPoints> Methods::rotatePoints(vtkSmartPointer<vtkPoints> points, vtkSmartPointer<vtkMatrix3x3> matrix) {
	for (vtkIdType i = 0; i < points->GetNumberOfPoints(); i++) {
		if (!::isnan(points->GetPoint(i)[0])) {
			double *currentPoint = points->GetPoint(i);
			matrix->MultiplyPoint(currentPoint, currentPoint);
			points->SetPoint(i, currentPoint);
		}
	}
	return points;
}

/*! Find a corresponding point in to lists with a tolerance of 0.001.
 \param data Pointer to the orginal mesh.
 \param pointList1 First list of point ids.
 \param pointList2 Second list of point ids.
 \return A list with tupel of corresponding points.
 */
unordered_map<vtkIdType, vtkIdType> Methods::getCommonPoint(DataFormat &data, vtkSmartPointer<vtkIdList> pointList1, vtkSmartPointer<vtkIdList> pointList2) {
	return getCommonPoint(data, pointList1, pointList2, 0.001);
}

/*! Find a corresponding point in to lists with a tolerance of 0.001.
 \param data Pointer to the orginal mesh.
 \param pointList1 First list of point ids.
 \param pointList2 Second list of point ids.
 \param tolerance Distance that would be accept between to points.
 \return A list with tupel of corresponding points.
 */
unordered_map<vtkIdType, vtkIdType> Methods::getCommonPoint(DataFormat &data, vtkSmartPointer<vtkIdList> pointList1, vtkSmartPointer<vtkIdList> pointList2, double tolerance) {
	unordered_map<vtkIdType, vtkIdType> PointTupels;

	for (vtkIdType i = 0; i < pointList1->GetNumberOfIds(); i++) {
		double point[3] = { data.getVtkData()->GetPoint(pointList1->GetId(i))[0], data.getVtkData()->GetPoint(pointList1->GetId(i))[1], data.getVtkData()->GetPoint(pointList1->GetId(i))[2] };
		double mindistance = std::numeric_limits<double>::max();
		vtkIdType k = -1;
		for (vtkIdType j = 0; j < pointList2->GetNumberOfIds(); j++) {
			double distance = sqrt(abs(vtkMath::Distance2BetweenPoints(point, data.getVtkData()->GetPoint(pointList2->GetId(j)))));
			if ((distance < mindistance) && (pointList1->GetId(i) != pointList2->GetId(j))) {
				mindistance = distance;
				k = j;
			}
		}
		if ((mindistance < tolerance) && (k != -1)) {
			PointTupels[pointList1->GetId(i)] = pointList2->GetId(k);
		}
	}

	return PointTupels;
}

/*! Find a corresponding cells in to list.
 \param data Pointer to the orginal mesh.
 \param pointList1 First list of point ids.
 \param pointList2 Second list of point ids.
 \return A list with tupel of corresponding points.
 */
unordered_map<vtkIdType, vtkIdType> Methods::getCommonPoints(DataFormat &data, vtkSmartPointer<vtkIdList> cellIdList1, vtkSmartPointer<vtkIdList> cellIdList2) {
	unordered_map<vtkIdType, vtkIdType> PointTupels;

	for (vtkIdType i = 0; i < cellIdList1->GetNumberOfIds(); i++) {
		vtkSmartPointer<vtkIdList> cellPoints1 = vtkSmartPointer<vtkIdList>::New();
		data.getVtkData()->GetCellPoints(cellIdList1->GetId(i), cellPoints1);

		for (vtkIdType j = 0; j < cellIdList2->GetNumberOfIds(); j++) {
			vtkSmartPointer<vtkIdList> cellPoints2 = vtkSmartPointer<vtkIdList>::New();
			data.getVtkData()->GetCellPoints(cellIdList2->GetId(j), cellPoints2);

			unordered_map<vtkIdType, vtkIdType> tempPointTupels = getCommonPoint(data, cellPoints1, cellPoints2);

			for (auto &tupel : tempPointTupels) {
				if (PointTupels.find(tupel.first) == PointTupels.end()) {
					PointTupels[tupel.first] = tupel.second;
				}
			}
		}
	}

	return PointTupels;
}

/*! Replace point ids of the cells with that point id of the common list.
 \param data Pointer to the orginal mesh.
 \param commonPoints List with replace ids.
 */
void Methods::replacePoints(DataFormat &data, unordered_map<vtkIdType, vtkIdType> commonPoints) {
	vtkSmartPointer<vtkUnstructuredGrid> dataUGrid;
	vtkSmartPointer<vtkPolyData> dataPoly;

	if (data.getInputType() == DataFormat::PossibleInputType::vtu) {
		dataUGrid = vtkUnstructuredGrid::SafeDownCast(data.getVtkData());
	}
	else if (data.getInputType() == DataFormat::PossibleInputType::vtp) {
		dataPoly = vtkPolyData::SafeDownCast(data.getVtkData());
	}

	for (auto &id : commonPoints) {
		vtkSmartPointer<vtkIdList> cellIdstoPoint1 = vtkSmartPointer<vtkIdList>::New();
		if (data.getInputType() == DataFormat::PossibleInputType::vtu) {
			dataUGrid->GetPointCells(id.first, cellIdstoPoint1);
		}
		else if (data.getInputType() == DataFormat::PossibleInputType::vtp) {
			dataPoly->GetPointCells(id.first, cellIdstoPoint1);
		}


		for (vtkIdType i = 0; i < cellIdstoPoint1->GetNumberOfIds(); i++) {
			vtkIdType currentCellId = cellIdstoPoint1->GetId(i);

			if (data.getInputType() == DataFormat::PossibleInputType::vtu) {
				vtkSmartPointer<vtkIdList> currentCellPointsIDs = dataUGrid->GetCell(currentCellId)->GetPointIds();

				std::vector<vtkIdType> pointIds((int)currentCellPointsIDs->GetNumberOfIds(), 0);

				for (vtkIdType j = 0; j < currentCellPointsIDs->GetNumberOfIds(); j++) {
					vtkIdType currentPointsId = currentCellPointsIDs->GetId(j);
					if (currentPointsId == id.first) {
						pointIds[j] = id.second;
					}
					else {
						pointIds[j] = currentPointsId;
					}
				}

				dataUGrid->ReplaceCell(currentCellId, (int)data.getVtkData()->GetCell(currentCellId)->GetNumberOfPoints(), &(pointIds[0]));
			}
			else if (data.getInputType() == DataFormat::PossibleInputType::vtp) {
				dataPoly->ReplaceCellPoint(currentCellId, id.first, id.second);
			}
		}
		if (data.getInputType() == DataFormat::PossibleInputType::vtu) {
			dataUGrid->BuildLinks();
			data.setVtkData(dataUGrid);
		}
		else if (data.getInputType() == DataFormat::PossibleInputType::vtp) {
			dataPoly->BuildCells();
			dataPoly->BuildLinks();
			data.setVtkData(dataPoly);
		}
	}
} // Methods::replacePoints

/*! Calculate the resolution of a voxel mesh
 \param data Pointer with the mesh data.
 \return The resolution as veector (x,y,Z).
 */
vector<double> Methods::getXYZResolution(DataFormat &data) {
	vector<double> XYZResolution = { 0, 0, 0 };
	double bounds[6];
	data.getVtkData()->GetCellBounds(0, bounds); // (xmin, xmax, ymin, ymax, zmin, zmax)

	XYZResolution[0] = std::abs(bounds[1] - bounds[0]);

	XYZResolution[1] = std::abs(bounds[3] - bounds[2]);

	XYZResolution[2] = std::abs(bounds[5] - bounds[4]);

	std::cout << "Resolution: " << XYZResolution[0] << " x " << XYZResolution[1] << " x " << XYZResolution[2] << "\n";

	return XYZResolution;
}

/*! Create a Box with voxel.
 \param data Pointer to the orginal mesh.
 \param centrePoint Center point of the voxel box.
 \param x_width x dimission
 \param y_width y dimission
 \param z_width z dimission
 \return A DataFormat that have included the voxel box.
 */

DataFormat Methods::createVoxelBox(DataFormat &data, vector<double> centrePoint, double x_width, double y_width, double z_width) {
	DataFormat box;

	vtkSmartPointer<vtkUnstructuredGrid> voxel = vtkSmartPointer<vtkUnstructuredGrid>::New();
	vtkSmartPointer<vtkUnstructuredGrid> centrePointVoxel = vtkSmartPointer<vtkUnstructuredGrid>::New();
	vtkSmartPointer<vtkPointLocator> points = vtkSmartPointer<vtkPointLocator>::New();
	vtkSmartPointer<vtkPoints> voxelPoints = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkPoints> centrePoints = vtkSmartPointer<vtkPoints>::New();

	vector<double> XYZResolution = Methods::getXYZResolution(data);

	double pointshiftX = XYZResolution.at(0);
	double pointshiftY = XYZResolution.at(1);
	double pointshiftZ = XYZResolution.at(2);


	vtkIdType centerPointId = data.getCentrePoints()->FindPoint(centrePoint.at(0), centrePoint.at(1), centrePoint.at(2));

	centrePoint.at(0) = data.getCentrePoints()->GetPoint(centerPointId)[0];
	centrePoint.at(1) = data.getCentrePoints()->GetPoint(centerPointId)[1];
	centrePoint.at(2) = data.getCentrePoints()->GetPoint(centerPointId)[2];

	double bounds[6] = { -x_width + centrePoint.at(0), x_width + centrePoint.at(0), -y_width + centrePoint.at(1), y_width + centrePoint.at(1), -z_width + centrePoint.at(2), z_width + centrePoint.at(2) };

	points->InitPointInsertion(voxelPoints, bounds);

	for (double x = -x_width; x <= x_width; x += pointshiftX) {
		for (double y = -y_width; y <= y_width; y += pointshiftY) {
			for (double z = -z_width; z <= z_width; z += pointshiftZ) {
				vector<double> currentCentrePoint = { centrePoint.at(0) + x, centrePoint.at(1) + y, centrePoint.at(2) + z };

				centrePoints->InsertNextPoint(currentCentrePoint.at(0), currentCentrePoint.at(1), currentCentrePoint.at(2));

				vtkIdType id1;
				double point1[3] = { currentCentrePoint.at(0) - pointshiftX / 2, currentCentrePoint.at(1) - pointshiftY / 2, currentCentrePoint.at(2) - pointshiftZ / 2 };
				points->InsertUniquePoint(point1, id1);

				vtkIdType id2;
				double point2[3] = { currentCentrePoint.at(0) + pointshiftX / 2, currentCentrePoint.at(1) - pointshiftY / 2, currentCentrePoint.at(2) - pointshiftZ / 2 };
				points->InsertUniquePoint(point2, id2);

				vtkIdType id3;
				double point3[3] = { currentCentrePoint.at(0) - pointshiftX / 2, currentCentrePoint.at(1) + pointshiftY / 2, currentCentrePoint.at(2) - pointshiftZ / 2 };
				points->InsertUniquePoint(point3, id3);

				vtkIdType id4;
				double point4[3] = { currentCentrePoint.at(0) + pointshiftX / 2, currentCentrePoint.at(1) + pointshiftY / 2, currentCentrePoint.at(2) - pointshiftZ / 2 };
				points->InsertUniquePoint(point4, id4);

				vtkIdType id5;
				double point5[3] = { currentCentrePoint.at(0) - pointshiftX / 2, currentCentrePoint.at(1) - pointshiftY / 2, currentCentrePoint.at(2) + pointshiftZ / 2 };
				points->InsertUniquePoint(point5, id5);

				vtkIdType id6;
				double point6[3] = { currentCentrePoint.at(0) + pointshiftX / 2, currentCentrePoint.at(1) - pointshiftY / 2, currentCentrePoint.at(2) + pointshiftZ / 2 };
				points->InsertUniquePoint(point6, id6);

				vtkIdType id7;
				double point7[3] = { currentCentrePoint.at(0) - pointshiftX / 2, currentCentrePoint.at(1) + pointshiftY / 2, currentCentrePoint.at(2) + pointshiftZ / 2 };
				points->InsertUniquePoint(point7, id7);

				vtkIdType id8;
				double point8[3] = { currentCentrePoint.at(0) + pointshiftX / 2, currentCentrePoint.at(1) + pointshiftY / 2, currentCentrePoint.at(2) + pointshiftZ / 2 };
				points->InsertUniquePoint(point8, id8);

				vtkIdType ptIds[8] = { id1, id2, id3, id4, id5, id6, id7, id8 };

				// insert points and cells into unstructured grid
				voxel->InsertNextCell(VTK_VOXEL, 8, ptIds);
			}
		}
	}

	centrePointVoxel->SetPoints(centrePoints);
	voxel->SetPoints(points->GetPoints());

	box.setCentrePoints(centrePointVoxel);
	box.setVtkData(voxel);

	return box;
} // Methods::createVoxelBox

/*! Calcuate the resolution of the mesh by consider 10% of the cells.
 \param data Pointer to the orginal mesh.
 \return The resolution in mm of the mesh.
 */

double Methods::getResolution(DataFormat &data) {
	double resolution = 0.0;
	int m = 0;

	set<vtkIdType> viewedCells;
	int numberOfCells = (int)data.getCentrePoints()->GetPoints()->GetNumberOfPoints();

	for (vtkIdType i = 0; i < numberOfCells / 10; i++) {
		srand((unsigned int)time(NULL));
		vtkIdType currentId = rand() % numberOfCells;

		while (viewedCells.find(currentId) != viewedCells.end()) {
			currentId = rand() % numberOfCells;
		}

		viewedCells.insert(currentId);

		double cellPoint[3] = { data.getCentrePoints()->GetPoint(currentId)[0], data.getCentrePoints()->GetPoint(currentId)[1], data.getCentrePoints()->GetPoint(currentId)[2] };

		vtkSmartPointer<vtkIdList> neighbours = vtkSmartPointer<vtkIdList>::New();
		if (data.getInputType() == DataFormat::PossibleInputType::vtu) {
			for (int j = 0; j < data.getVtkData()->GetCell(currentId)->GetNumberOfFaces(); j++) {
				vtkSmartPointer<vtkIdList> facePointIdList = data.getVtkData()->GetCell(currentId)->GetFace(j)->GetPointIds();

				vtkSmartPointer<vtkIdList> tempNeighborCellIds = vtkSmartPointer<vtkIdList>::New();
				data.getVtkData()->GetCellNeighbors(currentId, facePointIdList, tempNeighborCellIds);

				for (vtkIdType k = 0; k < tempNeighborCellIds->GetNumberOfIds(); k++) {
					neighbours->InsertUniqueId(tempNeighborCellIds->GetId(k));
				}
			}
		}
		else if (data.getInputType() == DataFormat::PossibleInputType::vtp) {
			for (int j = 0; j < data.getVtkData()->GetCell(currentId)->GetNumberOfEdges(); j++) {
				vtkSmartPointer<vtkIdList> edgePointIdList = data.getVtkData()->GetCell(currentId)->GetEdge(j)->GetPointIds();

				vtkSmartPointer<vtkIdList> tempNeighborCellIds = vtkSmartPointer<vtkIdList>::New();
				data.getVtkData()->GetCellNeighbors(currentId, edgePointIdList, tempNeighborCellIds);

				for (vtkIdType k = 0; k < tempNeighborCellIds->GetNumberOfIds(); k++) {
					neighbours->InsertUniqueId(tempNeighborCellIds->GetId(k));
				}
			}
		}

		for (vtkIdType l = 0; l < neighbours->GetNumberOfIds(); l++) {
			if (viewedCells.find(neighbours->GetId(l)) == viewedCells.end()) {
				double distance = round(sqrt(abs(vtkMath::Distance2BetweenPoints(cellPoint, data.getCentrePoints()->GetPoint(neighbours->GetId(l))))), 2);

				resolution += distance;
				m++;
			}
		}
	}

	resolution = resolution / m;
	resolution = round(resolution, 2);
	return resolution;
} // Methods::getResolution

/*! Create a voxels around a path in a define radius and resolution.
 \param data Pointer to the orginal mesh.
 \param path Path around the voxel are created.
 \param radius Radius in mm
 \param resolution Resolution of the mesh
 \return A DataFormat that have included the voxel around the path.
 */
DataFormat Methods::createVoxelAroundPath(DataFormat &data, vtkSmartPointer<vtkIdList> path, double radius) {
	vector<double> centroidandBounds = getCentroidandBounds(data, path);

	double boundingBoxOverAll[6] = {
	centroidandBounds.at(3) - radius, centroidandBounds.at(4) + radius, centroidandBounds.at(5) - radius, centroidandBounds.at(6) + radius, centroidandBounds.at(7) - radius, centroidandBounds.at(8) + radius
	};

	vector<double> centrepoint = { centroidandBounds.at(0), centroidandBounds.at(1), centroidandBounds.at(2) };

	DataFormat box;

	box = createVoxelBox(data, centrepoint, boundingBoxOverAll[1] - boundingBoxOverAll[0], boundingBoxOverAll[3] - boundingBoxOverAll[2], boundingBoxOverAll[5] - boundingBoxOverAll[4]);

	double startPointArray[3] = { data.getCentrePoints()->GetPoint(path->GetId(0))[0], data.getCentrePoints()->GetPoint(path->GetId(0))
	[1], data.getCentrePoints()->GetPoint(path->GetId(0))[2] };

	vtkIdType startPointIdBox = box.getCentrePoints()->FindPoint(startPointArray);
	double startPointBoxArray[3] = { box.getCentrePoints()->GetPoint(startPointIdBox)[0], box.getCentrePoints()->GetPoint(startPointIdBox)[1], box.getCentrePoints()->GetPoint(startPointIdBox)[2] };

	double transformVectorArray[3] = { 0, 0, 0 };
	vtkMath::Subtract(startPointArray, startPointBoxArray, transformVectorArray);
	vector<double> transformVector = { transformVectorArray[0], transformVectorArray[1], transformVectorArray[2] };
	box.getVtkData()->SetPoints(movePoints(box.getVtkData()->GetPoints(), transformVector));
	box.setCentrePoints(calcCellCentroids(box.getVtkData()));

	set<vtkIdType> visibleCells;
	vtkSmartPointer<vtkIdList> cylinderCells = vtkSmartPointer<vtkIdList>::New();

	for (vtkIdType i = 0; i < path->GetNumberOfIds(); i++) {
		vector<double> currentPoint = { data.getCentrePoints()->GetPoint(path->GetId(i))[0], data.getCentrePoints()->GetPoint(path->GetId(i))[1], data.getCentrePoints()->GetPoint(path->GetId(i))[2] };

		vtkSmartPointer<vtkIdList> possibleCells = getCellsinRadius(box, currentPoint, radius);
		for (vtkIdType i = 0; i < possibleCells->GetNumberOfIds(); i++) {
			vtkIdType currentID = possibleCells->GetId(i);
			if (visibleCells.find(currentID) == visibleCells.end()) {
				cylinderCells->InsertNextId(currentID);
			}
			visibleCells.insert(currentID);
		}
	}

	vtkSmartPointer<vtkUnstructuredGrid> cylinder = extractCellsVtu(box.getVtkData(), cylinderCells);

	box.setVtkData(cylinder);
	box.setCentrePoints(calcCellCentroids(cylinder));
	box.setInputType(data.getInputType());

	return box;
} // Methods::createVoxelAroundPath

/*! Create a voxel tube between two points in a define bridgewidth and resolution.
 \param data Pointer to the orginal mesh.
 \param startPoint Start point of the tube as vector (x,y,z).
 \param targetPoint Target point of the tube as vector (x,y,z).
 \param resolution Resolution of the mesh.
 \param bridgewidth Radius of the bridge in mm.
 \return A DataFormat that have included the voxel around the path.
 */

DataFormat Methods::createVoxelTube(DataFormat &data, vector<double> startPoint, vector<double> targetPoint, double resolution, double bridgewidth) {
	double startPointArray[3] = { startPoint.at(0), startPoint.at(1), startPoint.at(2) };
	double targetPointArray[3] = { targetPoint.at(0), targetPoint.at(1), targetPoint.at(2) };

	double boundingBoxstartPoint[6] = { startPointArray[0] - bridgewidth, startPointArray[0] + bridgewidth, startPointArray[1] - bridgewidth, startPointArray[1] + bridgewidth, startPointArray[2] - bridgewidth, startPointArray[2] + bridgewidth };

	double boundingBoxtargetPoint[6] = { targetPointArray[0] - bridgewidth, targetPointArray[0] + bridgewidth, targetPointArray[1] - bridgewidth, targetPointArray[1] + bridgewidth, targetPointArray[2] - bridgewidth, targetPointArray[2] + bridgewidth };

	double boundingBoxOverAll[6];

	for (int i = 0; i < 6; i += 2) {
		if (boundingBoxstartPoint[i] < boundingBoxtargetPoint[i])
			boundingBoxOverAll[i] = boundingBoxstartPoint[i];
		else
			boundingBoxOverAll[i] = boundingBoxtargetPoint[i];

		if (boundingBoxstartPoint[i + 1] > boundingBoxtargetPoint[i + 1])
			boundingBoxOverAll[i + 1] = boundingBoxstartPoint[i + 1];
		else
			boundingBoxOverAll[i + 1] = boundingBoxtargetPoint[i + 1];
	}

	vector<double> centrepoint = { (boundingBoxOverAll[0] + boundingBoxOverAll[1]) / 2, (boundingBoxOverAll[2] + boundingBoxOverAll[3]) / 2, (boundingBoxOverAll[4] + boundingBoxOverAll[5]) / 2 };

	DataFormat box;
	box = createVoxelBox(data, centrepoint, boundingBoxOverAll[1] - boundingBoxOverAll[0], boundingBoxOverAll[3] - boundingBoxOverAll[2], boundingBoxOverAll[5] - boundingBoxOverAll[4]);


	vtkIdType startPointIdBox = box.getCentrePoints()->FindPoint(startPointArray);

	double startPointBoxArray[3] = { box.getCentrePoints()->GetPoint(startPointIdBox)[0], box.getCentrePoints()->GetPoint(startPointIdBox)[1], box.getCentrePoints()->GetPoint(startPointIdBox)[2] };

	double transformVectorArray[3] = { 0, 0, 0 };
	vtkMath::Subtract(startPointArray, startPointBoxArray, transformVectorArray);
	vector<double> transformVector = { transformVectorArray[0], transformVectorArray[1], transformVectorArray[2] };
	box.getVtkData()->SetPoints(movePoints(box.getVtkData()->GetPoints(), transformVector));
	box.setCentrePoints(calcCellCentroids(box.getVtkData()));

	double normalVector[3];
	vtkMath::Subtract(targetPointArray, startPointArray, normalVector);
	vtkMath::Normalize(normalVector);

	double shiftVector[3] = { normalVector[0], normalVector[1], normalVector[2] };
	vtkMath::MultiplyScalar(shiftVector, resolution);

	double currentLinePointArray[3] = { startPoint.at(0), startPoint.at(1), startPoint.at(2) };

	set<vtkIdType> visibleCells;
	vtkSmartPointer<vtkIdList> cylinderCells = vtkSmartPointer<vtkIdList>::New();

	bool exit = false;

	while (!exit) {
		vector<double> currentLinePoint = { currentLinePointArray[0], currentLinePointArray[1], currentLinePointArray[2] };
		vtkSmartPointer<vtkIdList> possibleCells = getCellsinRadius(box, currentLinePoint, bridgewidth);
		for (vtkIdType i = 0; i < possibleCells->GetNumberOfIds(); i++) {
			vtkIdType currentID = possibleCells->GetId(i);
			if (visibleCells.find(currentID) == visibleCells.end()) {
				cylinderCells->InsertNextId(currentID);
			}
			visibleCells.insert(currentID);
		}

		vtkMath::Add(currentLinePointArray, shiftVector, currentLinePointArray);
		double evaluatePlaneTarget = vtkPlane::Evaluate(normalVector, targetPointArray, currentLinePointArray);
		if (evaluatePlaneTarget > 0)
			exit = true;
	}

	vtkSmartPointer<vtkUnstructuredGrid> cylinder = extractCellsVtu(box.getVtkData(), cylinderCells);


	box.setVtkData(cylinder);
	box.setCentrePoints(calcCellCentroids(cylinder));
	box.setInputType(data.getInputType());

	return box;
} // Methods::createVoxelTube

/*! Extract the cells from a vtu which are in the list.
 \param inputdata Mesh were the cells would be extract.
 \param cellIDs Cells-ids that would be extract form the mesh.
 \return The extracted cells.
 */
vtkSmartPointer<vtkUnstructuredGrid> Methods::extractCellsVtu(vtkSmartPointer<vtkPointSet> inputdata, vtkSmartPointer<vtkIdList> cellIDs) {
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkPointLocator> pointsLocator = vtkSmartPointer<vtkPointLocator>::New();

	pointsLocator->InitPointInsertion(points, inputdata->GetBounds());

	vtkSmartPointer<vtkUnstructuredGrid> outputdata = vtkSmartPointer<vtkUnstructuredGrid>::New();
	for (vtkIdType i = 0; i < cellIDs->GetNumberOfIds(); i++) {
		vtkSmartPointer<vtkIdList> cellPoints = vtkSmartPointer<vtkIdList>::New();
		inputdata->GetCellPoints(cellIDs->GetId(i), cellPoints);
		for (vtkIdType j = 0; j < cellPoints->GetNumberOfIds(); j++) {
			vtkIdType id = -1;
			double *point = inputdata->GetPoint(cellPoints->GetId(j));
			pointsLocator->InsertUniquePoint(point, id);
			cellPoints->SetId(j, id);
		}
		outputdata->InsertNextCell(inputdata->GetCellType(cellIDs->GetId(i)), cellPoints);
	}

	outputdata->SetPoints(pointsLocator->GetPoints());
	outputdata->BuildLinks();

	return outputdata;
}

/*! Extract the cells from a vtp which are in the list.
 \param inputdata Mesh were the cells would be extract.
 \param cellIDs Cells-ids that would be extract form the mesh.
 \return The extracted cells.
 */
vtkSmartPointer<vtkPolyData> Methods::extractCellsVtp(vtkSmartPointer<vtkPolyData> inputdata, vtkSmartPointer<vtkIdList> cellIDs) {
	vtkSmartPointer<vtkIdTypeArray> cells = vtkSmartPointer<vtkIdTypeArray>::New();

	cells->SetNumberOfComponents(1);

	for (vtkIdType i = 0; i < cellIDs->GetNumberOfIds(); i++) {
		cells->InsertNextValue(cellIDs->GetId(i));
	}


	vtkSmartPointer<vtkSelectionNode> selectionNode = vtkSmartPointer<vtkSelectionNode>::New();
	selectionNode->SetFieldType(vtkSelectionNode::CELL);
	selectionNode->SetContentType(vtkSelectionNode::INDICES);
	selectionNode->SetSelectionList(cells);

	vtkSmartPointer<vtkSelection> selection = vtkSmartPointer<vtkSelection>::New();
	selection->AddNode(selectionNode);

	vtkSmartPointer<vtkExtractSelection> extractSelection = vtkSmartPointer<vtkExtractSelection>::New();
	extractSelection->SetInputData(0, inputdata);
	extractSelection->SetInputData(1, selection);
	extractSelection->Update();

	vtkSmartPointer<vtkUnstructuredGrid> datavtu = vtkSmartPointer<vtkUnstructuredGrid>::New();
	datavtu->ShallowCopy(extractSelection->GetOutput());

	DataFormat data;
	data.setVtkData(datavtu);
	vtkSmartPointer<vtkPolyData> outputdata = extractSurface(data);
	outputdata = triangulateSurface(outputdata);

	return outputdata;
} // Methods::extractCellsVtp

/*! Calc the cell center of all cells.
 \param inputdata Mesh with the cells.
 \return Unstructured grid with all cell center points. The points have the same id as th cells.
 */
vtkSmartPointer<vtkUnstructuredGrid> Methods::calcCellCentroids(vtkSmartPointer<vtkPointSet> inputdata) {
	vtkSmartPointer<vtkUnstructuredGrid> centrepoints = vtkSmartPointer<vtkUnstructuredGrid>::New();
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

	for (vtkIdType i = 0; i < inputdata->GetNumberOfCells(); i++) {
		vtkSmartPointer<vtkPoints> pointsofCell = inputdata->GetCell(i)->GetPoints();
		long numofPoints = pointsofCell->GetNumberOfPoints();
		double point[3] = { 0, 0, 0 };
		for (int i = 0; i < numofPoints; i++) {
			vtkMath::Add(point, pointsofCell->GetPoint(i), point);
		}
		points->InsertNextPoint(point[0] / numofPoints, point[1] / numofPoints, point[2] / numofPoints);
	}
	centrepoints->SetPoints(points);
	return centrepoints;
}

/*! Return all cell-ids of a mesh.
 \param data Pointer with mesh data.
 \return A list with all cell-id of mesh.
 */
vtkSmartPointer<vtkIdList> Methods::getAllCellIds(DataFormat &data) {
	vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();

	for (vtkIdType i = 0; i < data.getVtkData()->GetNumberOfCells(); i++) {
		cellIds->InsertNextId(i);
	}

	return cellIds;
}

/*! Return the cell-ids which are not in both mesh.
 \param newdata Pointer to the new mesh data.
 \param olddata Pointer to the old mesh data.
 \return A list with the cell-id of the cell which are not equal in both meshs.
 */
vtkSmartPointer<vtkIdList> Methods::getCellsWereNotinOldData(DataFormat &newData, DataFormat &oldData) {
	vtkSmartPointer<vtkIdList> newCellIDs = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkPolyData> surface = extractSurface(oldData);

	vtkSmartPointer<vtkSelectEnclosedPoints> selectEnclosedPoints = vtkSmartPointer<vtkSelectEnclosedPoints>::New();

	selectEnclosedPoints->Initialize(surface);
	selectEnclosedPoints->SetTolerance(0);

	for (vtkIdType i = 0; i < newData.getCentrePoints()->GetNumberOfPoints(); i++) {
		double *point = newData.getCentrePoints()->GetPoint(i);
		int inside = selectEnclosedPoints->IsInsideSurface(point);
		if (inside != 1) {
			newCellIDs->InsertNextId(i);
		}
	}

	return newCellIDs;
}

/*! Return the cell-ids which are in both mesh.
 \param newdata Pointer to the new mesh data.
 \param olddata Pointer to the old mesh data.
 \return A list with the cell-id of the cell which are equal in both meshs.
 */
vtkSmartPointer<vtkIdList> Methods::getCellsWereInOldData(DataFormat &newData, DataFormat &oldData) {
	vtkSmartPointer<vtkPointLocator> centrePointLocator = vtkSmartPointer<vtkPointLocator>::New();

	centrePointLocator->SetDataSet(oldData.getCentrePoints());
	centrePointLocator->BuildLocator();

	vtkSmartPointer<vtkIdList> dupliactedCells = vtkSmartPointer<vtkIdList>::New();

	vtkSmartPointer<vtkIdList> cellIDsnewData = getAllCellIds(newData);

	for (vtkIdType i = 0; i < cellIDsnewData->GetNumberOfIds(); i++) {
		double point[3] = { newData.getCentrePoints()->GetPoint(cellIDsnewData->GetId(i))[0], newData.getCentrePoints()->GetPoint(cellIDsnewData->GetId(i))[1], newData.getCentrePoints()->GetPoint(cellIDsnewData->GetId(i))[2] };
		vtkIdType oldCellID = centrePointLocator->FindClosestPoint(point);
		if ((oldData.getCentrePoints()->GetPoint(oldCellID)[0] == point[0]) &&
			(oldData.getCentrePoints()->GetPoint(oldCellID)[1] == point[1]) &&
			(oldData.getCentrePoints()->GetPoint(oldCellID)[2] == point[2])) {
			dupliactedCells->InsertNextId(i);
		}
		else {
			double distance = sqrt(abs(vtkMath::Distance2BetweenPoints(point, oldData.getCentrePoints()->GetPoint(oldCellID))));
			if (distance <= 0.00001)
				dupliactedCells->InsertNextId(i);
		}
	}


	return dupliactedCells;
} // Methods::getCellsWereInOldData

/*! Return the point-ids which are in both mesh.
 \param newdata Pointer to the new mesh data.
 \param olddata Pointer to the old mesh data.
 \return A list with the point-id of the points which are equal in both meshs.
 */
vtkSmartPointer<vtkIdList> Methods::getPointsWereInOldData(DataFormat &newData, DataFormat &oldData) {
	vtkSmartPointer<vtkIdList> newCellIDs = vtkSmartPointer<vtkIdList>::New();

	vtkSmartPointer<vtkSelectEnclosedPoints> selectEnclosedPoints = vtkSmartPointer<vtkSelectEnclosedPoints>::New();

	selectEnclosedPoints->Initialize(extractSurface(oldData));
	selectEnclosedPoints->SetTolerance(0);

	for (vtkIdType i = 0; i < newData.getVtkData()->GetNumberOfPoints(); i++) {
		double *point = newData.getVtkData()->GetPoint(i);
		int inside = selectEnclosedPoints->IsInsideSurface(point);
		if (inside == 1) {
			newCellIDs->InsertNextId(i);
		}
	}

	return newCellIDs;
}

/*! Union to mesh together in one mesh.
 \param olddata Pointer to the old mesh data.
 \param newdata Pointer to the new mesh data.
 \return A list which are insert in the old mesh.
 */
vtkSmartPointer<vtkIdList> Methods::unionData(DataFormat &olddata, DataFormat &newData) {
	vtkSmartPointer<vtkIdList> insertCellIDs = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();


	vtkSmartPointer<vtkDoubleArray> material = vtkSmartPointer<vtkDoubleArray>::New();

	material->SetName("Material");
	material->SetNumberOfComponents(1);
	material->SetNumberOfTuples(olddata.getVtkData()->GetNumberOfCells() + newData.getVtkData()->GetNumberOfCells());
	material->FillComponent(0, 0);

	vtkSmartPointer<vtkDoubleArray> rawMaterial = vtkSmartPointer<vtkDoubleArray>::New();
	rawMaterial->SetName("RawMaterial");
	rawMaterial->SetNumberOfComponents(1);
	rawMaterial->SetNumberOfTuples(olddata.getVtkData()->GetNumberOfCells() + newData.getVtkData()->GetNumberOfCells());
	rawMaterial->FillComponent(0, 0);

	vtkSmartPointer<vtkDoubleArray> theta = vtkSmartPointer<vtkDoubleArray>::New();
	theta->SetName("Theta");
	theta->SetNumberOfComponents(1);
	theta->SetNumberOfTuples(olddata.getVtkData()->GetNumberOfCells() + newData.getVtkData()->GetNumberOfCells());
	theta->FillComponent(0, 0);

	vtkSmartPointer<vtkDoubleArray> phi = vtkSmartPointer<vtkDoubleArray>::New();
	phi->SetName("Phi");
	phi->SetNumberOfComponents(1);
	phi->SetNumberOfTuples(olddata.getVtkData()->GetNumberOfCells() + newData.getVtkData()->GetNumberOfCells());
	phi->FillComponent(0, 0);

	vtkSmartPointer<vtkDoubleArray> thetaPhi = vtkSmartPointer<vtkDoubleArray>::New();
	thetaPhi->SetName("ThetaPhi");
	thetaPhi->SetNumberOfComponents(2);
	thetaPhi->SetNumberOfTuples(olddata.getVtkData()->GetNumberOfCells() + newData.getVtkData()->GetNumberOfCells());
	thetaPhi->SetComponentName(0, "Theta");
	thetaPhi->SetComponentName(1, "Phi");
	thetaPhi->FillComponent(0, 0);
	thetaPhi->FillComponent(1, 0);

	vtkSmartPointer<vtkDoubleArray> diffVec = vtkSmartPointer<vtkDoubleArray>::New();
	diffVec->SetName("DifferenceVector");
	diffVec->SetNumberOfComponents(3);
	diffVec->SetNumberOfTuples(olddata.getVtkData()->GetNumberOfCells() + newData.getVtkData()->GetNumberOfCells());
	diffVec->FillComponent(0, 0);
	diffVec->FillComponent(1, 0);
	diffVec->FillComponent(2, 0);

	vtkSmartPointer<vtkDoubleArray> baseDiffVec = vtkSmartPointer<vtkDoubleArray>::New();
	baseDiffVec->SetName("BasedPathDifferenceVector");
	baseDiffVec->SetNumberOfComponents(3);
	baseDiffVec->SetNumberOfTuples(olddata.getVtkData()->GetNumberOfCells() + newData.getVtkData()->GetNumberOfCells());
	baseDiffVec->FillComponent(0, 0);
	baseDiffVec->FillComponent(1, 0);
	baseDiffVec->FillComponent(2, 0);

	vtkSmartPointer<vtkDoubleArray> status = vtkSmartPointer<vtkDoubleArray>::New();
	status->SetName("Status");
	status->SetNumberOfComponents(1);
	status->FillComponent(0, 0);


	if (olddata.getInputType() == DataFormat::PossibleInputType::vtu) {
		vtkSmartPointer<vtkPointLocator> pointLocator = vtkSmartPointer<vtkPointLocator>::New();

		pointLocator->InitPointInsertion(points, olddata.getVtkData()->GetBounds());

		for (vtkIdType i = 0; i < olddata.getVtkData()->GetPoints()->GetNumberOfPoints(); i++) {
			pointLocator->InsertNextPoint(olddata.getVtkData()->GetPoint(i));
		}

		vtkSmartPointer<vtkUnstructuredGrid> uGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

		uGrid->Allocate(olddata.getVtkData()->GetNumberOfCells() + newData.getVtkData()->GetNumberOfCells());

		uGrid->GetCellData()->AddArray(material);
		uGrid->GetCellData()->AddArray(rawMaterial);
		uGrid->GetCellData()->AddArray(theta);
		uGrid->GetCellData()->AddArray(phi);
		uGrid->GetCellData()->AddArray(thetaPhi);
		uGrid->GetCellData()->AddArray(diffVec);
		uGrid->GetCellData()->AddArray(baseDiffVec);
		uGrid->GetCellData()->AddArray(status);


		for (vtkIdType j = 0; j < olddata.getVtkData()->GetNumberOfCells(); j++) {
			vtkIdType newcellID = -1;
			newcellID = uGrid->InsertNextCell(olddata.getVtkData()->GetCellType(j), olddata.getVtkData()->GetCell(j)->GetPointIds());

			uGrid->GetCellData()->GetArray("Status")->InsertTuple(newcellID, olddata.getVtkData()->GetCellData()->GetArray("Status")->GetTuple(j));
			uGrid->GetCellData()->GetArray("Material")->InsertTuple(newcellID, olddata.getVtkData()->GetCellData()->GetArray("Material")->GetTuple(j));
			uGrid->GetCellData()->GetArray("RawMaterial")->InsertTuple(newcellID, olddata.getVtkData()->GetCellData()->GetArray("RawMaterial")->GetTuple(j));
			uGrid->GetCellData()->GetArray("Phi")->InsertTuple(newcellID, olddata.getVtkData()->GetCellData()->GetArray("Phi")->GetTuple(j));
			uGrid->GetCellData()->GetArray("Theta")->InsertTuple(newcellID, olddata.getVtkData()->GetCellData()->GetArray("Theta")->GetTuple(j));
			uGrid->GetCellData()->GetArray("ThetaPhi")->InsertTuple(newcellID, olddata.getVtkData()->GetCellData()->GetArray("ThetaPhi")->GetTuple(j));
			uGrid->GetCellData()->GetArray("DifferenceVector")->InsertTuple(newcellID, olddata.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetTuple(j));
			uGrid->GetCellData()->GetArray("BasedPathDifferenceVector")->InsertTuple(newcellID, olddata.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(j));
		}


		for (vtkIdType i = 0; i < newData.getVtkData()->GetNumberOfCells(); i++) {
			vtkSmartPointer<vtkIdList> cellPointsIds = vtkSmartPointer<vtkIdList>::New();
			newData.getVtkData()->GetCellPoints(i, cellPointsIds);

			vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();

			for (vtkIdType j = 0; j < cellPointsIds->GetNumberOfIds(); j++) {
				double *point = newData.getVtkData()->GetPoint(cellPointsIds->GetId(j));
				vtkIdType id;
				pointLocator->InsertUniquePoint(point, id);
				ptIds->InsertNextId(id);
			}

			vtkIdType newcellID = -1;
			int cellType = (int)newData.getVtkData()->GetCellType(i);
			newcellID = uGrid->InsertNextCell(cellType, ptIds);

			insertCellIDs->InsertNextId(newcellID);

			olddata.getCentrePoints()->GetPoints()->InsertPoint(newcellID, newData.getCentrePoints()->GetPoint(i));

			uGrid->GetCellData()->GetArray("Status")->InsertComponent(newcellID, 0, 0);
			uGrid->GetCellData()->GetArray("Material")->InsertComponent(newcellID, 0, 300);
			uGrid->GetCellData()->GetArray("RawMaterial")->InsertComponent(newcellID, 0, 300);
			uGrid->GetCellData()->GetArray("Phi")->InsertComponent(newcellID, 0, 0);
			uGrid->GetCellData()->GetArray("Theta")->InsertComponent(newcellID, 0, 0);
			uGrid->GetCellData()->GetArray("ThetaPhi")->InsertTuple2(newcellID, 0, 0);
			uGrid->GetCellData()->GetArray("DifferenceVector")->InsertTuple3(newcellID, 0, 0, 0);
			uGrid->GetCellData()->GetArray("BasedPathDifferenceVector")->InsertTuple3(newcellID, 0, 0, 0);
		}
		uGrid->SetPoints(pointLocator->GetPoints());
		uGrid->BuildLinks();
		olddata.setVtkData(uGrid);
	}
	else if (olddata.getInputType() == DataFormat::PossibleInputType::vtp) {
		vtkSmartPointer<vtkMergePoints> pointLocator = vtkSmartPointer<vtkMergePoints>::New();
		pointLocator->InitPointInsertion(points, olddata.getVtkData()->GetBounds());

		for (vtkIdType i = 0; i < olddata.getVtkData()->GetPoints()->GetNumberOfPoints(); i++) {
			pointLocator->InsertNextPoint(olddata.getVtkData()->GetPoint(i));
		}

		vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
		vtkSmartPointer<vtkCellArray> cellArray = vtkSmartPointer<vtkCellArray>::New();

		polydata->Allocate(olddata.getVtkData()->GetNumberOfCells() + newData.getVtkData()->GetNumberOfCells());

		polydata->GetCellData()->AddArray(material);
		polydata->GetCellData()->AddArray(rawMaterial);
		polydata->GetCellData()->AddArray(theta);
		polydata->GetCellData()->AddArray(phi);
		polydata->GetCellData()->AddArray(thetaPhi);
		polydata->GetCellData()->AddArray(diffVec);
		polydata->GetCellData()->AddArray(baseDiffVec);
		polydata->GetCellData()->AddArray(status);

		if (!olddata.getDoubleLayer()) {
			vtkSmartPointer<vtkDoubleArray> materialEndo = vtkSmartPointer<vtkDoubleArray>::New();
			materialEndo->SetName("MaterialEndo");
			materialEndo->SetNumberOfComponents(1);
			materialEndo->SetNumberOfTuples(olddata.getVtkData()->GetNumberOfCells() +
				newData.getVtkData()->GetNumberOfCells());
			materialEndo->FillComponent(0, 0);

			vtkSmartPointer<vtkDoubleArray> rawMaterialEndo = vtkSmartPointer<vtkDoubleArray>::New();
			rawMaterialEndo->SetName("RawMaterialEndo");
			rawMaterialEndo->SetNumberOfComponents(1);
			rawMaterialEndo->SetNumberOfTuples(olddata.getVtkData()->GetNumberOfCells() + newData.getVtkData()->GetNumberOfCells());
			rawMaterialEndo->FillComponent(0, 0);

			vtkSmartPointer<vtkDoubleArray> thetaEndo = vtkSmartPointer<vtkDoubleArray>::New();
			thetaEndo->SetName("ThetaEndo");
			thetaEndo->SetNumberOfComponents(1);
			thetaEndo->SetNumberOfTuples(olddata.getVtkData()->GetNumberOfCells() + newData.getVtkData()->GetNumberOfCells());
			thetaEndo->FillComponent(0, 0);

			vtkSmartPointer<vtkDoubleArray> phiEndo = vtkSmartPointer<vtkDoubleArray>::New();
			phiEndo->SetName("PhiEndo");
			phiEndo->SetNumberOfComponents(1);
			phiEndo->SetNumberOfTuples(olddata.getVtkData()->GetNumberOfCells() + newData.getVtkData()->GetNumberOfCells());
			phiEndo->FillComponent(0, 0);

			vtkSmartPointer<vtkDoubleArray> thetaPhiEndo = vtkSmartPointer<vtkDoubleArray>::New();
			thetaPhiEndo->SetName("ThetaPhiEndo");
			thetaPhiEndo->SetNumberOfComponents(2);
			thetaPhiEndo->SetNumberOfTuples(olddata.getVtkData()->GetNumberOfCells() +
				newData.getVtkData()->GetNumberOfCells());
			thetaPhiEndo->SetComponentName(0, "ThetaEndo");
			thetaPhiEndo->SetComponentName(1, "PhiEndo");
			thetaPhiEndo->FillComponent(0, 0);
			thetaPhiEndo->FillComponent(1, 0);

			vtkSmartPointer<vtkDoubleArray> diffVecEndo = vtkSmartPointer<vtkDoubleArray>::New();
			diffVecEndo->SetName("DifferenceVectorEndo");
			diffVecEndo->SetNumberOfComponents(3);
			diffVecEndo->SetNumberOfTuples(olddata.getVtkData()->GetNumberOfCells() +
				newData.getVtkData()->GetNumberOfCells());
			diffVecEndo->FillComponent(0, 0);
			diffVecEndo->FillComponent(1, 0);
			diffVecEndo->FillComponent(2, 0);

			vtkSmartPointer<vtkDoubleArray> statusEndo = vtkSmartPointer<vtkDoubleArray>::New();
			statusEndo->SetName("StatusEndo");
			statusEndo->SetNumberOfTuples(olddata.getVtkData()->GetNumberOfCells() + newData.getVtkData()->GetNumberOfCells());
			statusEndo->SetNumberOfComponents(1);
			statusEndo->FillComponent(0, 0);

			polydata->GetCellData()->AddArray(materialEndo);
			polydata->GetCellData()->AddArray(rawMaterialEndo);
			polydata->GetCellData()->AddArray(thetaEndo);
			polydata->GetCellData()->AddArray(phiEndo);
			polydata->GetCellData()->AddArray(thetaPhiEndo);
			polydata->GetCellData()->AddArray(diffVecEndo);
			polydata->GetCellData()->AddArray(statusEndo);
		}


		for (int k = 0; k < olddata.getVtkData()->GetNumberOfCells(); k++) {
			vtkIdType cellID = -1;
			cellID = cellArray->InsertNextCell(olddata.getVtkData()->GetCell(k));

			polydata->GetCellData()->GetArray("Status")->InsertTuple(cellID, olddata.getVtkData()->GetCellData()->GetArray("Status")->GetTuple(k));
			polydata->GetCellData()->GetArray("Material")->InsertTuple(cellID, olddata.getVtkData()->GetCellData()->GetArray("Material")->GetTuple(k));
			polydata->GetCellData()->GetArray("RawMaterial")->InsertTuple(cellID, olddata.getVtkData()->GetCellData()->GetArray("RawMaterial")->GetTuple(k));
			polydata->GetCellData()->GetArray("Phi")->InsertTuple(cellID, olddata.getVtkData()->GetCellData()->GetArray("Phi")->GetTuple(k));
			polydata->GetCellData()->GetArray("Theta")->InsertTuple(cellID, olddata.getVtkData()->GetCellData()->GetArray("Theta")->GetTuple(k));
			polydata->GetCellData()->GetArray("ThetaPhi")->InsertTuple(cellID, olddata.getVtkData()->GetCellData()->GetArray("ThetaPhi")->GetTuple(k));
			polydata->GetCellData()->GetArray("DifferenceVector")->InsertTuple(cellID, olddata.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetTuple(k));
			polydata->GetCellData()->GetArray("BasedPathDifferenceVector")->InsertTuple(cellID, olddata.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(k));


			if (!olddata.getDoubleLayer()) {
				polydata->GetCellData()->GetArray("StatusEndo")->InsertTuple(cellID, olddata.getVtkData()->GetCellData()->GetArray("StatusEndo")->GetTuple(k));
				polydata->GetCellData()->GetArray("MaterialEndo")->InsertTuple(cellID, olddata.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetTuple(k));
				polydata->GetCellData()->GetArray("RawMaterialEndo")->InsertTuple(cellID, olddata.getVtkData()->GetCellData()->GetArray("RawMaterialEndo")->GetTuple(k));
				polydata->GetCellData()->GetArray("PhiEndo")->InsertTuple(cellID, olddata.getVtkData()->GetCellData()->GetArray("PhiEndo")->GetTuple(k));
				polydata->GetCellData()->GetArray("ThetaEndo")->InsertTuple(cellID, olddata.getVtkData()->GetCellData()->GetArray("ThetaEndo")->GetTuple(k));
				polydata->GetCellData()->GetArray("ThetaPhiEndo")->InsertTuple(cellID, olddata.getVtkData()->GetCellData()->GetArray("ThetaPhiEndo")->GetTuple(k));
				polydata->GetCellData()->GetArray("DifferenceVectorEndo")->InsertTuple(cellID, olddata.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->GetTuple(k));
			}
		}

		for (vtkIdType i = 0; i < newData.getVtkData()->GetNumberOfCells(); i++) {
			vtkSmartPointer<vtkIdList> cellPointsIds = vtkSmartPointer<vtkIdList>::New();
			newData.getVtkData()->GetCellPoints(i, cellPointsIds);

			vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();

			for (vtkIdType j = 0; j < cellPointsIds->GetNumberOfIds(); j++) {
				double *point = newData.getVtkData()->GetPoint(cellPointsIds->GetId(j));
				vtkIdType id;
				pointLocator->InsertUniquePoint(point, id);
				ptIds->InsertNextId(id);
			}

			vtkIdType newcellID = -1;


			vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
			for (vtkIdType j = 0; j < ptIds->GetNumberOfIds(); j++) {
				triangle->GetPointIds()->SetId(j, ptIds->GetId(j));
			}
			newcellID = cellArray->InsertNextCell(triangle);

			insertCellIDs->InsertNextId(newcellID);

			olddata.getCentrePoints()->GetPoints()->InsertPoint(newcellID, newData.getCentrePoints()->GetPoint(i));

			polydata->GetCellData()->GetArray("Status")->InsertComponent(newcellID, 0, 0);
			polydata->GetCellData()->GetArray("Material")->InsertComponent(newcellID, 0, 300);
			polydata->GetCellData()->GetArray("RawMaterial")->InsertComponent(newcellID, 0, 300);
			polydata->GetCellData()->GetArray("Phi")->InsertComponent(newcellID, 0, 0);
			polydata->GetCellData()->GetArray("Theta")->InsertComponent(newcellID, 0, 0);
			polydata->GetCellData()->GetArray("ThetaPhi")->InsertTuple2(newcellID, 0, 0);
			polydata->GetCellData()->GetArray("DifferenceVector")->InsertTuple3(newcellID, 0, 0, 0);
			polydata->GetCellData()->GetArray("BasedPathDifferenceVector")->InsertTuple3(newcellID, 0, 0, 0);


			if (!olddata.getDoubleLayer()) {
				polydata->GetCellData()->GetArray("StatusEndo")->InsertComponent(newcellID, 0, 0);
				polydata->GetCellData()->GetArray("MaterialEndo")->InsertComponent(newcellID, 0, 300);
				polydata->GetCellData()->GetArray("RawMaterialEndo")->InsertComponent(newcellID, 0, 300);
				polydata->GetCellData()->GetArray("PhiEndo")->InsertComponent(newcellID, 0, 0);
				polydata->GetCellData()->GetArray("ThetaEndo")->InsertComponent(newcellID, 0, 0);
				polydata->GetCellData()->GetArray("ThetaPhiEndo")->InsertTuple2(newcellID, 0, 0);
				polydata->GetCellData()->GetArray("DifferenceVectorEndo")->InsertTuple3(newcellID, 0, 0, 0);
			}
		}

		polydata->SetPoints(pointLocator->GetPoints());
		polydata->SetPolys(cellArray);
		polydata->BuildCells();
		polydata->BuildLinks();
		olddata.setVtkData(polydata);
	}

	return insertCellIDs;
} // Methods::unionData

/*! Extract the surface of a mesh
 \param data Pointer to the mesh data.
 \return The surface of the mesh.
 */
vtkSmartPointer<vtkPolyData> Methods::extractSurface(DataFormat &data) {
	vtkSmartPointer<vtkPointSet> dataSurface = data.getVtkData();

	vtkSmartPointer<vtkGeometryFilter> surfaceFilter = vtkSmartPointer<vtkGeometryFilter>::New();

	surfaceFilter->SetInputData(dataSurface);
	surfaceFilter->Update();

	return surfaceFilter->GetOutput();
}

/*! Calc a rotation taranslation matrix.
 \param vector1 Frist Vector (x,y,z).
 \param vector2 Second Vecotr (x,y,z).
 \return A roation taranslation matrix.
 */
vtkSmartPointer<vtkMatrix3x3> Methods::getRotaionTranslationMatrix(vector<double> vector1, vector<double> vector2) {
	vtkSmartPointer<vtkMatrix3x3> matrix1 = vtkSmartPointer<vtkMatrix3x3>::New();
	vtkSmartPointer<vtkMatrix3x3> matrix1inv = vtkSmartPointer<vtkMatrix3x3>::New();
	vtkSmartPointer<vtkMatrix3x3> tempMatrix = vtkSmartPointer<vtkMatrix3x3>::New();

	vtkSmartPointer<vtkMatrix3x3> matrix2 = vtkSmartPointer<vtkMatrix3x3>::New();

	double vector1Array[3] = { vector1.at(0), vector1.at(1), vector1.at(2) };
	double vector2Array[3] = { vector2.at(0), vector2.at(1), vector2.at(2) };

	vtkMath::Normalize(vector1Array);
	vtkMath::Normalize(vector2Array);

	double vector3Array[3];

	vtkMath::Cross(vector1Array, vector2Array, vector3Array);
	vtkMath::Normalize(vector3Array);

	double vector4Array[3];
	vtkMath::Cross(vector3Array, vector1Array, vector4Array);

	matrix1->SetElement(0, 0, vector1Array[0]);
	matrix1->SetElement(0, 1, vector1Array[1]);
	matrix1->SetElement(0, 2, vector1Array[2]);

	matrix1->SetElement(1, 0, vector4Array[0]);
	matrix1->SetElement(1, 1, vector4Array[1]);
	matrix1->SetElement(1, 2, vector4Array[2]);

	matrix1->SetElement(2, 0, vector3Array[0]);
	matrix1->SetElement(2, 1, vector3Array[1]);
	matrix1->SetElement(2, 2, vector3Array[2]);

	double cos = vtkMath::Dot(vector2Array, vector1Array);
	double sin = vtkMath::Dot(vector2Array, vector4Array);

	matrix2->SetElement(0, 0, cos);
	matrix2->SetElement(0, 1, sin);
	matrix2->SetElement(0, 2, 0);

	matrix2->SetElement(1, 0, -sin);
	matrix2->SetElement(1, 1, cos);
	matrix2->SetElement(1, 2, 0);

	matrix2->SetElement(2, 0, 0);
	matrix2->SetElement(2, 1, 0);
	matrix2->SetElement(2, 2, 1);


	matrix1inv->DeepCopy(matrix1);
	matrix1inv->Invert();


	vtkMatrix3x3::Multiply3x3(matrix1inv, matrix2, tempMatrix);
	vtkMatrix3x3::Multiply3x3(tempMatrix, matrix1, tempMatrix);


	return tempMatrix;
} // Methods::getRotaionTranslationMatrix

/*! Delete all cells that are not connected to the cell in which the point lie.
 \param data Pointer with the mesh data.
 \param point Point to that the cells are must connected as vector (x,y,z).
 */
void Methods::cleanCase(DataFormat &data, vector<double> point) {
	vtkIdType pointID = data.getCentrePoints()->FindPoint(point.at(0), point.at(1), point.at(2));

	set<vtkIdType> cellIdList;

	cellIdList.insert(pointID);

	set<vtkIdType> viewedCells;

	bool end = false;

	while (!end) {
		set<vtkIdType> newcellIdList;
		for (set<vtkIdType>::iterator iter = cellIdList.begin(); iter != cellIdList.end(); iter++) {
			vtkIdType currentCellIds = *iter;
			vtkSmartPointer<vtkIdList> cellPointIds = vtkSmartPointer<vtkIdList>::New();

			data.getVtkData()->GetCellPoints(currentCellIds, cellPointIds);

			for (vtkIdType i = 0; i < cellPointIds->GetNumberOfIds(); i++) {
				vtkSmartPointer<vtkIdList> CellPointIdList = vtkSmartPointer<vtkIdList>::New();
				CellPointIdList->InsertNextId(cellPointIds->GetId(i));

				vtkSmartPointer<vtkIdList> tempNeighborCellIds = vtkSmartPointer<vtkIdList>::New();
				data.getVtkData()->GetCellNeighbors(currentCellIds, CellPointIdList, tempNeighborCellIds);

				for (vtkIdType j = 0; j < tempNeighborCellIds->GetNumberOfIds(); j++) {
					if (viewedCells.find(tempNeighborCellIds->GetId(j)) == viewedCells.end()) {
						newcellIdList.insert(tempNeighborCellIds->GetId(j));
						viewedCells.insert(tempNeighborCellIds->GetId(j));
					}
				}
			}
		}
		if (newcellIdList.size() == 0) {
			end = true;
		}
		cellIdList.clear();
		cellIdList = newcellIdList;
	}


	vtkSmartPointer<vtkUnstructuredGrid> dataUGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	vtkSmartPointer<vtkPolyData> dataPoly = vtkSmartPointer<vtkPolyData>::New();


	vtkSmartPointer<vtkDoubleArray> material = vtkSmartPointer<vtkDoubleArray>::New();
	material->SetName("Material");
	material->SetNumberOfComponents(1);


	vtkSmartPointer<vtkCellArray> cellArray = vtkSmartPointer<vtkCellArray>::New();

	for (vtkIdType i = 0; i < data.getVtkData()->GetNumberOfCells(); i++) {
		vtkIdType currentid = i;

		if (viewedCells.find(currentid) != viewedCells.end()) {
			if (data.getInputType() == DataFormat::PossibleInputType::vtu) {
				vtkSmartPointer<vtkIdList> pointIds = vtkSmartPointer<vtkIdList>::New();
				data.getVtkData()->GetCellPoints(i, pointIds);
				vtkIdType cellId = dataUGrid->InsertNextCell(data.getVtkData()->GetCellType(i), pointIds);
				material->InsertComponent(cellId, 0, data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(i, 0));
			}
			else if (data.getInputType() == DataFormat::PossibleInputType::vtp) {
				vtkSmartPointer<vtkIdList> pointIds = vtkSmartPointer<vtkIdList>::New();
				vtkSmartPointer<vtkCell> cell = data.getVtkData()->GetCell(i);
				vtkIdType cellId = cellArray->InsertNextCell(cell);
				material->InsertComponent(cellId, 0, data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(i, 0));
			}
		}
	}

	if (data.getInputType() == DataFormat::PossibleInputType::vtu) {
		dataUGrid->SetPoints(data.getVtkData()->GetPoints());
		dataUGrid->BuildLinks();
		dataUGrid->GetCellData()->AddArray(material);

		data.setVtkData(dataUGrid);
	}
	else if (data.getInputType() == DataFormat::PossibleInputType::vtp) {
		dataPoly->SetPoints(data.getVtkData()->GetPoints());
		dataPoly->SetPolys(cellArray);
		dataPoly->GetCellData()->AddArray(material);

		dataPoly->BuildCells();
		dataPoly->BuildLinks();
		data.setVtkData(dataPoly);
	}

	data.setCentrePoints(calcCellCentroids(data.getVtkData()));
} // Methods::cleanCase

/*! Delete all cells the are not connected in a surface to the cells in which one of the points lie.
 \param data Pointer with the mesh data.
 \param point1 First point to that the cells are must connected as vector (x,y,z).
 \param point2 Second point to that the cells are must connected as vector (x,y,z).
 */
void Methods::cleanCaseVTP(DataFormat &data, vector<double> point1, vector<double> point2) {
	vtkIdType pointID1 = data.getCentrePoints()->FindPoint(point1.at(0), point1.at(1), point1.at(2));

	set<vtkIdType> cellIdList;

	cellIdList.insert(pointID1);

	set<vtkIdType> viewedCells;

	bool end = false;

	while (!end) {
		set<vtkIdType> newcellIdList;
		for (set<vtkIdType>::iterator iter = cellIdList.begin(); iter != cellIdList.end(); iter++) {
			vtkIdType currentCellIds = *iter;
			vtkSmartPointer<vtkIdList> cellPointIds = vtkSmartPointer<vtkIdList>::New();

			data.getVtkData()->GetCellPoints(currentCellIds, cellPointIds);

			for (vtkIdType i = 0; i < cellPointIds->GetNumberOfIds(); i++) {
				vtkSmartPointer<vtkIdList> CellPointIdList = vtkSmartPointer<vtkIdList>::New();
				CellPointIdList->InsertNextId(cellPointIds->GetId(i));

				vtkSmartPointer<vtkIdList> tempNeighborCellIds = vtkSmartPointer<vtkIdList>::New();
				data.getVtkData()->GetCellNeighbors(currentCellIds, CellPointIdList, tempNeighborCellIds);

				for (vtkIdType j = 0; j < tempNeighborCellIds->GetNumberOfIds(); j++) {
					if (viewedCells.find(tempNeighborCellIds->GetId(j)) == viewedCells.end()) {
						newcellIdList.insert(tempNeighborCellIds->GetId(j));
						viewedCells.insert(tempNeighborCellIds->GetId(j));
					}
				}
			}
		}
		if (newcellIdList.size() == 0) {
			end = true;
		}
		cellIdList.clear();
		cellIdList = newcellIdList;
	}

	vtkIdType pointID2 = data.getCentrePoints()->FindPoint(point2.at(0), point2.at(1), point2.at(2));

	cellIdList.insert(pointID2);

	end = false;

	while (!end) {
		set<vtkIdType> newcellIdList;
		for (set<vtkIdType>::iterator iter = cellIdList.begin(); iter != cellIdList.end(); iter++) {
			vtkIdType currentCellIds = *iter;
			vtkSmartPointer<vtkIdList> cellPointIds = vtkSmartPointer<vtkIdList>::New();

			data.getVtkData()->GetCellPoints(currentCellIds, cellPointIds);

			for (vtkIdType i = 0; i < cellPointIds->GetNumberOfIds(); i++) {
				vtkSmartPointer<vtkIdList> CellPointIdList = vtkSmartPointer<vtkIdList>::New();
				CellPointIdList->InsertNextId(cellPointIds->GetId(i));

				vtkSmartPointer<vtkIdList> tempNeighborCellIds = vtkSmartPointer<vtkIdList>::New();
				data.getVtkData()->GetCellNeighbors(currentCellIds, CellPointIdList, tempNeighborCellIds);

				for (vtkIdType j = 0; j < tempNeighborCellIds->GetNumberOfIds(); j++) {
					if (viewedCells.find(tempNeighborCellIds->GetId(j)) == viewedCells.end()) {
						newcellIdList.insert(tempNeighborCellIds->GetId(j));
						viewedCells.insert(tempNeighborCellIds->GetId(j));
					}
				}
			}
		}
		if (newcellIdList.size() == 0) {
			end = true;
		}
		cellIdList.clear();
		cellIdList = newcellIdList;
	}

	vtkSmartPointer<vtkPolyData> dataPoly = vtkSmartPointer<vtkPolyData>::New();

	vtkSmartPointer<vtkDoubleArray> material = vtkSmartPointer<vtkDoubleArray>::New();
	material->SetName("Material");
	material->SetNumberOfComponents(1);


	vtkSmartPointer<vtkCellArray> cellArray = vtkSmartPointer<vtkCellArray>::New();

	for (vtkIdType i = 0; i < data.getVtkData()->GetNumberOfCells(); i++) {
		vtkIdType currentid = i;

		if (viewedCells.find(currentid) != viewedCells.end()) {
			vtkSmartPointer<vtkIdList> pointIds = vtkSmartPointer<vtkIdList>::New();
			vtkSmartPointer<vtkCell> cell = data.getVtkData()->GetCell(i);
			vtkIdType cellId = cellArray->InsertNextCell(cell);
			material->InsertComponent(cellId, 0, data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(i, 0));
		}
	}

	dataPoly->SetPoints(data.getVtkData()->GetPoints());
	dataPoly->SetPolys(cellArray);
	dataPoly->GetCellData()->AddArray(material);

	dataPoly->BuildCells();
	dataPoly->BuildLinks();
	data.setVtkData(dataPoly);


	data.setCentrePoints(calcCellCentroids(data.getVtkData()));
} // Methods::cleanCaseVTP

/*! Clean the mesh from all cells the are not connected to the cells in which one of the points lie.
 \param data Pointer with the mesh data.
 \param point Point to that the cells are must connected as vector (x,y,z).
 */
void Methods::endCleanCase(DataFormat &data, vector<double> point) {
	vtkIdType pointID = data.getCentrePoints()->FindPoint(point.at(0), point.at(1), point.at(2));

	vtkSmartPointer<vtkIdList> cellIDs = vtkSmartPointer<vtkIdList>::New();


	set<vtkIdType> cellIdList;

	cellIdList.insert(pointID);

	set<vtkIdType> viewedCells;

	bool end = false;

	while (!end) {
		set<vtkIdType> newcellIdList;
		for (set<vtkIdType>::iterator iter = cellIdList.begin(); iter != cellIdList.end(); iter++) {
			vtkIdType currentCellIds = *iter;
			vtkSmartPointer<vtkIdList> cellPointIds = vtkSmartPointer<vtkIdList>::New();

			data.getVtkData()->GetCellPoints(currentCellIds, cellPointIds);

			for (vtkIdType i = 0; i < cellPointIds->GetNumberOfIds(); i++) {
				vtkSmartPointer<vtkIdList> CellPointIdList = vtkSmartPointer<vtkIdList>::New();
				CellPointIdList->InsertNextId(cellPointIds->GetId(i));

				vtkSmartPointer<vtkIdList> tempNeighborCellIds = vtkSmartPointer<vtkIdList>::New();
				data.getVtkData()->GetCellNeighbors(currentCellIds, CellPointIdList, tempNeighborCellIds);

				for (vtkIdType j = 0; j < tempNeighborCellIds->GetNumberOfIds(); j++) {
					if (viewedCells.find(tempNeighborCellIds->GetId(j)) == viewedCells.end()) {
						newcellIdList.insert(tempNeighborCellIds->GetId(j));
						viewedCells.insert(tempNeighborCellIds->GetId(j));
						cellIDs->InsertNextId(tempNeighborCellIds->GetId(j));
					}
				}
			}
		}
		if (newcellIdList.size() == 0) {
			end = true;
		}
		cellIdList.clear();
		cellIdList = newcellIdList;
	}


	if (data.getVtkData()->GetNumberOfCells() != cellIDs->GetNumberOfIds()) {
		vtkSmartPointer<vtkDoubleArray> material = vtkSmartPointer<vtkDoubleArray>::New();
		material->SetName("Material");
		material->SetNumberOfComponents(1);
		material->SetNumberOfTuples(cellIDs->GetNumberOfIds());

		vtkSmartPointer<vtkDoubleArray> rawMaterial = vtkSmartPointer<vtkDoubleArray>::New();
		rawMaterial->SetName("RawMaterial");
		rawMaterial->SetNumberOfComponents(1);
		rawMaterial->SetNumberOfTuples(cellIDs->GetNumberOfIds());

		vtkSmartPointer<vtkDoubleArray> theta = vtkSmartPointer<vtkDoubleArray>::New();
		theta->SetName("Theta");
		theta->SetNumberOfComponents(1);
		theta->SetNumberOfTuples(cellIDs->GetNumberOfIds());
		theta->FillComponent(0, 0);

		vtkSmartPointer<vtkDoubleArray> phi = vtkSmartPointer<vtkDoubleArray>::New();
		phi->SetName("Phi");
		phi->SetNumberOfComponents(1);
		phi->SetNumberOfTuples(cellIDs->GetNumberOfIds());
		phi->FillComponent(0, 0);

		vtkSmartPointer<vtkDoubleArray> thetaPhi = vtkSmartPointer<vtkDoubleArray>::New();
		thetaPhi->SetName("ThetaPhi");
		thetaPhi->SetNumberOfComponents(2);
		thetaPhi->SetNumberOfTuples(cellIDs->GetNumberOfIds());
		thetaPhi->SetComponentName(0, "Theta");
		thetaPhi->SetComponentName(1, "Phi");
		thetaPhi->FillComponent(0, 0);
		thetaPhi->FillComponent(1, 0);

		vtkSmartPointer<vtkDoubleArray> diffVec = vtkSmartPointer<vtkDoubleArray>::New();
		diffVec->SetName("DifferenceVector");
		diffVec->SetNumberOfComponents(3);
		diffVec->SetNumberOfTuples(cellIDs->GetNumberOfIds());
		diffVec->FillComponent(0, 0);
		diffVec->FillComponent(1, 0);
		diffVec->FillComponent(2, 0);

		vtkSmartPointer<vtkDoubleArray> baseDiffVec = vtkSmartPointer<vtkDoubleArray>::New();
		baseDiffVec->SetName("BasedPathDifferenceVector");
		baseDiffVec->SetNumberOfComponents(3);
		baseDiffVec->SetNumberOfTuples(cellIDs->GetNumberOfIds());
		baseDiffVec->FillComponent(0, 0);
		baseDiffVec->FillComponent(1, 0);
		baseDiffVec->FillComponent(2, 0);

		vtkSmartPointer<vtkDoubleArray> status = vtkSmartPointer<vtkDoubleArray>::New();
		status->SetName("Status");
		status->SetNumberOfComponents(1);
		status->SetNumberOfTuples(cellIDs->GetNumberOfIds());
		status->FillComponent(0, 0);

		if (data.getInputType() == DataFormat::PossibleInputType::vtu) {
			vtkSmartPointer<vtkUnstructuredGrid> uGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
			uGrid->Allocate(data.getVtkData()->GetNumberOfCells());

			uGrid->GetCellData()->AddArray(material);
			uGrid->GetCellData()->AddArray(rawMaterial);
			uGrid->GetCellData()->AddArray(theta);
			uGrid->GetCellData()->AddArray(phi);
			uGrid->GetCellData()->AddArray(thetaPhi);
			uGrid->GetCellData()->AddArray(diffVec);
			uGrid->GetCellData()->AddArray(baseDiffVec);
			uGrid->GetCellData()->AddArray(status);

			uGrid->SetPoints(data.getVtkData()->GetPoints());

			for (vtkIdType j = 0; j < cellIDs->GetNumberOfIds(); j++) {
				vtkIdType newcellID = -1;
				newcellID = uGrid->InsertNextCell(data.getVtkData()->GetCellType(cellIDs->GetId(j)), data.getVtkData()->GetCell(cellIDs->GetId(j))->GetPointIds());

				uGrid->GetCellData()->GetArray("Status")->InsertTuple(newcellID, data.getVtkData()->GetCellData()->GetArray("Status")->GetTuple(cellIDs->GetId(j)));
				uGrid->GetCellData()->GetArray("Material")->InsertTuple(newcellID, data.getVtkData()->GetCellData()->GetArray("Material")->GetTuple(cellIDs->GetId(j)));
				uGrid->GetCellData()->GetArray("RawMaterial")->InsertTuple(newcellID, data.getVtkData()->GetCellData()->GetArray("RawMaterial")->GetTuple(cellIDs->GetId(j)));
				uGrid->GetCellData()->GetArray("Phi")->InsertTuple(newcellID, data.getVtkData()->GetCellData()->GetArray("Phi")->GetTuple(cellIDs->GetId(j)));
				uGrid->GetCellData()->GetArray("Theta")->InsertTuple(newcellID, data.getVtkData()->GetCellData()->GetArray("Theta")->GetTuple(cellIDs->GetId(j)));
				uGrid->GetCellData()->GetArray("ThetaPhi")->InsertTuple(newcellID, data.getVtkData()->GetCellData()->GetArray("ThetaPhi")->GetTuple(cellIDs->GetId(j)));
				uGrid->GetCellData()->GetArray("DifferenceVector")->InsertTuple(newcellID, data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetTuple(cellIDs->GetId(j)));
				uGrid->GetCellData()->GetArray("BasedPathDifferenceVector")->InsertTuple(newcellID, data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(cellIDs->GetId(j)));
			}

			uGrid->BuildLinks();
			data.setVtkData(uGrid);
			data.setCentrePoints(calcCellCentroids(data.getVtkData()));
		}
		else if (data.getInputType() == DataFormat::PossibleInputType::vtp) {
			vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
			vtkSmartPointer<vtkCellArray> cellArray = vtkSmartPointer<vtkCellArray>::New();
			polydata->Allocate(data.getVtkData()->GetNumberOfCells(), (int)data.getVtkData()->GetNumberOfCells());
			polydata->SetPoints(data.getVtkData()->GetPoints());

			polydata->GetCellData()->AddArray(material);
			polydata->GetCellData()->AddArray(rawMaterial);
			polydata->GetCellData()->AddArray(theta);
			polydata->GetCellData()->AddArray(phi);
			polydata->GetCellData()->AddArray(thetaPhi);
			polydata->GetCellData()->AddArray(diffVec);
			polydata->GetCellData()->AddArray(baseDiffVec);
			polydata->GetCellData()->AddArray(status);

			if (!data.getDoubleLayer()) {
				vtkSmartPointer<vtkDoubleArray> materialEndo = vtkSmartPointer<vtkDoubleArray>::New();
				materialEndo->SetName("MaterialEndo");
				materialEndo->SetNumberOfComponents(1);
				materialEndo->SetNumberOfTuples(cellIDs->GetNumberOfIds());

				vtkSmartPointer<vtkDoubleArray> rawMaterialEndo = vtkSmartPointer<vtkDoubleArray>::New();
				rawMaterialEndo->SetName("RawMaterialEndo");
				rawMaterialEndo->SetNumberOfComponents(1);
				rawMaterialEndo->SetNumberOfTuples(cellIDs->GetNumberOfIds());

				vtkSmartPointer<vtkDoubleArray> thetaEndo = vtkSmartPointer<vtkDoubleArray>::New();
				thetaEndo->SetName("ThetaEndo");
				thetaEndo->SetNumberOfComponents(1);
				thetaEndo->SetNumberOfTuples(cellIDs->GetNumberOfIds());
				thetaEndo->FillComponent(0, 0);

				vtkSmartPointer<vtkDoubleArray> phiEndo = vtkSmartPointer<vtkDoubleArray>::New();
				phiEndo->SetName("PhiEndo");
				phiEndo->SetNumberOfComponents(1);
				phiEndo->SetNumberOfTuples(cellIDs->GetNumberOfIds());
				phiEndo->FillComponent(0, 0);

				vtkSmartPointer<vtkDoubleArray> thetaPhiEndo = vtkSmartPointer<vtkDoubleArray>::New();
				thetaPhiEndo->SetName("ThetaPhiEndo");
				thetaPhiEndo->SetNumberOfComponents(2);
				thetaPhiEndo->SetNumberOfTuples(cellIDs->GetNumberOfIds());
				thetaPhiEndo->SetComponentName(0, "ThetaEndo");
				thetaPhiEndo->SetComponentName(1, "PhiEndo");
				thetaPhiEndo->FillComponent(0, 0);
				thetaPhiEndo->FillComponent(1, 0);

				vtkSmartPointer<vtkDoubleArray> diffVecEndo = vtkSmartPointer<vtkDoubleArray>::New();
				diffVecEndo->SetName("DifferenceVectorEndo");
				diffVecEndo->SetNumberOfComponents(3);
				diffVecEndo->SetNumberOfTuples(cellIDs->GetNumberOfIds());
				diffVecEndo->FillComponent(0, 0);
				diffVecEndo->FillComponent(1, 0);
				diffVecEndo->FillComponent(2, 0);

				vtkSmartPointer<vtkDoubleArray> statusEndo = vtkSmartPointer<vtkDoubleArray>::New();
				statusEndo->SetName("StatusEndo");
				statusEndo->SetNumberOfComponents(1);
				statusEndo->SetNumberOfTuples(cellIDs->GetNumberOfIds());
				statusEndo->FillComponent(0, 0);

				polydata->GetCellData()->AddArray(materialEndo);
				polydata->GetCellData()->AddArray(rawMaterialEndo);
				polydata->GetCellData()->AddArray(thetaEndo);
				polydata->GetCellData()->AddArray(phiEndo);
				polydata->GetCellData()->AddArray(thetaPhiEndo);
				polydata->GetCellData()->AddArray(diffVecEndo);
				polydata->GetCellData()->AddArray(statusEndo);
			}


			for (int k = 0; k < cellIDs->GetNumberOfIds(); k++) {
				vtkIdType newcellID = -1;
				newcellID = cellArray->InsertNextCell(data.getVtkData()->GetCell(cellIDs->GetId(k)));

				polydata->GetCellData()->GetArray("Status")->InsertTuple(newcellID, data.getVtkData()->GetCellData()->GetArray("Status")->GetTuple(cellIDs->GetId(k)));
				polydata->GetCellData()->GetArray("Material")->InsertTuple(newcellID, data.getVtkData()->GetCellData()->GetArray("Material")->GetTuple(cellIDs->GetId(k)));
				polydata->GetCellData()->GetArray("RawMaterial")->InsertTuple(newcellID, data.getVtkData()->GetCellData()->GetArray("RawMaterial")->GetTuple(cellIDs->GetId(k)));
				polydata->GetCellData()->GetArray("Phi")->InsertTuple(newcellID, data.getVtkData()->GetCellData()->GetArray("Phi")->GetTuple(cellIDs->GetId(k)));
				polydata->GetCellData()->GetArray("Theta")->InsertTuple(newcellID, data.getVtkData()->GetCellData()->GetArray("Theta")->GetTuple(cellIDs->GetId(k)));
				polydata->GetCellData()->GetArray("ThetaPhi")->InsertTuple(newcellID, data.getVtkData()->GetCellData()->GetArray("ThetaPhi")->GetTuple(cellIDs->GetId(k)));
				polydata->GetCellData()->GetArray("DifferenceVector")->InsertTuple(newcellID, data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetTuple(cellIDs->
					GetId(k)));
				polydata->GetCellData()->GetArray("BasedPathDifferenceVector")->InsertTuple(newcellID, data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector")->GetTuple(cellIDs->GetId(k)));

				if (!data.getDoubleLayer()) {
					polydata->GetCellData()->GetArray("StatusEndo")->InsertTuple(newcellID, data.getVtkData()->GetCellData()->GetArray("StatusEndo")->GetTuple(cellIDs->GetId(k)));
					polydata->GetCellData()->GetArray("MaterialEndo")->InsertTuple(newcellID, data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetTuple(cellIDs->GetId(k)));
					polydata->GetCellData()->GetArray("RawMaterialEndo")->InsertTuple(newcellID, data.getVtkData()->GetCellData()->GetArray("RawMaterialEndo")->GetTuple(cellIDs->
						GetId(k)));
					polydata->GetCellData()->GetArray("PhiEndo")->InsertTuple(newcellID, data.getVtkData()->GetCellData()->GetArray("PhiEndo")->GetTuple(cellIDs->GetId(k)));
					polydata->GetCellData()->GetArray("ThetaEndo")->InsertTuple(newcellID, data.getVtkData()->GetCellData()->GetArray("ThetaEndo")->GetTuple(cellIDs->GetId(k)));
					polydata->GetCellData()->GetArray("ThetaPhiEndo")->InsertTuple(newcellID, data.getVtkData()->GetCellData()->GetArray("ThetaPhiEndo")->GetTuple(cellIDs->GetId(k)));
					polydata->GetCellData()->GetArray("DifferenceVectorEndo")->InsertTuple(newcellID, data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->GetTuple(cellIDs->GetId(k)));
				}
			}

			polydata->SetPolys(cellArray);
			polydata->BuildCells();
			polydata->BuildLinks();

			vtkSmartPointer<vtkCleanPolyData> cleanFilter = vtkSmartPointer<vtkCleanPolyData>::New();
			cleanFilter->SetInputData(polydata);
			cleanFilter->Update();
			polydata = cleanFilter->GetOutput();

			data.setVtkData(polydata);
			data.setCentrePoints(calcCellCentroids(data.getVtkData()));
		}
	}
} // Methods::endCleanCase

/*! Construct a pill between two points with a define resolution and radius.
 \param point1 First point as vector (x,y,z).
 \param point2 Second point as vector (x,y,z).
 \param resolution Resolution of the mesh.
 \param radius Radius of the pill.
 \return The construct surface pill between the points.
 */

DataFormat Methods::getPill(vector<double> point1, vector<double> point2, double resolution, double radius) {
	double point1Array[3] = { point1.at(0), point1.at(1), point1.at(2) };
	double point2Array[3] = { point2.at(0), point2.at(1), point2.at(2) };
	double distancePoint1Point2 = sqrt(abs(vtkMath::Distance2BetweenPoints(point1Array, point2Array)));
	int resolution2 = ceil(2 * radius*M_PI / 3 / resolution);

	if (distancePoint1Point2 > 0) {
		vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
		sphereSource->SetCenter(0, 0, 0);
		sphereSource->SetRadius(radius);
		sphereSource->SetThetaResolution(resolution2);
		sphereSource->SetPhiResolution(resolution2);
		sphereSource->SetStartPhi(90);
		sphereSource->Update();

		vtkSmartPointer<vtkSphereSource> sphereSource2 = vtkSmartPointer<vtkSphereSource>::New();
		sphereSource2->SetCenter(0, 0, distancePoint1Point2);
		sphereSource2->SetRadius(radius);
		sphereSource2->SetThetaResolution(resolution2);
		sphereSource2->SetPhiResolution(resolution2);
		sphereSource2->SetEndPhi(90);
		sphereSource2->Update();

		vtkSmartPointer<vtkPolyData> hemisphere1 = sphereSource->GetOutput();
		vtkSmartPointer<vtkPolyData> hemisphere2 = sphereSource2->GetOutput();

		vtkSmartPointer<vtkAppendPolyData> hemispheres = vtkSmartPointer<vtkAppendPolyData>::New();
		hemispheres->AddInputData(hemisphere1);
		hemispheres->AddInputData(hemisphere2);
		hemispheres->Update();

		vtkSmartPointer<vtkCleanPolyData> cleanFilter = vtkSmartPointer<vtkCleanPolyData>::New();
		cleanFilter->SetInputData(hemispheres->GetOutput());
		cleanFilter->Update();

		vtkSmartPointer<vtkPolyData> hemispheresurfaces = cleanFilter->GetOutput();

		if (distancePoint1Point2 > resolution) {
			vector<double> e3 = { 0, 0, 1 };
			vector<double> e2 = { 0, 1, 0 };

			vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
			for (int i = 0; i*resolution < distancePoint1Point2; i += 2) {
				points->InsertNextPoint(0, i*resolution, 0);
			}
			points->InsertNextPoint(0, distancePoint1Point2, 0);

			vtkSmartPointer<vtkLineSource> lineSource = vtkSmartPointer<vtkLineSource>::New();
			lineSource->SetPoints(points);
			lineSource->Update();

			vtkSmartPointer<vtkTubeFilter> cylinder = vtkSmartPointer<vtkTubeFilter>::New();
			cylinder->SetInputConnection(lineSource->GetOutputPort());
			cylinder->SetRadius(radius);
			cylinder->SetNumberOfSides(resolution2);
			cylinder->Update();


			vtkSmartPointer<vtkTriangleFilter> triangleFilter = vtkSmartPointer<vtkTriangleFilter>::New();
			triangleFilter->SetInputConnection(cylinder->GetOutputPort());
			triangleFilter->Update();

			vtkSmartPointer<vtkPolyData> cylinderSurface = triangleFilter->GetOutput();
			cylinderSurface->SetPoints(rotatePoints(cylinderSurface->GetPoints(), getRotaionTranslationMatrix(e3, e2)));
			cylinderSurface->BuildLinks();

			vtkSmartPointer<vtkAppendPolyData> cylinderBridge = vtkSmartPointer<vtkAppendPolyData>::New();
			cylinderBridge->AddInputData(hemispheresurfaces);
			cylinderBridge->AddInputData(cylinderSurface);
			cylinderBridge->Update();

			vtkSmartPointer<vtkCleanPolyData> cleanFilterCylinderBridge = vtkSmartPointer<vtkCleanPolyData>::New();
			cleanFilterCylinderBridge->SetInputData(cylinderBridge->GetOutput());
			cleanFilterCylinderBridge->Update();

			vtkSmartPointer<vtkPolyData> pill = cleanFilterCylinderBridge->GetOutput();

			for (vtkIdType i = 0; i < cylinderSurface->GetNumberOfPoints(); i++) {
				double point[3] = { cylinderSurface->GetPoint(i)[0], cylinderSurface->GetPoint(i)[1], cylinderSurface->GetPoint(i)[2] };

				double mindistance = std::numeric_limits<double>::max();
				double newpoint[3];

				for (vtkIdType j = 0; j < hemispheresurfaces->GetNumberOfPoints(); j++) {
					double distance = sqrt(abs(vtkMath::Distance2BetweenPoints(point, hemispheresurfaces->GetPoint(j))));
					if (distance < mindistance) {
						mindistance = distance;
						newpoint[0] = hemispheresurfaces->GetPoint(j)[0];
						newpoint[1] = hemispheresurfaces->GetPoint(j)[1];
						newpoint[2] = hemispheresurfaces->GetPoint(j)[2];
					}
				}
				if (mindistance <= 1.51 * resolution) {
					vtkIdType oldpointId = pill->FindPoint(point);
					vtkIdType newpointId = pill->FindPoint(newpoint);

					vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
					pill->GetPointCells(oldpointId, cellIds);

					for (vtkIdType j = 0; j < cellIds->GetNumberOfIds(); j++) {
						pill->ReplaceCellPoint(cellIds->GetId(j), oldpointId, newpointId);
					}
					pill->BuildCells();
					pill->BuildLinks();
				}
			}


			double pillNormal[3] = { 0, 0, 0 };
			vtkMath::Subtract(point2Array, point1Array, pillNormal);
			vtkMath::Normalize(pillNormal);
			vector<double> pillNormalVector = { pillNormal[0], pillNormal[1], pillNormal[2] };

			pill->SetPoints(rotatePoints(pill->GetPoints(), getRotaionTranslationMatrix(pillNormalVector, e3)));
			pill->SetPoints(movePoints(pill->GetPoints(), point1));
			pill->BuildLinks();


			for (vtkIdType i = 0; i < pill->GetNumberOfCells(); i++) {
				vtkSmartPointer<vtkIdList> pointIds = pill->GetCell(i)->GetPointIds();
				if ((pointIds->GetId(0) == pointIds->GetId(1)) ||
					(pointIds->GetId(0) == pointIds->GetId(2)) ||
					(pointIds->GetId(1) == pointIds->GetId(2))) {
					pill->DeleteCell(i);
				}
			}
			pill->RemoveDeletedCells();
			pill->BuildLinks();

			vtkSmartPointer<vtkCleanPolyData> cleanFilter3 = vtkSmartPointer<vtkCleanPolyData>::New();
			cleanFilter3->SetInputData(pill);
			cleanFilter3->Update();

			pill = cleanFilter3->GetOutput();

			DataFormat pillData;
			pillData.setInputType(DataFormat::vtp);
			pillData.setVtkData(pill);
			pillData.setCentrePoints(Methods::calcCellCentroids(pill));

			return pillData;
		}
		else {
			vector<double> e3 = { 0, 0, 1 };

			for (vtkIdType i = 0; i < hemisphere1->GetNumberOfPoints(); i++) {
				double point[3] = { hemisphere1->GetPoint(i)[0], hemisphere1->GetPoint(i)[1], hemisphere1->GetPoint(i)[2] };
				vtkIdType pointid = hemisphere2->FindPoint(point);
				double distance = sqrt(abs(vtkMath::Distance2BetweenPoints(point, hemisphere2->GetPoint(pointid))));

				if (distance <= 1.5 * distancePoint1Point2) {
					vtkIdType oldpointId = hemispheresurfaces->FindPoint(point);
					vtkIdType newpointId = hemispheresurfaces->FindPoint(hemisphere2->GetPoint(pointid));

					vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
					hemispheresurfaces->GetPointCells(oldpointId, cellIds);
					for (vtkIdType j = 0; j < cellIds->GetNumberOfIds(); j++) {
						hemispheresurfaces->ReplaceCellPoint(cellIds->GetId(j), oldpointId, newpointId);
					}
					hemispheresurfaces->BuildLinks();
				}
			}

			double pillNormal[3] = { 0, 0, 0 };
			vtkMath::Subtract(point2Array, point1Array, pillNormal);
			vtkMath::Normalize(pillNormal);
			vector<double> pillNormalVector = { pillNormal[0], pillNormal[1], pillNormal[2] };

			hemispheresurfaces->SetPoints(rotatePoints(hemispheresurfaces->GetPoints(), getRotaionTranslationMatrix(pillNormalVector, e3)));
			hemispheresurfaces->SetPoints(movePoints(hemispheresurfaces->GetPoints(), point1));
			hemispheresurfaces->BuildLinks();


			for (vtkIdType i = 0; i < hemispheresurfaces->GetNumberOfCells(); i++) {
				vtkSmartPointer<vtkIdList> pointIds = hemispheresurfaces->GetCell(i)->GetPointIds();
				if ((pointIds->GetId(0) == pointIds->GetId(1)) ||
					(pointIds->GetId(0) == pointIds->GetId(2)) ||
					(pointIds->GetId(1) == pointIds->GetId(2))) {
					hemispheresurfaces->DeleteCell(i);
				}
			}
			hemispheresurfaces->RemoveDeletedCells();
			hemispheresurfaces->BuildLinks();

			vtkSmartPointer<vtkCleanPolyData> cleanFilter3 = vtkSmartPointer<vtkCleanPolyData>::New();
			cleanFilter3->SetInputData(hemispheresurfaces);
			cleanFilter3->Update();

			hemispheresurfaces = cleanFilter3->GetOutput();


			DataFormat pillData;
			pillData.setInputType(DataFormat::vtp);
			pillData.setVtkData(hemispheresurfaces);
			pillData.setCentrePoints(Methods::calcCellCentroids(hemispheresurfaces));

			return pillData;
		}
	}
	else {
		vtkSmartPointer<vtkSphereSource> sphere = vtkSmartPointer<vtkSphereSource>::New();
		sphere->SetCenter(0, 0, 0);
		sphere->SetRadius(radius);
		sphere->SetThetaResolution(resolution2);
		sphere->SetPhiResolution(resolution2);
		sphere->Update();

		vtkSmartPointer<vtkPolyData> sphereData = sphere->GetOutput();
		sphereData->SetPoints(movePoints(sphereData->GetPoints(), point1));
		sphereData->BuildLinks();

		DataFormat pillData;
		pillData.setInputType(DataFormat::vtp);
		pillData.setVtkData(sphereData);
		pillData.setCentrePoints(Methods::calcCellCentroids(sphereData));

		return pillData;
	}
} // Methods::getPill

/*! Union two surface mesh together.
 \param input1 First mesh
 \param input2 Second mesh
 \return The union surface.
 */
vtkSmartPointer<vtkPolyData> Methods::vtpUnion(vtkSmartPointer<vtkPolyData> input1, vtkSmartPointer<vtkPolyData> input2) {
	CorkTriMesh in1, in2, Un;

	in1.n_triangles = (int)input1->GetNumberOfCells();
	in1.n_vertices = (int)input1->GetNumberOfPoints();
	in1.triangles = new uint[(in1.n_triangles) * 3];
	in1.vertices = new float[(in1.n_vertices) * 3];

	for (uint i = 0; i < in1.n_triangles; i++) {
		(in1.triangles)[3 * i + 0] = (int)input1->GetCell(i)->GetPointId(0);
		(in1.triangles)[3 * i + 1] = (int)input1->GetCell(i)->GetPointId(1);
		(in1.triangles)[3 * i + 2] = (int)input1->GetCell(i)->GetPointId(2);
	}

	for (uint i = 0; i < in1.n_vertices; i++) {
		(in1.vertices)[3 * i + 0] = input1->GetPoint(i)[0];
		(in1.vertices)[3 * i + 1] = input1->GetPoint(i)[1];
		(in1.vertices)[3 * i + 2] = input1->GetPoint(i)[2];
	}

	in2.n_triangles = (int)input2->GetNumberOfCells();
	in2.n_vertices = (int)input2->GetNumberOfPoints();
	in2.triangles = new uint[(in2.n_triangles) * 3];
	in2.vertices = new float[(in2.n_vertices) * 3];

	for (uint i = 0; i < in2.n_triangles; i++) {
		(in2.triangles)[3 * i + 0] = (int)input2->GetCell(i)->GetPointId(0);
		(in2.triangles)[3 * i + 1] = (int)input2->GetCell(i)->GetPointId(1);
		(in2.triangles)[3 * i + 2] = (int)input2->GetCell(i)->GetPointId(2);
	}

	for (uint i = 0; i < in2.n_vertices; i++) {
		(in2.vertices)[3 * i + 0] = input2->GetPoint(i)[0];
		(in2.vertices)[3 * i + 1] = input2->GetPoint(i)[1];
		(in2.vertices)[3 * i + 2] = input2->GetPoint(i)[2];
	}
	computeUnion(in1, in2, &Un);

	// vtkSmartPointer<vtkBooleanOperationPolyDataFilter> booleanOperation = // vtkSmartPointer<vtkBooleanOperationPolyDataFilter>::New();
	// booleanOperation->SetOperationToUnion();
	// booleanOperation->SetInputData(0, input1);
	// booleanOperation->SetInputData(1, input2);
	// booleanOperation->SetTolerance(0);
	// booleanOperation->Update();

	vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
	polydata->Allocate((int)Un.n_triangles, (int)Un.n_triangles);
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

	for (int i = 0; i < Un.n_vertices; i++) {
		points->InsertNextPoint(Un.vertices[i * 3], Un.vertices[i * 3 + 1], Un.vertices[i * 3 + 2]);
	}

	polydata->SetPoints(points);

	for (int j = 0; j < Un.n_triangles; j++) {
		vtkSmartPointer<vtkIdList> pointsList = vtkSmartPointer<vtkIdList>::New();
		for (int k = 0; k < 3; k++) {
			int id = Un.triangles[j * 3 + k];
			pointsList->InsertNextId(id);
		}
		polydata->InsertNextCell(5, pointsList);
	}

	vtkSmartPointer<vtkTriangleFilter> triangleFilter = vtkSmartPointer<vtkTriangleFilter>::New();
	triangleFilter->SetInputData(polydata);
	triangleFilter->Update();
	vtkSmartPointer<vtkCleanPolyData> cleanFilter = vtkSmartPointer<vtkCleanPolyData>::New();
	cleanFilter->SetInputData(triangleFilter->GetOutput());
	cleanFilter->Update();
	vtkSmartPointer<vtkPolyData> output = cleanFilter->GetOutput();
	output->BuildCells();
	output->BuildLinks();
	return output;
} // Methods::vtpUnion

/*! Calc the different between two surface meshs.
 \param input1 First mesh
 \param input2 Second mesh
 \return The difference of two surface meshs.
 */
vtkSmartPointer<vtkPolyData> Methods::vtpDifference(vtkSmartPointer<vtkPolyData> input1, vtkSmartPointer<vtkPolyData> input2) {
	CorkTriMesh in1, in2, Un;

	in1.n_triangles = (int)input1->GetNumberOfCells();
	in1.n_vertices = (int)input1->GetNumberOfPoints();
	in1.triangles = new uint[(in1.n_triangles) * 3];
	in1.vertices = new float[(in1.n_vertices) * 3];

	for (uint i = 0; i < in1.n_triangles; i++) {
		(in1.triangles)[3 * i + 0] = (int)input1->GetCell(i)->GetPointId(0);
		(in1.triangles)[3 * i + 1] = (int)input1->GetCell(i)->GetPointId(1);
		(in1.triangles)[3 * i + 2] = (int)input1->GetCell(i)->GetPointId(2);
	}

	for (uint i = 0; i < in1.n_vertices; i++) {
		(in1.vertices)[3 * i + 0] = input1->GetPoint(i)[0];
		(in1.vertices)[3 * i + 1] = input1->GetPoint(i)[1];
		(in1.vertices)[3 * i + 2] = input1->GetPoint(i)[2];
	}

	in2.n_triangles = (int)input2->GetNumberOfCells();
	in2.n_vertices = (int)input2->GetNumberOfPoints();
	in2.triangles = new uint[(in2.n_triangles) * 3];
	in2.vertices = new float[(in2.n_vertices) * 3];

	for (uint i = 0; i < in2.n_triangles; i++) {
		(in2.triangles)[3 * i + 0] = (int)input2->GetCell(i)->GetPointId(0);
		(in2.triangles)[3 * i + 1] = (int)input2->GetCell(i)->GetPointId(1);
		(in2.triangles)[3 * i + 2] = (int)input2->GetCell(i)->GetPointId(2);
	}

	for (uint i = 0; i < in2.n_vertices; i++) {
		(in2.vertices)[3 * i + 0] = input2->GetPoint(i)[0];
		(in2.vertices)[3 * i + 1] = input2->GetPoint(i)[1];
		(in2.vertices)[3 * i + 2] = input2->GetPoint(i)[2];
	}

	computeDifference(in1, in2, &Un);

	// vtkSmartPointer<vtkBooleanOperationPolyDataFilter> booleanOperation = // vtkSmartPointer<vtkBooleanOperationPolyDataFilter>::New();
	// booleanOperation->SetOperationToDifference();
	// booleanOperation->SetInputData(0, input1);
	// booleanOperation->SetInputData(1, input2);
	// booleanOperation->SetTolerance(0);
	// booleanOperation->Update();

	vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
	polydata->Allocate((int)Un.n_triangles, (int)Un.n_triangles);
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

	for (int i = 0; i < Un.n_vertices; i++) {
		points->InsertNextPoint(Un.vertices[i * 3], Un.vertices[i * 3 + 1], Un.vertices[i * 3 + 2]);
	}

	polydata->SetPoints(points);

	for (int j = 0; j < Un.n_triangles; j++) {
		vtkSmartPointer<vtkIdList> pointsList = vtkSmartPointer<vtkIdList>::New();
		for (int k = 0; k < 3; k++) {
			int id = Un.triangles[(j * 3) + k];
			pointsList->InsertNextId(id);
		}
		polydata->InsertNextCell(5, pointsList);
	}

	vtkSmartPointer<vtkTriangleFilter> triangleFilter = vtkSmartPointer<vtkTriangleFilter>::New();
	triangleFilter->SetInputData(polydata);

	// triangleFilter->SetInputConnection(booleanOperation->GetOutputPort());
	triangleFilter->Update();

	vtkSmartPointer<vtkCleanPolyData> cleanFilter = vtkSmartPointer<vtkCleanPolyData>::New();
	cleanFilter->SetInputData(triangleFilter->GetOutput());
	cleanFilter->Update();

	vtkSmartPointer<vtkPolyData> output = cleanFilter->GetOutput();

	// vtkSmartPointer<vtkPolyData> output = triangleFilter->GetOutput();
	output->BuildCells();
	output->BuildLinks();

	return output;
} // Methods::vtpDifference

/*! Set the fiber orientation of cells a surface bridge between the two atrials.
 \param data Pointer with the mesh data.
 \param path The path between the two points.
 \param inMaterial Tissue class in that the fiber would be set.
 \param point1 First point of the bridge as vector (x,y,z).
 \param point2 Second point of the bridge as vector (x,y,z).
 */

void Methods::vtpBridgeMarker(DataFormat &data, vtkSmartPointer<vtkIdList> path, Material::Mat inMaterial, vector<double> point1, vector<double> point2) {
	double point1Array[3] = { point1.at(0), point1.at(1), point1.at(2) };
	double point2Array[3] = { point2.at(0), point2.at(1), point2.at(2) };
	double diffVector[3] = { 0, 0, 0 };

	vtkMath::Subtract(point1Array, point2Array, diffVector);
	vtkMath::Normalize(diffVector);

	for (vtkIdType i = 0; i < path->GetNumberOfIds(); i++) {
		if (!data.getDoubleLayer()) {
			int material = data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(path->GetId(i), 0);
			if ((data.getVtkData()->GetCellData()->GetArray("StatusEndo")->GetComponent(path->GetId(i), 0) == 0) &&
				(material == inMaterial)) {
				data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->InsertTuple(path->GetId(i), diffVector);
				data.getVtkData()->GetCellData()->GetArray("StatusEndo")->SetComponent(path->GetId(i), 0, 3);
			}
		}

		int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(path->GetId(i), 0);
		if ((data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(path->GetId(i), 0) == 0) && (material == inMaterial)) {
			data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->InsertTuple(path->GetId(i), diffVector);
			data.getVtkData()->GetCellData()->GetArray("Status")->SetComponent(path->GetId(i), 0, 3);
		}
	}
} // Methods::vtpBridgeMarker

/*! Construct a pipe around a path with a define radius and resolution.
 \param data Pointer with the mesh data.
 \param path The path between around the pipe will be construct.
 \param resolution Resolution of the mesh.
 \param radius Radius of the construct pipe.
 \return A pipe around a path.
 */

DataFormat Methods::getPillAroundPath(DataFormat &data, vtkSmartPointer<vtkIdList> path, double resolution, double radius) {
	double point1Array[3] = { data.getCentrePoints()->GetPoint(path->GetId(0))[0], data.getCentrePoints()->GetPoint(path->GetId(0))[1], data.getCentrePoints()->GetPoint(path->GetId(0))[2] };

	vector<double> point1Vector = { point1Array[0], point1Array[1], point1Array[2] };

	double point2Array[3] = { data.getCentrePoints()->GetPoint(path->GetId(path->GetNumberOfIds() - 1))[0], data.getCentrePoints()->GetPoint(path->GetId(path->GetNumberOfIds() - 1))[1], data.getCentrePoints()->GetPoint(path->GetId(path->GetNumberOfIds() - 1))[2] };

	vector<double> point2Vector = { point2Array[0], point2Array[1], point2Array[2] };

	int resolution2 = ceil(2 * radius*M_PI / 3 / resolution);

	if (path->GetNumberOfIds() < 1) {
		return getPill(point1Vector, point2Vector, radius, resolution);
	}

	vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
	sphereSource->SetCenter(0, 0, 0);
	sphereSource->SetRadius(radius);
	sphereSource->SetThetaResolution(resolution2);
	sphereSource->SetPhiResolution(resolution2 / 3);
	sphereSource->SetStartPhi(90);
	sphereSource->Update();

	vtkSmartPointer<vtkSphereSource> sphereSource2 = vtkSmartPointer<vtkSphereSource>::New();
	sphereSource2->SetCenter(0, 0, 0);
	sphereSource2->SetRadius(radius);
	sphereSource2->SetThetaResolution(resolution2);
	sphereSource2->SetPhiResolution(resolution2 / 3);
	sphereSource2->SetEndPhi(90);
	sphereSource2->Update();

	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

	for (vtkIdType i = 0; i < path->GetNumberOfIds(); i += 1) {
		double point[3] = { data.getCentrePoints()->GetPoint(path->GetId(i))[0], data.getCentrePoints()->GetPoint(path->GetId(i))[1], data.getCentrePoints()->GetPoint(path->GetId(i))[2] };
		points->InsertNextPoint(point);
	}

	vtkSmartPointer<vtkParametricSpline> spline = vtkSmartPointer<vtkParametricSpline>::New();
	spline->SetPoints(points);

	vtkSmartPointer<vtkParametricFunctionSource> functionSource2 = vtkSmartPointer<vtkParametricFunctionSource>::New();
	functionSource2->SetParametricFunction(spline);
	functionSource2->SetUResolution(ceil((int)points->GetNumberOfPoints() / 2)); // changed from original because of intersection of surfaces, it was /2 before, it was 2* for coarse mesh
	functionSource2->Update();

	// Create the tube
	vtkSmartPointer<vtkTubeFilter> tuber = vtkSmartPointer<vtkTubeFilter>::New();
	tuber->SetInputData(functionSource2->GetOutput());
	tuber->SetNumberOfSides(resolution2);
	tuber->SetRadius(radius);
	tuber->SidesShareVerticesOn();
	tuber->Update();

	vtkSmartPointer<vtkTriangleFilter> triangleFilter = vtkSmartPointer<vtkTriangleFilter>::New();
	triangleFilter->SetInputData(tuber->GetOutput());
	triangleFilter->Update();

	vtkSmartPointer<vtkPolyData> tube = triangleFilter->GetOutput();

	vtkSmartPointer<vtkPolyData> hemisphere1 = sphereSource->GetOutput();

	double sphere1Normal[3] = { 0.0, 0.0, 0.0 };

	double vec1[3] = { tube->GetPointData()->GetArray("TubeNormals")->GetComponent(0, 0), tube->GetPointData()->GetArray("TubeNormals")->GetComponent(0, 1), tube->GetPointData()->GetArray("TubeNormals")->GetComponent(0, 2) };

	double vec2[3] = { tube->GetPointData()->GetArray("TubeNormals")->GetComponent(1, 0), tube->GetPointData()->GetArray("TubeNormals")->GetComponent(1, 1), tube->GetPointData()->GetArray("TubeNormals")->GetComponent(1, 2) };

	vtkMath::Cross(vec2, vec1, sphere1Normal);
	vtkMath::Normalize(sphere1Normal);

	vector<double> sphere1NormalVector = { sphere1Normal[0], sphere1Normal[1], sphere1Normal[2] };
	vector<double> e3 = { 0, 0, 1 };

	hemisphere1->SetPoints(rotatePoints(hemisphere1->GetPoints(), getRotaionTranslationMatrix(sphere1NormalVector, e3)));
	movePoints(hemisphere1->GetPoints(), point1Vector);
	hemisphere1->BuildLinks();

	vtkSmartPointer<vtkPolyData> hemisphere2 = sphereSource2->GetOutput();

	double vec3[3] = { tube->GetPointData()->GetArray("TubeNormals")->GetComponent(tube->GetPointData()->GetNumberOfTuples() - 2, 0), tube->GetPointData()->GetArray("TubeNormals")->GetComponent(tube->GetPointData()->GetNumberOfTuples() - 2, 1), tube->GetPointData()->GetArray("TubeNormals")->GetComponent(tube->GetPointData()->GetNumberOfTuples() - 2, 2) };

	double vec4[3] = { tube->GetPointData()->GetArray("TubeNormals")->GetComponent(tube->GetPointData()->GetNumberOfTuples() - 1, 0), tube->GetPointData()->GetArray("TubeNormals")->GetComponent(tube->GetPointData()->GetNumberOfTuples() - 1, 1), tube->GetPointData()->GetArray("TubeNormals")->GetComponent(tube->GetPointData()->GetNumberOfTuples() - 1, 2) };

	double sphere2Normal[3] = { 0, 0, 0 };
	vtkMath::Cross(vec4, vec3, sphere2Normal);
	vtkMath::Normalize(sphere2Normal);
	vector<double> sphere2NormalVector = { sphere2Normal[0], sphere2Normal[1], sphere2Normal[2] };

	hemisphere2->SetPoints(rotatePoints(hemisphere2->GetPoints(), getRotaionTranslationMatrix(sphere2NormalVector, e3)));
	movePoints(hemisphere2->GetPoints(), point2Vector);
	hemisphere2->BuildLinks();

	vtkSmartPointer<vtkAppendPolyData> hemispheres = vtkSmartPointer<vtkAppendPolyData>::New();
	hemispheres->AddInputData(hemisphere1);
	hemispheres->AddInputData(hemisphere2);
	hemispheres->Update();

	vtkSmartPointer<vtkCleanPolyData> cleanFilter = vtkSmartPointer<vtkCleanPolyData>::New();
	cleanFilter->SetInputData(hemispheres->GetOutput());
	cleanFilter->Update();

	vtkSmartPointer<vtkPolyData> hemispheresurfaces = cleanFilter->GetOutput();

	tube = tubeIntersectionFilter(tube, functionSource2->GetOutput(), resolution2);

	vtkSmartPointer<vtkAppendPolyData> cylinderBridge = vtkSmartPointer<vtkAppendPolyData>::New();
	cylinderBridge->AddInputData(hemispheresurfaces);
	cylinderBridge->AddInputData(tube);
	cylinderBridge->Update();

	vtkSmartPointer<vtkCleanPolyData> cleanFilter2 = vtkSmartPointer<vtkCleanPolyData>::New();
	cleanFilter2->SetInputData(cylinderBridge->GetOutput());
	cleanFilter2->Update();

	vtkSmartPointer<vtkPolyData> pill = cleanFilter2->GetOutput();

	for (vtkIdType i = 0; i < tube->GetNumberOfPoints(); i++) {
		double point[3] = { tube->GetPoint(i)[0], tube->GetPoint(i)[1], tube->GetPoint(i)[2] };

		double mindistance = std::numeric_limits<double>::max();
		double newpoint[3];

		for (vtkIdType j = 0; j < hemispheresurfaces->GetNumberOfPoints(); j++) {
			double distance = sqrt(abs(vtkMath::Distance2BetweenPoints(point, hemispheresurfaces->GetPoint(j))));
			if (distance < mindistance) {
				mindistance = distance;
				newpoint[0] = hemispheresurfaces->GetPoint(j)[0];
				newpoint[1] = hemispheresurfaces->GetPoint(j)[1];
				newpoint[2] = hemispheresurfaces->GetPoint(j)[2];
			}
		}
		if (mindistance <= 1.51 * resolution) {
			vtkIdType oldpointId = pill->FindPoint(point);
			vtkIdType newpointId = pill->FindPoint(newpoint);


			vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
			pill->GetPointCells(oldpointId, cellIds);

			for (vtkIdType j = 0; j < cellIds->GetNumberOfIds(); j++) {
				pill->ReplaceCellPoint(cellIds->GetId(j), oldpointId, newpointId);
			}
			pill->BuildCells();
			pill->BuildLinks();
		}
	}

	for (vtkIdType i = 0; i < pill->GetNumberOfCells(); i++) {
		vtkSmartPointer<vtkIdList> pointIds = pill->GetCell(i)->GetPointIds();
		if ((pointIds->GetId(0) == pointIds->GetId(1)) ||
			(pointIds->GetId(0) == pointIds->GetId(2)) ||
			(pointIds->GetId(1) == pointIds->GetId(2))) {
			pill->DeleteCell(i);
		}
	}
	pill->RemoveDeletedCells();
	pill->BuildLinks();

	vtkSmartPointer<vtkCleanPolyData> cleanFilter3 = vtkSmartPointer<vtkCleanPolyData>::New();
	cleanFilter3->SetInputData(pill);
	cleanFilter3->Update();

	pill = cleanFilter3->GetOutput();

	DataFormat pillData;
	pillData.setInputType(DataFormat::vtp);
	pillData.setVtkData(pill);
	pillData.setCentrePoints(Methods::calcCellCentroids(pill));

	return pillData;
} // Methods::getPillAroundPath

/*! Set the fiber orientation of a pipe around a path on a surface.
 \param data Pointer with the mesh data.
 \param oldData Pointer with the mesh data without the insert pipe.
 \param path The path between around the pipe will be construct.
 \param oldPath The path in the old data.
 */
void Methods::vtpBridgeMarkerArroundPath(DataFormat &data, DataFormat &oldData, vtkSmartPointer<vtkIdList> path, vtkSmartPointer<vtkIdList> oldPath) {
	vtkSmartPointer<vtkUnstructuredGrid> pathSearchList = vtkSmartPointer<vtkUnstructuredGrid>::New();
	vtkSmartPointer<vtkPoints> pointsPath = vtkSmartPointer<vtkPoints>::New();

	for (vtkIdType i = 0; i < oldPath->GetNumberOfIds(); i++) {
		vtkIdType oldCellId = oldPath->GetId(i);
		double point[3] = { oldData.getCentrePoints()->GetPoint(oldCellId)[0], oldData.getCentrePoints()->GetPoint(oldCellId)[1], oldData.getCentrePoints()->GetPoint(oldCellId)[2] };
		pointsPath->InsertNextPoint(point);
	}
	pathSearchList->SetPoints(pointsPath);

	for (vtkIdType j = 0; j < path->GetNumberOfIds(); j++) {
		vtkIdType cellId = path->GetId(j);
		double point[3] = { data.getCentrePoints()->GetPoint(cellId)[0], data.getCentrePoints()->GetPoint(cellId)[1], data.getCentrePoints()->GetPoint(cellId)[2] };
		vtkIdType locatorId = pathSearchList->FindPoint(point);
		if (!data.getDoubleLayer()) {
			double *diffVector = oldData.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->GetTuple(oldPath->GetId(locatorId));
			data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->InsertTuple(cellId, diffVector);
			data.getVtkData()->GetCellData()->GetArray("StatusEndo")->SetComponent(path->GetId(j), 0, 3);
		}

		double *diffVector = oldData.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetTuple(oldPath->GetId(locatorId));
		data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->InsertTuple(cellId, diffVector);
		data.getVtkData()->GetCellData()->GetArray("Status")->SetComponent(path->GetId(j), 0, 3);
	}
} // Methods::vtpBridgeMarkerArroundPath

double Methods::DistanceBetweenEndoandEpiSurfaceAlongPath(DataFormat &data, vtkSmartPointer<vtkIdList> path, string distance) {
	vtkSmartPointer<vtkPolyData> polydata = Methods::extractSurface(data);

	polydata = Methods::triangulateSurface(polydata);

	polydata->GetCellData()->RemoveArray("Material");
	polydata->GetCellData()->GetArray("RawMaterial")->SetName("Material");

	DataFormat surface;
	surface.setVtkData(polydata);
	surface.setCentrePoints(Methods::calcCellCentroids(polydata));

	double maxDistance = -std::numeric_limits<double>::max();
	double minDistance = std::numeric_limits<double>::max();
	double avgDistance = 0;
	int count = 0;


	for (vtkIdType i = 0; i < path->GetNumberOfIds(); i++) {
		vector<double> pathPoint = { data.getCentrePoints()->GetPoint(path->GetId(i))[0], data.getCentrePoints()->GetPoint(path->GetId(i))[1], data.getCentrePoints()->GetPoint(path->GetId(i))[2] };

		vector<double> epiPointVec = { Methods::findClosedPointinMaterialInRightEpi(surface, pathPoint) };
		vector<double> endoPointVec = { Methods::findClosedPointinMaterialInRightEndo(surface, pathPoint) };
		if (!isnan(epiPointVec.at(0)) && !isnan(endoPointVec.at(0))) {
			double pointEpi[3] = { epiPointVec.at(0), epiPointVec.at(1), epiPointVec.at(2) };
			double pointEndo[3] = { endoPointVec.at(0), endoPointVec.at(1), endoPointVec.at(2) };
			double distanceBetweenEndoEpi = sqrt(abs(vtkMath::Distance2BetweenPoints(pointEndo, pointEpi)));

			if (distanceBetweenEndoEpi > maxDistance) {
				maxDistance = distanceBetweenEndoEpi;
			}
			if (distanceBetweenEndoEpi < minDistance) {
				minDistance = distanceBetweenEndoEpi;
			}

			avgDistance += distanceBetweenEndoEpi;
			count++;
		}
	}

	avgDistance /= count;
	if (distance.compare("max") == 0) {
		cout << "Max Distance: " << maxDistance << endl;
		return maxDistance;
	}
	else if (distance.compare("min") == 0) {
		cout << "Min Distance: " << minDistance << endl;
		return minDistance;
	}
	else if (distance.compare("avg") == 0) {
		cout << "Avg Distance: " << avgDistance << endl;
		return avgDistance;
	}
} // Methods::DistanceBetweenEndoandEpiSurfaceAlongPath

/*! Calculate the fiktive distance to the epicard in a surface mesh.
 \param data Pointer with the mesh data.
 \param cellId Current cellid
 \param closestPointId Closest point on the epicard.
 \param distanceEndoEpi Fictive distance between endo- and epicard.
 \return The fictive distance in mm from a point to the epicard.
 */
double Methods::calcDistancetoPointinEpi(DataFormat &data, vtkIdType cellId, vtkIdType closestPointId, double distanceEndoEpi) {
	double currentPoint[3] = { 0, 0, 0 };

	currentPoint[0] = data.getCentrePoints()->GetPoint(cellId)[0];
	currentPoint[1] = data.getCentrePoints()->GetPoint(cellId)[1];
	currentPoint[2] = data.getCentrePoints()->GetPoint(cellId)[2];

	double closestPoint[3] = { 0, 0, 0 };
	closestPoint[0] = data.getCentrePoints()->GetPoint(closestPointId)[0];
	closestPoint[1] = data.getCentrePoints()->GetPoint(closestPointId)[1];
	closestPoint[2] = data.getCentrePoints()->GetPoint(closestPointId)[2];

	vtkSmartPointer<vtkIdList> cellPoints = vtkSmartPointer<vtkIdList>::New();
	data.getVtkData()->GetCellPoints(cellId, cellPoints);
	double point1Array[3] = { data.getVtkData()->GetPoint(cellPoints->GetId(0))[0], data.getVtkData()->GetPoint(cellPoints->GetId(0))[1], data.getVtkData()->GetPoint(cellPoints->GetId(0))[2] };
	double point2Array[3] = { data.getVtkData()->GetPoint(cellPoints->GetId(1))[0], data.getVtkData()->GetPoint(cellPoints->GetId(1))[1], data.getVtkData()->GetPoint(cellPoints->GetId(1))[2] };
	double point3Array[3] = { data.getVtkData()->GetPoint(cellPoints->GetId(2))[0], data.getVtkData()->GetPoint(cellPoints->GetId(2))[1], data.getVtkData()->GetPoint(cellPoints->GetId(2))[2] };
	double vectorPoint1Point3[3] = { 0, 0, 0 };
	double vectorPoint1Point2[3] = { 0, 0, 0 };
	double cellNormalvector[3] = { 0, 0, 0 };

	vtkMath::Subtract(point2Array, point1Array, vectorPoint1Point2);
	vtkMath::Subtract(point3Array, point1Array, vectorPoint1Point3);
	vtkMath::Cross(vectorPoint1Point2, vectorPoint1Point3, cellNormalvector);
	vtkMath::Normalize(cellNormalvector);
	vtkMath::MultiplyScalar(cellNormalvector, distanceEndoEpi);

	double linePoint2[3] = { 0, 0, 0 };
	vtkMath::Add(currentPoint, cellNormalvector, linePoint2);

	double distancetoEndo = sqrt(abs(vtkMath::Distance2BetweenPoints(closestPoint, linePoint2)));

	return distancetoEndo;
} // Methods::calcDistancetoPointinEpi

/*! Calculate the fiktive distance to the endocard in a surface mesh.
 \param data Pointer with the mesh data.
 \param cellId Current cellid
 \param closestPointId Closest point on the endocard.
 \param distanceEndoEpi Fictive distance between endo- and epicard.
 \return The fictive distance in mm from a point to the endocard.
 */
double Methods::calcDistancetoPointinEndo(DataFormat &data, vtkIdType cellId, vtkIdType closestPointId, double distanceEndoEpi) {
	double currentPoint[3] = { 0, 0, 0 };

	currentPoint[0] = data.getCentrePoints()->GetPoint(cellId)[0];
	currentPoint[1] = data.getCentrePoints()->GetPoint(cellId)[1];
	currentPoint[2] = data.getCentrePoints()->GetPoint(cellId)[2];

	double closestPoint[3] = { 0, 0, 0 };
	closestPoint[0] = data.getCentrePoints()->GetPoint(closestPointId)[0];
	closestPoint[1] = data.getCentrePoints()->GetPoint(closestPointId)[1];
	closestPoint[2] = data.getCentrePoints()->GetPoint(closestPointId)[2];

	vtkSmartPointer<vtkIdList> cellPoints = vtkSmartPointer<vtkIdList>::New();
	data.getVtkData()->GetCellPoints(closestPointId, cellPoints);
	double point1Array[3] = { data.getVtkData()->GetPoint(cellPoints->GetId(0))[0], data.getVtkData()->GetPoint(cellPoints->GetId(0))[1], data.getVtkData()->GetPoint(cellPoints->GetId(0))[2] };
	double point2Array[3] = { data.getVtkData()->GetPoint(cellPoints->GetId(1))[0], data.getVtkData()->GetPoint(cellPoints->GetId(1))[1], data.getVtkData()->GetPoint(cellPoints->GetId(1))[2] };
	double point3Array[3] = { data.getVtkData()->GetPoint(cellPoints->GetId(2))[0], data.getVtkData()->GetPoint(cellPoints->GetId(2))[1], data.getVtkData()->GetPoint(cellPoints->GetId(2))[2] };
	double vectorPoint1Point3[3] = { 0, 0, 0 };
	double vectorPoint1Point2[3] = { 0, 0, 0 };
	double cellNormalvector[3] = { 0, 0, 0 };

	vtkMath::Subtract(point2Array, point1Array, vectorPoint1Point2);
	vtkMath::Subtract(point3Array, point1Array, vectorPoint1Point3);
	vtkMath::Cross(vectorPoint1Point2, vectorPoint1Point3, cellNormalvector);
	vtkMath::Normalize(cellNormalvector);
	vtkMath::MultiplyScalar(cellNormalvector, distanceEndoEpi);

	double linePoint2[3] = { 0, 0, 0 };
	vtkMath::Add(closestPoint, cellNormalvector, linePoint2);

	double distancetoEndo = sqrt(abs(vtkMath::Distance2BetweenPoints(currentPoint, linePoint2)));

	return distancetoEndo;
} // Methods::calcDistancetoPointinEndo

/*! Triangulate a surface mesh.
 \param polydata The surface that would be trianglate.
 \return The triangulate surface.
 */
vtkSmartPointer<vtkPolyData> Methods::triangulateSurface(vtkSmartPointer<vtkPolyData> polydata) {
	vtkSmartPointer<vtkTriangleFilter> triangleFilter = vtkSmartPointer<vtkTriangleFilter>::New();

	triangleFilter->SetInputData(polydata);
	triangleFilter->Update();

	vtkSmartPointer<vtkCleanPolyData> cleanFilter = vtkSmartPointer<vtkCleanPolyData>::New();
	cleanFilter->SetInputData(triangleFilter->GetOutput());
	cleanFilter->Update();

	return cleanFilter->GetOutput();
}

/*! Smooth a surface mesh.
 \param polydata The surface that would be smooth.
 \param smoothIterations Number of smoothing iteration.
 \return The smooth surface.
 */
vtkSmartPointer<vtkPolyData> Methods::smoothSurface(vtkSmartPointer<vtkPolyData> polydata, double smoothIterations) {
	vtkSmartPointer<vtkSmoothPolyDataFilter> smoother = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();

	smoother->SetInputData(polydata);
	smoother->SetNumberOfIterations(smoothIterations);
	smoother->Update();

	return smoother->GetOutput();
}

/*! Calculate the surface normals of a surface mesh.
 \param polydata The surface mesh
 \return The surface with the calculated normal vectors.
 */
vtkSmartPointer<vtkPolyData> Methods::getSurfaceCellNormals(vtkSmartPointer<vtkPolyData> polydata) {
	vtkSmartPointer<vtkDoubleArray> normals = vtkSmartPointer<vtkDoubleArray>::New();

	normals->SetName("Normals");
	normals->SetNumberOfComponents(3);
	normals->SetNumberOfTuples(polydata->GetNumberOfCells());
	normals->FillComponent(0, 0);
	normals->FillComponent(1, 0);
	normals->FillComponent(2, 0);

	for (vtkIdType i = 0; i < polydata->GetNumberOfCells(); i++) {
		vtkSmartPointer<vtkPoints> cellPoints = polydata->GetCell(i)->GetPoints();

		double point0[3] = { cellPoints->GetPoint(0)[0], cellPoints->GetPoint(0)[1], cellPoints->GetPoint(0)[2] };
		double point1[3] = { cellPoints->GetPoint(1)[0], cellPoints->GetPoint(1)[1], cellPoints->GetPoint(1)[2] };
		double point2[3] = { cellPoints->GetPoint(2)[0], cellPoints->GetPoint(2)[1], cellPoints->GetPoint(2)[2] };

		double vecPoint0Point1[3] = { 0, 0, 0 };
		double vecPoint0Point2[3] = { 0, 0, 0 };

		vtkMath::Subtract(point0, point1, vecPoint0Point1);
		vtkMath::Subtract(point0, point2, vecPoint0Point2);

		double cellNormal[3] = { 0, 0, 0 };
		vtkMath::Cross(vecPoint0Point1, vecPoint0Point2, cellNormal);
		vtkMath::Normalize(cellNormal);

		normals->SetTuple(i, cellNormal);
	}

	polydata->GetCellData()->AddArray(normals);

	return polydata;
} // Methods::getSurfaceCellNormals

/*! Test if the tube surface has self intersections.
 \param tube To test surface mesh
 \param line Ground line of the tube
 \param resolution The mesh resolution
 \return The surface with show outside of the mesh normal vectors.
 */
vtkSmartPointer<vtkPolyData> Methods::tubeIntersectionFilter(vtkSmartPointer<vtkPolyData> tube, vtkSmartPointer<vtkPolyData> line, double resolution) {
	std::vector<std::vector<double>> equalNormalsArray((int)resolution, std::vector<double>(3));

	vtkIdType sizeLine = line->GetNumberOfPoints();

	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

	vtkSmartPointer<vtkPointLocator> pointLocator = vtkSmartPointer<vtkPointLocator>::New();

	pointLocator->InitPointInsertion(points, tube->GetBounds());

	double equalNormalVec[3] = { 0.0, 0.0, 0.0 };
	double vec1[3] = { tube->GetPointData()->GetArray("TubeNormals")->GetComponent(floor((resolution - 1) / 4.0), 0), tube->GetPointData()->GetArray("TubeNormals")->GetComponent(floor((resolution - 1) / 4.0), 1), tube->GetPointData()->GetArray("TubeNormals")->GetComponent(floor((resolution - 1) / 4.0), 2) };

	double vec2[3] = { tube->GetPointData()->GetArray("TubeNormals")->GetComponent(resolution - 1, 0), tube->GetPointData()->GetArray("TubeNormals")->GetComponent(resolution - 1, 1), tube->GetPointData()->GetArray("TubeNormals")->GetComponent(resolution - 1, 2) };

	vtkMath::Cross(vec1, vec2, equalNormalVec);
	vtkMath::Normalize(equalNormalVec);

	for (vtkIdType j = 0; j < resolution; j++) {
		pointLocator->InsertNextPoint(tube->GetPoint(j));
		equalNormalsArray[j][0] = equalNormalVec[0];
		equalNormalsArray[j][1] = equalNormalVec[1];
		equalNormalsArray[j][2] = equalNormalVec[2];
	}

	for (vtkIdType i = 1; i < sizeLine - 1; i++) {
		vtkSmartPointer<vtkPoints> newPoints = vtkSmartPointer<vtkPoints>::New();

		for (vtkIdType l = i * (resolution); l < (i + 1)*(resolution); l++) {
			newPoints->InsertNextPoint(tube->GetPoint(l));
		}

		double vec1[3] = { tube->GetPointData()->GetArray("TubeNormals")->GetComponent((i*resolution + floor((resolution - 1) / 4.0)), 0), tube->GetPointData()->GetArray("TubeNormals")->GetComponent((i*resolution +
		floor((resolution - 1) / 4.0)), 1), tube->GetPointData()->GetArray("TubeNormals")->GetComponent((i*resolution +
		floor((resolution - 1) / 4.0)), 2) };

		double vec2[3] = { tube->GetPointData()->GetArray("TubeNormals")->GetComponent(((i + 1)*resolution) - 1, 0), tube->GetPointData()->GetArray("TubeNormals")->GetComponent(((i + 1)*resolution) -
		1, 1), tube->GetPointData()->GetArray("TubeNormals")->GetComponent(((i + 1)*resolution) - 1, 2) };

		vtkMath::Cross(vec1, vec2, equalNormalVec);
		vtkMath::Normalize(equalNormalVec);

		for (vtkIdType j = 0; j < newPoints->GetNumberOfPoints(); j++) {
			double newPointinView[3] = { newPoints->GetPoint(j)[0], newPoints->GetPoint(j)[1], newPoints->GetPoint(j)[2] };
			double equalVec[3] = { equalNormalsArray[j][0], equalNormalsArray[j][1], equalNormalsArray[j][2] };

			double eval = vtkPlane::Evaluate(equalVec, points->GetPoint(j), newPointinView);

			if (eval <= 0) {
				vtkIdType newPointIdinView = pointLocator->FindClosestInsertedPoint(newPointinView);

				double newPoint[3] = { points->GetPoint(newPointIdinView)[0], points->GetPoint(newPointIdinView)[1], points->GetPoint(newPointIdinView)[2] };
				vtkIdType newPointId = tube->FindPoint(newPoint);

				vtkIdType oldPointId = tube->FindPoint(newPointinView);

				vtkSmartPointer<vtkIdList> cells = vtkSmartPointer<vtkIdList>::New();
				tube->GetPointCells(oldPointId, cells);

				for (vtkIdType k = 0; k < cells->GetNumberOfIds(); k++) {
					tube->ReplaceCellPoint(cells->GetId(k), oldPointId, newPointId);
				}
			}
			if (eval > 0) {
				pointLocator->InsertPoint(j, newPoints->GetPoint(j));
				equalNormalsArray[j][0] = equalNormalVec[0];
				equalNormalsArray[j][1] = equalNormalVec[1];
				equalNormalsArray[j][2] = equalNormalVec[2];
			}
		}
	}
	tube->BuildCells();
	tube->BuildLinks();

	for (vtkIdType i = 0; i < tube->GetNumberOfCells(); i++) {
		vtkSmartPointer<vtkIdList> pointIds = tube->GetCell(i)->GetPointIds();
		if ((pointIds->GetId(0) == pointIds->GetId(1)) ||
			(pointIds->GetId(0) == pointIds->GetId(2)) ||
			(pointIds->GetId(1) == pointIds->GetId(2))) {
			tube->DeleteCell(i);
		}
	}
	tube->RemoveDeletedCells();
	tube->BuildLinks();

	vtkSmartPointer<vtkCleanPolyData> cleanFilter = vtkSmartPointer<vtkCleanPolyData>::New();
	cleanFilter->SetInputData(tube);
	cleanFilter->Update();

	return cleanFilter->GetOutput();
} // Methods::tubeIntersectionFilter

/*! Tetrahedralize a vtp surface with tetgen.
 \param inputData Th surface mesh
 \return The tetrahedralize mesh as unstructed grid.
 */
vtkSmartPointer<vtkUnstructuredGrid> Methods::tetrahedralizeVTPSurface(vtkSmartPointer<vtkPolyData> inputData) {
	vtkSmartPointer<vtkTriangleFilter> triangleFilter = vtkSmartPointer<vtkTriangleFilter>::New();

	triangleFilter->SetInputData(inputData);
	triangleFilter->Update();

	vtkSmartPointer<vtkCleanPolyData> cleanFilter = vtkSmartPointer<vtkCleanPolyData>::New();
	cleanFilter->SetInputData(triangleFilter->GetOutput());
	cleanFilter->Update();

	inputData = cleanFilter->GetOutput();

	tetgenio inData = ConvertTetgenio::convertVTPtoTetgenio(inputData);

	tetgenio outData = tetgenio();

	tetgenbehavior behavior = tetgenbehavior();

	char arg[] = "Yq1.0";
	behavior.parse_commandline(arg);

	tetrahedralize(&behavior, &inData, &outData);

	vtkSmartPointer<vtkUnstructuredGrid> ugrid = ConvertTetgenio::convertTetgeniotoVTU(outData);

	return ugrid;
} // Methods::tetrahedralizeVTPSurface

/*! Deleted all cell that has no neighbour on the edges.
 \param data The surface mesh
 \return The cleaned surface mesh.
 */
vtkSmartPointer<vtkPolyData> Methods::closeSurfaceFilter(vtkSmartPointer<vtkPolyData> data) {
	for (vtkIdType i = 0; i < data->GetNumberOfCells(); i++) {
		vtkSmartPointer<vtkIdList> pointIds = data->GetCell(i)->GetPointIds();
		if ((pointIds->GetId(0) == pointIds->GetId(1)) ||
			(pointIds->GetId(0) == pointIds->GetId(2)) ||
			(pointIds->GetId(1) == pointIds->GetId(2))) {
			data->DeleteCell(i);
		}
		if (pointIds->GetNumberOfIds() < 3) {
			data->DeleteCell(i);
		}
	}
	data->RemoveDeletedCells();
	data->BuildLinks();

	vtkSmartPointer<vtkTriangleFilter> triangleFilter = vtkSmartPointer<vtkTriangleFilter>::New();
	triangleFilter->SetInputData(data);
	triangleFilter->Update();

	data = triangleFilter->GetOutput();

	for (vtkIdType i = 0; i < data->GetNumberOfCells(); i++) {
		vtkSmartPointer<vtkIdList> cellPointsIds = vtkSmartPointer<vtkIdList>::New();
		cellPointsIds->DeepCopy(data->GetCell(i)->GetPointIds());
		for (vtkIdType j = 0; j < data->GetNumberOfCells(); j++) {
			if (j != i) {
				vtkSmartPointer<vtkIdList> equalCellPointsIds = vtkSmartPointer<vtkIdList>::New();
				equalCellPointsIds->DeepCopy(data->GetCell(j)->GetPointIds());
				vtkIdType count = 0;
				for (vtkIdType k = 0; k < equalCellPointsIds->GetNumberOfIds(); k++) {
					vtkIdType isId = cellPointsIds->IsId(equalCellPointsIds->GetId(k));
					if (isId > -1) {
						count++;
					}
				}
				if ((count == equalCellPointsIds->GetNumberOfIds()) && (equalCellPointsIds->GetNumberOfIds() > 0)) {
					data->DeleteCell(i);
				}
			}
		}
	}

	data->RemoveDeletedCells();

	vtkSmartPointer<vtkCleanPolyData> cleanFilter = vtkSmartPointer<vtkCleanPolyData>::New();
	cleanFilter->SetInputData(data);
	cleanFilter->Update();

	data = cleanFilter->GetOutput();


	data->BuildCells();
	data->BuildLinks();

	bool change = true;
	while (change) {
		change = false;
		for (vtkIdType i = 0; i < data->GetNumberOfCells(); i++) {
			vtkSmartPointer<vtkIdList> neighbours = vtkSmartPointer<vtkIdList>::New();

			for (int j = 0; j < data->GetCell(i)->GetNumberOfEdges(); j++) {
				vtkSmartPointer<vtkIdList> edgePointIdList = data->GetCell(i)->GetEdge(j)->GetPointIds();

				vtkSmartPointer<vtkIdList> tempNeighborCellIds = vtkSmartPointer<vtkIdList>::New();
				data->GetCellNeighbors(i, edgePointIdList, tempNeighborCellIds);

				for (vtkIdType k = 0; k < tempNeighborCellIds->GetNumberOfIds(); k++) {
					neighbours->InsertUniqueId(tempNeighborCellIds->GetId(k));
				}
			}
			if (neighbours->GetNumberOfIds() < 3) {
				data->DeleteCell(i);
				change = true;
			}
		}
		data->RemoveDeletedCells();
		data->BuildCells();
		data->BuildLinks();

		vtkSmartPointer<vtkCleanPolyData> cleanFilter = vtkSmartPointer<vtkCleanPolyData>::New();
		cleanFilter->SetInputData(data);
		cleanFilter->Update();

		data = cleanFilter->GetOutput();
	}
	return data;
} // Methods::closeSurfaceFilter

/*! Write out intermediate data.
 \param data Data with the mesh
 \param filename Stored filename.
 */
void Methods::writeIntermediateData(DataFormat data, string filename) {
	UniversalWriter writer;

	if (data.getInputType() == DataFormat::vtp) {
		filename = filename + ".vtp";
		writer.write(data, filename);
	}
	else if (data.getInputType() == DataFormat::vtu) {
		filename = filename + ".vtu";
		writer.write(data, filename);
	}
	else if (data.getInputType() == DataFormat::vtk) {
		filename = filename + ".vtk";
		writer.write(data, filename);
	}
}

/*! Dilatate the path transmural.
 \param data Pointer with the mesh data.
 \param path1 First path on the surface of the mesh.
 \param path2 Second path on the surface of the mesh.
 \param width Dilatation width
 \param atrium In which atrial the path would be dilatated.
 \return The cell ids of the dilatated path.
 */
vtkSmartPointer<vtkIdList> Methods::growPathTransmural(DataFormat &data, vtkSmartPointer<vtkIdList> path1, vtkSmartPointer<vtkIdList> path2, double width, string atrium) {
	DataFormat surface;

	surface.setVtkData(data.getSurfacewithNormals());
	surface.setCentrePoints(data.getCentrePointsSurface());

	vtkSmartPointer<vtkIdList> returnIdList = vtkSmartPointer<vtkIdList>::New();
	set<vtkIdType> viewedCells;

	vtkSmartPointer<vtkIdList> pathLarge = path1;
	vtkSmartPointer<vtkIdList> pathSmall = path2;

	for (vtkIdType i = 0; i < pathLarge->GetNumberOfIds(); i++) {
		vector<double> pillPoint1 = { data.getCentrePointsSurface()->GetPoint(pathLarge->GetId(i))[0], data.getCentrePointsSurface()->GetPoint(pathLarge->GetId(i))[1], data.getCentrePointsSurface()->GetPoint(pathLarge->GetId(i))[2] };

		vector<double> pillPoint2 = Methods::findClosedPointOnPath(surface, pillPoint1, pathSmall);

		DataFormat pillData;
		pillData = getPill(pillPoint1, pillPoint2, 0.2, width);
		vtkSmartPointer<vtkPolyData> pill = vtkPolyData::SafeDownCast(pillData.getVtkData());

		vtkSmartPointer<vtkSelectEnclosedPoints> selectEnclosedPoints = vtkSmartPointer<vtkSelectEnclosedPoints>::New();
		selectEnclosedPoints->Initialize(pill);
		selectEnclosedPoints->SetTolerance(0);

		vtkSmartPointer<vtkIdList> cells = getCellsinRadius(data, pillPoint1, 15);

		for (vtkIdType j = 0; j < cells->GetNumberOfIds(); j++) {
			double *point = data.getCentrePoints()->GetPoint(cells->GetId(j));
			int inside = selectEnclosedPoints->IsInsideSurface(point);

			if ((inside == 1) && (viewedCells.find(cells->GetId(j)) == viewedCells.end())) {
				int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(cells->GetId(j), 0);
				if (atrium.compare("right") == 0) {
					if (Material::isInRight(material)) {
						returnIdList->InsertNextId(cells->GetId(j));
					}
				}
				else if (atrium.compare("left") == 0) {
					if (Material::isInLeft(material)) {
						returnIdList->InsertNextId(cells->GetId(j));
					}
				}
				else if (Material::isInLeft(material) || Material::isInRight(material)) {
					returnIdList->InsertNextId(cells->GetId(j));
				}

				viewedCells.insert(cells->GetId(j));
			}
		}
	}

	pathLarge = path2;
	pathSmall = path1;

	for (vtkIdType i = 0; i < pathLarge->GetNumberOfIds(); i++) {
		vector<double> pillPoint1 = { data.getCentrePointsSurface()->GetPoint(pathLarge->GetId(i))[0], data.getCentrePointsSurface()->GetPoint(pathLarge->GetId(i))[1], data.getCentrePointsSurface()->GetPoint(pathLarge->GetId(i))[2] };
		vector<double> pillPoint2 = Methods::findClosedPointOnPath(surface, pillPoint1, pathSmall);

		DataFormat pillData;
		pillData = getPill(pillPoint1, pillPoint2, 0.2, width);
		vtkSmartPointer<vtkPolyData> pill = vtkPolyData::SafeDownCast(pillData.getVtkData());

		vtkSmartPointer<vtkSelectEnclosedPoints> selectEnclosedPoints = vtkSmartPointer<vtkSelectEnclosedPoints>::New();
		selectEnclosedPoints->Initialize(pill);
		selectEnclosedPoints->SetTolerance(0);

		vtkSmartPointer<vtkIdList> cells = getCellsinRadius(data, pillPoint1, 15);

		for (vtkIdType j = 0; j < cells->GetNumberOfIds(); j++) {
			double *point = data.getCentrePoints()->GetPoint(cells->GetId(j));
			int inside = selectEnclosedPoints->IsInsideSurface(point);

			if ((inside == 1) && (viewedCells.find(cells->GetId(j)) == viewedCells.end())) {
				int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(cells->GetId(j), 0);
				if (atrium.compare("right") == 0) {
					if (Material::isInRight(material)) {
						returnIdList->InsertNextId(cells->GetId(j));
					}
				}
				else if (atrium.compare("left") == 0) {
					if (Material::isInLeft(material)) {
						returnIdList->InsertNextId(cells->GetId(j));
					}
				}
				else if (Material::isInLeft(material) || Material::isInRight(material)) {
					returnIdList->InsertNextId(cells->GetId(j));
				}

				viewedCells.insert(cells->GetId(j));
			}
		}
	}

	return returnIdList;
} // Methods::growPathTransmural

// Methods::growPathTransmuralEndo

/*! Dilatate the path transmural only in one material.
 \param data Pointer with the mesh data.
 \param path Path which would be dilatated transmural.
 \param width Dilatation width
 \param inMaterial In which material the path would be dilatated.
 \param atrium In which atrial the path would be dilatated.
 \return The cell ids of the dilatated path.
 */
vtkSmartPointer<vtkIdList> Methods::growPathTransmural(DataFormat &data, vtkSmartPointer<vtkIdList> path, double width, Material::Mat inMaterial, string atrium) {
	return Methods::growPathTransmural(data, path, width, inMaterial, inMaterial, atrium);
}

vtkSmartPointer<vtkIdList> Methods::growPathTransmuralEndo(DataFormat &data, vtkSmartPointer<vtkIdList> path, double width, Material::Mat inMaterial) {
	return Methods::growPathTransmuralEndo(data, path, width, inMaterial, inMaterial);
}

// Methods::growPathTransmural
vtkSmartPointer<vtkIdList> Methods::growPathTransmuralEndo(DataFormat &data, vtkSmartPointer<vtkIdList> orgPath, double width, Material::Mat inMaterial1, Material::Mat inMaterial2) {

	DataFormat surface;

	surface.setVtkData(data.getSurfacewithNormals());
	surface.setCentrePoints(data.getCentrePointsSurface());

	vector<double> startPoint = pointAtPercentOfPath(data, orgPath, 0);
	vector<double> targetPoint = pointAtPercentOfPath(data, orgPath, 100);

	vector<double> pointEpiStart = Methods::findClosedPointinMaterial(surface, startPoint, Material::Vorhof_rechts);
	vector<double> pointEndoStart = Methods::findClosedPointinMaterial(surface, startPoint, Material::subendo_right_Atrial_Cardium);
	vector<double> pointEpiTarget = Methods::findClosedPointinMaterial(surface, targetPoint, Material::Vorhof_rechts);
	vector<double> pointEndoTarget = Methods::findClosedPointinMaterial(surface, targetPoint, Material::subendo_right_Atrial_Cardium);

	vtkSmartPointer<vtkIdList> epiPath = Methods::pathSearchNearToPath(data, pointEpiStart, pointEpiTarget, orgPath, "surfaceEpi");

	vtkSmartPointer<vtkIdList> endoPath = Methods::pathSearchNearToPath(data, pointEndoStart, pointEndoTarget, orgPath, "surfaceEndo");

	vtkSmartPointer<vtkIdList> returnIds = vtkSmartPointer<vtkIdList>::New();

	set<vtkIdType> changedId;

	vtkSmartPointer<vtkDoubleArray> tempDiffVec = vtkDoubleArray::SafeDownCast(data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector"));

	vtkSmartPointer<vtkDoubleArray> tempArray = vtkSmartPointer<vtkDoubleArray>::New();
	tempArray->SetNumberOfComponents(3);
	tempArray->SetNumberOfTuples(data.getCentrePoints()->GetNumberOfPoints());
	tempArray->FillComponent(0, 0);
	tempArray->FillComponent(1, 0);
	tempArray->FillComponent(2, 0);

	for (vtkIdType i = 0; i < epiPath->GetNumberOfIds(); i++) {
		vector<double> pillPoint1 = { data.getCentrePoints()->GetPoint(epiPath->GetId(i))[0], data.getCentrePoints()->GetPoint(epiPath->GetId(i))[1], data.getCentrePoints()->GetPoint(epiPath->GetId(i))[2] };

		vector<double> pillPoint2 = Methods::findClosedPointOnPath(data, pillPoint1, endoPath);

		DataFormat pillData;
		pillData = getPill(pillPoint1, pillPoint2, 0.2, width);

		vtkSmartPointer<vtkPolyData> pill = vtkPolyData::SafeDownCast(pillData.getVtkData());

		vtkSmartPointer<vtkSelectEnclosedPoints> selectEnclosedPoints = vtkSmartPointer<vtkSelectEnclosedPoints>::New();
		selectEnclosedPoints->Initialize(pill);
		selectEnclosedPoints->SetTolerance(0);

		vtkSmartPointer<vtkIdList> cells = getCellsinRadius(data, pillPoint1, 2 * ConfigFiberorientation::maxWitdthAtrialWall);

		vtkIdType parentId = findClosedPointIdOnPath(data, pillPoint1, orgPath);

		for (vtkIdType j = 0; j < cells->GetNumberOfIds(); j++) {
			vtkIdType currentId = cells->GetId(j);

			double *point = data.getCentrePoints()->GetPoint(currentId);
			int inside = selectEnclosedPoints->IsInsideSurface(point);

			if (inside == 1) {
				int Rawmaterial = data.getVtkData()->GetCellData()->GetArray("RawMaterial")->GetComponent(currentId, 0);
				int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(currentId, 0);
				if (Rawmaterial == inMaterial1) {
					if ((tempDiffVec->GetComponent(parentId, 0) != 0) ||
						(tempDiffVec->GetComponent(parentId, 1) != 0) || (tempDiffVec->GetComponent(parentId, 2) != 0)) {
						tempArray->InsertComponent(currentId, 0, (tempDiffVec->GetComponent(parentId, 0) + tempArray->GetComponent(currentId, 0)));
						tempArray->InsertComponent(currentId, 1, (tempDiffVec->GetComponent(parentId, 1) + tempArray->GetComponent(currentId, 1)));
						tempArray->InsertComponent(currentId, 2, (tempDiffVec->GetComponent(parentId, 2) + tempArray->GetComponent(currentId, 2)));
						changedId.insert(currentId);
					}
				}
			}
		}
	}

	for (vtkIdType i = 0; i < endoPath->GetNumberOfIds(); i++) {
		vector<double> pillPoint1 = { data.getCentrePoints()->GetPoint(endoPath->GetId(i))[0], data.getCentrePoints()->GetPoint(endoPath->GetId(i))[1], data.getCentrePoints()->GetPoint(endoPath->GetId(i))[2] };

		vector<double> pillPoint2 = Methods::findClosedPointOnPath(data, pillPoint1, epiPath);

		DataFormat pillData;
		pillData = getPill(pillPoint1, pillPoint2, 0.2, width);

		vtkSmartPointer<vtkPolyData> pill = vtkPolyData::SafeDownCast(pillData.getVtkData());

		vtkSmartPointer<vtkSelectEnclosedPoints> selectEnclosedPoints = vtkSmartPointer<vtkSelectEnclosedPoints>::New();
		selectEnclosedPoints->Initialize(pill);
		selectEnclosedPoints->SetTolerance(0);

		vtkSmartPointer<vtkIdList> cells = getCellsinRadius(data, pillPoint1, 2 * ConfigFiberorientation::maxWitdthAtrialWall);

		vtkIdType parentId = findClosedPointIdOnPath(data, pillPoint1, orgPath); // Einfügen

		for (vtkIdType j = 0; j < cells->GetNumberOfIds(); j++) {
			vtkIdType currentId = cells->GetId(j);

			double *point = data.getCentrePoints()->GetPoint(currentId);
			int inside = selectEnclosedPoints->IsInsideSurface(point);

			if (inside == 1) {
				int Rawmaterial = data.getVtkData()->GetCellData()->GetArray("RawMaterial")->GetComponent(currentId, 0);
				int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(currentId, 0);
				if (Rawmaterial == inMaterial1) {
					if ((tempDiffVec->GetComponent(parentId, 0) != 0) ||
						(tempDiffVec->GetComponent(parentId, 1) != 0) || (tempDiffVec->GetComponent(parentId, 2) != 0)) {
						tempArray->InsertComponent(currentId, 0, (tempDiffVec->GetComponent(parentId, 0) + tempArray->GetComponent(currentId, 0)));
						tempArray->InsertComponent(currentId, 1, (tempDiffVec->GetComponent(parentId, 1) + tempArray->GetComponent(currentId, 1)));
						tempArray->InsertComponent(currentId, 2, (tempDiffVec->GetComponent(parentId, 2) + tempArray->GetComponent(currentId, 2)));
						changedId.insert(currentId);
					}
				}
			}
		}
	}


	for (set<vtkIdType>::iterator iter = changedId.begin(); iter != changedId.end(); iter++) {
		int Rawmaterial = data.getVtkData()->GetCellData()->GetArray("RawMaterial")->GetComponent(*iter, 0);
		int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(*iter, 0);
		if (((data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(*iter, 0) == 0) &&
			(Rawmaterial == inMaterial1)) || ((Rawmaterial == inMaterial1) && (material == inMaterial2))) {
			double *tempPoint = tempArray->GetTuple(*iter);
			vtkMath::Normalize(tempPoint);
			data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->InsertTuple(*iter, tempPoint);
			data.getVtkData()->GetCellData()->GetArray("Status")->SetComponent(*iter, 0, 3);
			returnIds->InsertUniqueId(*iter);
		}
	}

	vtkSmartPointer<vtkIdList> pathInEndo = vtkSmartPointer<vtkIdList>::New();


	for (vtkIdType i = 0; i < orgPath->GetNumberOfIds(); i++) {
		int Rawmaterial = data.getVtkData()->GetCellData()->GetArray("RawMaterial")->GetComponent((orgPath->GetId(i)), 0);
		if (Rawmaterial == inMaterial1) {
			pathInEndo->InsertUniqueId((orgPath->GetId(i)));
		}
		else {
			data.getVtkData()->GetCellData()->GetArray("Status")->SetComponent((orgPath->GetId(i)), 0, 0);
		}
	}

	vtkSmartPointer<vtkIdList> returnList = Methods::add2Paths(returnIds, pathInEndo);


	return returnList;
} // Methods::growPathTransmuralEndo

/*! Dilatate the path transmural only in two material.
 \param data Pointer with the mesh data.
 \param path Path which would be dilatated transmural.
 \param width Dilatation width
 \param inMaterial1 First material in which the path would be dilatated.
 \param inMaterial1 Second material in which the path would be dilatated.
 \param atrium In which atrial the path would be dilatated.
 \return The cell ids of the dilatated path.
 */

vtkSmartPointer<vtkIdList> Methods::growPathfromEndo(DataFormat &data, vtkSmartPointer<vtkIdList> path, double width, Material::Mat inMaterial) {
	return Methods::growPathTransmuralEndo(data, path, width, inMaterial, inMaterial);
}

vtkSmartPointer<vtkIdList> Methods::growPathfromEndo(DataFormat &data, vtkSmartPointer<vtkIdList> orgPath, double width, Material::Mat inMaterial1, Material::Mat inMaterial2) {
	DataFormat surface;

	surface.setVtkData(data.getSurfacewithNormals());
	surface.setCentrePoints(data.getCentrePointsSurface());

	vector<double> startPoint = pointAtPercentOfPath(data, orgPath, 0);
	vector<double> targetPoint = pointAtPercentOfPath(data, orgPath, 100);

	vector<double> pointEndoStart = Methods::findClosedPointinMaterial(surface, startPoint, Material::subendo_right_Atrial_Cardium);
	vector<double> pointEndoTarget = Methods::findClosedPointinMaterial(surface, targetPoint, Material::subendo_right_Atrial_Cardium);

	vtkSmartPointer<vtkIdList> endoPath = Methods::pathSearchNearToPath(data, pointEndoStart, pointEndoTarget, orgPath, "surfaceEndo");

	vtkSmartPointer<vtkIdList> returnIds = vtkSmartPointer<vtkIdList>::New();

	set<vtkIdType> changedId;

	vtkSmartPointer<vtkDoubleArray> tempDiffVec = vtkDoubleArray::SafeDownCast(data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector"));

	vtkSmartPointer<vtkDoubleArray> tempArray = vtkSmartPointer<vtkDoubleArray>::New();
	tempArray->SetNumberOfComponents(3);
	tempArray->SetNumberOfTuples(data.getCentrePoints()->GetNumberOfPoints());
	tempArray->FillComponent(0, 0);
	tempArray->FillComponent(1, 0);
	tempArray->FillComponent(2, 0);


	for (vtkIdType i = 0; i < endoPath->GetNumberOfIds(); i++) {
		vector<double> currentPoint = { data.getCentrePoints()->GetPoint(endoPath->GetId(i))[0], data.getCentrePoints()->GetPoint(endoPath->GetId(i))[1], data.getCentrePoints()->GetPoint(endoPath->GetId(i))[2] };

		vtkIdType parentId = findClosedPointIdOnPath(data, currentPoint, orgPath);

		vtkSmartPointer<vtkIdList> cells = getCellsinRadius(data, currentPoint, width);

		for (vtkIdType j = 0; j < cells->GetNumberOfIds(); j++) {
			vtkIdType currentId = cells->GetId(j);

			int Rawmaterial = data.getVtkData()->GetCellData()->GetArray("RawMaterial")->GetComponent(currentId, 0);
			int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(currentId, 0);
			if (Rawmaterial == inMaterial1) {
				if ((tempDiffVec->GetComponent(parentId, 0) != 0) ||
					(tempDiffVec->GetComponent(parentId, 1) != 0) || (tempDiffVec->GetComponent(parentId, 2) != 0)) {
					tempArray->InsertComponent(currentId, 0, (tempDiffVec->GetComponent(parentId, 0) + tempArray->GetComponent(currentId, 0)));
					tempArray->InsertComponent(currentId, 1, (tempDiffVec->GetComponent(parentId, 1) + tempArray->GetComponent(currentId, 1)));
					tempArray->InsertComponent(currentId, 2, (tempDiffVec->GetComponent(parentId, 2) + tempArray->GetComponent(currentId, 2)));
					changedId.insert(currentId);
				}
			}
		}
	}


	for (set<vtkIdType>::iterator iter = changedId.begin(); iter != changedId.end(); iter++) {
		int Rawmaterial = data.getVtkData()->GetCellData()->GetArray("RawMaterial")->GetComponent(*iter, 0);
		int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(*iter, 0);
		if (((data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(*iter, 0) == 0) &&
			((Rawmaterial == inMaterial1))) || ((Rawmaterial == inMaterial1) && (material == inMaterial2))) {
			double *tempPoint = tempArray->GetTuple(*iter);
			vtkMath::Normalize(tempPoint);
			data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->InsertTuple(*iter, tempPoint);
			data.getVtkData()->GetCellData()->GetArray("Status")->SetComponent(*iter, 0, 3);
			returnIds->InsertUniqueId(*iter);
		}
	}

	vtkSmartPointer<vtkIdList> pathInEndo = vtkSmartPointer<vtkIdList>::New();


	for (vtkIdType i = 0; i < orgPath->GetNumberOfIds(); i++) {
		int Rawmaterial = data.getVtkData()->GetCellData()->GetArray("RawMaterial")->GetComponent((orgPath->GetId(i)), 0);
		if (Rawmaterial == inMaterial1) {
			pathInEndo->InsertUniqueId((orgPath->GetId(i)));
		}
		else {
			data.getVtkData()->GetCellData()->GetArray("Status")->SetComponent((orgPath->GetId(i)), 0, 0);
		}
	}

	vtkSmartPointer<vtkIdList> returnList = Methods::add2Paths(returnIds, pathInEndo);


	return returnList;
} // Methods::growPathfromEndo

/*! Dilatate the path transmural only in two material.
 \param data Pointer with the mesh data.
 \param path Path which would be dilatated transmural.
 \param width Dilatation width
 \param inMaterial1 First material in which the path would be dilatated.
 \param inMaterial1 Second material in which the path would be dilatated.
 \param atrium In which atrial the path would be dilatated.
 \return The cell ids of the dilatated path.
 */
vtkSmartPointer<vtkIdList> Methods::growPathTransmural(DataFormat &data, vtkSmartPointer<vtkIdList> path, double width, Material::Mat inMaterial1, Material::Mat inMaterial2, string atrium) {
	DataFormat surface;

	surface.setVtkData(data.getSurfacewithNormals());
	surface.setCentrePoints(data.getCentrePointsSurface());

	vtkSmartPointer<vtkIdList> returnIds = vtkSmartPointer<vtkIdList>::New();

	DataFormat surfaceData;
	surfaceData.setInputType(DataFormat::vtp);
	surfaceData.setVtkData(data.getSurfacewithNormals());
	surfaceData.setCentrePoints(data.getCentrePointsSurface());

	set<vtkIdType> changedId;

	vtkSmartPointer<vtkDoubleArray> tempDiffVec = vtkDoubleArray::SafeDownCast(data.getVtkData()->GetCellData()->GetArray("BasedPathDifferenceVector"));

	vtkSmartPointer<vtkDoubleArray> tempArray = vtkSmartPointer<vtkDoubleArray>::New();
	tempArray->SetNumberOfComponents(3);
	tempArray->SetNumberOfTuples(data.getCentrePoints()->GetNumberOfPoints());
	tempArray->FillComponent(0, 0);
	tempArray->FillComponent(1, 0);
	tempArray->FillComponent(2, 0);

	for (vtkIdType i = 0; i < path->GetNumberOfIds(); i++) {
		vtkIdType parentId = path->GetId(i);


		vector<double> pathPoint = { data.getCentrePoints()->GetPoint(parentId)[0], data.getCentrePoints()->GetPoint(parentId)[1], data.getCentrePoints()->GetPoint(parentId)[2] };
		vector<double> pillPoint1;
		vector<double> pillPoint2;


		if (atrium.compare("right") == 0) {
			pillPoint1 = { Methods::findClosedPointinMaterialInRightEpi(surface, pathPoint) };
			pillPoint2 = { Methods::findClosedPointinMaterialInRightEndo(surface, pathPoint) };
		}
		else if (atrium.compare("left") == 0) {
			pillPoint1 = { Methods::findClosedPointinMaterialInLeftEpi(surface, pathPoint) };
			pillPoint2 = { Methods::findClosedPointinMaterialInLeftEndo(surface, pathPoint) };
		}

		DataFormat pillData;
		pillData = getPill(pillPoint1, pillPoint2, 0.2, width);
		vtkSmartPointer<vtkPolyData> pill = vtkPolyData::SafeDownCast(pillData.getVtkData());

		vtkSmartPointer<vtkSelectEnclosedPoints> selectEnclosedPoints = vtkSmartPointer<vtkSelectEnclosedPoints>::New();
		selectEnclosedPoints->Initialize(pill);
		selectEnclosedPoints->SetTolerance(0);

		vtkSmartPointer<vtkIdList> cells = getCellsinRadius(data, pillPoint1, ConfigFiberorientation::maxWitdthAtrialWall);

		for (vtkIdType j = 0; j < cells->GetNumberOfIds(); j++) {
			vtkIdType currentId = cells->GetId(j);

			double *point = data.getCentrePoints()->GetPoint(currentId);
			int inside = selectEnclosedPoints->IsInsideSurface(point);

			if (inside == 1) {
				int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(currentId, 0);
				if ((material == inMaterial1) || (material == inMaterial2)) {
					if ((tempDiffVec->GetComponent(parentId, 0) != 0) ||
						(tempDiffVec->GetComponent(parentId, 1) != 0) || (tempDiffVec->GetComponent(parentId, 2) != 0)) {
						tempArray->InsertComponent(currentId, 0, (tempDiffVec->GetComponent(parentId, 0) + tempArray->GetComponent(currentId, 0)));
						tempArray->InsertComponent(currentId, 1, (tempDiffVec->GetComponent(parentId, 1) + tempArray->GetComponent(currentId, 1)));
						tempArray->InsertComponent(currentId, 2, (tempDiffVec->GetComponent(parentId, 2) + tempArray->GetComponent(currentId, 2)));
					}
					changedId.insert(currentId);
				}
			}
		}
	}

	for (set<vtkIdType>::iterator iter = changedId.begin(); iter != changedId.end(); iter++) {
		int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(*iter, 0);
		if ((data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(*iter, 0) == 0) &&
			((material == inMaterial1) || (material == inMaterial2))) {
			double *tempPoint = tempArray->GetTuple(*iter);
			vtkMath::Normalize(tempPoint);
			data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->InsertTuple(*iter, tempPoint);
			data.getVtkData()->GetCellData()->GetArray("Status")->SetComponent(*iter, 0, 3);
			returnIds->InsertUniqueId(*iter);
		}
	}

	vtkSmartPointer<vtkIdList> returnList = Methods::add2Paths(returnIds, path);
	return returnList;
} // Methods::growPathTransmural

/*! Delete temporray cell array that not more be need.
 \param data Pointer with the mesh data.
 */
void Methods::deleteArrays(DataFormat &data) {
	data.getVtkData()->GetCellData()->RemoveArray("Status");
	data.getVtkData()->GetCellData()->RemoveArray("BasedPathDifferenceVector");
	data.getVtkData()->GetCellData()->RemoveArray("basePath");
	if (!data.getDoubleLayer()) {
		data.getVtkData()->GetCellData()->RemoveArray("StatusEndo");
		data.getVtkData()->GetCellData()->RemoveArray("Normals");
		data.getVtkData()->GetCellData()->RemoveArray("basePathEndo");
	}
}

/*! Calculate the cell which are inside in another mesh.
 \param data Pointer with the mesh data.
 \param insideData Pointer with the other mesh data.
 \return The cell ids with are in the inside mseh form the data mesh.
 */
vtkSmartPointer<vtkIdList> Methods::getCellsWhichWereInside(DataFormat &data, DataFormat &insideData) {
	vtkSmartPointer<vtkPolyData> surface = vtkSmartPointer<vtkPolyData>::New();

	if (insideData.getInputType() == DataFormat::vtp) {
		surface = vtkPolyData::SafeDownCast(insideData.getVtkData());
	}
	else {
		surface = triangulateSurface(extractSurface(insideData));
	}
	vtkSmartPointer<vtkSelectEnclosedPoints> selectEnclosedPoints = vtkSmartPointer<vtkSelectEnclosedPoints>::New();
	selectEnclosedPoints->Initialize(surface);
	selectEnclosedPoints->SetTolerance(0);

	vtkSmartPointer<vtkIdList> returnIdList = vtkSmartPointer<vtkIdList>::New();

	for (vtkIdType j = 0; j < data.getCentrePoints()->GetNumberOfPoints(); j++) {
		double *point = data.getCentrePoints()->GetPoint(j);
		int inside = selectEnclosedPoints->IsInsideSurface(point);

		if (inside == 1) {
			returnIdList->InsertNextId(j);
		}
	}

	return returnIdList;
} // Methods::getCellsWhichWereInside

/*! Test if a surface has a hole.
 \param surface The test surface as vtk polydata.
 \return True if it has a hole or else false.
 */
bool Methods::hasHoles(vtkSmartPointer<vtkPolyData> surface) {
	for (vtkIdType i = 0; i < surface->GetNumberOfCells(); i++) {
		vtkSmartPointer<vtkIdList> neighbours = vtkSmartPointer<vtkIdList>::New();

		for (int j = 0; j < surface->GetCell(i)->GetNumberOfEdges(); j++) {
			vtkSmartPointer<vtkIdList> edgePointIdList = surface->GetCell(i)->GetEdge(j)->GetPointIds();

			vtkSmartPointer<vtkIdList> tempNeighborCellIds = vtkSmartPointer<vtkIdList>::New();
			surface->GetCellNeighbors(i, edgePointIdList, tempNeighborCellIds);

			for (vtkIdType k = 0; k < tempNeighborCellIds->GetNumberOfIds(); k++) {
				neighbours->InsertUniqueId(tempNeighborCellIds->GetId(k));
			}
		}
		if (neighbours->GetNumberOfIds() < 3) {
			return true;
		}
	}

	return false;
}

/*!
 \page Methods

 \section DESCRIPTION_Methods DESCRIPTION
 This is the methods container that contains the most methods for the atrials tools.

 \section SOURCE_Methods SOURCE

 Methods.cpp

 \section SEEALSO_mainFiberOrientation SEE ALSO
 \ref VcgTools \ref Reader \ref Writer \ref Config

 \section CHANGELOG_Methods CHANGELOG
 V1.0.0 - 01.09.2015 (Andreas Wachter): Starting with changelog\n
 */

/*! \file TestTetrahedralize.cpp

   \brief Test if an bridge surface is been tetralizied with tetgen for include it in SetFiberorientation.

   \version 1.0.0

   \date Andreas Wachter 01.09.15

   \author Andreas Wachter\n
   Institute of Biomedical Engineering\n
   Kit\n
   http://www.ibt.kit.de\n

   \sa Synopsis \ref TestTetrahedralize
 */

#include "TestTetrahedralize.h"

int main(int argc, char *const argv[]) {
	for (int i = 0; i < argc; i++) {
		if (strcmp(argv[i], "-help") == 0) {

			std::cerr << "\t" << "TestTetrahedralize: This programm can used for testing if a bridge can tetralizied with tetgen ." << "\n";
			std::cerr << "\t" << "<" << "path for input datafile *.vtu | vtp | vtk" << ">" << "\n";
			exit(0);
		}
	}
	try {
		if (argc > 1) {
			std::string fullLoadPath = argv[1];

			std::string path;
			std::string filename;
			std::string fileExtension;
			std::string preFileName = "";

			long extensionPos;
			long pathPos;

			pathPos = fullLoadPath.find_last_of("/");
			path = fullLoadPath.substr(0, pathPos);
			filename = fullLoadPath.substr(pathPos + 1, fullLoadPath.size());

			extensionPos = filename.find_last_of(".");
			fileExtension = filename.substr(extensionPos + 1, filename.size() - extensionPos);

			preFileName = path + "/" + filename.substr(0, extensionPos);

			vtkSmartPointer<vtkPolyData> surface = vtkSmartPointer<vtkPolyData>::New();

			DataFormat data;

			if ((fileExtension == "vtp") || (fileExtension == "vtk")) {
				UniversalReader reader;
				data = reader.read(fullLoadPath);

				surface = Methods::extractSurface(data);
				surface = Methods::triangulateSurface(surface);
			}
			else if (fileExtension == "stl") {
				std::cout << "read stl" << endl;

				vtkSmartPointer<vtkSTLReader> stlReader = vtkSmartPointer<vtkSTLReader>::New();
				stlReader->SetFileName(fullLoadPath.c_str());
				stlReader->Update();
				surface = stlReader->GetOutput();
			}


			vtkSmartPointer<vtkCleanPolyData> cleanFilter = vtkSmartPointer<vtkCleanPolyData>::New();
			cleanFilter->SetInputData(surface);
			cleanFilter->Update();

			surface = cleanFilter->GetOutput();

			surface = Methods::triangulateSurface(surface);

			for (vtkIdType i = 0; i < surface->GetNumberOfCells(); i++) {
				vtkSmartPointer<vtkIdList> pointIds = surface->GetCell(i)->GetPointIds();
				if ((pointIds->GetId(0) == pointIds->GetId(1)) ||
					(pointIds->GetId(0) == pointIds->GetId(2)) ||
					(pointIds->GetId(1) == pointIds->GetId(2))) {
					surface->DeleteCell(i);
				}
				if (pointIds->GetNumberOfIds() < 3) {
					surface->DeleteCell(i);
				}
			}
			surface->RemoveDeletedCells();
			surface->BuildLinks();


			cleanFilter = vtkSmartPointer<vtkCleanPolyData>::New();
			cleanFilter->SetInputData(surface);
			cleanFilter->Update();

			surface = cleanFilter->GetOutput();
			surface->BuildCells();
			surface->BuildLinks();

			cout << "cell before clean: " << surface->GetNumberOfCells() << endl;
			bool exit = false;

			while (!exit) {  // determinate if no cell without an nighbour are existing
				exit = true;

				for (vtkIdType i = 0; i < surface->GetNumberOfCells(); i++) {
					vtkSmartPointer<vtkIdList> neighbours = vtkSmartPointer<vtkIdList>::New();

					for (int j = 0; j < surface->GetCell(i)->GetNumberOfEdges(); j++) {
						vtkSmartPointer<vtkIdList> edgePointIdList = surface->GetCell(i)->GetEdge(j)->GetPointIds();

						vtkSmartPointer<vtkIdList> tempNeighborCellIds = vtkSmartPointer<vtkIdList>::New();
						surface->GetCellNeighbors(i, edgePointIdList, tempNeighborCellIds);

						if (tempNeighborCellIds->GetNumberOfIds() < 1) {
							cout << "cell kein Nachbar: " << i << endl;

							vtkSmartPointer<vtkPoints> pointsofCell = surface->GetCell(i)->GetPoints();
							long numofPoints = pointsofCell->GetNumberOfPoints();
							double point[3] = { 0, 0, 0 };
							for (int i = 0; i < numofPoints; i++) {
								vtkMath::Add(point, pointsofCell->GetPoint(i), point);
							}
							cout << point[0] / numofPoints << ", " << point[1] / numofPoints << ", " << point[2] / numofPoints << endl;

							surface->DeleteCell(i);
							exit = false;
						}
					}
				}
			}
			surface->RemoveDeletedCells();
			surface->BuildCells();
			surface->BuildLinks();


			cleanFilter = vtkSmartPointer<vtkCleanPolyData>::New();
			cleanFilter->SetInputData(surface);
			cleanFilter->Update();
			surface = cleanFilter->GetOutput();

			data.setVtkData(surface);

			cout << "cell after clean: " << surface->GetNumberOfCells() << endl;

			pathPos = fullLoadPath.find_last_of("/");
			path = fullLoadPath.substr(0, pathPos);
			filename = fullLoadPath.substr(pathPos + 1, fullLoadPath.size());

			extensionPos = filename.find_last_of(".");
			fileExtension = filename.substr(extensionPos + 1, filename.size() - extensionPos);

			preFileName = path + "/" + filename.substr(0, extensionPos);

			data.setVtkData(Methods::tetrahedralizeVTPSurface(surface));
		}
		else {

			std::cerr << "\t" << " TestTetrahedralize: This programm can used for testing if a bridge can tetralizied with tetgen ." << "\n\n";
			std::cerr << "\t" << "<" << "path for input datafile *.vtu | vtp | vtk" << ">" << "\n";

		}
	}
	catch (int e) {
		std::cerr << argv[0] << " Error: " << e << std::endl;
		exit(-1);
	}

	return 0;
}  // main

/*! \page TestTetrahedralize

   Test if a bridge surface is been tetralizied with tetgen for include it in SetFiberorientation.

   \section SYNOPSIS_TestTetrahedralize SYNOPSIS
   TestTetrahedralize \<path for input datafile\>\n

   \section OPTIONS_TestTetrahedralize OPTIONS
   \param "<path for input datafile>" The full path for the surface file of the bridge that would be tetrahedralized.
	  The surface can be saved as vtp|stl|vtk.

   \section DESCRIPTION_TestTetrahedralize DESCRIPTION
   To test if a bridge surface is been tetralizied. It clean the surface before it would be tatrahedralized with tetgen.
	  It remove all triangle that has no nighbour. So be carful that the mesh has no hole. The program do the same with
	  the bridge as SetFiberorientation.

   \section SOURCE_TestTetrahedralize SOURCE

   TestTetrahedralize.cpp

   \section SEEALSO_TestTetrahedralize SEE ALSO
   \ref Methods \ref DataFormat \ref Reader \ref Writer

   \section CHANGELOG_TestTetrahedralize CHANGELOG
   V1.0.0 - 01.09.2015 (Andreas Wachter): Starting with changelog\n
 */

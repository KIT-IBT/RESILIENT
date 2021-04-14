/*! \file FindAndMarkSeedPoints.cpp

   \brief Find and mark the seedpoint in the mesh by using the material class.
   \version 1.0.0

   \date Andreas Wachter 01.09.15

   \author Andreas Wachter\n
   Institute of Biomedical Engineering\n
   Kit\n
   http://www.ibt.kit.de\n

   \sa Synopsis \ref FindAndMarkSeedPoints
 */

#include "FindAndMarkSeedPoints.h"

int main(int argc, char *const argv[]) {
	for (int i = 0; i < argc; i++) {
		if (strcmp(argv[i], "-help") == 0) {
			std::cerr << "\t" << "SetAtrialFiberOrientation: This programm can used for setting semiautomatic atrial fieberorientation." << "\n";

			std::cerr << "\t" << "<" << "path for input datafile *.vtu | vtp | vtk" << ">" << "\n";
			std::cerr << "\t" << "<" << "path for output datafile *.vtu | vtp | vtk" << ">" << "\n";
			std::cerr << "\t" << "<" << "path for seedpoints file *.txt" << ">" << "\n";

			exit(0);
		}
	}

	try {
		if (argc > 3) {
			std::string fullLoadPath = argv[1];
			std::string fullStorePath = argv[2];
			std::string pointsPath = argv[3];

			DataFormat data;
			UniversalReader reader;
			data = reader.read(fullLoadPath);

			if (data.getInputType() == DataFormat::vtp) {
				data.setDoubleLayer(false);
			}
			else {
				data.setDoubleLayer(true);
			}

			vtkSmartPointer<vtkDoubleArray> status = vtkSmartPointer<vtkDoubleArray>::New();
			status->SetName("Status");
			status->SetNumberOfComponents(1);
			status->SetNumberOfTuples(data.getCentrePoints()->GetNumberOfPoints());
			status->FillComponent(0, 0);
			data.getVtkData()->GetCellData()->AddArray(status);

			if (!data.getDoubleLayer()) {
				vtkSmartPointer<vtkDoubleArray> materialEndo = vtkSmartPointer<vtkDoubleArray>::New();
				materialEndo->SetName("MaterialEndo");
				materialEndo->SetNumberOfComponents(1);
				materialEndo->SetNumberOfTuples(data.getCentrePoints()->GetNumberOfPoints());

				for (vtkIdType i = 0; i < data.getVtkData()->GetNumberOfCells(); i++) {
					if (data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(i, 0) == Material::Vorhof_links) {
						materialEndo->vtkDataArray::SetComponent(i, 0, Material::Vorhof_links_Endo);
					}
					else {
						materialEndo->vtkDataArray::SetComponent(i, 0, data.getVtkData()->GetCellData()->GetArray(
							"Material")->GetComponent(i, 0));
					}
				}

				vtkSmartPointer<vtkDoubleArray> statusEndo = vtkSmartPointer<vtkDoubleArray>::New();
				statusEndo->SetName("StatusEndo");
				statusEndo->SetNumberOfComponents(1);
				statusEndo->SetNumberOfTuples(data.getCentrePoints()->GetNumberOfPoints());
				statusEndo->FillComponent(0, 0);

				data.getVtkData()->GetCellData()->AddArray(statusEndo);
				data.getVtkData()->GetCellData()->AddArray(materialEndo);
			}

			std::string path;
			std::string filename;
			std::string fileExtension;
			long   extensionPos;
			long   pathPos;

			pathPos = pointsPath.find_last_of("/");
			path = pointsPath.substr(0, pathPos);
			filename = pointsPath.substr(pathPos + 1, pointsPath.size());

			extensionPos = filename.find_last_of(".");
			fileExtension = filename.substr(extensionPos + 1, filename.size() - extensionPos);

			chdir(path.c_str());


			std::string Names[] = { "SCV1", "SCV2", "SCV3", "SCV4", "SCV5", "SCV6", "SCV7", "SCV8", "SCV9" };
			std::string NamesL[] =
			{ "LV1", "LV2", "LV3", "LV4", "LV5", "LV6", "LV7", "LV8", "LV9", "LV10", "LV11", "LV12", "LV13" };
			std::string line = "\n";
			std::string buf = " ";
			std::string CoordsDef = "L";

			int i = 0;
			int j = 0;
			std::ifstream Coordinatefile(filename.c_str());
			double   tempPoint[4] = { 0.0, 0.0, 0.0, 1.0 };

			if (Coordinatefile.is_open()) {
				while (getline(Coordinatefile, line)) {
					if (line.size() > 0) {
						std::vector<std::string> CoordinatesReaded;
						std::stringstream   ss(line);
						while (ss >> buf) {  // reads string until next whitespace
							CoordinatesReaded.push_back(buf);
						}
						if (i < 9) {
							if (!Names[i].compare(CoordinatesReaded.at(0).c_str())) {
								std::cout << "read R" << i << " (" << std::stod(CoordinatesReaded.at(2).c_str()) << "," << std::stod(CoordinatesReaded.at(3).c_str()) << "," << std::stod(CoordinatesReaded.at(4).c_str()) << ") (global)" << std::endl;
								std::vector<double> point = { std::stod(CoordinatesReaded.at(2).c_str()), std::stod(CoordinatesReaded.at(3).c_str()), std::stod(CoordinatesReaded.at(4).c_str()) };
								point = Methods::findClosedPointinMaterial(data, point, Material::Vorhof_rechts);
								vtkSmartPointer<vtkIdList> path = Methods::getCellsinRadius(data, data.getCentrePoints()->FindPoint(point.at(0), point.at(1), point.at(2)), 1);
								Methods::setMaterial(data, path, Material::testMaterialRight);
							}
						}
						if (i >= 9) {
							if (!NamesL[j].compare(CoordinatesReaded.at(0).c_str())) {
								std::cout << "read L" << j << " (" << std::stod(CoordinatesReaded.at(2).c_str()) << "," << std::stod(CoordinatesReaded.at(3).c_str()) << "," << std::stod(CoordinatesReaded.at(4).c_str()) << ") (global)" << std::endl;
								std::vector<double> pointorg = { std::stod(CoordinatesReaded.at(2).c_str()), std::stod(CoordinatesReaded.at(3).c_str()), std::stod(CoordinatesReaded.at(4).c_str()) };
								std::vector<double> point = Methods::findClosedPointinMaterial(data, pointorg, Material::Vorhof_links);
								vtkSmartPointer<vtkIdList> path = Methods::getCellsinRadius(data, data.getCentrePoints()->FindPoint(point.at(0), point.at(1), point.at(2)), 1);
								Methods::setMaterial(data, path, Material::testMaterialLeft);

								point = Methods::findClosedPointinMaterial(data, pointorg, Material::Vorhof_links_Endo);
								vtkSmartPointer<vtkIdList> pathendo = Methods::getCellsinRadius(data, data.getCentrePoints()->FindPoint(point.at(0), point.at(1), point.at(2)), 1);
								Methods::setMaterial(data, pathendo, Material::testMaterialLeftEndo);
							}
						}
						j++;
					}
					i++;
				}

				Coordinatefile.close();

				Methods::deleteArrays(data);

				UniversalWriter writer;
				writer.write(data, fullStorePath);
			}
			else {
				std::cerr << "Unable to open file coordinationfile\n";
				exit(-1);
			}
		}
		else {
			std::cerr << "\t" << " FindAndMarkSeedPoints: Mark the localization of the seedpoint in the mesh by using the material class." << "\n";

			std::cerr << "\t" << "<" << "path for input datafile *.vtu | vtp | vtk" << ">" << "\n";
			std::cerr << "\t" << "<" << "path for output datafile *.vtu|vtp|vtk" << ">" << "\n";
			std::cerr << "\t" << "<" << "path for seedpoints file *.txt" << ">" << "\n";

		}
	}
	catch (int e) {
		std::cerr << argv[0] << " Error: " << e << std::endl;
		exit(-1);
	}

	return 0;
}  // main

/*! \page FindAndMarkSeedPoints

   Mark the seedpoint in the mesh by using the material class.
   \section SYNOPSIS_FindAndMarkSeedPoints SYNOPSIS
   FindAndMarkSeedPoints \<path for input datafile\> \<path for output datafile\> \<path for seedpoint datafile\> \n

   \section OPTIONS_FindAndMarkSeedPoints OPTIONS
   \param "<path for input datafile>" The full path for the input datafile. The filetype can be vtp|vtu|vtk.
   \param "<path for output datafile>" The full path for the output datafile. The filetype can be vtp|vtu|vtk.
   \param "<path for seedpoints datafile>" The full path for the seedpoint txt-datafile.

   \section DESCRIPTION_FindAndMarkSeedPoints DESCRIPTION
   Mark the seedpoint in the mesh by using the material class 500 in right atria, 501 in left epi atria and 502 in left
	  endo atria. The marked point radius is 1 mm.

   \section SOURCE_FindAndMarkSeedPoints SOURCE

   FindAndMarkSeedPoints.cpp

   \section SEEALSO_FindAndMarkSeedPoints SEE ALSO
   \ref Methods \ref DataFormat \ref Reader \ref Writer

   \section CHANGELOG_FindAndMarkSeedPoints CHANGELOG
   V1.0.0 - 01.09.2015 (Andreas Wachter): Starting with changelog\n
 */

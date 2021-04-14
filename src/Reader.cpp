/*! \file Reader.cpp

   \brief  For reading in the vtk-, vtu- and vtp- files.

   \version 1.0.0

   \date Andreas Wachter 01.09.15

   \author Andreas Wachter\n
   Institute of Biomedical Engineering\n
   Kit\n
   http://www.ibt.kit.de\n

   \sa Synopsis \ref Reader
 */

#include "Reader.h"

using namespace std;


Reader::Reader() {}

/*!Constructor for UniversalReader*/
UniversalReader::UniversalReader() {}

/*! Read in vtu, vtp or vtk data from a the path. It write out some informations.
   \param FullPath Path of the mesh that would be loaded.
   \return A DataFormat with the mesh data.
 */
DataFormat UniversalReader::read(string FullPath) {
	return read(FullPath, true);
}

/*! Read in vtu, vtp or vtk data from a the path. It write out some more informations during the prozess.
   \param FullPath Path of the mesh that would be loaded.
   \param verbose True to write some information out, False for not.
   \return A DataFormat with the mesh data.
 */
DataFormat UniversalReader::read(string FullPath, bool verbose) {
	string fullPath = FullPath;
	string path;
	string filename;
	string fileExtension;
	long   extensionPos;
	long   pathPos;
	int inputType = 0;

	pathPos = fullPath.find_last_of("/");
	path = fullPath.substr(0, pathPos);
	filename = fullPath.substr(pathPos + 1, fullPath.size());

	extensionPos = filename.find_last_of(".");
	fileExtension = filename.substr(extensionPos + 1, filename.size() - extensionPos);

	chdir(path.c_str());
	if (verbose) {
		cout << path << endl;
		cout << filename << endl;
	}
	if (fileExtension == "vtu")
		inputType = 2;
	if (fileExtension == "vtp")
		inputType = 3;
	if (fileExtension == "vtk")
		inputType = 4;

	DataFormat data;

	switch (inputType) {
	case 2: {
		if (verbose) {
			cout << "File extension: " + fileExtension << endl;
		}
		VtuReader vtuReader;
		data = vtuReader.read(filename);
		break;
	}
	case 3: {
		if (verbose) {
			cout << "File extension: " + fileExtension << endl;
		}
		VtpReader vtpReader;
		data = vtpReader.read(filename);
		break;
	}
	case 4: {
		if (verbose) {
			cout << "File extension: " + fileExtension << endl;
		}
		VtkReader vtkReader;
		data = vtkReader.read(filename);
		break;
	}
	default: {
		cerr << "File extension: " + fileExtension + " undefine" << endl;
		exit(1);
	}
	}  // switch

	return data;
}  // UniversalReader::read


// methods for VtuReader class

/*! Constructor for vtu file reader*/

VtuReader::VtuReader() {
	vtuReaderPtr = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
}

/*!For read in a vtu data file.
   \param filename Path of the vtu mesh that would be loaded.
   \return A DataFormat with the mesh data.
 */
DataFormat VtuReader::read(std::string filename) {
	DataFormat tempData;

	tempData.setInputType(DataFormat::PossibleInputType::vtu);

	vtuReaderPtr->SetFileName(filename.c_str());
	vtuReaderPtr->Update();
	tempData.setVtkData(vtuReaderPtr->GetOutput());

	vtkSmartPointer<vtkPoints> pointsofCell;

	// Centrepoints
	tempData.setCentrePoints(Methods::calcCellCentroids(tempData.getVtkData()));


	return tempData;
}

// methods for VtpReader class

/*! Constructor for vtp file reader*/

VtpReader::VtpReader() {
	vtpReaderPtr = vtkSmartPointer<vtkXMLPolyDataReader>::New();
}

/*!For read in a vtp data file.
   \param filename Path of the vtp mesh that would be loaded.
   \return A DataFormat with the mesh data.
 */
DataFormat VtpReader::read(std::string filename) {
	DataFormat tempData;

	tempData.setInputType(DataFormat::PossibleInputType::vtp);

	vtpReaderPtr->SetFileName(filename.c_str());
	vtpReaderPtr->Update();
	tempData.setVtkData(vtpReaderPtr->GetOutput());

	vtkSmartPointer<vtkPoints> pointsofCell;

	// Centrepoints
	tempData.setCentrePoints(Methods::calcCellCentroids(tempData.getVtkData()));

	return tempData;
}

// methods for VtkReader class

/*! Constructor for vtk file reader*/

VtkReader::VtkReader() {
	vtkReaderPtr = vtkSmartPointer<vtkDataSetReader>::New();
}

/*!For read in a vtk data file.
   \param filename Path of the vtk mesh that would be loaded.
   \return A DataFormat with the mesh data.
 */
DataFormat VtkReader::read(std::string filename) {
	vtkReaderPtr->SetFileName(filename.c_str());
	vtkReaderPtr->Update();

	DataFormat tempData;

	if (vtkReaderPtr->IsFilePolyData()) {
		tempData.setInputType(DataFormat::PossibleInputType::vtp);
		tempData.setVtkData(vtkReaderPtr->GetPolyDataOutput());
	}
	else if (vtkReaderPtr->IsFileUnstructuredGrid()) {
		tempData.setInputType(DataFormat::PossibleInputType::vtu);
		tempData.setVtkData(vtkReaderPtr->GetUnstructuredGridOutput());
	}
	else {
		cerr << "reader the vtk internal dataformat is not supported";
		exit(-1);
	}

	// Centrepoints
	tempData.setCentrePoints(Methods::calcCellCentroids(tempData.getVtkData()));

	return tempData;
}

/*!
   \page Reader

   \section DESCRIPTION_Reader DESCRIPTION
   For reading in the vtk-, vtu- and vtp- files.

   \section SOURCE_Reader SOURCE

   Reader.cpp

   \section SEEALSO_Reader SEE ALSO
   \ref DataFormat

   \section CHANGELOG_Reader CHANGELOG
   V1.0.0 - 01.09.2015 (Andreas Wachter): Starting with changelog\n
 */

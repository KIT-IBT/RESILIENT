/*! \file Writer.cpp

   \brief  For writing the data as a vtk-, vtu- or vtp- files.

   \version 1.0.0

   \date Andreas Wachter 01.09.15

   \author Andreas Wachter\n
   Institute of Biomedical Engineering\n
   Kit\n
   http://www.ibt.kit.de\n

   \sa Synopsis \ref Writer
 */

#include "Writer.h"

using namespace std;

Writer::Writer() {}

/*!Constructor for UniversalWriter*/
UniversalWriter::UniversalWriter() {}

/*! Write the data as vtu, vtp or vtk file at the path. It write out some more informations during the prozess.
   \param FullPath Path of the mesh that would be loaded.
   \param A Pointer with the mesh data.
 */
void UniversalWriter::write(DataFormat &data, string FullPath) {
	write(data, FullPath, true);
}

/*! Write the data as vtu, vtp or vtk file at the path. It write out some more informations during the prozess,
   if you want.
   \param data A Pointer with the mesh data.
   \param FullPath Path where the mesh would be stored.
   \param verbose True to write some information out, False for not.
 */
void UniversalWriter::write(DataFormat &data, string FullPath, bool verbose) {
	string fullPath = FullPath;
	string path;
	string filename;
	string fileExtension;
	long   extensionPos;
	long   pathPos;
	int outputType = 0;

	pathPos = fullPath.find_last_of("/");
	path = fullPath.substr(0, pathPos);
	filename = fullPath.substr(pathPos + 1, fullPath.size());

	extensionPos = filename.find_last_of(".");
	fileExtension = filename.substr(extensionPos + 1, filename.size() - extensionPos);

	chdir(path.c_str());
	if (verbose) {
		cout << filename << endl;
	}
	ifstream FileTest(filename.c_str());
	if (FileTest) {
		FileTest.close();
		if (verbose) {
			cerr << "File exists and will be deleted" << endl;
		}
		remove(filename.c_str());
	}


	if (fileExtension == "vtu")
		outputType = 2;
	if (fileExtension == "vtp")
		outputType = 3;
	if (fileExtension == "vtk")
		outputType = 4;

	switch (outputType) {
	case 2: {
		if (verbose) {
			cout << "File extension: " + fileExtension << endl;
		}
		VtuWriter vtuWriter;
		vtuWriter.write(data, filename);
		break;
	}
	case 3: {
		if (verbose) {
			cout << "File extension: " + fileExtension << endl;
		}
		VtpWriter vtpWriter;
		vtpWriter.write(data, filename);
		break;
	}
	case 4: {
		if (verbose) {
			cout << "File extension: " + fileExtension << endl;
		}
		VtkWriter vtkWriter;
		vtkWriter.write(data, filename);
		break;
	}

	default: {
		cerr << "File extension: " + fileExtension + " undefine" << endl;
		exit(1);
	}
	}  // switch
}  // UniversalWriter::write

// methods for VtuWriter class

/*! Constructor for VtuWriter*/
VtuWriter::VtuWriter() {
	vtuWriterPtr = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
}

/*! Write the mesh data in a vtu unstrucstured Grid file.
   \param writeData A Pointer with the mesh data.
   \param filename Filename where the mesh would be stored.
 */
void VtuWriter::write(DataFormat &writeData, string filename) {
	if (writeData.getInputType() == 2) {
		vtuWriterPtr->SetFileName(filename.c_str());
		vtuWriterPtr->SetInputData(writeData.getVtkData());
		vtuWriterPtr->Write();

		if (ConfigFiberorientation::writeViewFiberOrientation) {
			long extensionPos = filename.find_last_of(".");
			string filebase = filename.substr(0, extensionPos);

			vtkDataArray *diffVec = writeData.getVtkData()->GetCellData()->GetArray("DifferenceVector");
			writeData.getCentrePoints()->GetPointData()->AddArray(diffVec);
			vtkDataArray *mat = writeData.getVtkData()->GetCellData()->GetArray("Material");
			writeData.getCentrePoints()->GetPointData()->AddArray(mat);

			filename = filebase + "_orginalOrientation.vtu";
			vtkSmartPointer<vtkXMLUnstructuredGridWriter> vtuWriterPtr = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
			vtuWriterPtr->SetFileName(filename.c_str());
			vtuWriterPtr->SetInputData(writeData.getCentrePoints());

			vtuWriterPtr->Write();
		}
	}
	else {
		cerr << "Conversion from: " << writeData.getInputType() << " to: vtu are not supported" << endl;
		exit(1);
	}
}  // VtuWriter::write

// methods for VtpWriter class

/*! Constructor for VtpWriter*/
VtpWriter::VtpWriter() {
	vtpWriterPtr = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
}

/*! Write the mesh data in a vtp file.
   \param writeData A Pointer with the mesh data.
   \param filename Filename where the mesh would be stored.
 */
void VtpWriter::write(DataFormat &writeData, string filename) {
	if (writeData.getInputType() == 3) {
		vtpWriterPtr->SetFileName(filename.c_str());
		vtpWriterPtr->SetInputData(writeData.getVtkData());
		vtpWriterPtr->Write();

		if (ConfigFiberorientation::writeViewFiberOrientation) {
			long extensionPos = filename.find_last_of(".");
			string filebase = filename.substr(0, extensionPos);

			vtkDataArray *diffVec = writeData.getVtkData()->GetCellData()->GetArray("DifferenceVector");
			writeData.getCentrePoints()->GetPointData()->AddArray(diffVec);
			vtkDataArray *mat = writeData.getVtkData()->GetCellData()->GetArray("Material");
			writeData.getCentrePoints()->GetPointData()->AddArray(mat);
			if (!writeData.getDoubleLayer()) {
				vtkDataArray *diffVecEndo = writeData.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo");
				writeData.getCentrePoints()->GetPointData()->AddArray(diffVecEndo);
				vtkDataArray *matEndo = writeData.getVtkData()->GetCellData()->GetArray("MaterialEndo");
				writeData.getCentrePoints()->GetPointData()->AddArray(matEndo);
			}
			filename = filebase + "_orginalOrientation.vtu";
			vtkSmartPointer<vtkXMLUnstructuredGridWriter> vtuWriterPtr = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
			vtuWriterPtr->SetFileName(filename.c_str());
			vtuWriterPtr->SetInputData(writeData.getCentrePoints());
			vtuWriterPtr->Write();
		}
	}
	else {    // if it not a surface mesh then the mesh would be convert
		cerr << "Conversion from: " << writeData.getInputType() << " to: vtp" << endl;

		vtkSmartPointer<vtkPolyData> polydata = Methods::extractSurface(writeData);
		polydata = Methods::triangulateSurface(polydata);
		polydata = Methods::smoothSurface(polydata, 1000);

		vtpWriterPtr->SetFileName(filename.c_str());
		vtpWriterPtr->SetInputData(polydata);
		vtpWriterPtr->Write();
	}
}  // VtpWriter::write

// methods for VtkWriter class

/*! Constructor for VtkWriter*/
VtkWriter::VtkWriter() {
	vtkWriterPtr = vtkSmartPointer<vtkDataSetWriter>::New();
}

/*! Write the mesh data in a vtk file.
   \param writeData A Pointer with the mesh data.
   \param filename Filename where the mesh would be stored.

 */
void VtkWriter::write(DataFormat &writeData, string filename) {
	vtkWriterPtr->SetFileName(filename.c_str());
	vtkWriterPtr->SetInputData(writeData.getVtkData());
	//vtkWriterPtr->SetFileTypeToBinary();
	vtkWriterPtr->Write();

	if (ConfigFiberorientation::writeViewFiberOrientation) {
		long extensionPos = filename.find_last_of(".");
		string filebase = filename.substr(0, extensionPos);

		vtkDataArray *diffVec = writeData.getVtkData()->GetCellData()->GetArray("DifferenceVector");
		writeData.getCentrePoints()->GetPointData()->AddArray(diffVec);
		vtkDataArray *mat = writeData.getVtkData()->GetCellData()->GetArray("Material");
		writeData.getCentrePoints()->GetPointData()->AddArray(mat);
		if (!writeData.getDoubleLayer()) {
			vtkDataArray *diffVecEndo = writeData.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo");
			writeData.getCentrePoints()->GetPointData()->AddArray(diffVecEndo);
			vtkDataArray *matEndo = writeData.getVtkData()->GetCellData()->GetArray("MaterialEndo");
			writeData.getCentrePoints()->GetPointData()->AddArray(matEndo);
		}
		filename = filebase + "_orginalOrientation.vtu";
		vtkSmartPointer<vtkXMLUnstructuredGridWriter> vtuWriterPtr = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
		vtuWriterPtr->SetFileName(filename.c_str());
		vtuWriterPtr->SetInputData(writeData.getCentrePoints());
		vtuWriterPtr->Write();
	}
}  // VtkWriter::write

/*!
   \page Writer

   \section DESCRIPTION_Writer DESCRIPTION
   For writing the mesh data as a vtk-, vtu- or vtp- files.

   \section SOURCE_Writer SOURCE

   Writer.cpp

   \section SEEALSO_Writer SEE ALSO
   \ref DataFormat

   \section CHANGELOG_Writer CHANGELOG
   V1.0.0 - 01.09.2015 (Andreas Wachter): Starting with changelog\n
 */

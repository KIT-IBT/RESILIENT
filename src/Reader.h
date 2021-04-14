/*! For reading in the vtk-, vtu- and vtp- files.*/

//
//  Reader.h
//  RESILIENT
//
//  Created by Andreas Wachter on 01.09.15.
//  Copyright (c) 2015 IBT. All rights reserved.
//

#ifndef _RESILIENT_Reader_
#define _RESILIENT_Reader_

#include "Headers.h"
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkDataSetReader.h>
#include "Methods.h"


// abstract class for a Reader:
class Reader {
public:
	Reader();  // constructor
	/*!For abstract read function.
	   \param filename Path of the mesh that would be loaded.
	 */
	virtual DataFormat read(std::string) = 0;  // pure virtual function to read in the data; arguments, filename
};


// Universal Reader
class UniversalReader : public Reader {
public:
	UniversalReader();
	DataFormat            read(std::string, bool verbose);
	DataFormat            read(std::string);
};


// VTU Reader
class VtuReader : private Reader {
private:
	vtkSmartPointer<vtkXMLUnstructuredGridReader> vtuReaderPtr;

public:
	VtuReader();  // constructor
	DataFormat read(std::string);
};

// VTP Reader
class VtpReader : private Reader {
private:
	vtkSmartPointer<vtkXMLPolyDataReader> vtpReaderPtr;

public:
	VtpReader();  // constructor
	DataFormat read(std::string);
};

// VTK Reader
class VtkReader : private Reader {
private:
	vtkSmartPointer<vtkDataSetReader> vtkReaderPtr;

public:
	VtkReader();  // constructor
	DataFormat read(std::string);
};


#endif /* defined(_RESILIENT_Reader_) */

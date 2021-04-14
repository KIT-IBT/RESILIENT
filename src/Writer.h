//
//  Writer.h
//  RESILIENT
//
//  Created by Andreas Wachter on 08.08.14.
//  Copyright (c) 2014 IBT. All rights reserved.
//

#ifndef _RESILIENT_Writer_
#define _RESILIENT_Writer_

#include "Headers.h"
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkDataSetWriter.h>
#include "Methods.h"
#include <vtkTransform.h>
#include <vtkTransformFilter.h>


// abstract class for a Writer:
class Writer {
public:
	Writer();                                              // constructor
	virtual void write(DataFormat &data, std::string) = 0;  // pure virtual function to wirte in the data; arguments,
															// filename
};

// Universal Writer
class UniversalWriter : public Writer {
public:
	UniversalWriter();
	void write(DataFormat &data, std::string);
	void write(DataFormat &data, std::string, bool verbose);
};

// Vtu Writer
class VtuWriter : public Writer {
private:
	vtkSmartPointer<vtkXMLUnstructuredGridWriter> vtuWriterPtr;

public:
	VtuWriter();
	void write(DataFormat &data, std::string);
};

// Vtp Writer
class VtpWriter : public Writer {
private:
	vtkSmartPointer<vtkXMLPolyDataWriter> vtpWriterPtr;

public:
	VtpWriter();
	void write(DataFormat &data, std::string);
};

// Vtk Writer
class VtkWriter : public Writer {
private:
	vtkSmartPointer<vtkDataSetWriter> vtkWriterPtr;

public:
	VtkWriter();
	void write(DataFormat &data, std::string);
};

#endif /* defined(_RESILIENT_Writer_) */

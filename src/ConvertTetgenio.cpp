/*! \file ConvertTetgenio

   \brief For converting an vtk polydata in a tetgen format and a tetgen mesh in an vtk unstrucerd grid.
   \version 1.0.0

   \date Andreas Wachter 01.09.15

   \author Andreas Wachter\n
   Institute of Biomedical Engineering\n
   Kit\n
   http://www.ibt.kit.de\n

   \sa Synopsis \ref ConvertTetgenio
 */

#include "ConvertTetgenio.h"

 /*!To convert a polydata mesh in the tetgen format.
	\param poly Pointer to the polydata mesh.
	\return The mesh in tegen format.
  */
tetgenio ConvertTetgenio::convertVTPtoTetgenio(vtkSmartPointer<vtkPolyData> &poly) {
	tetgenio tetgenData = tetgenio();

	tetgenio::facet *f;
	tetgenio::polygon *p;

	// All indices start from 0.
	tetgenData.firstnumber = 0;

	tetgenData.numberofpoints = (int)poly->GetNumberOfPoints();
	tetgenData.pointlist = new REAL[tetgenData.numberofpoints * 3];

	for (int i = 0; i < tetgenData.numberofpoints; i++) {
		tetgenData.pointlist[i * 3] = poly->GetPoint(i)[0];
		tetgenData.pointlist[i * 3 + 1] = poly->GetPoint(i)[1];
		tetgenData.pointlist[i * 3 + 2] = poly->GetPoint(i)[2];
	}

	tetgenData.numberoffacets = (int)poly->GetNumberOfCells();

	tetgenData.facetlist = new tetgenio::facet[tetgenData.numberoffacets];
	tetgenData.facetmarkerlist = new int[tetgenData.numberoffacets];


	tetgenData.numberoffacets = (int)poly->GetNumberOfCells();

	tetgenData.facetlist = new tetgenio::facet[tetgenData.numberoffacets];
	tetgenData.facetmarkerlist = new int[tetgenData.numberoffacets];

	for (int j = 0; j < tetgenData.numberoffacets; j++) {
		f = &tetgenData.facetlist[j];
		f->numberofpolygons = 1;
		f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
		f->numberofholes = 0;
		f->holelist = NULL;
		p = &f->polygonlist[0];
		p->numberofvertices = 3;
		p->vertexlist = new int[p->numberofvertices];
		p->vertexlist[0] = (int)poly->GetCell(j)->GetPointId(0);
		p->vertexlist[1] = (int)poly->GetCell(j)->GetPointId(1);
		p->vertexlist[2] = (int)poly->GetCell(j)->GetPointId(2);
	}

	return tetgenData;
}  // ConvertTetgenio::convertVTPtoTetgenio

/*!To convert a tetgen mesh in a vtk unstructured grid format.
   \param tetgenData Pointer to the tegen mesh.
   \return The mesh in vtkUnstructuredGrid format.
 */
vtkSmartPointer<vtkUnstructuredGrid> ConvertTetgenio::convertTetgeniotoVTU(tetgenio &tetgenData) {
	vtkSmartPointer<vtkUnstructuredGrid> ugrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();


	for (int i = 0; i < tetgenData.numberofpoints; i++) {
		points->InsertNextPoint(tetgenData.pointlist[i * 3], tetgenData.pointlist[i * 3 + 1], tetgenData.pointlist[i * 3 + 2]);
	}

	ugrid->SetPoints(points);

	for (int j = 0; j < tetgenData.numberoftetrahedra; j++) {
		vtkSmartPointer<vtkIdList> pointsList = vtkSmartPointer<vtkIdList>::New();
		for (int k = 0; k < 4; k++) {
			int id = tetgenData.tetrahedronlist[(j * 4) + k];
			pointsList->InsertNextId(id);
		}
		ugrid->InsertNextCell(10, pointsList);
	}

	return ugrid;
}


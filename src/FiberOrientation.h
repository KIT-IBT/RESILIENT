/*!Definition file for the annotation of the fiber orientation and insertion of the automatic bridges.*/

//
// FiberOrientation.h
// RESILIENT
//
// Created by Andreas Wachter on 01.09.14.
// Copyright (c) 2015 IBT. All rights reserved.
//

#ifndef _RESILIENT_FiberOrientation_
#define _RESILIENT_FiberOrientation_

#include "Headers.h"
#include "Methods.h"



class FiberOrientation {
 private:
 vtkSmartPointer<vtkPoints> BPoints;

 vtkSmartPointer<vtkPoints> RPoints;
 vtkSmartPointer<vtkPoints> LPointsEpi;
 vtkSmartPointer<vtkPoints> LPointsEndo;

 vtkSmartPointer<vtkPoints> FreeBridgePoints;

 vtkSmartPointer<vtkPoints> PectiStartPoints;
 vtkSmartPointer<vtkPoints> PectiTargetPoints;

 std::vector<double> PectinateCentre;
 std::map<std::string, vtkSmartPointer<vtkIdList>> Path;

 double Resolution;

 public:
 std::vector<int> freeBridgeMaterial;
 std::vector<double> freeBridgeWidth;

 // constructor
 FiberOrientation();
 static void calculateFiberOriented(DataFormat &data, std::string loadpath, std::string filename, std::string freeBridgePath);
 void readOutFiberPoints(std::string filename);
 int readOutFreeBridgePointsPoints(std::string fullPath);

 void setPectStartPoint(vtkIdType position, std::vector<double> point);
 std::vector<double> getPectStartPoint(vtkIdType position);
 void setPectTargetPoint(vtkIdType position, std::vector<double> point);
 std::vector<double> getPectTargetPoint(vtkIdType position);
 void setPath(std::string pathName, vtkSmartPointer<vtkIdList> path);
 vtkSmartPointer<vtkIdList> getPath(std::string pathName);
 void setRPoint(vtkIdType position, double x, double y, double z);
 void setRPoint(vtkIdType position, std::vector<double> point);
 std::vector<double> getRPoint(vtkIdType position);
 void setRPoints(vtkSmartPointer<vtkPoints> rPoints);
 vtkSmartPointer<vtkPoints> getRPoints();
 void setBPoint(vtkIdType position, double x, double y, double z);
 void setBPoint(vtkIdType position, std::vector<double> point);
 std::vector<double> getBPoint(vtkIdType position);
 void setLPointEpi(vtkIdType position, double x, double y, double z);
 void setLPointEpi(vtkIdType position, std::vector<double> point);
 std::vector<double> getLPointEpi(vtkIdType position);
 void setLPointsEpi(vtkSmartPointer<vtkPoints> lPointsEpi);
 vtkSmartPointer<vtkPoints> getLPointsEpi();
 void setLPointEndo(vtkIdType position, double x, double y, double z);
 void setLPointEndo(vtkIdType position, std::vector<double> point);
 std::vector<double> getLPointEndo(vtkIdType position);
 void setLPointsEndo(vtkSmartPointer<vtkPoints> lPointsEndo);
 vtkSmartPointer<vtkPoints> getLPointsEndo();
 void setFreeBridgePoint(vtkIdType position, double x, double y, double z);
 void setFreeBridgePoint(vtkIdType position, std::vector<double> point);
 std::vector<double> getFreeBridgePoint(vtkIdType position);
 void setFiberOrientationRightAtrium(DataFormat &data);
 void setFiberOrientationLeftEndoAtrium(DataFormat &data);
 void setFiberOrientationLeftEpiAtrium(DataFormat &data);
 void setFiberOrientationLeftAppendes(DataFormat &data);
 void setFiberOrientationLeftVeins(DataFormat &data);
 void setBridges(DataFormat &data, std::string freeBridgePath);
 void setSinusNode(DataFormat &data);
 void closeOrientation(DataFormat &data);
 void setBridge(DataFormat &data, std::vector<double> pointLeft, std::vector<double> pointRight, int bridgenumber, double bridgewidth, double bridgeLengthinRight, double bridgeLengthinLeft, Material::Mat bridgeMaterialLeft, Material::Mat bridgeMaterialRight, double searchRadius);
 void setBachmannBridgelong(DataFormat &data, std::vector<double> pointLeft, std::vector<double> pointRight, int bridgenumber, double bridgewidth, double searchRadius, Material::Mat bridgeMaterialLeft, Material::Mat bridgeMaterialRight);
 void initPoints(vtkSmartPointer<vtkPoints> &points);
 static vtkSmartPointer<vtkPoints> testPointsandMove(DataFormat &data, vtkSmartPointer<vtkPoints> points,  Material::Mat inMaterial);
 DataFormat update(DataFormat data);
 void setResolution(double resolution);
 double getResolution();
 void smoothOrientationToSurface(DataFormat &data, DataFormat &dataForSmooth);
 void setMaterialtoOldMaterial(DataFormat &data);
 void setFreeBridge(DataFormat &data, std::vector<double> pointLeft, std::vector<double> pointRight,  int bridgeNumber, double bridgewidth, int bridgeMaterialLeft);
 void closeOrientationFreeBridge(DataFormat &data, int freeBridgeMaterial);
}; // class FiberOrientation


#endif /* defined(_RESILIENT_FiberOrientation_) */

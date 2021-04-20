/*! \file FiberOrientation.cpp

 \brief Definition file for the annotation of the fiber orientation and insertion of the automatic bridges.
 \version 1.0.0

 \date Andreas Wachter 01.09.15

 \author Andreas Wachter\n
 Institute of Biomedical Engineering\n
 Kit\n
 http://www.ibt.kit.de\n

 \sa Synopsis \ref FiberOrientation
 */


#include "FiberOrientation.h"
#include "Reader.h"

using namespace std;

/*! Constructor for the FiberOrientation. It contains all informations (orginal seed point, additional point, brigde
 points, free bridge points, pectinate points, mesh resolution and fiber paths). */

FiberOrientation::FiberOrientation() {
	RPoints = vtkSmartPointer<vtkPoints>::New();
	RPoints->SetNumberOfPoints(45);
	initPoints(RPoints);

	BPoints = vtkSmartPointer<vtkPoints>::New();
	BPoints->SetNumberOfPoints(35);
	initPoints(BPoints);

	LPointsEndo = vtkSmartPointer<vtkPoints>::New();
	LPointsEndo->SetNumberOfPoints(60);
	initPoints(LPointsEndo);

	LPointsEpi = vtkSmartPointer<vtkPoints>::New();
	LPointsEpi->SetNumberOfPoints(60);
	initPoints(LPointsEpi);

	FreeBridgePoints = vtkSmartPointer<vtkPoints>::New();
	FreeBridgePoints->SetNumberOfPoints(30);
	initPoints(FreeBridgePoints);

	PectiStartPoints = vtkSmartPointer<vtkPoints>::New();
	PectiStartPoints->SetNumberOfPoints(15);
	initPoints(PectiStartPoints);

	PectiTargetPoints = vtkSmartPointer<vtkPoints>::New();
	PectiTargetPoints->SetNumberOfPoints(15);
	initPoints(PectiTargetPoints);

	Path = map<string, vtkSmartPointer<vtkIdList>>();

	// paths
	PectinateCentre = vector<double>(3);

	Resolution = -1;
}

/*! Initalize all points with not a number (NAN).
 \param points pointer for the vtkPoints
 */
void FiberOrientation::initPoints(vtkSmartPointer<vtkPoints> &points) {
	for (vtkIdType i = 0; i < points->GetNumberOfPoints(); i++) {
		points->SetPoint(i, NAN, NAN, NAN);
	}
}

/*! Set function for path.
 \param pathName unique pathname for a calculated fiber path
 \param path vtkIdList with the path point ids
 */
void FiberOrientation::setPath(string pathName, vtkSmartPointer<vtkIdList> path) {
	Path[pathName] = path;
}

/*! Get function for path.
 \param pathName unique pathname for a calculated fiber path
 */
vtkSmartPointer<vtkIdList> FiberOrientation::getPath(string pathName) {
	return Path[pathName];
}

/*! Get function for a pectinate start point.
 \param position of the pectinate start point
 \return the pectinate start point at position as vector with the coordinates (x,y,z). If it exist, if not then it
 write out "Point don't exixts!".
 */
vector<double> FiberOrientation::getPectStartPoint(vtkIdType position) {
	vector<double> point(3);

	point.at(0) = PectiStartPoints->GetPoint(position)[0];
	point.at(1) = PectiStartPoints->GetPoint(position)[1];
	point.at(2) = PectiStartPoints->GetPoint(position)[2];

	if (::isnan(point.at(0)) && ::isnan(point.at(1)) && ::isnan(point.at(2)))
		cerr << "Point don't exixts!" << endl;

	return point;
}

/*! Get function for a pectinate target point.
 \param number of the pectinate target point
 \return the pectinate target point at number as vector with the coordinates (x,y,z). If it exist, if not then it
 write out "Point don't exixts!".
 */
vector<double> FiberOrientation::getPectTargetPoint(vtkIdType position) {
	vector<double> point(3);

	point.at(0) = PectiTargetPoints->GetPoint(position)[0];
	point.at(1) = PectiTargetPoints->GetPoint(position)[1];
	point.at(2) = PectiTargetPoints->GetPoint(position)[2];

	if (::isnan(point.at(0)) && ::isnan(point.at(1)) && ::isnan(point.at(2)))
		cerr << "Point don't exixts!" << endl;

	return point;
}

/*! Set function for a pectinate start point.
 \param position of the pectinate start point
 \param point vector with the coordinates (x,y,z) of the pectinate start point
 */
void FiberOrientation::setPectStartPoint(vtkIdType position, vector<double> point) {
	double pointArray[3] = { point.at(0), point.at(1), point.at(2) };

	PectiStartPoints->SetPoint(position, pointArray);
}

/*! Set function for a pectinate target point.
 \param position of the pectinate target point
 \param point vector with the coordinates (x,y,z) of the pectinate target point
 */
void FiberOrientation::setPectTargetPoint(vtkIdType position, vector<double> point) {
	double pointArray[3] = { point.at(0), point.at(1), point.at(2) };

	PectiTargetPoints->SetPoint(position, pointArray);
}

/*! Set function for all the right seed points.
 \param rPoints the right seed points in vtkPoints format
 */
void FiberOrientation::setRPoints(vtkSmartPointer<vtkPoints> rPoints) {
	RPoints = rPoints;
}

/*! Set function for a right seed point.
 \param position of the right seed point
 \param x x coordinate of the right seed point
 \param y y coordinate of the right seed point
 \param z z coordinate of the right seed point
 */
void FiberOrientation::setRPoint(vtkIdType position, double x, double y, double z) {
	RPoints->SetPoint(position, x, y, z);
}

/*! Set function for a right seed point.
 \param position of the right seed point
 \param point vector with the coordinates (x,y,z) of the right seed point
 */
void FiberOrientation::setRPoint(vtkIdType position, vector<double> point) {
	double x = point.at(0);
	double y = point.at(1);
	double z = point.at(2);

	setRPoint(position, x, y, z);
}

/*! Get function for a right seed point.
 \param position of the right seed point
 \return the right seed point at position as vector with the coordinates (x,y,z). If it exist, if not then it write
 out "Point don't exixts!".
 */
vector<double> FiberOrientation::getRPoint(vtkIdType position) {
	vector<double> point(3);

	point.at(0) = RPoints->GetPoint(position)[0];
	point.at(1) = RPoints->GetPoint(position)[1];
	point.at(2) = RPoints->GetPoint(position)[2];

	if (::isnan(point.at(0)) && ::isnan(point.at(1)) && ::isnan(point.at(2)))
		cerr << "Point don't exixts!" << endl;

	return point;
}

/*! Get function for all right seed points.
 \return all right seed points as vtkPoints
 */
vtkSmartPointer<vtkPoints> FiberOrientation::getRPoints() {
	return RPoints;
}

/*! Set function for a left epicaridal seed point.
 \param position of the left epicaridal seed point
 \param point vector with the coordinates (x,y,z) of the left epicaridal seed point
 */
void FiberOrientation::setLPointEpi(vtkIdType position, vector<double> point) {
	double x = point.at(0);
	double y = point.at(1);
	double z = point.at(2);

	setLPointEpi(position, x, y, z);
}

/*! Set function for a left epicaridal seed point.
 \param position of the left epicaridal seed point
 \param x x coordinate of the left epicaridal seed point
 \param y y coordinate of the left epicaridal seed point
 \param z z coordinate of the left epicaridal seed point
 */
void FiberOrientation::setLPointEpi(vtkIdType position, double x, double y, double z) {
	LPointsEpi->SetPoint(position, x, y, z);
}

/*! Get function for a left epicaridal seed point.
 \param position of the left epicaridal seed point
 \return the left epicaridal seed point at position as vector with the coordinates (x,y,z). If it exist, if not then
 it write out "Point don't exixts!".
 */
vector<double> FiberOrientation::getLPointEpi(vtkIdType position) {
	vector<double> point(3);

	point.at(0) = LPointsEpi->GetPoint(position)[0];
	point.at(1) = LPointsEpi->GetPoint(position)[1];
	point.at(2) = LPointsEpi->GetPoint(position)[2];

	if (::isnan(point.at(0)) && ::isnan(point.at(1)) && ::isnan(point.at(2)))
		cerr << "Point don't exixts!" << endl;

	return point;
}

/*! Set function for all the left epicaridal seed points.
 \param lPoints the left epicaridal seed points in vtkPoints format
 */
void FiberOrientation::setLPointsEpi(vtkSmartPointer<vtkPoints> lPoints) {
	LPointsEpi = lPoints;
}

/*! Get function for all left epicaridal seed points.
 \return all left epicaridal seed points as vtkPoints
 */
vtkSmartPointer<vtkPoints> FiberOrientation::getLPointsEpi() {
	return LPointsEpi;
}

/*! Set function for a left endocaridal seed point.
 \param position of the left endocaridal seed point
 \param point vector with the coordinates (x,y,z) of the left endocaridal seed point
 */
void FiberOrientation::setLPointEndo(vtkIdType position, vector<double> point) {
	double x = point.at(0);
	double y = point.at(1);
	double z = point.at(2);

	setLPointEndo(position, x, y, z);
}

/*! Set function for a left endocaridal seed point.
 \param position of the left endocaridal seed point
 \param x x coordinate of the left endocaridal seed point
 \param y y coordinate of the left endocaridal seed point
 \param z z coordinate of the left endocaridal seed point
 */
void FiberOrientation::setLPointEndo(vtkIdType position, double x, double y, double z) {
	LPointsEndo->SetPoint(position, x, y, z);
}

/*! Get function for a left endocaridal seed point.
 \param position of the left endocaridal seed point
 \return the left endocaridal seed point at position as vector with the coordinates (x,y,z). If it exist, if not then
 it write out "Point don't exixts!".
 */
vector<double> FiberOrientation::getLPointEndo(vtkIdType position) {
	vector<double> point(3);

	point.at(0) = LPointsEndo->GetPoint(position)[0];
	point.at(1) = LPointsEndo->GetPoint(position)[1];
	point.at(2) = LPointsEndo->GetPoint(position)[2];

	if (::isnan(point.at(0)) && ::isnan(point.at(1)) && ::isnan(point.at(2)))
		cerr << "Point don't exixts!" << endl;

	return point;
}

/*! Set function for all the left endocaridal seed points.
 \param lPoints the left endocaridal seed points in vtkPoints format
 */
void FiberOrientation::setLPointsEndo(vtkSmartPointer<vtkPoints> lPoints) {
	LPointsEndo = lPoints;
}

/*! Get function for all left endocaridal seed points.
 \return all left endocaridal seed points as vtkPoints
 */
vtkSmartPointer<vtkPoints> FiberOrientation::getLPointsEndo() {
	return LPointsEndo;
}

/*! Set function for a free bridge point.
 \param position of the free bridge point
 \param point vector with the coordinates (x,y,z) of the free bridge point
 */
void FiberOrientation::setFreeBridgePoint(vtkIdType position, vector<double> point) {
	double x = point.at(0);
	double y = point.at(1);
	double z = point.at(2);

	setFreeBridgePoint(position, x, y, z);
}

/*! Set function for a free bridge point.
 \param position of the free bridge point
 \param x x coordinate of the free bridge point
 \param y y coordinate of the free bridge point
 \param z z coordinate of the free bridge point
 */
void FiberOrientation::setFreeBridgePoint(vtkIdType position, double x, double y, double z) {
	FreeBridgePoints->SetPoint(position, x, y, z);
}

/*! Get function for a free bridge point.
 \param position of the free bridge point
 \return the free bridge point at position as vector with the coordinates (x,y,z). If it exist, if not then it write
 out "Point don't exixts!".
 */
vector<double> FiberOrientation::getFreeBridgePoint(vtkIdType position) {
	vector<double> point(3);

	point.at(0) = FreeBridgePoints->GetPoint(position)[0];
	point.at(1) = FreeBridgePoints->GetPoint(position)[1];
	point.at(2) = FreeBridgePoints->GetPoint(position)[2];

	if (::isnan(point.at(0)) && ::isnan(point.at(1)) && ::isnan(point.at(2)))
		cerr << "Point don't exixts!" << endl;

	return point;
}

/*! Set function for a bridge point.
 \param position of the bridge point
 \param x x coordinate of the bridge point
 \param y y coordinate of the bridge point
 \param z z coordinate of the bridge point
 */
void FiberOrientation::setBPoint(vtkIdType position, double x, double y, double z) {
	BPoints->SetPoint(position, x, y, z);
}

/*! Set function for a bridge point.
 \param position of the bridge point
 \param point vector with the coordinates (x,y,z) of the bridge point
 */
void FiberOrientation::setBPoint(vtkIdType position, vector<double> point) {
	double x = point.at(0);
	double y = point.at(1);
	double z = point.at(2);

	setBPoint(position, x, y, z);
}

/*! Get function for a bridge point.
 \param position of the bridge point
 \return the bridge point at position as vector with the coordinates (x,y,z). If it exist, if not then it write out
 "Point don't exixts!".
 */
vector<double> FiberOrientation::getBPoint(vtkIdType position) {
	vector<double> point(3);

	point.at(0) = BPoints->GetPoint(position)[0];
	point.at(1) = BPoints->GetPoint(position)[1];
	point.at(2) = BPoints->GetPoint(position)[2];

	if (::isnan(point.at(0)) && ::isnan(point.at(1)) && ::isnan(point.at(2)))
		cerr << "Point don't exixts!" << endl;

	return point;
}

/*! Set function for the mesh resolution.
 \param resolution the resulution of the mesh
 */
void FiberOrientation::setResolution(double resolution) {
	Resolution = resolution;
}

/*! Get function for the mesh resolution.
 \return the resulution of the mesh
 */
double FiberOrientation::getResolution() {
	return Resolution;
}

/*! Read out all seed point from txt-file and move them in the corresponded material.
 \param data Pointer to the orginal mesh with marked insite and outsite
 \param fullPath path to the txt-file with the seed points
 */
void FiberOrientation::readOutFiberPoints(string FullPath) {
	string fullPath = FullPath;
	string path;
	string filename;
	long pathPos;

	pathPos = fullPath.find_last_of("/");
	path = fullPath.substr(0, pathPos);
	filename = fullPath.substr(pathPos + 1, fullPath.size());

	chdir(path.c_str());

	string Names[] = { "SCV1", "SCV2", "SCV3", "SCV4", "SCV5", "SCV6", "SCV7", "SCV8", "SCV9" };
	string NamesL[] = { "LV1", "LV2", "LV3", "LV4", "LV5", "LV6", "LV7", "LV8", "LV9", "LV10", "LV11", "LV12", "LV13" };
	string line = "\n";
	string buf = " ";

	int i = 0;
	int j = 0;
	ifstream Coordinatefile(filename.c_str());

	if (Coordinatefile.is_open()) {
		while (getline(Coordinatefile, line)) {
			if (line.size() > 0) {
				vector<string> CoordinatesReaded;
				stringstream ss(line);
				while (ss >> buf) { // reads string until next whitespace
					CoordinatesReaded.push_back(buf);
				}

				if (ConfigFiberorientation::rightAtrium) {
					if (i < 9) {
						if (!Names[i].compare(CoordinatesReaded.at(0).c_str())) {
							setRPoint(i, stod(CoordinatesReaded.at(2).c_str()), stod(CoordinatesReaded.at(3).c_str()), stod(CoordinatesReaded.at(4).c_str()));
							cout << "read R" << i << " (" << getRPoint(i).at(0) << ", " << getRPoint(i).at(1) << ", " <<
								getRPoint(i).at(2) << ") (global)" << endl;
						}
					}
					if (i >= 9) {
						if (!NamesL[j].compare(CoordinatesReaded.at(0).c_str())) {
							setLPointEpi(j, stod(CoordinatesReaded.at(2).c_str()), stod(CoordinatesReaded.at(3).c_str()), stod(CoordinatesReaded.at(4).c_str()));
							cout << "read L" << j << " (" << getLPointEpi(j).at(0) << ", " << getLPointEpi(j).at(1) << ", " <<
								getLPointEpi(j).at(2) << ") (global)" << endl;
						}
						j++;
					}
					i++;
				}
				else {
					if (!NamesL[i].compare(CoordinatesReaded.at(0).c_str())) {
						setLPointEpi(i, stod(CoordinatesReaded.at(2).c_str()), stod(CoordinatesReaded.at(3).c_str()), stod(CoordinatesReaded.at(4).c_str()));
						cout << "read L" << i << " (" << getLPointEpi(i).at(0) << ", " << getLPointEpi(i).at(1) << ", " <<
							getLPointEpi(i).at(2) << ") (global)" << endl;
					}
					i++;
				}
			}
		}
	}
	else {
		cerr << "Unable to open file coordinationfile";
		exit(-1);
	}
	Coordinatefile.close();
} // FiberOrientation::readOutFiberPoints

/*! Read out all free bridge point, bridge radius, bridge material from txt-file.
 \param fullPath path to the txt-file
 */
int FiberOrientation::readOutFreeBridgePointsPoints(string fullPath) {
	string path;
	string filename;
	long pathPos;

	pathPos = fullPath.find_last_of("/");
	path = fullPath.substr(0, pathPos);
	filename = fullPath.substr(pathPos + 1, fullPath.size());

	chdir(path.c_str());

	string line = "\n";
	string buf = " ";

	int i = 0;
	int numberOfLines = 0;
	ifstream Coordinatefile(filename.c_str());

	if (Coordinatefile.is_open()) {
		while (getline(Coordinatefile, line)) {
			if (line.size() > 0) {
				numberOfLines++;

				vector<string> CoordinatesReaded;
				stringstream ss(line);
				while (ss >> buf) { // reads string until next whitespace
					CoordinatesReaded.push_back(buf);
				}

				setFreeBridgePoint(i, Methods::round(stod(CoordinatesReaded.at(0).c_str()), 2), Methods::round(stod(CoordinatesReaded.at(1).c_str()), 2), Methods::round(stod(CoordinatesReaded.at(2).c_str()), 2));
				cout << "read left FreeBridgePoint: " << i << " (" << CoordinatesReaded[0].c_str() << "," <<
					CoordinatesReaded[1].c_str() << "," << CoordinatesReaded[2].c_str() << ")" << endl;

				setFreeBridgePoint(i + 1, Methods::round(stod(CoordinatesReaded.at(3).c_str()), 2), Methods::round(stod(CoordinatesReaded.at(4).c_str()), 2), Methods::round(stod(CoordinatesReaded.at(5).c_str()), 2));
				cout << "read right FreeBridgePoint: " << i + 1 << " (" << CoordinatesReaded[3].c_str() << "," <<
					CoordinatesReaded[4].c_str() << "," << CoordinatesReaded[5].c_str() << ")" << endl;

				cout << "read material class for bridge: " << CoordinatesReaded.at(6).c_str() << endl;
				freeBridgeMaterial.push_back(Methods::round(stod(CoordinatesReaded.at(6).c_str()), 3));

				cout << "read width for bridge: " << CoordinatesReaded.at(7).c_str() << endl;
				freeBridgeWidth.push_back(stoi(CoordinatesReaded.at(7).c_str()));

				i = i + 2;
			}
		}
		Coordinatefile.close();
	}
	else {
		cerr << "Unable to open point file for free bridgepoints ";
		exit(-1);
	}


	return numberOfLines;
} // FiberOrientation::readOutFreeBridgePointsPoints

/*! The definition part and calculation function for the fiber orientiation in the right atrium.
 \param data Pointer to the orginal mesh
 */
void FiberOrientation::setFiberOrientationRightAtrium(DataFormat &data) {
	cout << "calculate pathR1R2" << endl;
	setPath("pathR1R2", Methods::pathSearch(data, getRPoint(0), getRPoint(1), Material::Vorhof_rechts, "right"));
	cout << "calculate pathR2R3" << endl;
	setPath("pathR2R3", Methods::pathSearch(data, getRPoint(1), getRPoint(2), Material::Vorhof_rechts, "right"));
	cout << "calculate pathR3R1" << endl;
	setPath("pathR3R1", Methods::pathSearch(data, getRPoint(2), getRPoint(0), Material::Vorhof_rechts, "right"));

	setPath("scvRingComplete", Methods::add3Paths(getPath("pathR1R2"), getPath("pathR2R3"), getPath("pathR3R1")));
	Methods::pathMarkerRing(data, getPath("scvRingComplete"), Material::Vorhof_rechts);
	if (ConfigFiberorientation::growPathTransmural) {
		setPath("scvRingCompleteGrownPath", Methods::growPathTransmural(data, getPath("scvRingComplete"), 4.29, Material::Vorhof_rechts, "right"));
	}
	else {
		setPath("scvRingCompleteGrownPath", Methods::growPath(data, getPath("scvRingComplete"), 4.29, Material::Vorhof_rechts, Material::Vorhof_rechts));
	}
	Methods::smoothingRingPath(data, getPath("scvRingComplete"));
	Methods::setMaterial(data, getPath("scvRingCompleteGrownPath"), Material::SVC, Material::Vorhof_rechts);
	cout << "SVC finished" << endl;

	if (ConfigGeneral::debug) {
		Methods::writeIntermediateData(data, "SVC_right_debug");
	}

	cout << "calculate pathR8R7" << endl;
	setPath("pathR8R7", Methods::pathSearch(data, getRPoint(7), getRPoint(6), Material::Vorhof_rechts, "right"));
	cout << "calculate pathR7R9" << endl;
	setPath("pathR7R9", Methods::pathSearch(data, getRPoint(6), getRPoint(8), Material::Vorhof_rechts, "right"));
	cout << "calculate pathR9R8" << endl;
	setPath("pathR9R8", Methods::pathSearch(data, getRPoint(8), getRPoint(7), Material::Vorhof_rechts, "right"));

	setPath("tricusRC", Methods::add3Paths(getPath("pathR9R8"), getPath("pathR8R7"), getPath("pathR7R9")));
	Methods::pathMarkerRing(data, getPath("tricusRC"), Material::Vorhof_rechts);
	if (ConfigFiberorientation::growPathTransmural) {
		setPath("tricusRCGrownPath", Methods::growPathTransmural(data, getPath("tricusRC"), 6.27, Material::Vorhof_rechts, "right"));
	}
	else {
		setPath("tricusRCGrownPath", Methods::growPath(data, getPath("tricusRC"), 6.27, Material::Vorhof_rechts, Material::Vorhof_rechts));
	}
	Methods::smoothingRingPath(data, getPath("tricusRC"));
	Methods::setMaterial(data, getPath("tricusRCGrownPath"), Material::Tricuspid_Valve_Ring, Material::Vorhof_rechts);
	cout << "tricusring finished" << endl;

	if (ConfigGeneral::debug) {
		Methods::writeIntermediateData(data, "tricusring_right_debug");
	}

	setRPoint(39, Methods::pointAtPercentOfPath(data, getPath("pathR7R9"), 25));
	setRPoint(11, Methods::pointAtPercentOfPath(data, getPath("pathR3R1"), 3));


	cout << "calculate pathR12R7" << endl;
	setPath("pathR12R7", Methods::pathSearch(data, getRPoint(11), getRPoint(6), Material::Vorhof_rechts, Material::Tricuspid_Valve_Ring, Material::SVC, "right"));

	Methods::pathMarker(data, getPath("pathR12R7"), Material::Vorhof_rechts);
	if (ConfigFiberorientation::growPathTransmural) {
		setPath("tricuslineGrownPath", Methods::growPathTransmural(data, getPath("pathR12R7"), 3.96, Material::Vorhof_rechts, "right"));
	}
	else {
		setPath("tricuslineGrownPath", Methods::growPath(data, getPath("pathR12R7"), 3.96, Material::Vorhof_rechts, Material::Vorhof_rechts));
	}
	Methods::smoothingPath(data, getPath("pathR12R7"));
	Methods::setMaterial(data, getPath("tricuslineGrownPath"), Material::Tricusline, Material::Vorhof_rechts);
	cout << "tricusline finished" << endl;

	if (ConfigGeneral::debug) {
		Methods::writeIntermediateData(data, "tricusline_right_debug");
	}

	setRPoint(29, Methods::pointAtPercentOfPath(data, getPath("pathR12R7"), 40));

	cout << "calculate pathR4R5" << endl;
	setPath("pathR4R5", Methods::pathSearch(data, getRPoint(3), getRPoint(4), Material::Vorhof_rechts, "right"));

	setRPoint(13, Methods::pointAtPercentOfPath(data, getPath("pathR4R5"), 60));

	setRPoint(40, Methods::pointAtPercentOfPath(data, getPath("pathR4R5"), 30));

	cout << "calculate pathR3R12" << endl;
	setPath("pathR3R12", Methods::pathSearch(data, getRPoint(2), getRPoint(11), Material::Vorhof_rechts, Material::SVC, "right"));

	cout << "calculate pathR12R4" << endl;
	setPath("pathR12R4", Methods::pathSearchOverPlane(data, getRPoint(11), getRPoint(3), getRPoint(13)));

	setRPoint(28, Methods::pointAtPercentOfPath(data, getPath("pathR12R4"), 80));
	setRPoint(30, Methods::pointAtPercentOfPath(data, getPath("pathR12R4"), 65));


	cout << "calculate pathR4R14" << endl;
	setPath("pathR4R14", Methods::getPercentOfPath(data, getPath("pathR4R5"), 60));

	setRPoint(14, Methods::pointAtPercentOfPath(data, getPath("pathR12R7"), 20));
	cout << "calculate pathR14R15" << endl;
	setPath("pathR14R15", Methods::pathSearch(data, getRPoint(13), getRPoint(14), Material::Vorhof_rechts, Material::Tricusline, Material::SVC, "right"));

	setRPoint(15, Methods::pointAtPercentOfPath(data, getPath("pathR14R15"), 25));

	cout << "calculate pathR14R16 " << endl;
	setPath("pathR14R16", Methods::getPercentOfPath(data, getPath("pathR14R15"), 25));
	cout << "Crista terminalis calculation" << endl;


	setPath("christaTerminalis1", Methods::add2Paths(getPath("pathR3R12"), getPath("pathR12R4")));
	setPath("christaTerminalis2", Methods::add2Paths(getPath("pathR4R14"), getPath("pathR14R16")));

	setPath("christaTerminalis", Methods::add2Paths(getPath("christaTerminalis1"), getPath("christaTerminalis2")));

	Methods::pathMarker(data, getPath("christaTerminalis"));

	if (ConfigFiberorientation::growPathTransmural) {
		setPath("christaTerminalis1GrownPath", Methods::growCTPathTransmural(data, getPath("christaTerminalis1"), Material::Vorhof_rechts, Material::SVC, Material::Tricusline, "right"));
		setPath("christaTerminalis2GrownPath", Methods::growPathTransmural(data, getPath("christaTerminalis2"), 2.64, Material::Vorhof_rechts, Material::Tricusline, "right"));
	}
	else if (ConfigFiberorientation::growCristaDepOnWallthick) {
		double dist_along_crista = Methods::DistanceBetweenEndoandEpiSurfaceAlongPath(data, getPath("pathR12R4"), "avg");

		setPath("christaTerminalis1GrownPath", Methods::growCTPath(data, getPath("christaTerminalis1"), dist_along_crista, Material::Vorhof_rechts, Material::SVC, Material::Tricusline));
		setPath("christaTerminalis2GrownPath", Methods::growPath(data, getPath("christaTerminalis2"), dist_along_crista, Material::Vorhof_rechts, Material::Vorhof_rechts));
	}
	else {
		setPath("christaTerminalis1GrownPath", Methods::growCTPath(data, getPath("christaTerminalis1"), 2.64, Material::Vorhof_rechts, Material::SVC, Material::Tricusline));
		setPath("christaTerminalis2GrownPath", Methods::growPath(data, getPath("christaTerminalis2"), 2.64, Material::Vorhof_rechts, Material::Vorhof_rechts));
	}
	Methods::smoothingPath(data, getPath("christaTerminalis"));

	setPath("christaTerminalisGrownPath", Methods::add2Paths(getPath("christaTerminalis1GrownPath"), getPath("christaTerminalis2GrownPath")));
	Methods::setMaterial(data, getPath("christaTerminalisGrownPath"), Material::Crista_Terminalis_2);

	cout << "Crista Terminalis finished" << endl;

	if (ConfigGeneral::debug) {
		Methods::writeIntermediateData(data, "Crista_Terminalis_right_debug");
	}

	cout << "calculate pathR16R15 " << endl;
	setPath("pathR16R15", Methods::subPathFromPath(getPath("pathR14R15"), getPath("pathR14R16")));

	Methods::pathMarker(data, getPath("pathR16R15"), Material::Vorhof_rechts); // CTverlang
	if (ConfigFiberorientation::growPathTransmural) {
		setPath("cTVerlangGrownPath", Methods::growPathTransmural(data, getPath("pathR16R15"), 4.62, Material::Vorhof_rechts, "right"));
	}
	else {
		setPath("cTVerlangGrownPath", Methods::growPath(data, getPath("pathR16R15"), 4.62, Material::Vorhof_rechts, Material::Vorhof_rechts));
	}
	Methods::smoothingPath(data, getPath("pathR16R15"));
	Methods::setMaterial(data, getPath("cTVerlangGrownPath"), Material::otherwise_right_or_all, Material::Vorhof_rechts);
	cout << "CT verlang finished" << endl;

	if (ConfigGeneral::debug) {
		Methods::writeIntermediateData(data, "Crista_Terminalis_extended_right_debug");
	}

	cout << "calculate Pectinate muscles " << endl;
	setPath("helpPath1", Methods::getPercentOfPath(data, getPath("pathR12R4"), 3));
	setPath("helpPath2", Methods::subPathFromPath(getPath("pathR12R4"), getPath("helpPath1")));
	setPath("chrisTermA", Methods::add2Paths(getPath("helpPath2"), getPath("pathR4R14")));
	setPath("helpPath3", Methods::getPercentOfPath(data, getPath("pathR8R7"), 90));
	setPath("tricusRAa", Methods::inversePath(getPath("helpPath3")));

	PectinateCentre.at(0) = (getRPoint(13).at(0) + getRPoint(14).at(0)) / 2;
	PectinateCentre.at(1) = (getRPoint(13).at(1) + getRPoint(14).at(1)) / 2;
	PectinateCentre.at(2) = (getRPoint(13).at(2) + getRPoint(14).at(2)) / 2;

	setPectStartPoint(0, Methods::pointAtPercentOfPath(data, getPath("tricusRAa"), 0));
	setPectTargetPoint(0, Methods::pointAtPercentOfPath(data, getPath("chrisTermA"), 0));
	setPath("pathPect1", Methods::pathSearchOverPlane(data, getPectStartPoint(0), getPectTargetPoint(0), PectinateCentre));

	double lengthSeptumSpurium = 4.95;
	double lengthPathPecti1 = Methods::pathLength(data, getPath("pathPect1"));

	if (lengthPathPecti1 < lengthSeptumSpurium) {
		lengthSeptumSpurium = lengthPathPecti1;
	}

	setPath("SeptumSpurium", Methods::getMMOfPath(data, getPath("pathPect1"), lengthSeptumSpurium));

	Methods::pathMarker(data, getPath("SeptumSpurium"), Material::Vorhof_rechts);

	if (ConfigFiberorientation::growPathTransmural) {
		setPath("SeptumSpuriumGrownPath", Methods::growPathTransmural(data, getPath("SeptumSpurium"), 1.32, Material::Vorhof_rechts, "right"));
	}
	else if (ConfigFiberorientation::growPathTransmuralPectEndo) {
		setPath("SeptumSpuriumGrownPath", Methods::growPathTransmuralEndo(data, getPath("SeptumSpurium"), 1.32, Material::subendo_right_Atrial_Cardium, Material::SVC));
	}
	else if (ConfigFiberorientation::growPathPectfromEndo) {
		setPath("SeptumSpuriumGrownPath", Methods::growPathfromEndo(data, getPath("SeptumSpurium"), 1.32, Material::subendo_right_Atrial_Cardium, Material::SVC));
	}
	else {
		setPath("SeptumSpuriumGrownPath", Methods::growPath(data, getPath("SeptumSpurium"), 1.32, Material::Vorhof_rechts, Material::Vorhof_rechts));
	}

	Methods::smoothingPath(data, getPath("SeptumSpurium"));
	Methods::setMaterial(data, getPath(
		"SeptumSpuriumGrownPath"), Material::Pektinatmuskeln, Material::Vorhof_rechts, Material::Tricusline);

	cout << "SeptumSpurium calculated" << endl;

	if (ConfigGeneral::debug) {
		Methods::writeIntermediateData(data, "Septum_Spurium_right_debug");
	}

	Methods::pathMarker(data, getPath("pathPect1"), Material::Vorhof_rechts);
	if (ConfigFiberorientation::growPathTransmural) {
		setPath("pathPect1GrownPath", Methods::growPathTransmural(data, getPath("pathPect1"), 0.66, Material::Vorhof_rechts, "right"));
	}
	else if (ConfigFiberorientation::growPathTransmuralPectEndo) {
		setPath("pathPect1GrownPath", Methods::growPathTransmuralEndo(data, getPath("pathPect1"), 0.66, Material::subendo_right_Atrial_Cardium, Material::Tricusline));
	}
	else if (ConfigFiberorientation::growPathPectfromEndo) {
		setPath("pathPect1GrownPath", Methods::growPathfromEndo(data, getPath("pathPect1"), 0.66, Material::subendo_right_Atrial_Cardium, Material::Tricusline));
	}
	else {
		setPath("pathPect1GrownPath", Methods::growPath(data, getPath("pathPect1"), 0.66, Material::Vorhof_rechts, Material::Vorhof_rechts));
	}
	Methods::smoothingPath(data, getPath("pathPect1"));
	Methods::setMaterial(data, getPath(
		"pathPect1GrownPath"), Material::Pektinatmuskeln, Material::Vorhof_rechts, Material::Tricusline);

	if (ConfigGeneral::debug) {
		Methods::writeIntermediateData(data, "pathPect1_right_debug");
	}

	cout << "pathPect1 calculated" << endl;

	for (int i = 1; i < ConfigFiberorientation::numbersofPectinates; i++) {
		setPectStartPoint(i, Methods::pointAtPercentOfPath(data, getPath("tricusRAa"), (100.0 / (ConfigFiberorientation::numbersofPectinates - 1))*i));
		setPectTargetPoint(i, Methods::pointAtPercentOfPath(data, getPath("chrisTermA"), (100.0 / (ConfigFiberorientation::numbersofPectinates - 1))*i));
		setPath("pathPect" + to_string(i + 1), Methods::pathSearchOverPlane(data, getPectStartPoint(i), getPectTargetPoint(i), PectinateCentre));

		Methods::pathMarker(data, getPath("pathPect" + to_string(i + 1)), Material::Vorhof_rechts);
		if (ConfigFiberorientation::growPathTransmural) {
			setPath("pathPect" + to_string(i + 1) + "GrownPath", Methods::growPathTransmural(data, getPath("pathPect" + to_string(i + 1)), 0.66, Material::Vorhof_rechts, "right"));
		}
		else if (ConfigFiberorientation::growPathTransmuralPectEndo) {
			setPath("pathPect" + to_string(i + 1) + "GrownPath", Methods::growPathTransmuralEndo(data, getPath("pathPect" + to_string(i + 1)), 0.66, Material::subendo_right_Atrial_Cardium, Material::Tricusline));
		}
		else if (ConfigFiberorientation::growPathPectfromEndo) {
			setPath("pathPect" + to_string(i + 1) + "GrownPath", Methods::growPathfromEndo(data, getPath("pathPect" + to_string(i + 1)), 0.66, Material::subendo_right_Atrial_Cardium, Material::Tricusline));
		}
		else {
			setPath("pathPect" + to_string(i + 1) + "GrownPath", Methods::growPath(data, getPath("pathPect" + to_string(i + 1)), 0.66, Material::Vorhof_rechts, Material::Vorhof_rechts));
		}
		Methods::smoothingPath(data, getPath("pathPect" + to_string(i + 1)));
		Methods::setMaterial(data, getPath("pathPect" + to_string(
			i + 1) + "GrownPath"), Material::Pektinatmuskeln, Material::Vorhof_rechts, Material::Tricusline);
		cout << "pathPect" + to_string(i + 1) + " calculated" << endl;

		setPath("regionBetweenpathPect" + to_string(i) + "andpathPect" + to_string(i + 1), Methods::growBetweenPecti(data, getPath("pathPect" + to_string(i)), getPath("pathPect" + to_string(i + 1)), PectinateCentre, Material::Vorhof_rechts));
		Methods::setMaterial(data, getPath("regionBetweenpathPect" + to_string(i) + "andpathPect" + to_string(
			i + 1)), Material::otherwise_right_or_all, Material::Vorhof_rechts);
	}

	if (ConfigGeneral::debug) {
		Methods::writeIntermediateData(data, "Pectinate_muscles_right_debug");
	}


	cout << "calculate right appendages" << endl;
	cout << "calculate pathR13R6" << endl;
	setRPoint(12, Methods::pointAtPercentOfPath(data, getPath("pathR2R3"), 85));
	setPath("pathR13R6", Methods::pathSearch(data, getRPoint(12), getRPoint(5), Material::Vorhof_rechts, Material::Tricusline, Material::SVC, Material::Crista_Terminalis_2, "right"));

	cout << "calculate pathR13R17" << endl;
	setRPoint(16, Methods::pointAtMMOfPath(data, getPath("pathR13R6"), 6.6));
	setPath("pathR13R17", Methods::pathSearch(data, getRPoint(12), getRPoint(16), Material::Vorhof_rechts, Material::Tricusline, Material::SVC, Material::Crista_Terminalis_2, "right"));

	Methods::pathMarker(data, getPath("pathR13R17"), Material::Vorhof_rechts);
	if (ConfigFiberorientation::growPathTransmural) {
		setPath("AppendaGrownPath", Methods::growPathTransmural(data, getPath("pathR13R17"), 1.32, Material::Vorhof_rechts, "right"));
	}
	else if (ConfigFiberorientation::growPathTransmuralPectEndo) {
		setPath("AppendaGrownPath", Methods::growPathTransmuralEndo(data, getPath("pathR13R17"), 1.32, Material::subendo_right_Atrial_Cardium, Material::Vorhof_rechts));
	}
	else if (ConfigFiberorientation::growPathPectfromEndo) {
		setPath("AppendaGrownPath", Methods::growPathfromEndo(data, getPath("pathR13R17"), 1.32, Material::subendo_right_Atrial_Cardium, Material::Vorhof_rechts));
	}
	else {
		setPath("AppendaGrownPath", Methods::growPath(data, getPath("pathR13R17"), 1.32, Material::Vorhof_rechts, Material::Vorhof_rechts));
	}
	Methods::smoothingPath(data, getPath("pathR13R17"));
	Methods::setMaterial(data, getPath("AppendaGrownPath"), Material::Pektinatmuskeln, Material::Vorhof_rechts);

	if (ConfigGeneral::debug) {
		Methods::writeIntermediateData(data, "AppendaGrownPath_right_debug");
	}

	setRPoint(20, Methods::pointAtPercentOfPath(data, getPath("pathR12R7"), 55));
	setRPoint(21, Methods::pointAtPercentOfPath(data, getPath("pathPect1"), 59));
	setRPoint(22, Methods::pointAtPercentOfPath(data, getPath("pathPect1"), 43));
	setRPoint(25, Methods::pointAtPercentOfPath(data, getPath("pathPect1"), 20));

	// for marked the area between Pecti1 and AppendagesRing
	setPath("pathpathPect1Inverse", Methods::inversePath(getPath("pathPect1")));
	setPath("pathPectiTargetPoint1R22", Methods::getPercentOfPath(data, getPath("pathpathPect1Inverse"), 41));

	// for marked the area between Pecti1, AppendagesRing and Tricusline
	setPath("pathR7R26", Methods::getPercentOfPath(data, getPath("pathPect1"), 20));
	setPath("pathR31R21", Methods::getPercentOfPath(data, getPath("pathR12R7"), 55));
	setPath("pathR21R7", Methods::subPathFromPath(getPath("pathR12R7"), getPath("pathR31R21")));

	cout << "calculate pathR17R6" << endl;
	setPath("pathR17R6", Methods::pathSearch(data, getRPoint(16), getRPoint(5), Material::Vorhof_rechts, Material::Pektinatmuskeln, Material::Crista_Terminalis_2, Material::Tricusline, Material::SVC, "right"));

	cout << "calculate pathR22R6" << endl;
	setPath("pathR22R6", Methods::pathSearch(data, getRPoint(21), getRPoint(5), Material::Vorhof_rechts, Material::Pektinatmuskeln, "right"));

	cout << "calculate pathR23R6" << endl;
	setPath("pathR23R6", Methods::pathSearch(data, getRPoint(22), getRPoint(5), Material::Vorhof_rechts, Material::Pektinatmuskeln, "right"));

	cout << "calculate pathR26R6" << endl;
	setPath("pathR26R6", Methods::pathSearch(data, getRPoint(25), getRPoint(5), Material::Vorhof_rechts, Material::Pektinatmuskeln, Material::Tricusline, Material::Tricuspid_Valve_Ring, "right"));

	cout << "calculate pathR21R6" << endl;
	setPath("pathR21R6", Methods::pathSearch(data, getRPoint(20), getRPoint(
		5), Material::Vorhof_rechts, Material::Tricusline, "right"));

	setRPoint(33, Methods::getPointAlongPathOfPointInMaterial(data, getPath("pathR17R6"), Material::Vorhof_rechts, 0.66));
	setRPoint(34, Methods::getPointAlongPathOfPointInMaterial(data, getPath("pathR21R6"), Material::Vorhof_rechts, 1.32));
	setRPoint(35, Methods::getPointAlongPathOfPointInMaterial(data, getPath("pathR22R6"), Material::Vorhof_rechts, 0.66));
	setRPoint(36, Methods::getPointAlongPathOfPointInMaterial(data, getPath("pathR23R6"), Material::Vorhof_rechts, 0.66));
	setRPoint(37, Methods::getPointAlongPathOfPointInMaterial(data, getPath("pathR26R6"), Material::Vorhof_rechts, 1.32));

	cout << "calculate pathR17R22" << endl;
	setPath("pathR17R22", Methods::pathSearch(data, getRPoint(33), getRPoint(35), Material::Vorhof_rechts, "right"));

	cout << "calculate pathR22R23" << endl;
	setPath("pathR22R23", Methods::pathSearch(data, getRPoint(35), getRPoint(36), Material::Vorhof_rechts, "right"));

	cout << "calculate pathR23R26" << endl;
	setPath("pathR23R26", Methods::pathSearch(data, getRPoint(36), getRPoint(37), Material::Vorhof_rechts, "right"));

	cout << "calculate pathR26R21" << endl;
	setPath("pathR26R21", Methods::pathSearch(data, getRPoint(37), getRPoint(34), Material::Vorhof_rechts, "right"));

	cout << "calculate pathR21R17" << endl;
	setPath("pathR21R17", Methods::pathSearch(data, getRPoint(34), getRPoint(33), Material::Vorhof_rechts, "right"));

	setPath("appendageRing", Methods::add5Paths(getPath("pathR17R22"), getPath("pathR22R23"), getPath("pathR23R26"), getPath("pathR26R21"), getPath("pathR21R17")));
	Methods::pathMarkerRing(data, getPath("appendageRing"), Material::Vorhof_rechts);
	Methods::smoothingRingPath(data, getPath("appendageRing"));
	setPath("borderAppendagesRing", Methods::boundaryLayerAppend(data, getPath("appendageRing"), getRPoint(5)));

	setPath("growBorderAppendagesRing", Methods::growPathInRegion(data, getPath("borderAppendagesRing"), 2.64, Material::Vorhof_rechts, Material::Vorhof_rechts, getRPoint(5)));
	setPath("rightAppendages", Methods::growPathtoPoint(data, getPath("growBorderAppendagesRing"), getRPoint(5), Material::Vorhof_rechts));
	Methods::setMaterial(data, getPath("rightAppendages"), Material::Right_Atrial_Appendage, Material::Vorhof_rechts);

	if (ConfigGeneral::debug) {
		Methods::writeIntermediateData(data, "right_appendages_right_debug");
	}

	cout << "right appendages calculated" << endl;
	if (ConfigFiberorientation::growPectToRAppendage) {
		Methods::pathMarker(data, getPath("pathR17R6"), Material::Right_Atrial_Appendage);
		if (ConfigFiberorientation::growPathTransmural) {
			setPath("PectAppendaGrownPath", Methods::growPathTransmural(data, getPath("pathR17R6"), 0.66, Material::Right_Atrial_Appendage, "right"));
		}
		else if (ConfigFiberorientation::growPathTransmuralPectEndo) {
			setPath("PectAppendaGrownPath", Methods::growPathTransmuralEndo(data, getPath("pathR17R6"), 0.66, Material::subendo_right_Atrial_Cardium, Material::Right_Atrial_Appendage));
		}
		else if (ConfigFiberorientation::growPathPectfromEndo) {
			setPath("PectAppendaGrownPath", Methods::growPathfromEndo(data, getPath("pathR17R6"), 0.66, Material::subendo_right_Atrial_Cardium, Material::Right_Atrial_Appendage));
		}
		else {
			vtkSmartPointer<vtkIdList> tempPathAppendage = Methods::getCellsArroundPath(data, getPath(
				"pathR17R6"), 0.66, Material::Right_Atrial_Appendage);
			Methods::setStatus(data, tempPathAppendage, 0, Material::Right_Atrial_Appendage);
			setPath("PectAppendaGrownPath", Methods::growPath(data, getPath("pathR17R6"), 0.66, Material::Right_Atrial_Appendage, Material::Right_Atrial_Appendage));
		}
		Methods::smoothingPath(data, getPath("pathR17R6"));
		Methods::setMaterial(data, getPath(
			"PectAppendaGrownPath"), Material::Pektinatmuskeln, Material::Right_Atrial_Appendage);
	}
	if (ConfigGeneral::debug) {
		Methods::writeIntermediateData(data, "right_appendage_pathto_tip_right_debug");
	}

	cout << "calculated pectinate in appendage" << endl;

	cout << "mark area between Pecti1 and AppendagesRing" << endl;
	vector<double> seedPoint1 = Methods::getSeedPoint(data, getPath("pathPectiTargetPoint1R22"), getPath(
		"pathR17R22"), Material::Vorhof_rechts, Material::Pektinatmuskeln, Material::Right_Atrial_Appendage, 30, "right");
	setPath("areaBetweenPecti1Appendages", Methods::regionInterpolateBetween2PathsInDirection(data, getPath("pathPectiTargetPoint1R22"), getPath("pathR17R22"), Material::Vorhof_rechts, seedPoint1, getRPoint(12), getRPoint(11), getRPoint(14), getRPoint(21)));
	Methods::setMaterial(data, getPath(
		"areaBetweenPecti1Appendages"), Material::otherwise_right_or_all, Material::Vorhof_rechts);

	if (ConfigGeneral::debug) {
		Methods::writeIntermediateData(data, "areaBetweenPecti1Appendages_right_debug");
	}
	cout << "calculate Inferior Isthmus and under intercavale" << endl;
	setRPoint(31, Methods::pointAtPercentOfPath(data, getPath("pathR9R8"), 40));
	setPath("pathR32R14", Methods::pathSearch(data, getRPoint(31), getRPoint(13), Material::Vorhof_rechts, Material::Crista_Terminalis_2, Material::Tricuspid_Valve_Ring, "right"));
	setRPoint(32, Methods::pointAtPercentOfPath(data, getPath("pathR32R14"), 50));

	setPath("rigthHelpPath1", Methods::add3Paths(getPath("pathR12R7"), getPath("pathR7R9"), getPath("pathR9R8")));

	Methods::pathMarker(data, getPath("rigthHelpPath1"), Material::Tricusline, Material::Tricuspid_Valve_Ring);
	Methods::smoothingPath(data, getPath("rigthHelpPath1"));

	setPath("rigthHelpPath2", Methods::add2Paths(getPath("pathPect15"), getPath("pathR14R15")));
	Methods::pathMarker(data, getPath(
		"rigthHelpPath2"), Material::Crista_Terminalis_2, Material::Pektinatmuskeln, Material::otherwise_right_or_all);
	Methods::smoothingPath(data, getPath("rigthHelpPath2"));

	setPath("regionNotIsthmusUnderIntercavale", Methods::regionGrowBetween2PathwithDifferentOrientationInDirectionAndNotInDirection(
		data, getPath(
			"rigthHelpPath1"), getPath("rigthHelpPath2"), Material::Vorhof_rechts, getRPoint(32), getRPoint(13), getRPoint(16), getRPoint(3), getRPoint(7), getRPoint(13), getRPoint(7), PectinateCentre, getRPoint(3)));

	setRPoint(17, Methods::pointAtPercentOfPath(data, getPath("pathR12R4"), 80));

	Methods::setMaterialNotinDirection(data, getPath(
		"regionNotIsthmusUnderIntercavale"), Material::Inferior_Isthmus_right_atrium, Material::Vorhof_rechts, getRPoint(
			17), getRPoint(15), getRPoint(31), getRPoint(5));

	if (ConfigGeneral::debug) {
		Methods::writeIntermediateData(data, "isthmus_and_under_intercavale_right_debug");
	}

	cout << "isthmus and under intercavale calculated" << endl;

	cout << "calculate intercaval bundel" << endl;

	setRPoint(19, Methods::pointAtPercentOfPath(data, getPath("pathR12R4"), 20));

	setPath("pathR1R20", Methods::pathSearch(data, getRPoint(0), getRPoint(19), Material::Vorhof_rechts, Material::Crista_Terminalis_2, Material::SVC, "right"));
	Methods::pathMarker(data, getPath("pathR1R20"), Material::Vorhof_rechts);
	if (ConfigFiberorientation::growPathTransmural) {
		setPath("pathR1R20GrownPath", Methods::growPathTransmural(data, getPath("pathR1R20"), 2.64, Material::Vorhof_rechts, "right"));
	}
	else {
		setPath("pathR1R20GrownPath", Methods::growPath(data, getPath("pathR1R20"), 2.64, Material::Vorhof_rechts, Material::Vorhof_rechts));
	}
	Methods::smoothingPath(data, getPath("pathR1R20"));

	Methods::setMaterial(data, getPath("pathR1R20GrownPath"), Material::IntercavanalBundel, Material::Vorhof_rechts);

	if (ConfigGeneral::debug) {
		Methods::writeIntermediateData(data, "pathR1R20_right_debug");
	}

	setRPoint(23, Methods::pointAtPercentOfPath(data, getPath("pathR1R2"), 50));
	setRPoint(27, Methods::pointAtPercentOfPath(data, getPath("pathR12R4"), 90));

	setPath("pathR2R24", Methods::pathSearch(data, getRPoint(1), getRPoint(23), Material::Vorhof_rechts, Material::Crista_Terminalis_2, Material::SVC, Material::IntercavanalBundel, "right"));
	setPath("pathR24R28", Methods::pathSearch(data, getRPoint(23), getRPoint(27), Material::Vorhof_rechts, Material::Crista_Terminalis_2, Material::SVC, Material::IntercavanalBundel, "right"));
	setPath("pathR2R28", Methods::add2Paths(getPath("pathR2R24"), getPath("pathR24R28")));

	Methods::pathMarker(data, getPath("pathR2R28"), Material::Vorhof_rechts);
	if (ConfigFiberorientation::growPathTransmural) {
		setPath("pathR2R28GrownPath", Methods::growPathTransmural(data, getPath("pathR2R28"), 4.8, Material::Vorhof_rechts, "right"));
	}
	else {
		setPath("pathR2R28GrownPath", Methods::growPath(data, getPath("pathR2R28"), 4.8, Material::Vorhof_rechts, Material::Vorhof_rechts));
	}
	Methods::smoothingPath(data, getPath("pathR2R28"));
	Methods::setMaterial(data, getPath("pathR2R28GrownPath"), Material::IntercavanalBundel, Material::Vorhof_rechts);

	if (ConfigGeneral::debug) {
		Methods::writeIntermediateData(data, "pathR2R28_right_debug");
	}

	if (ConfigFiberorientation::markCoronary) {
		setPath("pathR5R9", Methods::pathSearch(data, getRPoint(4), getRPoint(8), Material::Crista_Terminalis_2, Material::Tricuspid_Valve_Ring, Material::Vorhof_rechts, "right"));


		DataFormat surface;
		surface.setVtkData(data.getSurfacewithNormals());
		surface.setCentrePoints(data.getCentrePointsSurface());

		setRPoint(43, Methods::findClosedPointinMaterialInRightEpi(surface, getRPoint(40)));
		setRPoint(44, Methods::findClosedPointinMaterialInRightEpi(surface, getRPoint(4)));
		setRPoint(45, Methods::findClosedPointinMaterialInRightEpi(surface, getRPoint(8)));


		setPath("pathR5R9EpiSurface", Methods::pathSearchNearToPath(data, getRPoint(44), getRPoint(45), getPath("pathR5R9"), "surfaceEpi"));
		setPath("pathR41R5EpiSurface", Methods::pathSearchNearToPath(data, getRPoint(43), getRPoint(44), getPath("pathR4R5"), "surfaceEpi"));
		setPath("pathR41R5R9EpiSurface", Methods::add2Paths(getPath("pathR41R5EpiSurface"), getPath("pathR5R9EpiSurface")));


		setPath("pathR41R5R9Epi", Methods::mapPathFromSurfaceToVolume(data, getPath("pathR41R5R9EpiSurface")));

		// Methods::setMaterial(data, getPath("pathR41R5R9Epi"), Material::testMaterialRight, Material::Vorhof_rechts);
		if (ConfigGeneral::debug) {
			Methods::writeIntermediateData(data, "preprecoronaryright_debug");
		}
		double dist_along_coronary = Methods::DistanceBetweenEndoandEpiSurfaceAlongPath(data, getPath(
			"pathR41R5R9Epi"), "max");

		setPath("pathR41R5R9EpiGrownPath", Methods::getCellsArroundPathwithStatusZero(data, getPath("pathR41R5R9Epi"), dist_along_coronary, Material::Vorhof_rechts));

		Methods::setMaterial(data, getPath("pathR41R5R9EpiGrownPath"), Material::testMaterialRight, Material::Vorhof_rechts);
		Methods::writeIntermediateData(data, "precoronary");

		// mark Coronary Sinus

		setRPoint(41, Methods::pointAtPercentOfPath(data, getPath("pathR4R5"), 75));
		setRPoint(42, Methods::pointAtPercentOfPath(data, getPath("pathR8R7"), 90));

		vector<double> seedPoint5 = Methods::findClosedPointinMaterialinDirectionPoint(data, getRPoint(4), getRPoint(8), 10, Material::Vorhof_rechts);

		Methods::setMaterialinInMaterialNotinDirectionandDirection(data, Material::Coronary_Sinus_ostium_tissue, Material::Vorhof_rechts, getRPoint(13), getRPoint(16), getRPoint(3), getRPoint(7), getRPoint(5), getRPoint(41), getRPoint(42), getRPoint(
			3), seedPoint5);
	}


	setRPoint(38, Methods::pointAtPercentOfPath(data, getPath("pathR2R28"), 50));
	vector<double> seedPoint2 = Methods::findFurthestPointinMaterialtoPoint(data, getRPoint(38), getRPoint(
		19), Material::Vorhof_rechts, 5);

	// Commented section to avoid having no fiber orientation in the area between sup & inf vena cava
	//
	Methods::setMaterialinInMaterialNotinDirectionandDirection(data, Material::not_right_Isthmus, Material::Vorhof_rechts, getRPoint(
		13), getRPoint(16), getRPoint(3), getRPoint(
			7), getRPoint(27), getRPoint(1), getRPoint(
				38), getRPoint(19), seedPoint2);
	if (ConfigGeneral::debug) {
		Methods::writeIntermediateData(data, "notrightIsthmus_right_debug");
	}

	if (ConfigFiberorientation::markCoronary) {
		Methods::setMaterial(data, getPath("pathR41R5R9EpiGrownPath"), Material::Vorhof_rechts, Material::testMaterialRight);
		Methods::setStatus(data, getPath("pathR41R5R9EpiGrownPath"), 0, Material::Vorhof_rechts);

		vector<double> seedPoint6 = Methods::findClosedPointinMaterial(data, getRPoint(4), Material::Vorhof_rechts);
		Methods::regionGrowMaterial(data, seedPoint6, Material::not_right_Isthmus, Material::Coronary_Sinus_ostium_tissue, Material::Vorhof_rechts);

		if (ConfigGeneral::debug) {
			Methods::writeIntermediateData(data, "coronary");
		}
	}


	cout << "coronary sinus calculated calculated" << endl;
	cout << "intercaval bundel calculated" << endl;

	vector<double> seedPoint3 = Methods::findClosedPointinMaterialinRadius(data, getRPoint(
		13), Material::Vorhof_rechts, Material::Crista_Terminalis_2, 5);

	// Commented section to avoid having no fiber orientation in the area between sup & inf vena cava
	//
	if (!isnan(seedPoint3.at(0))) {
		Methods::setMaterialinInMaterialNotinDirectionandDirection(data, Material::not_right_Isthmus, Material::Vorhof_rechts, getRPoint(13), getRPoint(
			16), getRPoint(3), getRPoint(7), getRPoint(
				27), getRPoint(1), getRPoint(38), getRPoint(
					19), seedPoint3);

	}

	if (ConfigGeneral::debug) {
		Methods::writeIntermediateData(data, "notrightIsthmus_right_debug");
	}

	cout << "marked intercaval bundel" << endl;
	vector<double> seedPoint4 = Methods::getSeedPoint(data, getPath("pathR1R20"), getPath(
		"pathR2R28"), Material::Vorhof_rechts, Material::IntercavanalBundel, Material::SVC, Material::Crista_Terminalis_2, 50, "right");
	setPath("pathR2R20", Methods::add2Paths(Methods::inversePath(getPath("pathR1R2")), getPath("pathR1R20")));

	Methods::pathMarker(data, getPath("pathR2R20"), Material::IntercavanalBundel, Material::SVC);
	Methods::smoothingPath(data, getPath("pathR2R20"));

	setPath("intercavalRegion", Methods::regionInterpolateBetween2PathwithSeed(data, getPath("pathR2R20"), getPath("pathR2R28"), Material::Vorhof_rechts, seedPoint4));
	Methods::setMaterial(data, getPath("intercavalRegion"), Material::IntercavanalBundel, Material::Vorhof_rechts);

	if (ConfigGeneral::debug) {
		Methods::writeIntermediateData(data, "intercavalRegion_right_debug");
	}

	cout << "right atrium calculated" << endl;
} // FiberOrientation::setFiberOrientationRightAtrium

/*! The definition part and calculation function for the fiber orientiation in the left enocardial atrium.
 \param data Pointer to the orginal mesh
 */
void FiberOrientation::setFiberOrientationLeftEndoAtrium(DataFormat &data) {
	cout << "calculate pathL3L1" << endl;
	setPath("pathL3L1", Methods::pathSearch(data, getLPointEndo(2), getLPointEndo(0), Material::Vorhof_links_Endo, "left"));
	cout << "calculate pathL1L2" << endl;
	setPath("pathL1L2", Methods::pathSearch(data, getLPointEndo(0), getLPointEndo(1), Material::Vorhof_links_Endo, "left"));
	cout << "calculate pathL2L3" << endl;
	setPath("pathL2L3", Methods::pathSearch(data, getLPointEndo(1), getLPointEndo(2), Material::Vorhof_links_Endo, "left"));

	setPath("MVRing", Methods::add3Paths(getPath("pathL3L1"), getPath("pathL1L2"), getPath("pathL2L3")));
	Methods::pathMarkerRing(data, getPath("MVRing"), Material::Vorhof_links_Endo);
	if (ConfigFiberorientation::growPathTransmural) {
		setPath("MVRingGrownPath", Methods::growPathTransmural(data, getPath("MVRing"), 4.5 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links_Endo, "left"));
	}
	else {
		setPath("MVRingGrownPath", Methods::growPath(data, getPath("MVRing"), 4.5 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links_Endo, Material::Vorhof_links_Endo));
	}
	Methods::smoothingRingPath(data, getPath("MVRing"));
	Methods::setMaterial(data, getPath(
		"MVRingGrownPath"), Material::subendo_left_Atrial_Cardium, Material::Vorhof_links_Endo);
	cout << "MVRing finished" << endl;

	if (ConfigGeneral::debug) {
		Methods::writeIntermediateData(data, "MVRing_left_endo_debug");
	}

	cout << "calculate pathL4L6" << endl;
	setPath("pathL4L6", Methods::pathSearch(data, getLPointEndo(3), getLPointEndo(5), Material::subendo_left_Atrial_Cardium, Material::Vorhof_links_Endo, "left"));

	cout << "calculate pathL6L8" << endl;
	setPath("pathL6L8", Methods::pathSearch(data, getLPointEndo(5), getLPointEndo(7), Material::subendo_left_Atrial_Cardium, Material::Vorhof_links_Endo, "left"));

	cout << "calculate pathL8L13" << endl;
	setPath("pathL8L13", Methods::pathSearch(data, getLPointEndo(7), getLPointEndo(12), Material::subendo_left_Atrial_Cardium, Material::Vorhof_links_Endo, "left"));

	cout << "calculate pathL13L6" << endl;
	setPath("pathL13L6", Methods::pathSearch(data, getLPointEndo(12), getLPointEndo(5), Material::subendo_left_Atrial_Cardium, Material::Vorhof_links_Endo, "left"));

	setPath("pathL4L6L8L13L6", Methods::add4Paths(getPath("pathL4L6"), getPath("pathL6L8"), getPath("pathL8L13"), getPath("pathL13L6")));

	Methods::pathMarker(data, getPath("pathL4L6L8L13L6"), Material::Vorhof_links_Endo);
	if (ConfigFiberorientation::growPathTransmural) {
		setPath("pathL4L6L8L13L6GrownPath", Methods::growPathTransmural(data, getPath("pathL4L6L8L13L6"), 2.64 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links_Endo, "left"));
	}
	else {
		setPath("pathL4L6L8L13L6GrownPath", Methods::growPath(data, getPath("pathL4L6L8L13L6"), 2.64 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links_Endo, Material::Vorhof_links_Endo));
	}
	Methods::smoothingRingPath(data, getPath("pathL4L6L8L13L6"));
	Methods::setMaterial(data, getPath(
		"pathL4L6L8L13L6GrownPath"), Material::subendo_left_Atrial_Cardium, Material::Vorhof_links_Endo);

	if (ConfigGeneral::debug) {
		Methods::writeIntermediateData(data, "pathL4L6L8L13L6_left_endo_debug");
	}

	cout << "calculate pathL5L7" << endl;
	setPath("pathL5L7", Methods::pathSearch(data, getLPointEndo(4), getLPointEndo(6), Material::subendo_left_Atrial_Cardium, Material::Vorhof_links_Endo, "left"));

	cout << "calculate pathL7L9" << endl;
	setPath("pathL7L9", Methods::pathSearch(data, getLPointEndo(6), getLPointEndo(8), Material::subendo_left_Atrial_Cardium, Material::Vorhof_links_Endo, "left"));

	cout << "calculate pathL9L1" << endl;
	setPath("pathL9L1", Methods::pathSearch(data, getLPointEndo(8), getLPointEndo(0), Material::subendo_left_Atrial_Cardium, Material::Vorhof_links_Endo, "left"));
	setPath("pathL9L31", Methods::getPercentOfPath(data, getPath("pathL9L1"), 30));
	setLPointEndo(30, Methods::pointAtPercentOfPath(data, getPath("pathL9L1"), 30));

	setPath("pathL1L9", Methods::inversePath(getPath("pathL9L1")));
	setPath("pathL1L32", Methods::getPercentOfPath(data, getPath("pathL1L9"), 30));
	setLPointEndo(31, Methods::pointAtPercentOfPath(data, getPath("pathL1L9"), 30));

	cout << "calculate pathL31L12" << endl;
	setPath("pathL31L12", Methods::pathSearch(data, getLPointEndo(30), getLPointEndo(11), Material::subendo_left_Atrial_Cardium, Material::Vorhof_links_Endo, "left"));

	cout << "calculate pathL12L5" << endl;
	setPath("pathL12L5", Methods::pathSearch(data, getLPointEndo(11), getLPointEndo(4), Material::subendo_left_Atrial_Cardium, Material::Vorhof_links_Endo, "left"));

	setPath("pathL5L7L9L31L12L5", Methods::add5Paths(getPath("pathL5L7"), getPath("pathL7L9"), getPath("pathL9L31"), getPath("pathL31L12"), getPath("pathL12L5")));
	Methods::pathMarkerRing(data, getPath("pathL5L7L9L31L12L5"), Material::Vorhof_links_Endo);
	if (ConfigFiberorientation::growPathTransmural) {
		setPath("pathL5L7L9L31L12L5GrownPath", Methods::growPathTransmural(data, getPath("pathL5L7L9L31L12L5"), 2.64 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links_Endo, "left"));
	}
	else {
		setPath("pathL5L7L9L31L12L5GrownPath", Methods::growPath(data, getPath("pathL5L7L9L31L12L5"), 2.64 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links_Endo, Material::Vorhof_links_Endo));
	}
	Methods::smoothingRingPath(data, getPath("pathL5L7L9L31L12L5"));
	Methods::setMaterial(data, getPath(
		"pathL5L7L9L31L12L5GrownPath"), Material::subendo_left_Atrial_Cardium, Material::Vorhof_links_Endo);

	if (ConfigGeneral::debug) {
		Methods::writeIntermediateData(data, "pathL5L7L9L31L12L5_left_endo_debug");
	}

	cout << "calculate pathL32L12" << endl;
	setPath("pathL32L12", Methods::pathSearch(data, getLPointEndo(31), getLPointEndo(11), Material::subendo_left_Atrial_Cardium, Material::Vorhof_links_Endo, "left"));

	cout << "calculate pathL12L7" << endl;
	setPath("pathL12L7", Methods::pathSearch(data, getLPointEndo(11), getLPointEndo(6), Material::subendo_left_Atrial_Cardium, Material::Vorhof_links_Endo, "left"));

	setPath("pathL1overL32L12L7", Methods::add3Paths(getPath("pathL1L32"), getPath("pathL32L12"), getPath("pathL12L7")));
	Methods::pathMarker(data, getPath("pathL1overL32L12L7"), Material::Vorhof_links_Endo);
	if (ConfigFiberorientation::growPathTransmural) {
		setPath("pathL1overL32L12L7GrownPath", Methods::growPathTransmural(data, getPath("pathL1overL32L12L7"), 2.64 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links_Endo, "left"));
	}
	else {
		setPath("pathL1overL32L12L7GrownPath", Methods::growPath(data, getPath("pathL1overL32L12L7"), 2.64 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links_Endo, Material::Vorhof_links_Endo));
	}
	Methods::smoothingPath(data, getPath("pathL1overL32L12L7"));
	Methods::setMaterial(data, getPath(
		"pathL1overL32L12L7GrownPath"), Material::subendo_left_Atrial_Cardium, Material::Vorhof_links_Endo);

	if (ConfigGeneral::debug) {
		Methods::writeIntermediateData(data, "pathL1overL32L12L7_left_endo_debug");
	}

	cout << "calculate " << endl;
	setPath("pathL5L4", Methods::pathSearch(data, getLPointEndo(4), getLPointEndo(3), Material::subendo_left_Atrial_Cardium, Material::Vorhof_links_Endo, "left"));
	setLPointEndo(17, Methods::pointAtPercentOfPath(data, getPath("pathL5L4"), 50));
	setPath("pathL7L6", Methods::pathSearch(data, getLPointEndo(6), getLPointEndo(5), Material::subendo_left_Atrial_Cardium, Material::Vorhof_links_Endo, "left"));
	setLPointEndo(18, Methods::pointAtPercentOfPath(data, getPath("pathL7L6"), 50));
	setPath("pathL9L8", Methods::pathSearch(data, getLPointEndo(8), getLPointEndo(7), Material::subendo_left_Atrial_Cardium, Material::Vorhof_links_Endo, "left"));
	setLPointEndo(19, Methods::pointAtPercentOfPath(data, getPath("pathL9L8"), 50));
	setPath("pathL19L20", Methods::pathSearch(data, getLPointEndo(18), getLPointEndo(19), Material::subendo_left_Atrial_Cardium, Material::Vorhof_links_Endo, "left"));
	setLPointEndo(20, Methods::pointAtPercentOfPath(data, getPath("pathL19L20"), 50));
	setPath("pathL19L21", Methods::pathSearch(data, getLPointEndo(18), getLPointEndo(20), Material::subendo_left_Atrial_Cardium, Material::Vorhof_links_Endo, "left"));
	setPath("pathL18L19", Methods::pathSearch(data, getLPointEndo(17), getLPointEndo(18), Material::subendo_left_Atrial_Cardium, Material::Vorhof_links_Endo, "left"));
	setPath("", Methods::add2Paths(getPath("pathL18L19"), getPath("pathL19L21")));

	Methods::pathMarker(data, getPath(""), Material::Vorhof_links_Endo);
	if (ConfigFiberorientation::growPathTransmural) {
		setPath("GrownPath", Methods::growPathTransmural(data, getPath(""), 2.64 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links_Endo, "left"));
	}
	else {
		setPath("GrownPath", Methods::growPath(data, getPath(""), 2.64 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links_Endo, Material::Vorhof_links_Endo));
	}
	Methods::smoothingPath(data, getPath(""));
	Methods::setMaterial(data, getPath(
		"GrownPath"), Material::subendo_left_Atrial_Cardium, Material::Vorhof_links_Endo);

	if (ConfigGeneral::debug) {
		Methods::writeIntermediateData(data, "_left_endo_debug");
	}

	cout << " finished" << endl;

	cout << "calculate pathL21L9" << endl;
	setPath("pathL21L9", Methods::pathSearch(data, getLPointEndo(20), getLPointEndo(8), Material::subendo_left_Atrial_Cardium, Material::Vorhof_links_Endo, "left"));

	Methods::pathMarker(data, getPath("pathL21L9"), Material::Vorhof_links_Endo);
	if (ConfigFiberorientation::growPathTransmural) {
		setPath("pathL21L9GrownPath", Methods::growPathTransmural(data, getPath("pathL21L9"), 2.64 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links_Endo, "left"));
	}
	else {
		setPath("pathL21L9GrownPath", Methods::growPath(data, getPath("pathL21L9"), 2.64 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links_Endo, Material::Vorhof_links_Endo));
	}
	Methods::smoothingPath(data, getPath("pathL21L9"));
	Methods::setMaterial(data, getPath(
		"pathL21L9GrownPath"), Material::subendo_left_Atrial_Cardium, Material::Vorhof_links_Endo);

	if (ConfigGeneral::debug) {
		Methods::writeIntermediateData(data, "pathL21L9_left_endo_debug");
	}

	cout << "pathL21L9 finished" << endl;

	cout << "calculate pathL21L8" << endl;
	setPath("pathL21L8", Methods::pathSearch(data, getLPointEndo(20), getLPointEndo(7), Material::subendo_left_Atrial_Cardium, Material::Vorhof_links_Endo, "left"));

	Methods::pathMarker(data, getPath("pathL21L8"), Material::Vorhof_links_Endo);
	if (ConfigFiberorientation::growPathTransmural) {
		setPath("pathL21L8GrownPath", Methods::growPathTransmural(data, getPath("pathL21L8"), 2.64 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links_Endo, "left"));
	}
	else {
		setPath("pathL21L8GrownPath", Methods::growPath(data, getPath("pathL21L8"), 2.64 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links_Endo, Material::Vorhof_links_Endo));
	}
	Methods::smoothingPath(data, getPath("pathL21L8"));
	Methods::setMaterial(data, getPath(
		"pathL21L8GrownPath"), Material::subendo_left_Atrial_Cardium, Material::Vorhof_links_Endo);

	if (ConfigGeneral::debug) {
		Methods::writeIntermediateData(data, "pathL21L8_left_endo_debug");
	}

	cout << "pathL21L8 finished" << endl;

	cout << "marked Area L5 L9 L21 L8 L4" << endl;
	Methods::setStatus(data, getPath("pathL5L4"), 1, Material::Vorhof_links_Endo);
	if (ConfigFiberorientation::growPathTransmural) {
		setPath("pathL5L4GrownPath", Methods::growPathTransmural(data, getPath("pathL5L4"), 2.64 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links_Endo, "left"));
	}
	else {
		setPath("pathL5L4GrownPath", Methods::getCellsArroundPath(data, getPath("pathL5L4"), 2.64 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links_Endo));
	}
	Methods::setMaterial(data, getPath("pathL5L4GrownPath"), Material::otherwise_left, Material::Vorhof_links_Endo);


	setPath("pathL5overL7toL9", Methods::add2Paths(getPath("pathL5L7"), getPath("pathL7L9")));
	setPath("pathandA", Methods::add2Paths(getPath(""), getPath("pathL21L9")));

	vector<double> seedPoint1 = Methods::getSeedPoint(data, getPath("pathL5overL7toL9"), getPath(
		"pathandA"), Material::Vorhof_links_Endo, Material::subendo_left_Atrial_Cardium, 50, "left");
	setPath("areaBetweenMiddleLineandLeftPaths", Methods::regionInterpolateBetween2PathwithSeed(data, getPath("pathL5overL7toL9"), getPath("pathandA"),
		Material::Vorhof_links_Endo, seedPoint1));
	Methods::setMaterial(data, getPath(
		"areaBetweenMiddleLineandLeftPaths"), Material::subendo_left_Atrial_Cardium, Material::Vorhof_links_Endo);

	setPath("pathL4overL6toL8", Methods::add2Paths(getPath("pathL4L6"), getPath("pathL6L8")));
	setPath("pathandB", Methods::add2Paths(getPath(""), getPath("pathL21L8")));

	vector<double> seedPoint2 = Methods::getSeedPoint(data, getPath("pathL4overL6toL8"), getPath("pathandB"), Material::Vorhof_links_Endo, Material::subendo_left_Atrial_Cardium, 50, "left");
	setPath("areaBetweenMiddleLineandRightPaths", Methods::regionInterpolateBetween2PathwithSeed(data, getPath("pathL4overL6toL8"), getPath("pathandB"), Material::Vorhof_links_Endo, seedPoint2));
	Methods::setMaterial(data, getPath(
		"areaBetweenMiddleLineandRightPaths"), Material::subendo_left_Atrial_Cardium, Material::Vorhof_links_Endo);

	if (ConfigGeneral::debug) {
		Methods::writeIntermediateData(data, "areaBetweenMiddleLineandRight_left_endo_debug");
	}

	cout << "Area finished" << endl;


	cout << "marked Area L14L4 and L1L15" << endl;
	Methods::setMaterial(data, getPath("pathL5L4GrownPath"), Material::Vorhof_links_Endo, Material::otherwise_left);
	Methods::setStatus(data, getPath("pathL5L4GrownPath"), 0, Material::Vorhof_links_Endo);


	cout << "calculate pathL1L4" << endl;
	setPath("pathL1L4", Methods::pathSearch(data, getLPointEndo(0), getLPointEndo(3), Material::subendo_left_Atrial_Cardium, Material::Vorhof_links_Endo, "left"));
	setPath("pathL1L16", Methods::getPercentOfPath(data, getPath("pathL1L4"), 50));

	setLPointEndo(15, Methods::pointAtPercentOfPath(data, getPath("pathL1L4"), 50));

	cout << "calculate pathL16L5" << endl;
	setPath("pathL16L5", Methods::pathSearch(data, getLPointEndo(15), getLPointEndo(4), Material::subendo_left_Atrial_Cardium, Material::Vorhof_links_Endo, "left"));

	setPath("pathL1L16L5", Methods::add3Paths(getPath("pathL1L16"), getPath("pathL16L5"), getPath("pathL5L7")));
	setPath("pathL1L16L5", Methods::inversePath(getPath("pathL1L16L5")));

	Methods::pathMarker(data, getPath("pathL1L16L5"), Material::Vorhof_links_Endo);
	if (ConfigFiberorientation::growPathTransmural) {
		setPath("pathL1L16L5GrownPath", Methods::growPathTransmural(data, getPath("pathL1L16L5"), 3.3 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links_Endo, "left"));
	}
	else {
		setPath("pathL1L16L5GrownPath", Methods::growPath(data, getPath("pathL1L16L5"), 3.3 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links_Endo, Material::Vorhof_links_Endo));
	}
	Methods::smoothingPath(data, getPath("pathL1L16L5"));
	Methods::setMaterial(data, getPath(
		"pathL1L16L5GrownPath"), Material::subendo_left_Atrial_Cardium, Material::Vorhof_links_Endo);

	if (ConfigGeneral::debug) {
		Methods::writeIntermediateData(data, "pathL1L16L5_left_endo_debug");
	}

	setLPointEndo(13, Methods::pointAtPercentOfPath(data, getPath("pathL2L3"), 80));

	cout << "calculate pathL14L4" << endl;
	setPath("pathL14L4", Methods::pathSearch(data, getLPointEndo(13), getLPointEndo(3), Material::subendo_left_Atrial_Cardium, Material::Vorhof_links_Endo, "left"));

	setPath("pathL14L4L6", Methods::add2Paths(getPath("pathL14L4"), getPath("pathL4L6")));
	setPath("pathL14L4L6", Methods::inversePath(getPath("pathL14L4L6")));

	Methods::pathMarker(data, getPath("pathL14L4L6"), Material::Vorhof_links_Endo);
	if (ConfigFiberorientation::growPathTransmural) {
		setPath("pathL14L4L6GrownPath", Methods::growPathTransmural(data, getPath("pathL14L4L6"), 2.64 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links_Endo, "left"));
	}
	else {
		setPath("pathL14L4L6GrownPath", Methods::growPath(data, getPath("pathL14L4L6"), 2.64 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links_Endo, Material::Vorhof_links_Endo));
	}
	Methods::smoothingPath(data, getPath("pathL14L4L6"));
	Methods::setMaterial(data, getPath(
		"pathL14L4L6GrownPath"), Material::subendo_left_Atrial_Cardium, Material::Vorhof_links_Endo);

	if (ConfigGeneral::debug) {
		Methods::writeIntermediateData(data, "pathL14L4L6_left_endo_debug");
	}

	vector<double> seedPoint3 = Methods::findClosedPointinMaterialwithoutOrientation(data, getLPointEndo(
		17), Material::Vorhof_links_Endo);

	setPath("areaBettweenpathL14L6andL1L7", Methods::regionInterpolateBetween2PathwithSeed(data, getPath("pathL14L4L6"), getPath("pathL1L16L5"), Material::Vorhof_links_Endo, seedPoint3));
	Methods::setMaterial(data, getPath(
		"areaBettweenpathL14L6andL1L7"), Material::subendo_left_Atrial_Cardium, Material::Vorhof_links_Endo);

	if (ConfigGeneral::debug) {
		Methods::writeIntermediateData(data, "areaBettweenpathL14L6andL1L7");
	}

	cout << "Area finished" << endl;

	Methods::setMaterial(data, getPath(
		"MVRingGrownPath"), Material::Mitral_Valve_Ring, Material::subendo_left_Atrial_Cardium);


	cout << "left endo-atrium calculated" << endl;
} // FiberOrientation::setFiberOrientationLeftEndoAtrium

/*! The definition part and calculation function for the fiber orientiation in the left epicardial atrium.
 \param data Pointer to the orginal mesh
 */
void FiberOrientation::setFiberOrientationLeftEpiAtrium(DataFormat &data) {
	cout << "calculate pathL3L1" << endl;
	setPath("pathL3L1", Methods::pathSearch(data, getLPointEpi(2), getLPointEpi(0), Material::Vorhof_links, "left"));
	cout << "calculate pathL1L2" << endl;
	setPath("pathL1L2", Methods::pathSearch(data, getLPointEpi(0), getLPointEpi(1), Material::Vorhof_links, "left"));
	cout << "calculate pathL2L3" << endl;
	setPath("pathL2L3", Methods::pathSearch(data, getLPointEpi(1), getLPointEpi(2), Material::Vorhof_links, "left"));

	setPath("MVRing", Methods::add3Paths(getPath("pathL3L1"), getPath("pathL1L2"), getPath("pathL2L3")));
	Methods::pathMarkerRing(data, getPath("MVRing"), Material::Vorhof_links);
	if (ConfigFiberorientation::growPathTransmural) {
		setPath("MVRingGrownPath", Methods::growPathTransmural(data, getPath("MVRing"), 5.5 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links, "left"));
	}
	else {
		setPath("MVRingGrownPath", Methods::growPath(data, getPath("MVRing"), 5.5 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links, Material::Vorhof_links));
	}

	Methods::smoothingRingPath(data, getPath("MVRing"));
	Methods::setMaterial(data, getPath("MVRingGrownPath"), Material::subepi_Atrial_Cardium, Material::Vorhof_links);
	cout << "MVRing finished" << endl;

	if (ConfigGeneral::debug) {
		Methods::writeIntermediateData(data, "MVRing_left_epi_debug");
	}

	setLPointEpi(13, Methods::pointAtPercentOfPath(data, getPath("pathL2L3"), 80));
	setLPointEpi(56, Methods::pointAtPercentOfPath(data, getPath("pathL3L1"), 25));

	cout << "calculate pathL4L14" << endl;
	setPath("pathL4L14", Methods::pathSearch(data, getLPointEpi(3), getLPointEpi(13), Material::subepi_Atrial_Cardium, Material::Vorhof_links, "left"));
	setLPointEpi(21, Methods::pointAtPercentOfPath(data, getPath("pathL4L14"), 10));
	setLPointEpi(22, Methods::pointAtPercentOfPath(data, getPath("pathL4L14"), 80));

	cout << "calculate pathL8L2" << endl;
	setPath("pathL8L2", Methods::pathSearch(data, getLPointEpi(7), getLPointEpi(1), Material::subepi_Atrial_Cardium, Material::Vorhof_links, "left"));
	setLPointEpi(23, Methods::pointAtPercentOfPath(data, getPath("pathL8L2"), 10));
	setLPointEpi(24, Methods::pointAtPercentOfPath(data, getPath("pathL8L2"), 80));

	cout << "calculate pathL9L1" << endl;
	setPath("pathL9L1", Methods::pathSearch(data, getLPointEpi(8), getLPointEpi(0), Material::subepi_Atrial_Cardium, Material::Vorhof_links, "left"));
	setLPointEpi(25, Methods::pointAtPercentOfPath(data, getPath("pathL9L1"), 20));
	setLPointEpi(26, Methods::pointAtPercentOfPath(data, getPath("pathL9L1"), 75));

	cout << "calculate pathL1L4" << endl;
	setPath("pathL1L4", Methods::pathSearch(data, getLPointEpi(0), getLPointEpi(3), Material::subepi_Atrial_Cardium, Material::Vorhof_links, "left"));
	setLPointEpi(15, Methods::pointAtPercentOfPath(data, getPath("pathL1L4"), 50));

	cout << "calculate pathL4L6" << endl;
	setPath("pathL4L6", Methods::pathSearch(data, getLPointEpi(3), getLPointEpi(5), Material::subepi_Atrial_Cardium, Material::Vorhof_links, "left"));

	cout << "calculate pathL6L8" << endl;
	setPath("pathL6L8", Methods::pathSearch(data, getLPointEpi(5), getLPointEpi(7), Material::subepi_Atrial_Cardium, Material::Vorhof_links, "left"));

	cout << "calculate pathL8L13" << endl;
	setPath("pathL8L13", Methods::pathSearch(data, getLPointEpi(7), getLPointEpi(12), Material::subepi_Atrial_Cardium, Material::Vorhof_links, "left"));

	cout << "calculate pathL13L6" << endl;
	setPath("pathL13L6", Methods::pathSearch(data, getLPointEpi(12), getLPointEpi(5), Material::subepi_Atrial_Cardium, Material::Vorhof_links, "left"));

	setPath("pathL6L8L13L6", Methods::add3Paths(getPath("pathL6L8"), getPath("pathL8L13"), getPath("pathL13L6")));

	Methods::pathMarkerRing(data, getPath("pathL6L8L13L6"), Material::Vorhof_links);

	if (ConfigFiberorientation::growPathTransmural) {
		setPath("pathL6L8L13L6GrownPathinEpi", Methods::growPathTransmural(data, getPath("pathL6L8L13L6"), 3.3 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links, Material::Vorhof_links_Endo, "left"));
		setPath("pathL6L8L13L6GrownPathinEndo", getPath("pathL6L8L13L6GrownPathinEpi"));
	}
	else {
		setPath("pathL6L8L13L6GrownPathinEpi", Methods::growPath(data, getPath("pathL6L8L13L6"), 3.3 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links, Material::Vorhof_links));
		setPath("pathL6L8L13L6GrownPathinEndo", Methods::growPath(data, getPath("pathL6L8L13L6"), 4.29 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links, Material::Vorhof_links_Endo));
	}

	Methods::smoothingRingPath(data, getPath("pathL6L8L13L6"));
	Methods::setMaterial(data, getPath(
		"pathL6L8L13L6GrownPathinEpi"), Material::subepi_Atrial_Cardium, Material::Vorhof_links);
	Methods::setMaterial(data, getPath(
		"pathL6L8L13L6GrownPathinEndo"), Material::subepi_Atrial_Cardium, Material::Vorhof_links_Endo);

	if (ConfigGeneral::debug) {
		Methods::writeIntermediateData(data, "pathL6L8L13L6_left_epi_debug");
	}

	cout << "calculate pathL26L24" << endl;
	setPath("pathL26L24", Methods::pathSearch(data, getLPointEpi(25), getLPointEpi(23), Material::subepi_Atrial_Cardium, Material::Vorhof_links, "left"));

	cout << "calculate pathL24L13" << endl;
	setPath("pathL24L13", Methods::pathSearch(data, getLPointEpi(23), getLPointEpi(12), Material::subepi_Atrial_Cardium, Material::Vorhof_links, "left"));

	cout << "calculate pathL13L22" << endl;
	setPath("pathL13L22", Methods::pathSearch(data, getLPointEpi(12), getLPointEpi(21), Material::subepi_Atrial_Cardium, Material::Vorhof_links, "left"));
	setLPointEpi(57, Methods::pointAtPercentOfPath(data, getPath("pathL13L22"), 10));
	setLPointEpi(55, Methods::pointAtPercentOfPath(data, getPath("pathL13L22"), 45));

	cout << "calculate pathL22L4" << endl;
	setPath("pathL22L4", Methods::pathSearch(data, getLPointEpi(21), getLPointEpi(3), Material::subepi_Atrial_Cardium, Material::Vorhof_links, "left"));

	setPath("pathL26L24L13L22L4L6", Methods::add5Paths(getPath("pathL26L24"), getPath("pathL24L13"), getPath("pathL13L22"), getPath("pathL22L4"), getPath("pathL4L6")));

	Methods::pathMarker(data, getPath("pathL26L24L13L22L4L6"), Material::Vorhof_links);

	if (ConfigFiberorientation::growPathTransmural) {
		setPath("pathL26L24L13L22L4L6GrownPathinEpi", Methods::growPathTransmural(data, getPath("pathL26L24L13L22L4L6"), 3.3 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links, Material::Vorhof_links_Endo, "left"));
		setPath("pathL26L24L13L22L4L6GrownPathinEndo", getPath("pathL26L24L13L22L4L6GrownPathinEpi"));
	}
	else {
		setPath("pathL26L24L13L22L4L6GrownPathinEpi", Methods::growPath(data, getPath("pathL26L24L13L22L4L6"), 3.3 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links, Material::Vorhof_links));
		setPath("pathL26L24L13L22L4L6GrownPathinEndo", Methods::growPath(data, getPath("pathL26L24L13L22L4L6"), 4.29 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links, Material::Vorhof_links_Endo));
	}

	Methods::smoothingPath(data, getPath("pathL26L24L13L22L4L6"));
	Methods::setMaterial(data, getPath(
		"pathL26L24L13L22L4L6GrownPathinEpi"), Material::subepi_Atrial_Cardium, Material::Vorhof_links);
	Methods::setMaterial(data, getPath(
		"pathL26L24L13L22L4L6GrownPathinEndo"), Material::subepi_Atrial_Cardium, Material::Vorhof_links_Endo);

	if (ConfigGeneral::debug) {
		Methods::writeIntermediateData(data, "pathL26L24L13L22L4L6_left_epi_debug");
	}

	cout << "calculate pathL7L5" << endl;
	setPath("pathL7L5", Methods::pathSearch(data, getLPointEpi(6), getLPointEpi(4), Material::subepi_Atrial_Cardium, Material::Vorhof_links, "left"));

	cout << "calculate pathL9L7" << endl;
	setPath("pathL9L7", Methods::pathSearch(data, getLPointEpi(8), getLPointEpi(6), Material::subepi_Atrial_Cardium, Material::Vorhof_links, "left"));

	cout << "calculate pathL26L9" << endl;
	setPath("pathL26L9", Methods::pathSearch(data, getLPointEpi(25), getLPointEpi(8), Material::subepi_Atrial_Cardium, Material::Vorhof_links, "left"));

	cout << "calculate pathL12L26" << endl;
	setPath("pathL12L26", Methods::pathSearch(data, getLPointEpi(11), getLPointEpi(25), Material::subepi_Atrial_Cardium, Material::Vorhof_links, "left"));

	cout << "calculate pathL5L12" << endl;
	setPath("pathL5L12", Methods::pathSearch(data, getLPointEpi(4), getLPointEpi(11), Material::subepi_Atrial_Cardium, Material::Vorhof_links, "left"));


	setPath("pathL5L12L26L9L7L5", Methods::add5Paths(getPath("pathL5L12"), getPath("pathL12L26"), getPath("pathL26L9"), getPath("pathL9L7"), getPath("pathL7L5")));
	Methods::pathMarkerRing(data, getPath("pathL5L12L26L9L7L5"), Material::Vorhof_links);

	if (ConfigFiberorientation::growPathTransmural) {
		setPath("pathL5L12L26L9L7L5GrownPathinEpi", Methods::growPathTransmural(data, getPath("pathL5L12L26L9L7L5"), 3.3 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links, Material::Vorhof_links_Endo, "left"));
		setPath("pathL5L12L26L9L7L5GrownPathinEndo", getPath("pathL5L12L26L9L7L5GrownPathinEpi"));
	}
	else {
		setPath("pathL5L12L26L9L7L5GrownPathinEpi", Methods::growPath(data, getPath("pathL5L12L26L9L7L5"), 3.3 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links, Material::Vorhof_links));
		setPath("pathL5L12L26L9L7L5GrownPathinEndo", Methods::growPath(data, getPath("pathL5L12L26L9L7L5"), 4.29 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links, Material::Vorhof_links_Endo));
	}


	Methods::smoothingRingPath(data, getPath("pathL5L12L26L9L7L5"));

	Methods::setMaterial(data, getPath(
		"pathL5L12L26L9L7L5GrownPathinEpi"), Material::subepi_Atrial_Cardium, Material::Vorhof_links);
	Methods::setMaterial(data, getPath(
		"pathL5L12L26L9L7L5GrownPathinEndo"), Material::subepi_Atrial_Cardium, Material::Vorhof_links_Endo);

	if (ConfigGeneral::debug) {
		Methods::writeIntermediateData(data, "pathL5L12L26L9L7L5_left_epi_debug");
	}

	cout << "calculate pathL5L16" << endl;
	setPath("pathL5L16", Methods::pathSearch(data, getLPointEpi(4), getLPointEpi(15), Material::subepi_Atrial_Cardium, Material::Vorhof_links, "left"));

	cout << "calculate pathL16L23" << endl;
	setPath("pathL16L23", Methods::pathSearch(data, getLPointEpi(15), getLPointEpi(22), Material::subepi_Atrial_Cardium, Material::Vorhof_links, "left"));

	cout << "calculate pathL23L25" << endl;
	setPath("pathL23L25", Methods::pathSearch(data, getLPointEpi(22), getLPointEpi(24), Material::subepi_Atrial_Cardium, Material::Vorhof_links, "left"));

	cout << "calculate pathL25L27" << endl;
	setPath("pathL25L27", Methods::pathSearch(data, getLPointEpi(24), getLPointEpi(26), Material::subepi_Atrial_Cardium, Material::Vorhof_links, "left"));

	cout << "calculate pathL27L12" << endl;
	setPath("pathL27L12", Methods::pathSearch(data, getLPointEpi(26), getLPointEpi(11), Material::subepi_Atrial_Cardium, Material::Vorhof_links, "left"));

	setPath("pathL5L16L23L25L27L12", Methods::add5Paths(getPath("pathL5L16"), getPath("pathL16L23"), getPath("pathL23L25"), getPath("pathL25L27"), getPath("pathL27L12")));

	Methods::pathMarker(data, getPath("pathL5L16L23L25L27L12"), Material::Vorhof_links);

	if (ConfigFiberorientation::growPathTransmural) {
		setPath("pathL5L16L23L25L27L12GrownPathinEpi", Methods::growPathTransmural(data, getPath("pathL5L16L23L25L27L12"), 3.3 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links, Material::Vorhof_links_Endo, "left"));
		setPath("pathL5L16L23L25L27L12GrownPathinEndo", getPath("pathL5L16L23L25L27L12GrownPathinEpi"));
	}
	else {
		setPath("pathL5L16L23L25L27L12GrownPathinEpi", Methods::growPath(data, getPath("pathL5L16L23L25L27L12"), 3.3 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links, Material::Vorhof_links));
		setPath("pathL5L16L23L25L27L12GrownPathinEndo", Methods::growPath(data, getPath("pathL5L16L23L25L27L12"), 4.29 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links, Material::Vorhof_links_Endo));
	}


	Methods::smoothingPath(data, getPath("pathL5L16L23L25L27L12"));
	Methods::setMaterial(data, getPath(
		"pathL5L16L23L25L27L12GrownPathinEpi"), Material::subepi_Atrial_Cardium, Material::Vorhof_links);
	Methods::setMaterial(data, getPath(
		"pathL5L16L23L25L27L12GrownPathinEndo"), Material::subepi_Atrial_Cardium, Material::Vorhof_links_Endo);

	if (ConfigGeneral::debug) {
		Methods::writeIntermediateData(data, "pathL5L16L23L25L27L12_left_epi_debug");
	}

	cout << "calculate pathL26L1" << endl;
	setPath("pathL26L1", Methods::pathSearch(data, getLPointEpi(0), getLPointEpi(25), Material::subepi_Atrial_Cardium, Material::Vorhof_links, "left"));
	setLPointEpi(32, Methods::pointAtPercentOfPath(data, getPath("pathL26L1"), 5));
	setPath("pathL26L33", Methods::getPercentOfPath(data, getPath("pathL26L1"), 5));

	cout << "calculate pathL33L12" << endl;
	setPath("pathL33L12", Methods::pathSearch(data, getLPointEpi(32), getLPointEpi(11), Material::subepi_Atrial_Cardium, Material::Vorhof_links, "left"));

	cout << "calculate pathL12L7" << endl;
	setPath("pathL12L7", Methods::pathSearch(data, getLPointEpi(11), getLPointEpi(6), Material::subepi_Atrial_Cardium, Material::Vorhof_links, "left"));


	setPath("pathL26L33L12L7", Methods::add3Paths(getPath("pathL26L33"), getPath("pathL33L12"), getPath("pathL12L7")));
	Methods::pathMarker(data, getPath("pathL26L33L12L7"), Material::Vorhof_links);

	if (ConfigFiberorientation::growPathTransmural) {
		setPath("pathL26L33L12L7GrownPathinEpi", Methods::growPathTransmural(data, getPath("pathL26L33L12L7"), 3.3 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links, Material::Vorhof_links_Endo, "left"));
		setPath("pathL26L33L12L7GrownPathinEndo", getPath("pathL26L33L12L7GrownPathinEpi"));
	}
	else {
		setPath("pathL26L33L12L7GrownPathinEpi", Methods::growPath(data, getPath("pathL26L33L12L7"), 3.3 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links, Material::Vorhof_links));
		setPath("pathL26L33L12L7GrownPathinEndo", Methods::growPath(data, getPath("pathL26L33L12L7"), 4.29 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links, Material::Vorhof_links_Endo));
	}


	Methods::smoothingPath(data, getPath("pathL26L33L12L7"));
	Methods::setMaterial(data, getPath(
		"pathL26L33L12L7GrownPathinEpi"), Material::subepi_Atrial_Cardium, Material::Vorhof_links);
	Methods::setMaterial(data, getPath(
		"pathL26L33L12L7GrownPathinEndo"), Material::subepi_Atrial_Cardium, Material::Vorhof_links_Endo);

	if (ConfigGeneral::debug) {
		Methods::writeIntermediateData(data, "pathL26L33L12L7_left_epi_debug");
	}
	cout << "calculate pathL8L24" << endl;
	setPath("pathL8L24", Methods::pathSearch(data, getLPointEpi(7), getLPointEpi(23), Material::subepi_Atrial_Cardium, Material::Vorhof_links, "left"));
	Methods::pathMarker(data, getPath("pathL8L24"), Material::Vorhof_links);

	if (ConfigFiberorientation::growPathTransmural) {
		setPath("pathL8L24GrownPathinEpi", Methods::growPathTransmural(data, getPath("pathL8L24"), 3.3 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links, Material::Vorhof_links_Endo, "left"));
		setPath("pathL8L24GrownPathinEndo", getPath("pathL8L24GrownPathinEpi"));
	}
	else {
		setPath("pathL8L24GrownPathinEpi", Methods::growPath(data, getPath("pathL8L24"), 3.3 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links, Material::Vorhof_links));
		setPath("pathL8L24GrownPathinEndo", Methods::growPath(data, getPath("pathL8L24"), 4.29 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links, Material::Vorhof_links_Endo));
	}

	Methods::smoothingPath(data, getPath("pathL8L24"));

	Methods::setMaterial(data, getPath("pathL8L24GrownPathinEpi"), Material::subepi_Atrial_Cardium, Material::Vorhof_links);
	Methods::setMaterial(data, getPath(
		"pathL8L24GrownPathinEndo"), Material::subepi_Atrial_Cardium, Material::Vorhof_links_Endo);

	Methods::setMaterial(data, getPath("MVRingGrownPath"), Material::Mitral_Valve_Ring, Material::subepi_Atrial_Cardium);

	cout << "calculate pathL5L4" << endl;
	setPath("pathL5L4", Methods::pathSearch(data, getLPointEpi(4), getLPointEpi(3), Material::subepi_Atrial_Cardium, Material::Vorhof_links, "left"));
	Methods::setStatus(data, getPath("pathL5L4"), 1, Material::Vorhof_links);

	if (ConfigFiberorientation::growPathTransmural) {
		setPath("pathL5L4GrownPathEpi", Methods::growPathTransmural(data, getPath("pathL5L4"), 3.63 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links, Material::Vorhof_links_Endo, "left"));
		setPath("pathL5L4GrownPathEndo", getPath("pathL5L4GrownPathEpi"));
	}
	else {
		setPath("pathL5L4GrownPathEpi", Methods::getCellsArroundPath(data, getPath("pathL5L4"), 3.63 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links));
		setPath("pathL5L4GrownPathEndo", Methods::getCellsArroundPath(data, getPath("pathL5L4"), 3.63 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links_Endo));
	}

	Methods::setMaterial(data, getPath("pathL5L4GrownPathEpi"), Material::otherwise_left, Material::Vorhof_links);
	Methods::setMaterial(data, getPath(
		"pathL5L4GrownPathEndo"), Material::otherwise_right_or_all, Material::Vorhof_links_Endo);

	if (ConfigGeneral::debug) {
		Methods::writeIntermediateData(data, "pathL5L4_left_epi_debug");
	}

	cout << "calculate area between pathL22L4L6L8L24 and pathL26L9L7L5L16" << endl;

	setPath("pathL22L4L6L8L24", Methods::add4Paths(getPath("pathL22L4"), getPath("pathL4L6"), getPath("pathL6L8"), getPath("pathL8L24")));
	setPath("pathL26L9L7L5L16", Methods::add4Paths(getPath("pathL26L9"), getPath("pathL9L7"), getPath("pathL7L5"), getPath("pathL5L16")));

	if (ConfigFiberorientation::growPathTransmural) {
		setPath("pathL22L4L6L8L24GrownPathinEpi", Methods::growPathTransmural(data, getPath("pathL22L4L6L8L24"), 4.29 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links, Material::Vorhof_links_Endo, "left"));
		setPath("pathL26L9L7L5L16GrownPathinEpi", Methods::growPathTransmural(data, getPath("pathL26L9L7L5L16"), 4.29 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links, Material::Vorhof_links_Endo, "left"));
	}
	else {
		setPath("pathL22L4L6L8L24GrownPathinEpi", Methods::growPath(data, getPath("pathL22L4L6L8L24"), 4.29 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links, Material::Vorhof_links));
		setPath("pathL26L9L7L5L16GrownPathinEpi", Methods::growPath(data, getPath("pathL26L9L7L5L16"), 4.29 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links, Material::Vorhof_links));
	}


	Methods::setMaterial(data, getPath(
		"pathL22L4L6L8L24GrownPathinEpi"), Material::subepi_Atrial_Cardium, Material::Vorhof_links);
	Methods::setMaterial(data, getPath(
		"pathL26L9L7L5L16GrownPathinEpi"), Material::subepi_Atrial_Cardium, Material::Vorhof_links);
	vector<double> seedPoint1 = Methods::getSeedPoint(data, getPath("pathL22L4L6L8L24"), getPath(
		"pathL26L9L7L5L16"), Material::Vorhof_links, Material::subepi_Atrial_Cardium, 50, "left");

	Methods::pathMarkerRing(data, getPath("pathL5L12L26L9L7L5"), Material::subepi_Atrial_Cardium);
	Methods::smoothingRingPath(data, getPath("pathL5L12L26L9L7L5"));

	setPath("areaBetweenL5L26andL4L24", Methods::regionInterpolateBetween2PathwithSeed(data, getPath("pathL22L4L6L8L24"), getPath("pathL26L9L7L5L16"), Material::Vorhof_links_Endo, Material::Vorhof_links, seedPoint1));

	Methods::setMaterial(data, getPath(
		"areaBetweenL5L26andL4L24"), Material::subepi_Atrial_Cardium, Material::Vorhof_links);
	Methods::setMaterial(data, getPath(
		"areaBetweenL5L26andL4L24"), Material::subepi_Atrial_Cardium, Material::Vorhof_links_Endo);

	Methods::setMaterial(data, getPath("pathL5L4GrownPathEpi"), Material::Vorhof_links, Material::otherwise_left);
	Methods::setMaterial(data, getPath(
		"pathL5L4GrownPathEndo"), Material::Vorhof_links_Endo, Material::otherwise_right_or_all);

	Methods::setStatus(data, getPath("pathL5L4GrownPathEpi"), 0, Material::Vorhof_links);
	Methods::setStatus(data, getPath("pathL5L4GrownPathEndo"), 0, Material::Vorhof_links_Endo);

	if (ConfigGeneral::debug) {
		Methods::writeIntermediateData(data, "areaBetweenL5L26andL4L24_left_epi_debug");
	}

	cout << "calculate area between pathL7L5L16L23L25L27L12 and pathL26L24L13L22L4L6" << endl;

	setPath("pathL7L5L16L23L25L27L12", Methods::add2Paths(getPath("pathL7L5"), getPath("pathL5L16L23L25L27L12")));
	Methods::pathMarker(data, getPath("pathL7L5L16L23L25L27L12"), Material::subepi_Atrial_Cardium);
	Methods::smoothingPath(data, getPath("pathL7L5L16L23L25L27L12"));

	vector<double> seedPoint2 = Methods::getSeedPoint(data, getPath("pathL23L25"), getPath(
		"pathL26L24"), Material::Mitral_Valve_Ring, Material::Vorhof_links, Material::subepi_Atrial_Cardium, 70, "left");

	setPath("areaBetweenL7L25L12andL26L6", Methods::regionInterpolateBetween2PathwithSeed(data, getPath("pathL7L5L16L23L25L27L12"), getPath("pathL26L24L13L22L4L6"), Material::Vorhof_links_Endo, Material::Vorhof_links, seedPoint2));

	Methods::setMaterial(data, getPath(
		"areaBetweenL7L25L12andL26L6"), Material::subepi_Atrial_Cardium, Material::Vorhof_links);
	Methods::setMaterial(data, getPath(
		"areaBetweenL7L25L12andL26L6"), Material::subepi_Atrial_Cardium, Material::Vorhof_links_Endo);
} // FiberOrientation::setFiberOrientationLeftEpiAtrium

/*! The definition part and calculation function for the fiber orientiation in the left pulmonary veins.
 \param data Pointer to the orginal mesh
 */
void FiberOrientation::setFiberOrientationLeftVeins(DataFormat &data) {
	cout << "calculate RIPV" << endl;
	cout << "calculate centroid RIPV" << endl;
	setPath("path6L8", Methods::pathSearch(data, getLPointEpi(5), getLPointEpi(7), Material::subepi_Atrial_Cardium, Material::Vorhof_links, "left"));
	setPath("pathL8L13", Methods::pathSearch(data, getLPointEpi(7), getLPointEpi(12), Material::subepi_Atrial_Cardium, Material::Vorhof_links, "left"));
	setPath("pathL13L6", Methods::pathSearch(data, getLPointEpi(12), getLPointEpi(5), Material::subepi_Atrial_Cardium, Material::Vorhof_links, "left"));
	setPath("pathL6L8L13L6", Methods::add3Paths(getPath("path6L8"), getPath("pathL8L13"), getPath("pathL13L6")));
	setLPointEpi(38, Methods::getCentroid(data, getPath("pathL6L8L13L6")));
	double maxDistanzRIPVRing = Methods::getMaxDistance(data, getPath("pathL6L8L13L6"), getLPointEpi(38));

	bool isRIPVVein = Methods::isVein(data, getLPointEpi(5), getLPointEpi(7), getLPointEpi(12), getLPointEpi(
		38), getLPointEpi(
			11), maxDistanzRIPVRing, Material::Vorhof_links, Material::Vorhof_links_Endo);

	if (isRIPVVein) {
		vector<double> averageRIPVVeinPoint = { 0, 0, 0 };
		averageRIPVVeinPoint.at(0) = (getLPointEpi(5).at(0) + getLPointEpi(7).at(0) + getLPointEpi(12).at(0)) / 3;
		averageRIPVVeinPoint.at(1) = (getLPointEpi(5).at(1) + getLPointEpi(7).at(1) + getLPointEpi(12).at(1)) / 3;
		averageRIPVVeinPoint.at(2) = (getLPointEpi(5).at(2) + getLPointEpi(7).at(2) + getLPointEpi(12).at(2)) / 3;

		maxDistanzRIPVRing = Methods::getMaxDistance(data, getPath("pathL6L8L13L6"), averageRIPVVeinPoint);
		vtkSmartPointer<vtkIdList> region = Methods::getCellsinRadius(data, averageRIPVVeinPoint, maxDistanzRIPVRing);

		setLPointEpi(39, Methods::getOuterVeinCentrepoint(getLPointEpi(5), getLPointEpi(7), getLPointEpi(12), getLPointEpi(38), 15));

		if (!data.getDoubleLayer()) {
			setPath("RIPV", Methods::growVeinInRegionOneSurface(data, Material::Vorhof_links_Endo, Material::Vorhof_links, Material::RIPV, getLPointEpi(39), 2.2, region));
		}
		else {
			setPath("RIPV", Methods::growVeinInRegion(data, Material::Vorhof_links_Endo, Material::Vorhof_links, Material::RIPV, getLPointEpi(39), region));
		}
	}
	else {
		cout << "RIPV could not be calculated" << endl;
	}

	cout << "calculate RSPV" << endl;
	cout << "calculate centroid RSPV" << endl;

	setPath("pathL4L6", Methods::pathSearch(data, getLPointEpi(3), getLPointEpi(5), Material::subepi_Atrial_Cardium, Material::Vorhof_links, "left"));
	setPath("pathL6L13", Methods::pathSearch(data, getLPointEpi(5), getLPointEpi(12), Material::subepi_Atrial_Cardium, Material::Vorhof_links, "left"));
	setPath("pathL13L22", Methods::pathSearch(data, getLPointEpi(12), getLPointEpi(21), Material::subepi_Atrial_Cardium, Material::Vorhof_links, "left"));
	setPath("pathL22L4", Methods::pathSearch(data, getLPointEpi(21), getLPointEpi(3), Material::subepi_Atrial_Cardium, Material::Vorhof_links, "left"));

	setPath("pathL4L6L13L22L4", Methods::add4Paths(getPath("pathL4L6"), getPath("pathL6L13"), getPath("pathL13L22"), getPath("pathL22L4")));

	setLPointEpi(43, Methods::getCentroid(data, getPath("pathL4L6L13L22L4")));
	double maxDistanzRSPVRing = Methods::getMaxDistance(data, getPath("pathL4L6L13L22L4"), getLPointEpi(43));

	bool isRSPVVein = Methods::isVein(data, getLPointEpi(3), getLPointEpi(5), getLPointEpi(12), getLPointEpi(
		43), getLPointEpi(
			11), maxDistanzRSPVRing, Material::Vorhof_links, Material::Vorhof_links_Endo);


	if (isRSPVVein) {
		vector<double> averageRSPVVeinPoint = { 0, 0, 0 };
		averageRSPVVeinPoint.at(0) = (getLPointEpi(3).at(0) + getLPointEpi(5).at(0) + getLPointEpi(12).at(0) + getLPointEpi(21).at(0)) / 4;
		averageRSPVVeinPoint.at(1) = (getLPointEpi(3).at(1) + getLPointEpi(5).at(1) + getLPointEpi(12).at(1) + getLPointEpi(21).at(1)) / 4;
		averageRSPVVeinPoint.at(2) = (getLPointEpi(3).at(2) + getLPointEpi(5).at(2) + getLPointEpi(12).at(2) + getLPointEpi(21).at(2)) / 4;

		maxDistanzRSPVRing = Methods::getMaxDistance(data, getPath("pathL4L6L13L22L4"), averageRSPVVeinPoint);
		vtkSmartPointer<vtkIdList> region = Methods::getCellsinRadius(data, averageRSPVVeinPoint, maxDistanzRSPVRing);

		Methods::pathMarkerRing(data, getPath("pathL4L6L13L22L4"), Material::Vorhof_links);
		if (ConfigFiberorientation::growPathTransmural) {
			setPath("pathL4L6L13L22L4GrownPathinEpi", Methods::growPathTransmural(data, getPath("pathL4L6L13L22L4"), 5.3 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links, Material::Vorhof_links_Endo, "left"));
			setPath("pathL4L6L13L22L4GrownPathinEndo", getPath("pathL4L6L13L22L4GrownPathinEpi"));
		}
		else {
			setPath("pathL4L6L13L22L4GrownPathinEndo", Methods::growPath(data, getPath("pathL4L6L13L22L4"), 6.3 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links, Material::Vorhof_links_Endo));
			setPath("pathL4L6L13L22L4GrownPathinEpi", Methods::growPath(data, getPath("pathL4L6L13L22L4"), 5.3 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links, Material::Vorhof_links));
		}

		Methods::setMaterial(data, getPath("pathL4L6L13L22L4GrownPathinEndo"), Material::RSPV, Material::Vorhof_links_Endo);
		Methods::setMaterial(data, getPath("pathL4L6L13L22L4GrownPathinEpi"), Material::RSPV, Material::Vorhof_links);

		setLPointEpi(44, Methods::getOuterVeinCentrepoint(getLPointEpi(3), getLPointEpi(5), getLPointEpi(12), getLPointEpi(43), 15));

		if (!data.getDoubleLayer()) {
			setPath("RSPV", Methods::growVeinInRegionOneSurface(data, Material::Vorhof_links_Endo, Material::Vorhof_links, Material::RSPV, getLPointEpi(44), 2.2, region));
		}
		else {
			setPath("RSPV", Methods::growVeinInRegion(data, Material::Vorhof_links_Endo, Material::Vorhof_links, Material::RSPV, getLPointEpi(44), region));
		}
	}
	else {
		cout << "RSPV could not be calculated" << endl;
	}

	cout << "calculate LSPV" << endl;

	cout << "calculate centroid LSPV" << endl;
	setPath("pathL5L7", Methods::pathSearch(data, getLPointEpi(4), getLPointEpi(6), Material::subepi_Atrial_Cardium, Material::Vorhof_links, "left"));
	setPath("pathL7L12", Methods::pathSearch(data, getLPointEpi(6), getLPointEpi(11), Material::subepi_Atrial_Cardium, Material::Vorhof_links, "left"));
	setPath("pathL12L5", Methods::pathSearch(data, getLPointEpi(11), getLPointEpi(4), Material::subepi_Atrial_Cardium, Material::Vorhof_links, "left"));
	setPath("pathL5L7L12L5", Methods::add3Paths(getPath("pathL5L7"), getPath("pathL7L12"), getPath("pathL12L5")));

	setLPointEpi(48, Methods::getCentroid(data, getPath("pathL5L7L12L5")));
	double maxDistanzLSPVRing = Methods::getMaxDistance(data, getPath("pathL5L7L12L5"), getLPointEpi(48));

	bool isLSPVVein = Methods::isVein(data, getLPointEpi(4), getLPointEpi(6), getLPointEpi(11), getLPointEpi(
		48), getLPointEpi(
			5), maxDistanzLSPVRing, Material::Vorhof_links, Material::Vorhof_links_Endo);

	if (isLSPVVein) {
		vector<double> averageLSPVVeinPoint = { 0, 0, 0 };
		averageLSPVVeinPoint.at(0) = (getLPointEpi(4).at(0) + getLPointEpi(6).at(0) + getLPointEpi(11).at(0)) / 3;
		averageLSPVVeinPoint.at(1) = (getLPointEpi(4).at(1) + getLPointEpi(6).at(1) + getLPointEpi(11).at(1)) / 3;
		averageLSPVVeinPoint.at(2) = (getLPointEpi(4).at(2) + getLPointEpi(6).at(2) + getLPointEpi(11).at(2)) / 3;

		maxDistanzLSPVRing = Methods::getMaxDistance(data, getPath("pathL5L7L12L5"), averageLSPVVeinPoint);
		vtkSmartPointer<vtkIdList> region = Methods::getCellsinRadius(data, averageLSPVVeinPoint, maxDistanzLSPVRing);


		setLPointEpi(49, Methods::getOuterVeinCentrepoint(getLPointEpi(11), getLPointEpi(6), getLPointEpi(4), getLPointEpi(48), 15));

		if (!data.getDoubleLayer()) {
			setPath("LSPV", Methods::growVeinInRegionOneSurface(data, Material::Vorhof_links_Endo, Material::Vorhof_links, Material::LSPV, getLPointEpi(49), 2.2, region));
		}
		else {
			setPath("LSPV", Methods::growVeinInRegion(data, Material::Vorhof_links_Endo, Material::Vorhof_links, Material::LSPV, getLPointEpi(49), region));
		}
	}
	else {
		cout << "LSPV could not be calculated" << endl;
	}


	cout << "calculate LIPV" << endl;
	cout << "calculate centroid LIPV" << endl;


	setPath("pathL7L9", Methods::pathSearch(data, getLPointEpi(6), getLPointEpi(8), Material::subepi_Atrial_Cardium, Material::Vorhof_links, Material::Mitral_Valve_Ring, "left"));
	setPath("pathL9L26", Methods::pathSearch(data, getLPointEpi(8), getLPointEpi(25), Material::subepi_Atrial_Cardium, Material::Vorhof_links, Material::Mitral_Valve_Ring, "left"));
	setPath("pathL26L12", Methods::pathSearch(data, getLPointEpi(25), getLPointEpi(11), Material::subepi_Atrial_Cardium, Material::Vorhof_links, Material::Mitral_Valve_Ring, "left"));
	setPath("pathL12L7", Methods::pathSearch(data, getLPointEpi(11), getLPointEpi(6), Material::subepi_Atrial_Cardium, Material::Vorhof_links, "left"));
	setPath("pathL7L9L26L12L7", Methods::add4Paths(getPath("pathL7L9"), getPath("pathL9L26"), getPath("pathL26L12"), getPath("pathL12L7")));
	setLPointEpi(53, Methods::getCentroid(data, getPath("pathL7L9L26L12L7")));

	double maxDistanzLIPVRing = Methods::getMaxDistance(data, getPath("pathL7L9L26L12L7"), getLPointEpi(53));

	bool isLIPVVein = Methods::isVein(data, getLPointEpi(6), getLPointEpi(8), getLPointEpi(11), getLPointEpi(
		53), getLPointEpi(
			5), maxDistanzLIPVRing, Material::Vorhof_links, Material::Vorhof_links_Endo);

	if (isLIPVVein) {
		vector<double> averageLIPVVeinPoint = { 0, 0, 0 };
		averageLIPVVeinPoint.at(0) = (getLPointEpi(6).at(0) + getLPointEpi(8).at(0) + getLPointEpi(25).at(0) + getLPointEpi(11).at(0)) / 4;
		averageLIPVVeinPoint.at(1) = (getLPointEpi(6).at(1) + getLPointEpi(8).at(1) + getLPointEpi(25).at(1) + getLPointEpi(11).at(1)) / 4;
		averageLIPVVeinPoint.at(2) = (getLPointEpi(6).at(2) + getLPointEpi(8).at(2) + getLPointEpi(25).at(2) + getLPointEpi(11).at(2)) / 4;

		maxDistanzLIPVRing = Methods::getMaxDistance(data, getPath("pathL7L9L26L12L7"), averageLIPVVeinPoint);
		vtkSmartPointer<vtkIdList> region = Methods::getCellsinRadius(data, averageLIPVVeinPoint, maxDistanzLIPVRing);

		setLPointEpi(54, Methods::getOuterVeinCentrepoint(getLPointEpi(8), getLPointEpi(6), getLPointEpi(11), getLPointEpi(53), 15));

		if (!data.getDoubleLayer()) {
			setPath("LIPV", Methods::growVeinInRegionOneSurface(data, Material::Vorhof_links_Endo, Material::Vorhof_links, Material::LIPV, getLPointEpi(54), 2.2, region));
		}
		else {
			setPath("LIPV", Methods::growVeinInRegion(data, Material::Vorhof_links_Endo, Material::Vorhof_links, Material::LIPV, getLPointEpi(54), region));
		}
	}
	else {
		cout << "LIPV could not be calculated" << endl;
	}
} // FiberOrientation::setFiberOrientationLeftVeins

/*! The definition part and calculation function for the fiber orientiation in the left atrial appandage.
 \param data Pointer to the orginal mesh
 */
void FiberOrientation::setFiberOrientationLeftAppendes(DataFormat &data) {
	cout << "calculate left appendages" << endl;

	cout << "calculate pathL1L11" << endl;
	setPath("pathL1L11", Methods::pathSearch(data, getLPointEpi(0), getLPointEpi(10), Material::Mitral_Valve_Ring, Material::subepi_Atrial_Cardium, Material::Vorhof_links, Material::Vorhof_links_Endo, "left"));

	cout << "calculate pathL5L11" << endl;
	setPath("pathL5L11", Methods::pathSearch(data, getLPointEpi(4), getLPointEpi(10), Material::subepi_Atrial_Cardium, Material::Vorhof_links, "left"));

	cout << "calculate pathL10L11" << endl;
	setPath("pathL10L11", Methods::pathSearch(data, getLPointEpi(9), getLPointEpi(10), Material::subepi_Atrial_Cardium, Material::Vorhof_links, "left"));


	setLPointEpi(33, Methods::pointAtMMOfPath(data, getPath("pathL1L11"), 6.93));
	setPath("pathL34L11", Methods::subPathFromPath(getPath("pathL1L11"), Methods::getMMOfPath(data, getPath("pathL1L11"), 6.93)));
	setLPointEpi(33, Methods::getPointAlongPathOfPointInMaterial(data, getPath("pathL34L11"), Material::Vorhof_links, 1.65));
	setLPointEpi(34, Methods::getPointAlongPathOfPointInMaterial(data, getPath(
		"pathL5L11"), Material::Vorhof_links, 1.65));
	setLPointEpi(35, Methods::getPointAlongPathOfPointInMaterial(data, getPath(
		"pathL10L11"), Material::Vorhof_links, 1.65));
	setLPointEpi(36, Methods::findClosedPointinMaterialinDirectionPoint(data, getLPointEpi(11), getLPointEpi(33), 6, Material::Vorhof_links));

	cout << "calculate pathL34L37" << endl;
	setPath("pathL34L37", Methods::pathSearch(data, getLPointEpi(33), getLPointEpi(36), Material::Vorhof_links, Material::Vorhof_links_Endo, Material::LIPV, Material::subepi_Atrial_Cardium, "left"));

	cout << "calculate pathL37L35" << endl;
	setPath("pathL37L35", Methods::pathSearch(data, getLPointEpi(36), getLPointEpi(34), Material::Vorhof_links, Material::Vorhof_links_Endo, Material::LIPV, Material::subepi_Atrial_Cardium, "left"));

	cout << "calculate pathL35L36" << endl;
	setPath("pathL35L36", Methods::pathSearch(data, getLPointEpi(34), getLPointEpi(35), Material::Vorhof_links, Material::Vorhof_links_Endo, "left"));

	cout << "calculate pathL36L34" << endl;
	setPath("pathL36L34", Methods::pathSearch(data, getLPointEpi(35), getLPointEpi(33), Material::Vorhof_links, Material::Vorhof_links_Endo, "left"));

	setPath("appendageLeftRing", Methods::add4Paths(getPath("pathL34L37"), getPath("pathL37L35"), getPath("pathL35L36"), getPath("pathL36L34")));
	Methods::pathMarkerRing(data, getPath("appendageLeftRing"), Material::Vorhof_links, Material::Vorhof_links_Endo);
	Methods::smoothingRingPath(data, getPath("appendageLeftRing"));

	setPath("borderAppendagesLeftRing", Methods::boundaryLayerAppend(data, getPath("appendageLeftRing"), getLPointEpi(10)));
	Methods::setStatus(data, getPath("borderAppendagesLeftRing"), 1, Material::Vorhof_links, Material::Vorhof_links_Endo);

	Methods::pathMarkerRing(data, getPath("pathL4L6L13L22L4"), Material::Vorhof_links);

	if (ConfigFiberorientation::growPathTransmural) {
		setPath("growBorderLeftAppendagesRingEpi", Methods::growPathTransmural(data, getPath("borderAppendagesLeftRing"), 2.6 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links, Material::Vorhof_links_Endo, "left"));
		setPath("growBorderLeftAppendagesRingEndo", getPath("growBorderLeftAppendagesRingEpi"));
	}
	else {
		setPath("growBorderLeftAppendagesRingEpi", Methods::growPath(data, getPath("borderAppendagesLeftRing"), 2.6 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links, Material::Vorhof_links));
		setPath("growBorderLeftAppendagesRingEndo", Methods::growPath(data, getPath("borderAppendagesLeftRing"), 2.6 + ConfigFiberorientation::growPathOffsetLeft, Material::Vorhof_links_Endo, Material::Vorhof_links_Endo));
	}


	setLPointEpi(10, Methods::findClosedPointinMaterialwithoutOrientation(data, getLPointEpi(10), Material::Vorhof_links));

	setPath("growBorderLeftAppendagesRing", Methods::add2Paths(getPath("growBorderLeftAppendagesRingEndo"), getPath("growBorderLeftAppendagesRingEpi")));
	setPath("leftAppendages", Methods::growPathtoPoint(data, getPath("growBorderLeftAppendagesRing"), getLPointEpi(10), Material::Vorhof_links, Material::Vorhof_links_Endo));
	Methods::setMaterial(data, getPath("leftAppendages"), Material::Left_Atrial_Appendage, Material::Vorhof_links);
	Methods::setMaterial(data, getPath("leftAppendages"), Material::Left_Atrial_Appendage, Material::Vorhof_links_Endo);

	cout << "left appendages calculated" << endl;
} // FiberOrientation::setFiberOrientationLeftAppendes

/*! The definition part for calculation and insertion of the interatrial bridges.
 \param data Pointer to the orginal mesh
 \param freeBridgePath path of the txt-file with the defintion of the free bridge points
 */
void FiberOrientation::setBridges(DataFormat &data, string freeBridgePath) {
	// Bachmann

	if (ConfigFiberorientation::bachmannBundle) {
		cout << "calculate BB" << endl;

		setPath("bridgeR2R29", Methods::pathSearch(data, getRPoint(2), getRPoint(29), "right"));
		setPath("bridgeR29R7", Methods::pathSearch(data, getRPoint(29), getRPoint(6), "right"));

		setPath("bridgeR2R29R7", Methods::add2Paths(getPath("bridgeR2R29"), getPath("bridgeR29R7")));

		Methods::pathMarker(data, getPath("bridgeR2R29R7"));
		setPath("bridgeR2R29R7Grown", Methods::growPathInRight(data, getPath("bridgeR2R29R7"), 2));
		Methods::smoothingPath(data, getPath("bridgeR2R29R7"));
		Methods::setMaterialinRight(data, getPath("bridgeR2R29R7Grown"), Material::Bachmann_Bundle_right_atrium_or_all);

		setPath("bridgeL1L10", Methods::pathSearch(data, getLPointEpi(0), getLPointEpi(9), "left"));
		setPath("bridgeL10L5", Methods::pathSearch(data, getLPointEpi(9), getLPointEpi(4), "left"));

		setPath("bridgeL1L10L5", Methods::add2Paths(getPath("bridgeL1L10"), getPath("bridgeL10L5")));

		Methods::pathMarker(data, getPath("bridgeL1L10L5"));
		setPath("bridgeL1L10L5Grown", Methods::growPathInLeft(data, getPath("bridgeL1L10L5"), 2));
		Methods::smoothingPath(data, getPath("bridgeL1L10L5"));
		Methods::setMaterial(data, getPath("bridgeL1L10L5Grown"), Material::Bachmann_Bundle_left_atrium);

		setBachmannBridgelong(data, getLPointEndo(9), getRPoint(
			29), 0, 2.31, 8, Material::Bachmann_Bundle_left_atrium, Material::Bachmann_Bundle_right_atrium_or_all);

		data = update(data);
		closeOrientation(data);
    
        if (ConfigGeneral::debug) {
            Methods::writeIntermediateData(data, "BachmannBundle_debug");
        }
    
    }
    
    
	double searchradius = 2 * getResolution();
	if (searchradius < 1) {
		searchradius = 1;
	}

	if (ConfigFiberorientation::upperPosteriorBridge) {
		cout << "calculate upper posterior bridge" << endl;

		setBridge(data, getLPointEpi(55), getRPoint(
			30), 1, 1.65, 1.65, 1.65, Material::Upper_Posterior_Bridge_left, Material::Upper_Posterior_Bridge_right, 2 * searchradius);
		data = update(data);
		closeOrientation(data);
	
        if (ConfigGeneral::debug) {
            Methods::writeIntermediateData(data, "upperPosteriorBridge_debug");
        }
    }

	if (ConfigFiberorientation::coronarySinusBridge) {
		cout << "calculate coronary Sinus bridge" << endl;

		setBridge(data, getLPointEpi(1), getRPoint(
			13), 2, 1.65, 1.65, 1.65, Material::Coronary_Sinus_Bridge_left, Material::Coronary_Sinus_Bridge_right, 3 * searchradius);
		data = update(data);
		closeOrientation(data);
        
        if (ConfigGeneral::debug) {
            Methods::writeIntermediateData(data, "coronarySinusBridge_debug");
        }
	}

	if (ConfigFiberorientation::middlePosteriorBridge) {
		cout << "calculate middle posterior bridge" << endl;

		setBridge(data, getLPointEpi(57), getRPoint(
			28), 3, 1.65, 1.65, 1.65, Material::Middle_Posterior_Bridge_left, Material::Middle_Posterior_Bridge_right, 2 * searchradius);
		data = update(data);
		closeOrientation(data);
        
        if (ConfigGeneral::debug) {
            Methods::writeIntermediateData(data, "middlePosteriorBridge_debug");
        }
	}

	if (ConfigFiberorientation::lowerAnteriorBridge) {
		cout << "calculate lower anterior bridge" << endl;

		setBridge(data, getLPointEpi(56), getRPoint(
			39), 4, 1.65, 1.65, 1.65, Material::Lower_Anterior_Bridge_left, Material::Lower_Anterior_Bridge_right, 2 * searchradius);
		data = update(data);
		closeOrientation(data);
        
        if (ConfigGeneral::debug) {
            Methods::writeIntermediateData(data, "lowerAnteriorBridge_debug");
        }
	}

	if (ConfigFiberorientation::upperAnteriorBridge) {
		cout << "calculate upper anterior bridge" << endl;

		setBridge(data, getLPointEpi(21), getRPoint(
			1), 5, 1.65, 1.65, 1.65, Material::Upper_Anterior_Bridge_left, Material::Upper_Anterior_Bridge_right, 2 * searchradius);
		data = update(data);
		closeOrientation(data);
        
        if (ConfigGeneral::debug) {
            Methods::writeIntermediateData(data, "upperAnteriorBridge_debug");
        }
        
	}
	if (ConfigFiberorientation::freeBridge) {
		int numberOfFreeBridges = readOutFreeBridgePointsPoints(freeBridgePath);

		for (int i = 0; i < numberOfFreeBridges; i++) {
			setFreeBridge(data, getFreeBridgePoint(i * 2), getFreeBridgePoint(i * 2 + 1), 6 + i, freeBridgeWidth.at(
				i), freeBridgeMaterial.at(i));
			data = update(data);
			closeOrientationFreeBridge(data, freeBridgeMaterial.at(i));
            
            if (ConfigGeneral::debug) {
                Methods::writeIntermediateData(data, "freeBridge_"+std::to_string(i)+"_debug");
            }
		}
	}
} // FiberOrientation::setBridges

/*! Function for insertion a interatrial bridge into a mesh.
 \param data Pointer to the orginal mesh
 \param pointLeft left point of the bridge
 \param pointRight right point of the bridge
 \param bridgeNumber number of the bridge
 \param bridgewidth radius of the inserted bridge
 \param bridgeLengthinRight length of the bridge in the mesh of the right atrium
 \param bridgeLengthinLeft length of the bridge in the mesh of the left atrium
 \param bridgeMaterialLeft tissue class of the left and inserted bridge
 \param bridgeMaterialRight tissue class of the right bridge
 \param searchRadius search radius for the origin point in the orginal atrial mesh
 */
void FiberOrientation::setBridge(DataFormat &data, vector<double> pointLeft, vector<double> pointRight, int bridgeNumber, double bridgewidth, double bridgeLengthinRight, double bridgeLengthinLeft, Material::Mat bridgeMaterialLeft, Material::Mat bridgeMaterialRight, double searchRadius) {
	double resolution = getResolution();

	vector<double> bridgeTempPoint1 = Methods::findClosedPointinMaterialinRadiustoPointLeft(data, pointLeft, pointRight, searchRadius, resolution);
	vector<double> bridgeTempPoint2 = Methods::findClosedPointinMaterialinRadiustoPointRight(data, pointRight, pointLeft, searchRadius, resolution);


	vector<double> bridgeTempPoint3 = { (bridgeTempPoint1.at(0) + bridgeTempPoint2.at(0)) / 2, (bridgeTempPoint1.at(1) + bridgeTempPoint2.at(1)) / 2,  (bridgeTempPoint1.at(2) + bridgeTempPoint2.at(2)) / 2 };
	int pointPos = bridgeNumber * 2;

	DataFormat pillForRegion = Methods::getPill(pointLeft, pointRight, resolution, searchRadius);
	vtkSmartPointer<vtkIdList> region = Methods::getCellsWhichWereInside(data, pillForRegion);

	setBPoint(pointPos, Methods::findClosedPointinMaterialInLeftInRegion(data, bridgeTempPoint3, region));
	setBPoint(pointPos + 1, Methods::findClosedPointinMaterialInRightInRegion(data, bridgeTempPoint3, region));

	setPath("tempbridge" + to_string(bridgeNumber) + ".1", Methods::pathSearch(data, pointLeft, getBPoint(
		pointPos), "left"));

	setPath("bridge" + to_string(bridgeNumber) + ".1", Methods::getMMInverseOfPath(data, getPath("tempbridge" + to_string(bridgeNumber) + ".1"), bridgeLengthinLeft));

	setPath("tempbridge" + to_string(bridgeNumber) + ".3", Methods::pathSearch(data, getBPoint(pointPos + 1), pointRight, "right"));

	setPath("bridge" + to_string(bridgeNumber) + ".3", Methods::getMMOfPath(data, getPath("tempbridge" + to_string(bridgeNumber) + ".3"), bridgeLengthinRight));


	if (data.getInputType() == DataFormat::PossibleInputType::vtu) {
		vtkSmartPointer<vtkIdList> cellsInRadiusPointLeft = Methods::getCellsinRadius(data, data.getCentrePoints()->FindPoint(
			getBPoint(pointPos + 1).at(0), getBPoint(pointPos + 1).at(1), getBPoint(pointPos + 1).at(
				2)), bridgewidth);
		vtkSmartPointer<vtkIdList> cellsInRadiusPointRight = Methods::getCellsinRadius(data, data.getCentrePoints()->FindPoint(
			getBPoint(pointPos).at(0), getBPoint(pointPos).at(1), getBPoint(pointPos).at(
				2)), bridgewidth);

		vtkSmartPointer<vtkIdList> Vorhof_linksIds = vtkSmartPointer<vtkIdList>::New();

		for (vtkIdType i = 0; i < cellsInRadiusPointLeft->GetNumberOfIds(); i++) {
			int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(cellsInRadiusPointLeft->GetId(i), 0);
			if (Material::isInLeft(material)) {
				Vorhof_linksIds->InsertNextId(cellsInRadiusPointLeft->GetId(i));
			}
		}
		vtkSmartPointer<vtkIdList> Vorhof_rechtsIds = vtkSmartPointer<vtkIdList>::New();

		for (vtkIdType i = 0; i < cellsInRadiusPointRight->GetNumberOfIds(); i++) {
			int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(cellsInRadiusPointRight->GetId(
				i), 0);
			if (Material::isInRight(material)) {
				Vorhof_rechtsIds->InsertNextId(cellsInRadiusPointRight->GetId(i));
			}
		}

		unordered_map<vtkIdType, vtkIdType> commonPoints = Methods::getCommonPoints(data, Vorhof_linksIds, Vorhof_rechtsIds);
		cout << "direct connection points: " << commonPoints.size() << endl;
		if (commonPoints.size() > 0) {
			Methods::replacePoints(data, commonPoints);
		}
		if (data.getVtkData()->GetCellType(data.getVtkData()->GetNumberOfCells() / 2) != 10) {
			DataFormat cylinder;
			cylinder = Methods::createVoxelTube(data, getBPoint(pointPos), getBPoint(pointPos + 1), resolution, bridgewidth);
			Methods::init(cylinder);

			vtkSmartPointer<vtkIdList> diffcellIDs = Methods::getCellsWereNotinOldData(cylinder, data);
			if (diffcellIDs->GetNumberOfIds() > 0) {
				vtkSmartPointer<vtkUnstructuredGrid> diffCylinder = Methods::extractCellsVtu(cylinder.getVtkData(), diffcellIDs);
				DataFormat newVoxel;
				newVoxel.setVtkData(diffCylinder);
				newVoxel.setCentrePoints(Methods::calcCellCentroids(diffCylinder));

				diffcellIDs = Methods::getAllCellNeighbours(newVoxel, newVoxel.getCentrePoints()->FindPoint(bridgeTempPoint3.at(0), bridgeTempPoint3.at(1), bridgeTempPoint3.at(2)));
				diffCylinder = Methods::extractCellsVtu(newVoxel.getVtkData(), diffcellIDs);

				newVoxel.setVtkData(diffCylinder);
				newVoxel.setCentrePoints(Methods::calcCellCentroids(diffCylinder));

				Methods::unionData(data, newVoxel);
			}
		}
		else if (data.getVtkData()->GetCellType(data.getVtkData()->GetNumberOfCells() / 2) == 10) {
			DataFormat pill;
			pill = Methods::getPill(getBPoint(pointPos), getBPoint(pointPos + 1), resolution, bridgewidth);

			vtkSmartPointer<vtkPolyData> surface = data.getSurfaceForBridges();

			vtkSmartPointer<vtkPolyData> pillvtp = Methods::vtpDifference(vtkPolyData::SafeDownCast(
				pill.getVtkData()), surface);

			pill.setVtkData(pillvtp);
			pill.setCentrePoints(Methods::calcCellCentroids(pillvtp));
			pill.setInputType(DataFormat::vtp);

			if (ConfigGeneral::debug) {
				Methods::writeIntermediateData(pill, "bridge_difference_debug");
			}

			vtkSmartPointer<vtkIdList> diffcellIDs = Methods::getAllCellNeighbours(pill, pill.getCentrePoints()->FindPoint(bridgeTempPoint3.at(0), bridgeTempPoint3.at(1), bridgeTempPoint3.at(2)));
			pillvtp = Methods::extractCellsVtp(pillvtp, diffcellIDs);

			/*bool hasSelfIntersect = vcgTools::VcgTools::hasSelfIntersections(pillvtp);
			 if (hasSelfIntersect) {
			 cout << "The Mesh has self- intersections" << endl;
			 } Error mesh.isct.tpp line 1063:Ran out of tries to perturb the mesh*/

			bool hasHoles = Methods::hasHoles(pillvtp);
			if (hasHoles) {
				cout << "The Mesh has holes" << endl;
				cout << "write stl" << endl;
				vtkSmartPointer<vtkSTLWriter> stlWriter = vtkSmartPointer<vtkSTLWriter>::New();
				stlWriter->SetInputData(pillvtp);
				stlWriter->SetFileTypeToBinary();
				stlWriter->SetFileName("Bridge.stl");
				stlWriter->Write();

				string fixedMeshPath;
				cout << "Check if the mesh in the file Bridge.stl is closed and has no self-intersections" << endl;
				cout << "Please enter the path where the fixed mesh is: " << endl;

				bool exit = false;
				while (!exit) {
					cin >> fixedMeshPath;

					ifstream FileTest(fixedMeshPath.c_str());

					if (FileTest) {
						exit = true;
					}
					else {
						cout << "write stl" << endl;
						vtkSmartPointer<vtkSTLWriter> stlWriter = vtkSmartPointer<vtkSTLWriter>::New();
						stlWriter->SetInputData(pillvtp);
						stlWriter->SetFileTypeToBinary();
						stlWriter->SetFileName("Bridge.stl");
						stlWriter->Write();

						cout << "try again:" << endl;
					}

					FileTest.close();
				}

				cout << "read stl" << endl;

				vtkSmartPointer<vtkSTLReader> stlReader = vtkSmartPointer<vtkSTLReader>::New();
				stlReader->SetFileName(fixedMeshPath.c_str());
				stlReader->Update();

				pillvtp = stlReader->GetOutput();
			}

			if (pillvtp->GetNumberOfCells() > 0) {
				pillvtp = Methods::closeSurfaceFilter(pillvtp);
				vtkSmartPointer<vtkUnstructuredGrid> pillvtu = Methods::tetrahedralizeVTPSurface(pillvtp);

				pill.setVtkData(pillvtu);
				pill.setCentrePoints(Methods::calcCellCentroids(pillvtu));
				pill.setInputType(DataFormat::vtu);

				Methods::unionData(data, pill);
			}
		}
		else {
			cerr << "Bridge: celltype undefined" << endl;
		}

		setPath("bridge" + to_string(bridgeNumber) + ".2", Methods::pathSearch(data, getBPoint(pointPos), getBPoint(pointPos + 1), "right/left"));

		setPath("bridge" + to_string(bridgeNumber), Methods::add3Paths(getPath("bridge" + to_string(bridgeNumber) + ".1"), getPath("bridge" + to_string(bridgeNumber) + ".2"), getPath("bridge" + to_string(bridgeNumber) + ".3")));

		Methods::pathMarker(data, getPath("bridge" + to_string(bridgeNumber)));

		setPath("growbridge" + to_string(bridgeNumber), Methods::growPath(data, getPath("bridge" + to_string(bridgeNumber)), bridgewidth));
		Methods::smoothingPath(data, getPath("bridge" + to_string(bridgeNumber)));

		Methods::setMaterialinLeft(data, getPath("growbridge" + to_string(bridgeNumber)), bridgeMaterialLeft);
		Methods::setMaterialinRight(data, getPath("growbridge" + to_string(bridgeNumber)), bridgeMaterialRight);
		Methods::setMaterial(data, getPath("growbridge" + to_string(
			bridgeNumber)), bridgeMaterialLeft, Material::insert_new_Material);
	}
	else if (data.getInputType() == DataFormat::PossibleInputType::vtp) {
		setPath("bridge" + to_string(bridgeNumber), Methods::add2Paths(getPath("bridge" + to_string(bridgeNumber) + ".1"), getPath("bridge" + to_string(bridgeNumber) + ".3")));

		if (getPath("bridge" + to_string(bridgeNumber))->GetNumberOfIds() > 0) {
			Methods::pathMarker(data, getPath("bridge" + to_string(bridgeNumber)), Material::insert_new_Material);

			setPath("growbridge" + to_string(bridgeNumber), Methods::growPath(data, getPath("bridge" + to_string(bridgeNumber)), bridgewidth, Material::insert_new_Material, Material::insert_new_Material));
			Methods::smoothingPath(data, getPath("bridge" + to_string(bridgeNumber)));

			Methods::setMaterial(data, getPath("growbridge" + to_string(
				bridgeNumber)), bridgeMaterialLeft, Material::insert_new_Material);
		}

		DataFormat pill;
		pill = Methods::getPill(getBPoint(pointPos), getBPoint(pointPos + 1), resolution, bridgewidth);

        if (ConfigGeneral::debug) {
            Methods::writeIntermediateData(pill, "pillBridge_debug");
        }
        
        DataFormat helpBridge;
		helpBridge.setInputType(DataFormat::vtp);
		helpBridge.setVtkData(Methods::vtpUnion(vtkPolyData::SafeDownCast(data.getVtkData()), vtkPolyData::SafeDownCast(pill.getVtkData())));
		helpBridge.setCentrePoints(Methods::calcCellCentroids(helpBridge.getVtkData()));

        if (ConfigGeneral::debug) {
            Methods::writeIntermediateData(helpBridge, "helpBridgeBridge1_debug");
        }
        
		vtkSmartPointer<vtkIdList> removeCells = Methods::getCellsWereInOldData(helpBridge, data); // , centrepoint, // radius);
		vtkSmartPointer<vtkPolyData> polyDataBridge = vtkPolyData::SafeDownCast(helpBridge.getVtkData());
		for (vtkIdType i = 0; i < removeCells->GetNumberOfIds(); i++) {
			polyDataBridge->DeleteCell(removeCells->GetId(i));
		}
		polyDataBridge->RemoveDeletedCells();
		polyDataBridge->BuildLinks();

		vtkSmartPointer<vtkCleanPolyData> cleanFilterPolyDataBridge = vtkSmartPointer<vtkCleanPolyData>::New();
		cleanFilterPolyDataBridge->SetInputData(polyDataBridge);
		cleanFilterPolyDataBridge->Update();
		helpBridge.setVtkData(cleanFilterPolyDataBridge->GetOutput());
        helpBridge.setCentrePoints(Methods::calcCellCentroids(helpBridge.getVtkData()));

        vtkSmartPointer<vtkIdList> diffcellIDs = Methods::getAllCellNeighbours(helpBridge, helpBridge.getCentrePoints()->FindPoint(bridgeTempPoint3.at(0), bridgeTempPoint3.at(1), bridgeTempPoint3.at(2)));

        helpBridge.setVtkData(Methods::extractCellsVtp(cleanFilterPolyDataBridge->GetOutput(), diffcellIDs));
        helpBridge.setCentrePoints(Methods::calcCellCentroids(helpBridge.getVtkData()));
        
		if (ConfigFiberorientation::openVTPBridgeInAtrial) {
			set<vtkIdType> removedCells;

			vtkSmartPointer<vtkIdList> removeDataCells = Methods::getPointsWereInOldData(data, pill);
			vtkSmartPointer<vtkPolyData> polyData = vtkPolyData::SafeDownCast(data.getVtkData());
			for (vtkIdType i = 0; i < removeDataCells->GetNumberOfIds(); i++) {
				vtkSmartPointer<vtkIdList> cells = vtkSmartPointer<vtkIdList>::New();
				polyData->GetPointCells(removeDataCells->GetId(i), cells);
				for (vtkIdType j = 0; j < cells->GetNumberOfIds(); j++) {
					if (removedCells.find(cells->GetId(j)) == removedCells.end()) {
						polyData->DeleteCell(cells->GetId(j));
						data.getCentrePoints()->GetPoints()->SetPoint(cells->GetId(j), NAN, NAN, NAN);
						removedCells.insert(cells->GetId(j));
					}
				}
			}

			polyData->RemoveDeletedCells();
			polyData->BuildLinks();

			vtkSmartPointer<vtkCleanPolyData> cleanFilter = vtkSmartPointer<vtkCleanPolyData>::New();
			cleanFilter->SetInputData(polyData);
			cleanFilter->Update();

			data.setVtkData(cleanFilter->GetOutput());
            
            if (ConfigGeneral::debug) {
                Methods::writeIntermediateData(helpBridge, "helpBridgeBridge2_debug");
            }
		}

		setPath("bridge" + to_string(bridgeNumber) + ".2", Methods::unionData(data, helpBridge));

		Methods::vtpBridgeMarker(data, getPath("bridge" + to_string(
			bridgeNumber) + ".2"), Material::insert_new_Material, getBPoint(
				pointPos), getBPoint(pointPos + 1));
		Methods::setMaterial(data, getPath("bridge" + to_string(
			bridgeNumber) + ".2"), bridgeMaterialLeft, Material::insert_new_Material);
	}
} // FiberOrientation::setBridge

/*! Function for insertion a free interatrial bridge into a mesh.
 \param data Pointer to the orginal mesh
 \param point1 point 1 of the bridge
 \param point2 point 2 of the bridge
 \param bridgeNumber number of the bridge
 \param bridgewidth radius of the inserted bridge
 \param bridgeMaterial tissue class of the inserted bridge
 */
void FiberOrientation::setFreeBridge(DataFormat &data, vector<double> point1, vector<double> point2, int bridgeNumber, double bridgewidth, int bridgeMaterial) {
	double resolution = getResolution();

	int pointPos = bridgeNumber * 2;

	if (data.getInputType() == DataFormat::PossibleInputType::vtu) {
		double point1Array[3] = { point1.at(0), point1.at(1), point1.at(2) };
		double point2Array[3] = { point2.at(0), point2.at(1), point2.at(2) };

		double distance = sqrt(abs(vtkMath::Distance2BetweenPoints(point1Array, point2Array)));

		bool calcCommonPoints = false;
		if (distance < 2 * bridgewidth) {
			vtkIdType point1Id = data.getCentrePoints()->FindPoint(point1Array);
			int material1 = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(point1Id, 0);

			vtkIdType point2Id = data.getCentrePoints()->FindPoint(point2Array);
			int material2 = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(point2Id, 0);

			if ((Material::isInLeft(material1) && Material::isInRight(material2)) ||
				(Material::isInLeft(material2) && Material::isInRight(material1))) {
				if (Material::isInLeft(material1)) {
					setBPoint(pointPos, point1);
				}
				else if (Material::isInRight(material1)) {
					setBPoint(pointPos + 1, point1);
				}

				if (Material::isInLeft(material2)) {
					setBPoint(pointPos, point2);
				}
				else if (Material::isInRight(material2)) {
					setBPoint(pointPos + 1, point2);
				}

				calcCommonPoints = true;
			}
		}
		else {
			vtkIdType point1Id = data.getCentrePoints()->FindPoint(point1Array);
			setBPoint(pointPos, data.getCentrePoints()->GetPoint(point1Id)[0], data.getCentrePoints()->GetPoint(
				point1Id)[1], data.getCentrePoints()->GetPoint(point1Id)[2]);
			setBPoint(pointPos + 1, point2);
		}

		if (calcCommonPoints) {
			vtkSmartPointer<vtkIdList> cellsInRadiusPointLeft = Methods::getCellsinRadius(data, data.getCentrePoints()->FindPoint(
				getBPoint(pointPos + 1).at(0), getBPoint(pointPos + 1).at(1), getBPoint(pointPos + 1).at(
					2)), bridgewidth);
			vtkSmartPointer<vtkIdList> cellsInRadiusPointRight = Methods::getCellsinRadius(data, data.getCentrePoints()->FindPoint(
				getBPoint(pointPos).at(0), getBPoint(pointPos).at(1), getBPoint(pointPos).at(
					2)), bridgewidth);

			vtkSmartPointer<vtkIdList> Vorhof_linksIds = vtkSmartPointer<vtkIdList>::New();

			for (vtkIdType i = 0; i < cellsInRadiusPointLeft->GetNumberOfIds(); i++) {
				int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(cellsInRadiusPointLeft->GetId(
					i), 0);
				if (Material::isInLeft(material)) {
					Vorhof_linksIds->InsertNextId(cellsInRadiusPointLeft->GetId(i));
				}
			}
			vtkSmartPointer<vtkIdList> Vorhof_rechtsIds = vtkSmartPointer<vtkIdList>::New();

			for (vtkIdType i = 0; i < cellsInRadiusPointRight->GetNumberOfIds(); i++) {
				int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(cellsInRadiusPointRight->GetId(
					i), 0);
				if (Material::isInRight(material)) {
					Vorhof_rechtsIds->InsertNextId(cellsInRadiusPointRight->GetId(i));
				}
			}

			unordered_map<vtkIdType, vtkIdType> commonPoints = Methods::getCommonPoints(data, Vorhof_linksIds, Vorhof_rechtsIds);
			cout << "direct conection points: " << commonPoints.size() << endl;
			if (commonPoints.size() > 0) {
				Methods::replacePoints(data, commonPoints);
			}
		}
		if (data.getVtkData()->GetCellType(data.getVtkData()->GetNumberOfCells() / 2) != 10) {
			DataFormat cylinder;
			cylinder = Methods::createVoxelTube(data, getBPoint(pointPos), getBPoint(pointPos + 1), resolution, bridgewidth);
			Methods::init(cylinder);

			vtkSmartPointer<vtkIdList> diffcellIDs = Methods::getCellsWereNotinOldData(cylinder, data);
			if (diffcellIDs->GetNumberOfIds() > 0) {
				vtkSmartPointer<vtkUnstructuredGrid> diffCylinder = Methods::extractCellsVtu(cylinder.getVtkData(), diffcellIDs);
				DataFormat newVoxel;
				newVoxel.setVtkData(diffCylinder);
				newVoxel.setCentrePoints(Methods::calcCellCentroids(diffCylinder));

				Methods::unionData(data, newVoxel);
			}
		}
		else if (data.getVtkData()->GetCellType(data.getVtkData()->GetNumberOfCells() / 2) == 10) {
			DataFormat pill;
			pill = Methods::getPill(getBPoint(pointPos), getBPoint(pointPos + 1), resolution, bridgewidth);

			vtkSmartPointer<vtkPolyData> surface = Methods::extractSurface(data);

			vtkSmartPointer<vtkPolyData> pillvtp = Methods::vtpDifference(vtkPolyData::SafeDownCast(
				pill.getVtkData()), surface);
			vtkSmartPointer<vtkUnstructuredGrid> pillvtu = vtkSmartPointer<vtkUnstructuredGrid>::New();

			/*bool hasSelfIntersect = vcgTools::VcgTools::hasSelfIntersections(pillvtp);
			 if (hasSelfIntersect) {
			 cout << "The Mesh has self- intersections" << endl;
			 }*/
			bool hasHoles = Methods::hasHoles(pillvtp);
			if (hasHoles) {
				cout << "The Mesh has holes" << endl;
				cout << "write stl" << endl;
				vtkSmartPointer<vtkSTLWriter> stlWriter = vtkSmartPointer<vtkSTLWriter>::New();
				stlWriter->SetInputData(pillvtp);
				stlWriter->SetFileTypeToBinary();
				stlWriter->SetFileName("Bridge.stl");
				stlWriter->Write();

				string fixedMeshPath;
				cout << "Check if the mesh in the file Bridge.stl is closed and has no self-intersections" << endl;
				cout << "Please enter the path where the fixed mesh is: " << endl;

				bool exit = false;
				while (!exit) {
					cin >> fixedMeshPath;

					ifstream FileTest(fixedMeshPath.c_str());

					if (FileTest) {
						exit = true;
					}
					else {
						cout << "write stl" << endl;
						vtkSmartPointer<vtkSTLWriter> stlWriter = vtkSmartPointer<vtkSTLWriter>::New();
						stlWriter->SetInputData(pillvtp);
						stlWriter->SetFileTypeToBinary();
						stlWriter->SetFileName("Bridge.stl");
						stlWriter->Write();

						cout << "try again:" << endl;
					}

					FileTest.close();
				}

				cout << "read stl" << endl;

				vtkSmartPointer<vtkSTLReader> stlReader = vtkSmartPointer<vtkSTLReader>::New();
				stlReader->SetFileName(fixedMeshPath.c_str());
				stlReader->Update();

				if (pillvtp->GetNumberOfCells() > 0) {
					pillvtp = Methods::closeSurfaceFilter(pillvtp);

					pillvtu = Methods::tetrahedralizeVTPSurface(pillvtp);
				}
			}

			if (pillvtp->GetNumberOfCells() > 0) {
				pill.setVtkData(pillvtu);
				pill.setCentrePoints(Methods::calcCellCentroids(pillvtu));
				pill.setInputType(DataFormat::vtu);

				Methods::unionData(data, pill);
			}
		}
		else {
			cerr << "Bridge: celltype undefined" << endl;
		}

		setPath("bridge" + to_string(bridgeNumber), Methods::pathSearch(data, getBPoint(pointPos), getBPoint(pointPos + 1), "right/left"));
		Methods::pathMarker(data, getPath("bridge" + to_string(bridgeNumber)));

		setPath("growbridge" + to_string(bridgeNumber), Methods::growPath(data, getPath("bridge" + to_string(bridgeNumber)), bridgewidth));
		Methods::smoothingPath(data, getPath("bridge" + to_string(bridgeNumber)));

		Methods::setMaterial(data, getPath("growbridge" + to_string(
			bridgeNumber)), bridgeMaterial, Material::insert_new_Material);
	}
	else if (data.getInputType() == DataFormat::PossibleInputType::vtp) {
		setBPoint(pointPos, point1);
		setBPoint(pointPos + 1, point2);

		DataFormat pill;
		pill = Methods::getPill(getBPoint(pointPos), getBPoint(pointPos + 1), resolution, bridgewidth);

		DataFormat helpBridge;
		helpBridge.setInputType(DataFormat::vtp);
		helpBridge.setVtkData(Methods::vtpUnion(vtkPolyData::SafeDownCast(data.getVtkData()), vtkPolyData::SafeDownCast(pill.getVtkData())));
		helpBridge.setCentrePoints(Methods::calcCellCentroids(helpBridge.getVtkData()));

		vtkSmartPointer<vtkIdList> removeCells = Methods::getCellsWereInOldData(helpBridge, data); // , centrepoint, // radius);
		vtkSmartPointer<vtkPolyData> polyDataBridge = vtkPolyData::SafeDownCast(helpBridge.getVtkData());
		for (vtkIdType i = 0; i < removeCells->GetNumberOfIds(); i++) {
			polyDataBridge->DeleteCell(removeCells->GetId(i));
		}
		polyDataBridge->RemoveDeletedCells();
		polyDataBridge->BuildLinks();

		vtkSmartPointer<vtkCleanPolyData> cleanFilterPolyDataBridge = vtkSmartPointer<vtkCleanPolyData>::New();
		cleanFilterPolyDataBridge->SetInputData(polyDataBridge);
		cleanFilterPolyDataBridge->Update();
		helpBridge.setVtkData(cleanFilterPolyDataBridge->GetOutput());

		if (ConfigFiberorientation::openVTPBridgeInAtrial) {
			set<vtkIdType> removedCells;

			vtkSmartPointer<vtkIdList> removeDataCells = Methods::getPointsWereInOldData(data, pill);
			vtkSmartPointer<vtkPolyData> polyData = vtkPolyData::SafeDownCast(data.getVtkData());
			for (vtkIdType i = 0; i < removeDataCells->GetNumberOfIds(); i++) {
				vtkSmartPointer<vtkIdList> cells = vtkSmartPointer<vtkIdList>::New();
				polyData->GetPointCells(removeDataCells->GetId(i), cells);
				for (vtkIdType j = 0; j < cells->GetNumberOfIds(); j++) {
					if (removedCells.find(cells->GetId(j)) == removedCells.end()) {
						polyData->DeleteCell(cells->GetId(j));
						data.getCentrePoints()->GetPoints()->SetPoint(cells->GetId(j), NAN, NAN, NAN);
						removedCells.insert(cells->GetId(j));
					}
				}
			}

			polyData->RemoveDeletedCells();
			polyData->BuildLinks();

			vtkSmartPointer<vtkCleanPolyData> cleanFilter = vtkSmartPointer<vtkCleanPolyData>::New();
			cleanFilter->SetInputData(polyData);
			cleanFilter->Update();

			data.setVtkData(cleanFilter->GetOutput());
		}

		setPath("bridge" + to_string(bridgeNumber), Methods::unionData(data, helpBridge));

		Methods::vtpBridgeMarker(data, getPath("bridge" + to_string(
			bridgeNumber)), Material::insert_new_Material, getBPoint(
				pointPos), getBPoint(pointPos + 1));
		Methods::setMaterial(data, getPath("bridge" + to_string(
			bridgeNumber)), bridgeMaterial, Material::insert_new_Material);
	}
} // FiberOrientation::setFreeBridge

/*! Function for insertion the bachmann bundle into a mesh.
 \param data Pointer to the orginal mesh
 \param pointLeft left point of the bachmann bundle
 \param pointRight right point of the bachmann bundle
 \param bridgeNumber number of the bridge
 \param bridgewidth radius of the inserted bachmann bundle
 \param searchRadius search radius for the origin point in the orginal atrial mesh
 \param bridgeMaterialLeft tissue class of the left and inserted bridge
 \param bridgeMaterialRight tissue class of the right bridge
 */
void FiberOrientation::setBachmannBridgelong(DataFormat &data, vector<double> pointLeft, vector<double> pointRight, int bridgeNumber, double bridgewidth, double searchRadius, Material::Mat bridgeMaterialLeft, Material::Mat bridgeMaterialRight) {
	vector<double> bridge1TempPoint1 = Methods::findClosedPointinMaterialinRadiustoPointLeft(data, pointLeft, pointRight, searchRadius, getResolution());
	vector<double> bridge1TempPoint2 = Methods::findClosedPointinMaterialinRadiustoPointRight(data, pointRight, pointLeft, searchRadius, getResolution());

	vector<double> bridge1TempPoint3 = { (bridge1TempPoint1.at(0) + bridge1TempPoint2.at(0)) / 2, (bridge1TempPoint1.at(1) + bridge1TempPoint2.at(1)) / 2,  (bridge1TempPoint1.at(2) + bridge1TempPoint2.at(2)) / 2 };
	int pointPos = bridgeNumber * 2;

	DataFormat pillForRegion = Methods::getPill(pointLeft, pointRight, getResolution(), searchRadius);
	vtkSmartPointer<vtkIdList> region = Methods::getCellsWhichWereInside(data, pillForRegion);

	setBPoint(pointPos, Methods::findClosedPointinMaterialInLeftInRegion(data, bridge1TempPoint3, region));
	setBPoint(pointPos + 1, Methods::findClosedPointinMaterialInRightInRegion(data, bridge1TempPoint3, region));

	vector<double> pathPoint = Methods::findClosedPointOnPath(data, getBPoint(pointPos + 1), getPath("bridgeR2R29R7"));


	setPath("bridge" + to_string(bridgeNumber) + ".1", Methods::pathSearch(data, pointLeft, getBPoint(pointPos), "left"));
	setPath("bridge" + to_string(bridgeNumber) + ".3", Methods::pathSearch(data, pathPoint, getBPoint(pointPos + 1), "right"));


	if (ConfigGeneral::debug) {
		cout << "Bachmann wide Point in left atria:" << pointLeft.at(0) << ", " << pointLeft.at(1) << ", " << pointLeft.at(
			2) << endl;
		cout << "Bachmann bridge point in left atria" << getBPoint(pointPos).at(0) << ", " << getBPoint(pointPos).at(1) <<
			", " << getBPoint(pointPos).at(2) << endl;


		cout << "Bachmann wide Point in right atria:" << pathPoint.at(0) << ", " << pathPoint.at(1) << ", " << pathPoint.at(
			2) << endl;
		cout << "Bachmann bridge point in right atria" << getBPoint(pointPos + 1).at(0) << ", " <<
			getBPoint(pointPos + 1).at(1) << ", " << getBPoint(pointPos + 1).at(2) << endl;
	}


	double resolution = getResolution();


	// volume mesh

	if (data.getInputType() == DataFormat::PossibleInputType::vtu) {
		vtkSmartPointer<vtkIdList> cellsInRadiusPointLeft = Methods::getCellsinRadius(data, data.getCentrePoints()->FindPoint(
			getBPoint(pointPos + 1).at(0), getBPoint(pointPos + 1).at(1), getBPoint(
				pointPos + 1).at(2)), bridgewidth);
		vtkSmartPointer<vtkIdList> cellsInRadiusPointRight = Methods::getCellsinRadius(data, data.getCentrePoints()->FindPoint(
			getBPoint(pointPos).at(0), getBPoint(pointPos).at(1), getBPoint(pointPos).at(
				2)), bridgewidth);

		vtkSmartPointer<vtkIdList> Vorhof_linksIds = vtkSmartPointer<vtkIdList>::New();

		for (vtkIdType i = 0; i < cellsInRadiusPointLeft->GetNumberOfIds(); i++) {
			int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(cellsInRadiusPointLeft->GetId(i), 0);
			if (Material::isInLeft(material)) {
				Vorhof_linksIds->InsertNextId(cellsInRadiusPointLeft->GetId(i));
			}
		}
		vtkSmartPointer<vtkIdList> Vorhof_rechtsIds = vtkSmartPointer<vtkIdList>::New();

		for (vtkIdType i = 0; i < cellsInRadiusPointRight->GetNumberOfIds(); i++) {
			int material = data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(cellsInRadiusPointRight->GetId(
				i), 0);
			if (Material::isInRight(material)) {
				Vorhof_rechtsIds->InsertNextId(cellsInRadiusPointLeft->GetId(i));
			}
		}
		unordered_map<vtkIdType, vtkIdType> commonPoints = Methods::getCommonPoints(data, Vorhof_linksIds, Vorhof_rechtsIds);
		cout << "direct conection points: " << commonPoints.size() << endl;
		if (commonPoints.size() > 0) {
			Methods::replacePoints(data, commonPoints);
		}

		// voxel mesh
		if (data.getVtkData()->GetCellType(data.getVtkData()->GetNumberOfCells() / 2) != 10) {
			// wide Bachmann bundle
			if (ConfigFiberorientation::wideBachmann) {
				DataFormat bachmannLong1;
				bachmannLong1 = Methods::createVoxelAroundPath(data, getPath("bridge" + to_string(
					bridgeNumber) + ".1"), bridgewidth);

				vtkSmartPointer<vtkIdList> diffcellIDs1 = Methods::getCellsWereNotinOldData(bachmannLong1, data);
				if (diffcellIDs1->GetNumberOfIds() > 0) {
					vtkSmartPointer<vtkUnstructuredGrid> diffCylinderPath1 = Methods::extractCellsVtu(
						bachmannLong1.getVtkData(), diffcellIDs1);
					bachmannLong1.setVtkData(diffCylinderPath1);
					bachmannLong1.setCentrePoints(Methods::calcCellCentroids(diffCylinderPath1));

					diffcellIDs1 = Methods::getAllCellNeighbours(bachmannLong1, Methods::findClosedCellIdinDirectionPoint(bachmannLong1, getBPoint(pointPos), getRPoint(5), 15));
					if (diffcellIDs1->GetNumberOfIds() > 0) {
						diffCylinderPath1 = Methods::extractCellsVtu(bachmannLong1.getVtkData(), diffcellIDs1);
						bachmannLong1.setVtkData(diffCylinderPath1);
						bachmannLong1.setCentrePoints(Methods::calcCellCentroids(diffCylinderPath1));
						Methods::unionData(data, bachmannLong1);
					}
				}


				DataFormat bachmannLong3;
				bachmannLong3 = Methods::createVoxelAroundPath(data, getPath("bridge" + to_string(
					bridgeNumber) + ".3"), bridgewidth);

				vtkSmartPointer<vtkIdList> diffcellIDs3 = Methods::getCellsWereNotinOldData(bachmannLong3, data);
				if (diffcellIDs3->GetNumberOfIds() > 0) {
					vtkSmartPointer<vtkUnstructuredGrid> diffCylinderPath3 = Methods::extractCellsVtu(
						bachmannLong3.getVtkData(), diffcellIDs3);
					bachmannLong3.setVtkData(diffCylinderPath3);
					bachmannLong3.setCentrePoints(Methods::calcCellCentroids(diffCylinderPath3));

					diffcellIDs3 = Methods::getAllCellNeighbours(bachmannLong3, Methods::findClosedCellIdinDirectionPoint(bachmannLong3, getBPoint(pointPos + 1), getLPointEpi(10), 15));
					if (diffcellIDs3->GetNumberOfIds() > 0) {
						diffCylinderPath3 = Methods::extractCellsVtu(bachmannLong3.getVtkData(), diffcellIDs3);
						bachmannLong3.setVtkData(diffCylinderPath3);
						bachmannLong3.setCentrePoints(Methods::calcCellCentroids(diffCylinderPath3));
						Methods::unionData(data, bachmannLong3);
					}
				}
			}

			DataFormat bachmannLong2;
			bachmannLong2 = Methods::createVoxelTube(data, getBPoint(pointPos), getBPoint(pointPos + 1), resolution, bridgewidth);

			vtkSmartPointer<vtkIdList> diffcellIDs2 = Methods::getCellsWereNotinOldData(bachmannLong2, data);
			if (diffcellIDs2->GetNumberOfIds() > 0) {
				vtkSmartPointer<vtkUnstructuredGrid> diffCylinder = Methods::extractCellsVtu(
					bachmannLong2.getVtkData(), diffcellIDs2);
				bachmannLong2.setVtkData(diffCylinder);
				bachmannLong2.setCentrePoints(Methods::calcCellCentroids(diffCylinder));

				diffcellIDs2 = Methods::getAllCellNeighbours(bachmannLong2, Methods::findClosedCellIdinDirectionPoint(bachmannLong2, getBPoint(pointPos), getRPoint(5), 15));
				if (diffcellIDs2->GetNumberOfIds() > 0) {
					diffCylinder = Methods::extractCellsVtu(bachmannLong2.getVtkData(), diffcellIDs2);
					bachmannLong2.setVtkData(diffCylinder);
					bachmannLong2.setCentrePoints(Methods::calcCellCentroids(diffCylinder));
					Methods::unionData(data, bachmannLong2);
				}
				else {
					double bPointArray[3] = { getBPoint(pointPos).at(0), getBPoint(pointPos).at(1), getBPoint(pointPos).at(2) };
					diffcellIDs2 = Methods::getAllCellNeighbours(bachmannLong2, bachmannLong2.getCentrePoints()->FindPoint(bPointArray));
					if (diffcellIDs2->GetNumberOfIds() > 0) {
						diffCylinder = Methods::extractCellsVtu(bachmannLong2.getVtkData(), diffcellIDs2);
						bachmannLong2.setVtkData(diffCylinder);
						bachmannLong2.setCentrePoints(Methods::calcCellCentroids(diffCylinder));
						Methods::unionData(data, bachmannLong2);
					}
				}
			}
		}

		// tetrahedron mesh and no wide Bachmann bundle
		else if ((data.getVtkData()->GetCellType(data.getVtkData()->GetNumberOfCells() / 2) == 10) &&
			!ConfigFiberorientation::wideBachmann) {
			DataFormat bachmannLong2;
			bachmannLong2 = Methods::getPill(getBPoint(pointPos), getBPoint(pointPos + 1), resolution, bridgewidth);

			vtkSmartPointer<vtkPolyData> surface = Methods::extractSurface(data);
			vtkSmartPointer<vtkPolyData> pill2vtp = Methods::vtpDifference(vtkPolyData::SafeDownCast(bachmannLong2.getVtkData()), surface);


			bachmannLong2.setVtkData(pill2vtp);
			bachmannLong2.setCentrePoints(Methods::calcCellCentroids(pill2vtp));

			vtkSmartPointer<vtkIdList> diffcellIDs = Methods::getAllCellNeighbours(bachmannLong2, bachmannLong2.getCentrePoints()->FindPoint(bridge1TempPoint3.at(0), bridge1TempPoint3.at(1), bridge1TempPoint3.at(2)));
			pill2vtp = Methods::extractCellsVtp(pill2vtp, diffcellIDs);

			vtkSmartPointer<vtkUnstructuredGrid> pill2vtu = vtkSmartPointer<vtkUnstructuredGrid>::New();

			/*bool hasSelfIntersect = vcgTools::VcgTools::hasSelfIntersections(pill2vtp);
			 if (hasSelfIntersect) {
			 cout << "The Mesh has self- intersections" << endl;
			 }*/
			bool hasHoles = Methods::hasHoles(pill2vtp);

			if (hasHoles) {
				cout << "The Mesh has holes" << endl;
				cout << "write stl" << endl;
				vtkSmartPointer<vtkSTLWriter> stlWriter = vtkSmartPointer<vtkSTLWriter>::New();
				stlWriter->SetInputData(pill2vtp);
				stlWriter->SetFileTypeToBinary();
				stlWriter->SetFileName("BachmannBridge.stl");
				stlWriter->Write();

				string fixedMeshPath;
				cout << "Check is the mesh in the file BachmannBridge.stl closed and has no self-intersections" << endl;
				cout << "Please enter the path where the fixed mesh is: " << endl;

				bool exit = false;
				while (!exit) {
					cin >> fixedMeshPath;

					ifstream FileTest(fixedMeshPath.c_str());

					if (FileTest) {
						exit = true;
					}
					else {
						cout << "write stl" << endl;
						vtkSmartPointer<vtkSTLWriter> stlWriter = vtkSmartPointer<vtkSTLWriter>::New();
						stlWriter->SetInputData(pill2vtp);
						stlWriter->SetFileTypeToBinary();
						stlWriter->SetFileName("BachmannBridge.stl");
						stlWriter->Write();

						cout << "try again:" << endl;
					}

					FileTest.close();
				}

				cout << "read stl" << endl;


				vtkSmartPointer<vtkSTLReader> stlReader = vtkSmartPointer<vtkSTLReader>::New();
				stlReader->SetFileName(fixedMeshPath.c_str());
				stlReader->Update();

				pill2vtp = stlReader->GetOutput();
				pill2vtp = Methods::closeSurfaceFilter(pill2vtp);

				bachmannLong2.setVtkData(pill2vtp);
				bachmannLong2.setCentrePoints(Methods::calcCellCentroids(pill2vtp));

				diffcellIDs = Methods::getAllCellNeighbours(bachmannLong2, bachmannLong2.getCentrePoints()->FindPoint(bridge1TempPoint3.at(0), bridge1TempPoint3.at(1), bridge1TempPoint3.at(2)));
				pill2vtp = Methods::extractCellsVtp(pill2vtp, diffcellIDs);
				pill2vtu = Methods::tetrahedralizeVTPSurface(pill2vtp);
			}

			bachmannLong2.setVtkData(pill2vtu);
			bachmannLong2.setCentrePoints(Methods::calcCellCentroids(pill2vtu));
			bachmannLong2.setInputType(DataFormat::vtu);

			Methods::unionData(data, bachmannLong2);
		}

		// tetrahedron mesh and wide Bachmann bundle
		else if ((data.getVtkData()->GetCellType(data.getVtkData()->GetNumberOfCells() / 2) == 10) &&
			ConfigFiberorientation::wideBachmann) {
			cout << "Calc Bachmann Bundel Part1" << endl;
			DataFormat bachmannLong1;
			bachmannLong1 = Methods::getPillAroundPath(data, getPath("bridge" + to_string(
				bridgeNumber) + ".1"), resolution, bridgewidth);
			if (ConfigGeneral::debug) {
				Methods::writeIntermediateData(bachmannLong1, "BundelPart1");
			}

			vtkSmartPointer<vtkPolyData> surface = data.getSurfaceForBridges();

			cout << "Calc Bachmann Bundel Part2" << endl;
			DataFormat bachmannLong2;
			bachmannLong2 = Methods::getPill(getBPoint(pointPos), getBPoint(pointPos + 1), resolution, bridgewidth);

			if (ConfigGeneral::debug) {
				Methods::writeIntermediateData(bachmannLong2, "BundelPart2");
			}

			cout << "Calc Bachmann Bundel Part3" << endl;
			DataFormat bachmannLong3;
			bachmannLong3 = Methods::getPillAroundPath(data, getPath("bridge" + to_string(
				bridgeNumber) + ".3"), resolution, bridgewidth);


			if (ConfigGeneral::debug) {
				Methods::writeIntermediateData(bachmannLong3, "BundelPart3");
			}

			vtkSmartPointer<vtkPolyData> allBachmann = Methods::vtpUnion(vtkPolyData::SafeDownCast(
				bachmannLong1.getVtkData()), vtkPolyData::SafeDownCast(bachmannLong2.getVtkData()));

			allBachmann = Methods::vtpUnion(allBachmann, vtkPolyData::SafeDownCast(bachmannLong3.getVtkData()));

			DataFormat Bachmann_union;
			Bachmann_union.setVtkData(allBachmann);
			Bachmann_union.setInputType(DataFormat::vtp);
			Bachmann_union.setCentrePoints(Methods::calcCellCentroids(allBachmann));

			if (ConfigGeneral::debug) {
				Methods::writeIntermediateData(Bachmann_union, "Bachmann_union");
			}


			allBachmann = Methods::vtpDifference(allBachmann, surface);

			DataFormat bachmannAll;
			bachmannAll.setVtkData(allBachmann);
			bachmannAll.setInputType(DataFormat::vtp);
			bachmannAll.setCentrePoints(Methods::calcCellCentroids(allBachmann));

			if (ConfigGeneral::debug) {
				Methods::writeIntermediateData(bachmannAll, "bachmannAll_differenz");
			}

			vtkSmartPointer<vtkIdList> diffcellIDs = Methods::getAllCellNeighbours(bachmannAll, bachmannAll.getCentrePoints()->FindPoint(bridge1TempPoint3.at(0), bridge1TempPoint3.at(1), bridge1TempPoint3.at(2)));

			allBachmann = Methods::extractCellsVtp(allBachmann, diffcellIDs);

			bachmannAll.setVtkData(allBachmann);
			bachmannAll.setInputType(DataFormat::vtp);
			bachmannAll.setCentrePoints(Methods::calcCellCentroids(allBachmann));

			// remeshing Bachmann

			DataFormat bachmannAll_with_elemTag;
			bachmannAll_with_elemTag.setVtkData(allBachmann);
			bachmannAll_with_elemTag.setInputType(DataFormat::vtk);
			bachmannAll_with_elemTag.setCentrePoints(Methods::calcCellCentroids(allBachmann));

			vtkSmartPointer<vtkIntArray> elemTag = vtkSmartPointer<vtkIntArray>::New();
			elemTag->SetName("elemTag");
			elemTag->SetNumberOfComponents(1);
			elemTag->SetNumberOfTuples(bachmannAll_with_elemTag.getCentrePoints()->GetNumberOfPoints());
			elemTag->FillComponent(0, 0);
			bachmannAll_with_elemTag.getVtkData()->GetCellData()->AddArray(elemTag);

			DataFormat surface_data;
			surface_data.setVtkData(data.getSurfacewithNormals());
			surface_data.setCentrePoints(data.getCentrePointsSurface());

			vtkSmartPointer<vtkIdList> cellPointsOnSurface = Methods::getCellsWhichWereInside(surface_data, Bachmann_union);

			vtkSmartPointer<vtkIdList> cellPointWhonotChange = vtkSmartPointer<vtkIdList>::New();
			for (vtkIdType i = 0; i < cellPointsOnSurface->GetNumberOfIds(); i++) {
				vtkIdType cellIdInBridge = bachmannAll_with_elemTag.getCentrePoints()->FindPoint(
					surface_data.getCentrePoints()->GetPoint(cellPointsOnSurface->GetId(i)));

				if ((surface_data.getCentrePoints()->GetPoint(cellPointsOnSurface->GetId(i))[0] == bachmannAll_with_elemTag.getCentrePoints()->GetPoint(cellIdInBridge)[0]) &&
					(surface_data.getCentrePoints()->GetPoint(cellPointsOnSurface->GetId(i))[1] == bachmannAll_with_elemTag.getCentrePoints()->GetPoint(cellIdInBridge)[1]) &&
					(surface_data.getCentrePoints()->GetPoint(cellPointsOnSurface->GetId(i))[2] == bachmannAll_with_elemTag.getCentrePoints()->GetPoint(cellIdInBridge)[2])) {
					vtkIdType cellIdInMesh = data.getCentrePoints()->FindPoint(surface_data.getCentrePoints()->GetPoint(cellPointsOnSurface->GetId(i)));

					int material = data.getVtkData()->GetCellData()->GetArray("RawMaterial")->GetComponent(cellIdInMesh, 0);
					if (Material::isInEpi(material)) {
						cellPointWhonotChange->InsertNextId(cellIdInBridge);
						bachmannAll_with_elemTag.getVtkData()->GetCellData()->GetArray("elemTag")->SetComponent(cellIdInBridge, 0, 600);
					}
				}
			}
			if (ConfigGeneral::debug) {
				Methods::writeIntermediateData(bachmannAll_with_elemTag, "bachmannAll_with_elemTag");
			}

			// string cmdv;
			// string out;
			// cmdv = "meshtool resample surfmesh -msh=bachmannAll_with_elemTag -ifmt=vtk -tags=0 -outmsh=bachmannAll_resample -ofmt=vtk -min=" + to_string(resolution)+" -max=" + to_string(resolution*3);
			// int n = cmdv.length();
			// char char_array[500];
			// char char_array[]=cmdv;
			// strcpy(char_array, cmdv.c_str());

			// out = Methods::exec(char_array);
			// stringstream stream(cmdv);
			// string oneWord;
			// unsigned int count = 0;

			// while(stream >> oneWord) { ++count;}

			// resample_mode(count, cmdv);

			// }
			// end remeshing


			/*bool hasSelfIntersect = vcgTools::VcgTools::hasSelfIntersections(allBachmann);
			 if (hasSelfIntersect) {
			 cout << "The Mesh has self- intersections" << endl;
			 }*/
			bool hasHoles = Methods::hasHoles(allBachmann);
			if (hasHoles) {
				cout << "The Mesh has holes" << endl;
				cout << "write stl" << endl;
				vtkSmartPointer<vtkSTLWriter> stlWriter = vtkSmartPointer<vtkSTLWriter>::New();
				stlWriter->SetInputData(allBachmann);
				stlWriter->SetFileTypeToBinary();
				stlWriter->SetFileName("BachmannBridge.stl");
				stlWriter->Write();

				string fixedMeshPath;
				cout << "Check is the mesh in the file BachmannBridge.stl closed and has no self-intersections" << endl;
				cout << "Please enter the path where the fixed mesh is: " << endl;

				bool exit = false;
				while (!exit) {
					cin >> fixedMeshPath;

					ifstream FileTest(fixedMeshPath.c_str());

					if (FileTest) {
						exit = true;
					}
					else {
						cout << "write stl" << endl;
						vtkSmartPointer<vtkSTLWriter> stlWriter = vtkSmartPointer<vtkSTLWriter>::New();
						stlWriter->SetInputData(allBachmann);
						stlWriter->SetFileTypeToBinary();
						stlWriter->SetFileName("BachmannBridge.stl");
						stlWriter->Write();

						cout << "try again:" << endl;
					}

					FileTest.close();
				}

				cout << "read fixed BB" << endl;


				// fixedMeshPath = "bachmannAll_resample";
				// UniversalReader reader;
				// bachmannAll = reader.read(fixedMeshPath+".stl");

				vtkSmartPointer<vtkSTLReader> stlReader = vtkSmartPointer<vtkSTLReader>::New();
				stlReader->SetFileName(fixedMeshPath.c_str());
				stlReader->Update();

				vtkSmartPointer<vtkPolyData> bachmannVtp = stlReader->GetOutput();

				DataFormat bachmannAll;
				bachmannAll.setVtkData(bachmannVtp);
				bachmannAll.setInputType(DataFormat::vtp);
				bachmannAll.setCentrePoints(Methods::calcCellCentroids(allBachmann));
			}
			vtkSmartPointer<vtkUnstructuredGrid> allBachmannVtu = Methods::tetrahedralizeVTPSurface(Methods::extractSurface(bachmannAll));


			DataFormat allBachmannData;
			allBachmannData.setVtkData(allBachmannVtu);
			allBachmannData.setCentrePoints(Methods::calcCellCentroids(allBachmannVtu));
			allBachmannData.setInputType(DataFormat::vtu);

			if (ConfigGeneral::debug) {
				Methods::writeIntermediateData(allBachmannData, "wideBachmannBridge");
			}

			Methods::unionData(data, allBachmannData);
		}
		else {
			cerr << "Bridge: celltype undefined" << endl;
		}

		setPath("bridge" + to_string(bridgeNumber) + ".2", Methods::pathSearch(data, getBPoint(pointPos), getBPoint(pointPos + 1), "right/left"));

		setPath("bridge" + to_string(bridgeNumber), Methods::add3Paths(getPath("bridge" + to_string(bridgeNumber) + ".1"), getPath("bridge" + to_string(bridgeNumber) + ".2"), getPath("bridge" + to_string(bridgeNumber) + ".3")));

		Methods::pathMarker(data, getPath("bridge" + to_string(bridgeNumber)));

		setPath("growbridge" + to_string(bridgeNumber), Methods::growPath(data, getPath("bridge" + to_string(bridgeNumber)), bridgewidth));
		Methods::smoothingPath(data, getPath("bridge" + to_string(bridgeNumber)));

		Methods::setMaterialinLeft(data, getPath("growbridge" + to_string(bridgeNumber)), bridgeMaterialLeft);
		Methods::setMaterialinRight(data, getPath("growbridge" + to_string(bridgeNumber)), bridgeMaterialRight);
		Methods::setMaterial(data, getPath("growbridge" + to_string(
			bridgeNumber)), bridgeMaterialLeft, Material::insert_new_Material); // Matreialklassen
			// anpassen
	}

	// surface mesh
	else if (data.getInputType() == DataFormat::PossibleInputType::vtp) {
		Methods::pathMarker(data, getPath("bridge" + to_string(bridgeNumber) + ".1"));
		setPath("grownbridge" + to_string(bridgeNumber) + ".1", Methods::growPathInRight(data, getPath("bridge" + to_string(bridgeNumber) + ".1"), bridgewidth));
		Methods::smoothingPath(data, getPath("bridge" + to_string(bridgeNumber) + ".1"));
		Methods::setMaterial(data, getPath("grownbridge" + to_string(bridgeNumber) + ".1"), bridgeMaterialLeft);

		if (ConfigFiberorientation::wideBachmann) {
			DataFormat pill1;
			pill1 = Methods::getPillAroundPath(data, getPath("bridge" + to_string(
				bridgeNumber) + ".1"), resolution, bridgewidth);

			DataFormat bachmannLong1;
			bachmannLong1.setInputType(DataFormat::vtp);
			bachmannLong1.setVtkData(Methods::vtpUnion(vtkPolyData::SafeDownCast(data.getVtkData()), vtkPolyData::SafeDownCast(pill1.getVtkData())));
			bachmannLong1.setCentrePoints(Methods::calcCellCentroids(bachmannLong1.getVtkData()));

			vtkSmartPointer<vtkIdList> removeCells1 = Methods::getCellsWereInOldData(bachmannLong1, data);
			vtkSmartPointer<vtkPolyData> polyDataBridge1 = vtkPolyData::SafeDownCast(bachmannLong1.getVtkData());
			for (vtkIdType i = 0; i < removeCells1->GetNumberOfIds(); i++) {
				polyDataBridge1->DeleteCell(removeCells1->GetId(i));
			}
			polyDataBridge1->RemoveDeletedCells();
			polyDataBridge1->BuildLinks();

			vtkSmartPointer<vtkCleanPolyData> cleanFilterPolyDataBridge1 = vtkSmartPointer<vtkCleanPolyData>::New();
			cleanFilterPolyDataBridge1->SetInputData(polyDataBridge1);
			cleanFilterPolyDataBridge1->Update();
			bachmannLong1.setVtkData(cleanFilterPolyDataBridge1->GetOutput());

			set<vtkIdType> removedCells1;

			vtkSmartPointer<vtkIdList> removeDataCells1 = Methods::getPointsWereInOldData(data, pill1);
			vtkSmartPointer<vtkPolyData> polyData1 = vtkPolyData::SafeDownCast(data.getVtkData());
			vtkSmartPointer<vtkPolyData> forOldDatpolyData = vtkSmartPointer<vtkPolyData>::New();
			forOldDatpolyData->DeepCopy(polyData1);
			vtkSmartPointer<vtkUnstructuredGrid> forOldCentrePoints = vtkSmartPointer<vtkUnstructuredGrid>::New();
			forOldCentrePoints->DeepCopy(data.getCentrePoints());

			DataFormat oldData1;
			oldData1.setVtkData(forOldDatpolyData);
			oldData1.setCentrePoints(forOldCentrePoints);

			for (vtkIdType i = 0; i < removeDataCells1->GetNumberOfIds(); i++) {
				vtkSmartPointer<vtkIdList> cells = vtkSmartPointer<vtkIdList>::New();
				polyData1->GetPointCells(removeDataCells1->GetId(i), cells);
				for (vtkIdType j = 0; j < cells->GetNumberOfIds(); j++) {
					if (removedCells1.find(cells->GetId(j)) == removedCells1.end()) {
						polyData1->DeleteCell(cells->GetId(j));
						data.getCentrePoints()->GetPoints()->SetPoint(cells->GetId(j), NAN, NAN, NAN);
						removedCells1.insert(cells->GetId(j));
					}
				}
			}

			polyData1->RemoveDeletedCells();
			polyData1->BuildLinks();

			vtkSmartPointer<vtkCleanPolyData> cleanFilter1 = vtkSmartPointer<vtkCleanPolyData>::New();
			cleanFilter1->SetInputData(polyData1);
			cleanFilter1->Update();

			data.setVtkData(cleanFilter1->GetOutput());

			setPath("grownbridge" + to_string(bridgeNumber) + ".1", Methods::unionData(data, bachmannLong1));

			data.setCentrePoints(Methods::calcCellCentroids(data.getVtkData()));

			Methods::vtpBridgeMarkerArroundPath(data, oldData1, getPath("grownbridge" + to_string(
				bridgeNumber) + ".1"), getPath("bridge" + to_string(bridgeNumber) + ".1"));

			Methods::setMaterial(data, getPath("grownbridge" + to_string(bridgeNumber) + ".1"), bridgeMaterialLeft);

			data = update(data);
        
            if (ConfigGeneral::debug) {
                Methods::writeIntermediateData(data, "BachmannBundle_1_debug");
            }
        }

		setPath("bridge" + to_string(bridgeNumber) + ".3", Methods::pathSearch(data, getBPoint(pointPos + 1), pathPoint, "right"));

		Methods::pathMarker(data, getPath("bridge" + to_string(bridgeNumber) + ".3"));
		setPath("grownbridge" + to_string(bridgeNumber) + ".3", Methods::growPath(data, getPath("bridge" + to_string(bridgeNumber) + ".3"), bridgewidth));
		Methods::smoothingPath(data, getPath("bridge" + to_string(bridgeNumber) + ".3"));
		Methods::setMaterial(data, getPath("grownbridge" + to_string(bridgeNumber) + ".3"), bridgeMaterialRight);


		// wide Bachmann bundle
		if (ConfigFiberorientation::wideBachmann) {
			DataFormat pill3;
			pill3 = Methods::getPillAroundPath(data, getPath("bridge" + to_string(
				bridgeNumber) + ".3"), resolution, bridgewidth);

			DataFormat bachmannLong3;
			bachmannLong3.setInputType(DataFormat::vtp);
			bachmannLong3.setVtkData(Methods::vtpUnion(vtkPolyData::SafeDownCast(data.getVtkData()), vtkPolyData::SafeDownCast(pill3.getVtkData())));
			bachmannLong3.setCentrePoints(Methods::calcCellCentroids(bachmannLong3.getVtkData()));

			vtkSmartPointer<vtkIdList> removeCells3 = Methods::getCellsWereInOldData(bachmannLong3, data);
			vtkSmartPointer<vtkPolyData> polyDataBridge3 = vtkPolyData::SafeDownCast(bachmannLong3.getVtkData());
			for (vtkIdType i = 0; i < removeCells3->GetNumberOfIds(); i++) {
				polyDataBridge3->DeleteCell(removeCells3->GetId(i));
			}
			polyDataBridge3->RemoveDeletedCells();
			polyDataBridge3->BuildLinks();

			vtkSmartPointer<vtkCleanPolyData> cleanFilterPolyDataBridge3 = vtkSmartPointer<vtkCleanPolyData>::New();
			cleanFilterPolyDataBridge3->SetInputData(polyDataBridge3);
			cleanFilterPolyDataBridge3->Update();
			bachmannLong3.setVtkData(cleanFilterPolyDataBridge3->GetOutput());

			set<vtkIdType> removedCells3;

			vtkSmartPointer<vtkIdList> removeDataCells3 = Methods::getPointsWereInOldData(data, pill3);
			vtkSmartPointer<vtkPolyData> polyData3 = vtkPolyData::SafeDownCast(data.getVtkData());

			vtkSmartPointer<vtkPolyData> forOldDatpolyData3 = vtkSmartPointer<vtkPolyData>::New();
			forOldDatpolyData3->DeepCopy(polyData3);
			vtkSmartPointer<vtkUnstructuredGrid> forOldCentrePoints3 = vtkSmartPointer<vtkUnstructuredGrid>::New();
			forOldCentrePoints3->DeepCopy(data.getCentrePoints());

			DataFormat oldData3;
			oldData3.setVtkData(forOldDatpolyData3);
			oldData3.setCentrePoints(forOldCentrePoints3);

			for (vtkIdType i = 0; i < removeDataCells3->GetNumberOfIds(); i++) {
				vtkSmartPointer<vtkIdList> cells3 = vtkSmartPointer<vtkIdList>::New();
				polyData3->GetPointCells(removeDataCells3->GetId(i), cells3);
				for (vtkIdType j = 0; j < cells3->GetNumberOfIds(); j++) {
					if (removedCells3.find(cells3->GetId(j)) == removedCells3.end()) {
						polyData3->DeleteCell(cells3->GetId(j));
						data.getCentrePoints()->GetPoints()->SetPoint(cells3->GetId(j), NAN, NAN, NAN);
						removedCells3.insert(cells3->GetId(j));
					}
				}
			}

			polyData3->RemoveDeletedCells();
			polyData3->BuildLinks();

			vtkSmartPointer<vtkCleanPolyData> cleanFilter3 = vtkSmartPointer<vtkCleanPolyData>::New();
			cleanFilter3->SetInputData(polyData3);
			cleanFilter3->Update();

			data.setVtkData(cleanFilter3->GetOutput());

			setPath("grownbridge" + to_string(bridgeNumber) + ".3", Methods::unionData(data, bachmannLong3));
			data.setCentrePoints(Methods::calcCellCentroids(data.getVtkData()));

			Methods::vtpBridgeMarkerArroundPath(data, oldData3, getPath("grownbridge" + to_string(
				bridgeNumber) + ".3"), getPath("bridge" + to_string(bridgeNumber) + ".3"));

			Methods::setMaterial(data, getPath("grownbridge" + to_string(bridgeNumber) + ".3"), bridgeMaterialRight);

			data = update(data);
		
            if (ConfigGeneral::debug) {
                Methods::writeIntermediateData(data, "BachmannBundle_3_debug");
            }
        }
        
		DataFormat pill2;
		pill2 = Methods::getPill(getBPoint(pointPos), getBPoint(pointPos + 1), resolution, bridgewidth);

		DataFormat bachmannLong2;
		bachmannLong2.setInputType(DataFormat::vtp);
		bachmannLong2.setVtkData(Methods::vtpUnion(vtkPolyData::SafeDownCast(data.getVtkData()), vtkPolyData::SafeDownCast(pill2.getVtkData())));
		bachmannLong2.setCentrePoints(Methods::calcCellCentroids(bachmannLong2.getVtkData()));

        if (ConfigGeneral::debug) {
            Methods::writeIntermediateData(bachmannLong2, "BachmannBundle_2_debug");
        }
        
		vtkSmartPointer<vtkIdList> removeCells2 = Methods::getCellsWereInOldData(bachmannLong2, data);
		vtkSmartPointer<vtkPolyData> polyDataBridge2 = vtkPolyData::SafeDownCast(bachmannLong2.getVtkData());
		for (vtkIdType i = 0; i < removeCells2->GetNumberOfIds(); i++) {
			polyDataBridge2->DeleteCell(removeCells2->GetId(i));
		}
		polyDataBridge2->RemoveDeletedCells();
		polyDataBridge2->BuildLinks();

		vtkSmartPointer<vtkCleanPolyData> cleanFilterPolyDataBridge2 = vtkSmartPointer<vtkCleanPolyData>::New();
		cleanFilterPolyDataBridge2->SetInputData(polyDataBridge2);
		cleanFilterPolyDataBridge2->Update();
		bachmannLong2.setVtkData(cleanFilterPolyDataBridge2->GetOutput());
		bachmannLong2.setCentrePoints(Methods::calcCellCentroids(bachmannLong2.getVtkData()));

        vtkSmartPointer<vtkIdList> diffcellIDs = Methods::getAllCellNeighbours(bachmannLong2, bachmannLong2.getCentrePoints()->FindPoint(bridge1TempPoint3.at(0), bridge1TempPoint3.at(1), bridge1TempPoint3.at(2)));

        bachmannLong2.setVtkData(Methods::extractCellsVtp(cleanFilterPolyDataBridge2->GetOutput(), diffcellIDs));
        bachmannLong2.setCentrePoints(Methods::calcCellCentroids(bachmannLong2.getVtkData()));

       if (ConfigGeneral::debug) {
          Methods::writeIntermediateData(bachmannLong2, "BachmannBundle_2.1_debug");
       }
		set<vtkIdType> removedCells2;

		vtkSmartPointer<vtkIdList> removeDataCells2 = Methods::getPointsWereInOldData(data, pill2);
		vtkSmartPointer<vtkPolyData> polyData2 = vtkPolyData::SafeDownCast(data.getVtkData());
		for (vtkIdType i = 0; i < removeDataCells2->GetNumberOfIds(); i++) {
			vtkSmartPointer<vtkIdList> cells = vtkSmartPointer<vtkIdList>::New();
			polyData2->GetPointCells(removeDataCells2->GetId(i), cells);
			for (vtkIdType j = 0; j < cells->GetNumberOfIds(); j++) {
				if (removedCells2.find(cells->GetId(j)) == removedCells2.end()) {
					polyData2->DeleteCell(cells->GetId(j));
					data.getCentrePoints()->GetPoints()->SetPoint(cells->GetId(j), NAN, NAN, NAN);
					removedCells2.insert(cells->GetId(j));
				}
			}
		}

		polyData2->RemoveDeletedCells();
		polyData2->BuildLinks();

		vtkSmartPointer<vtkCleanPolyData> cleanFilter2 = vtkSmartPointer<vtkCleanPolyData>::New();
		cleanFilter2->SetInputData(polyData2);
		cleanFilter2->Update();

		data.setVtkData(cleanFilter2->GetOutput());
        
        if (ConfigGeneral::debug) {
            Methods::writeIntermediateData(data, "BachmannBundle_2.2_debug");
        }
        
		setPath("grownbridge" + to_string(bridgeNumber) + ".2", Methods::unionData(data, bachmannLong2));
		Methods::vtpBridgeMarker(data, getPath("grownbridge" + to_string(
			bridgeNumber) + ".2"), Material::insert_new_Material, getBPoint(
				pointPos), getBPoint(pointPos + 1));
		Methods::setMaterial(data, getPath("grownbridge" + to_string(
			bridgeNumber) + ".2"), bridgeMaterialLeft, Material::insert_new_Material);
	}
} // FiberOrientation::setBachmannBridgelong

/*! The definition part of sinus node in the right atrial.
 \param data Pointer to the orginal mesh
 */
void FiberOrientation::setSinusNode(DataFormat &data) {
	if (Path.find("pathR4R5") == Path.end()) {
		cout << "calculate pathR4R5" << endl;
		setPath("pathR4R5", Methods::pathSearch(data, getRPoint(3), getRPoint(4), "right"));

		setRPoint(13, Methods::pointAtPercentOfPath(data, getPath("pathR4R5"), 60));
	}

	if (Path.find("pathR3R1") == Path.end()) {
		cout << "calculate pathR3R1" << endl;
		setPath("pathR3R1", Methods::pathSearch(data, getRPoint(2), getRPoint(0), "right"));

		setRPoint(11, Methods::pointAtPercentOfPath(data, getPath("pathR3R1"), 3));
	}

	if (Path.find("pathR12R4") == Path.end()) {
		cout << "calculate pathR12R4" << endl;
		setPath("pathR12R4", Methods::pathSearchOverPlane(data, getRPoint(11), getRPoint(3), getRPoint(13)));

		setRPoint(19, Methods::pointAtPercentOfPath(data, getPath("pathR12R4"), 20));
	}

	vtkSmartPointer<vtkParametricFunctionSource> parametricFunctionSource = vtkSmartPointer<vtkParametricFunctionSource>::New();

	parametricFunctionSource->SetUResolution(64);
	parametricFunctionSource->SetVResolution(64);
	parametricFunctionSource->SetWResolution(64);
	parametricFunctionSource->SetScalarModeToNone();
	parametricFunctionSource->GenerateTextureCoordinatesOff();

	double resolution = getResolution();

	vtkSmartPointer<vtkParametricEllipsoid> ellipsoid = vtkSmartPointer<vtkParametricEllipsoid>::New();
	ellipsoid->SetXRadius(2 + resolution);
	ellipsoid->SetYRadius(2 + resolution);
	ellipsoid->SetZRadius(4 + resolution);

	parametricFunctionSource->SetParametricFunction(ellipsoid);
	parametricFunctionSource->Update();


	vtkSmartPointer<vtkPolyData> polyellipsoid = parametricFunctionSource->GetOutput();

	vector<double> vectorR3R20 = { getRPoint(19).at(0) - getRPoint(2).at(0), getRPoint(19).at(1) - getRPoint(2).at(1), getRPoint(19).at(2) - getRPoint(
	2).at(2) };
	vector<double> vectorEllipsoid = { 0, 0, 1 };
	vtkSmartPointer<vtkMatrix3x3> rotationsMatrix = Methods::getRotaionTranslationMatrix(vectorR3R20, vectorEllipsoid);


	double *rPoint3 = data.getCentrePoints()->GetPoint(data.getCentrePoints()->FindPoint(getRPoint(2).at(0), getRPoint(2).at(1), getRPoint(2).at(2)));

	vector<double> translation = { rPoint3[0], rPoint3[1], rPoint3[2] };

	polyellipsoid->SetPoints(Methods::rotatePoints(polyellipsoid->GetPoints(), rotationsMatrix));
	polyellipsoid->BuildLinks();

	polyellipsoid->SetPoints(Methods::movePoints(polyellipsoid->GetPoints(), translation));
	polyellipsoid->BuildLinks();

	vtkSmartPointer<vtkSelectEnclosedPoints> selectEnclosedPoints = vtkSmartPointer<vtkSelectEnclosedPoints>::New();
	selectEnclosedPoints->Initialize(polyellipsoid);
	selectEnclosedPoints->SetTolerance(0.0001);

	vtkSmartPointer<vtkIdList> cells = Methods::getCellsinRadius(data, data.getCentrePoints()->FindPoint(getRPoint(2).at(0), getRPoint(2).at(1), getRPoint(2).at(2)), ellipsoid->GetZRadius());


	for (vtkIdType i = 0; i < cells->GetNumberOfIds(); i++) {
		double *point = data.getCentrePoints()->GetPoint(cells->GetId(i));
		int inside = selectEnclosedPoints->IsInsideSurface(point);
		if (inside == 1) {
			data.getVtkData()->GetCellData()->GetArray("Material")->SetComponent(cells->GetId(
				i), 0, (double)Material::Sinus_Knoten);
			if (!data.getDoubleLayer()) {
				data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->SetComponent(cells->GetId(
					i), 0, (double)Material::Sinus_Knoten);
			}
		}
	}
	cout << "finished SinusNode" << endl;
} // FiberOrientation::setSinusNode

/*! The definition part to close the fiber orientation of the atrials.
 \param data Pointer to the orginal mesh
 */
void FiberOrientation::closeOrientation(DataFormat &data) {
	for (vtkIdType i = 0; i < data.getVtkData()->GetNumberOfCells(); i++) {
		if ((data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(i, 0) == 0) &&
			(Material::isInLeft(data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(i, 0)) ||
				Material::isInRight(data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(i, 0)))) {
			double *point = data.getCentrePoints()->GetPoint(i);
			vector<double> seedPoint = { point[0], point[1], point[2] };
			if (!data.getDoubleLayer()) {
				vtkSmartPointer<vtkIdList> cells = Methods::regionGrowOneSurface(data, seedPoint, 2.2);
			}
			else {
				vtkSmartPointer<vtkIdList> cells = Methods::regionGrow(data, seedPoint);
			}
		}
	}
}

/*! To close the fiber orientation of the free definiton interatrial bridges.
 \param data Pointer to the orginal mesh
 \param freeBridgeMaterial material of the free bridge

 */
void FiberOrientation::closeOrientationFreeBridge(DataFormat &data, int freeBridgeMaterial) {
	for (vtkIdType i = 0; i < data.getVtkData()->GetNumberOfCells(); i++) {
		if ((data.getVtkData()->GetCellData()->GetArray("Status")->GetComponent(i, 0) == 0) &&
			(Material::isInLeft(data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(i, 0)) ||
				Material::isInRight(data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(i, 0)))) {
			double *point = data.getCentrePoints()->GetPoint(i);
			vector<double> seedPoint = { point[0], point[1], point[2] };
			if (!data.getDoubleLayer()) {
				vtkSmartPointer<vtkIdList> cells = Methods::regionGrowOneSurfaceFreeBridge(data, seedPoint, 2.2, freeBridgeMaterial);
			}
			else {
				vtkSmartPointer<vtkIdList> cells = Methods::regionGrowFreeBridge(data, seedPoint, freeBridgeMaterial);
			}
		}
	}
}

/*! To turn the fiber orientation inmto the atrial wall.
 \param data Pointer to the orginal mesh
 \param data Pointer to the orginal mesh with not seperated rigth an left atrials
 */
void FiberOrientation::smoothOrientationToSurface(DataFormat &data, DataFormat &smoothData) {
	data.setCentrePoints(Methods::calcCellCentroids(data.getVtkData()));

	if (data.getInputType() == DataFormat::vtu) {
		vtkSmartPointer<vtkIdList> newCellIds = Methods::getCellsWereNotinOldData(data, smoothData);
		vtkSmartPointer<vtkUnstructuredGrid> newCellsVtu = Methods::extractCellsVtu(data.getVtkData(), newCellIds);

		DataFormat newCells;
		newCells.setVtkData(newCellsVtu);
		newCells.setCentrePoints(Methods::calcCellCentroids(newCellsVtu));

		Methods::unionData(smoothData, newCells);

		vtkSmartPointer<vtkPolyData> polydata = Methods::extractSurface(smoothData);
		polydata = Methods::triangulateSurface(polydata);
		polydata = Methods::smoothSurface(polydata, 1000);

		DataFormat surface;
		surface.setVtkData(polydata);
		surface.setCentrePoints(Methods::calcCellCentroids(polydata));

		vtkSmartPointer<vtkPolyData> polydataWithNormals = Methods::getSurfaceCellNormals(polydata);

		for (vtkIdType i = 0; i < data.getVtkData()->GetNumberOfCells(); i++) {
			double *currentFiberorientation = data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetTuple(i);

			if ((currentFiberorientation[0] + currentFiberorientation[1] + currentFiberorientation[2]) != 0) {
				vtkIdType surfacePointId = surface.getCentrePoints()->FindPoint(data.getCentrePoints()->GetPoint(i));
				double *currentCellNormal = polydataWithNormals->GetCellData()->GetArray("Normals")->GetTuple(surfacePointId);

				vtkMath::Normalize(currentFiberorientation);
				vtkMath::Normalize(currentCellNormal);

				double newOrientation[3] = { 0, 0, 0 };
				vtkPlane::ProjectVector(currentFiberorientation, data.getCentrePoints()->GetPoint(
					i), currentCellNormal, newOrientation);

				data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->SetTuple(i, newOrientation);
			}
		}
	}
	else if (data.getInputType() == DataFormat::vtp) {
		vtkSmartPointer<vtkPolyData> polydataWithNormals = Methods::getSurfaceCellNormals(vtkPolyData::SafeDownCast(data.getVtkData()));

		for (vtkIdType i = 0; i < data.getVtkData()->GetNumberOfCells(); i++) {
			double *currentFiberorientation = data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->GetTuple(i);
			if ((currentFiberorientation[0] + currentFiberorientation[1] + currentFiberorientation[2]) != 0) {
				double *currentCellNormal = polydataWithNormals->GetCellData()->GetArray("Normals")->GetTuple(i);

				vtkMath::Normalize(currentFiberorientation);
				vtkMath::Normalize(currentCellNormal);

				double newOrientation[3] = { 0, 0, 0 };
				vtkPlane::ProjectVector(currentFiberorientation, data.getCentrePoints()->GetPoint(
					i), currentCellNormal, newOrientation);
				vtkMath::Normalize(newOrientation);

				data.getVtkData()->GetCellData()->GetArray("DifferenceVector")->SetTuple(i, newOrientation);
			}
			if (!data.getDoubleLayer()) {
				double *currentFiberorientation = data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->GetTuple(i);
				if ((currentFiberorientation[0] + currentFiberorientation[1] + currentFiberorientation[2]) != 0) {
					double *currentCellNormal = polydataWithNormals->GetCellData()->GetArray("Normals")->GetTuple(i);


					vtkMath::Normalize(currentFiberorientation);
					vtkMath::Normalize(currentCellNormal);

					double newOrientation[3] = { 0, 0, 0 };

					vtkPlane::ProjectVector(currentFiberorientation, data.getCentrePoints()->GetPoint(
						i), currentCellNormal, newOrientation);
					vtkMath::Normalize(newOrientation);

					data.getVtkData()->GetCellData()->GetArray("DifferenceVectorEndo")->SetTuple(i, newOrientation);
				}
			}
		}
	}
} // FiberOrientation::smoothOrientationToSurface

/*! Substituted some Materials in the mesh for the Simulation.
 \param data Pointer to the orginal mesh
 */
void FiberOrientation::setMaterialtoOldMaterial(DataFormat &data) {
	for (vtkIdType i = 0; i < data.getVtkData()->GetNumberOfCells(); i++) {
		// in the right atrial
		switch ((int)data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(i, 0)) {
		case 159:
			data.getVtkData()->GetCellData()->GetArray("Material")->SetComponent(i, 0, 32);
			break;
		case 175:
			data.getVtkData()->GetCellData()->GetArray("Material")->SetComponent(i, 0, 32);
			break;
		case 131:
			data.getVtkData()->GetCellData()->GetArray("Material")->SetComponent(i, 0, 32);
			break;
		case 187:
			data.getVtkData()->GetCellData()->GetArray("Material")->SetComponent(i, 0, 32);
			break;
		case 176:
			data.getVtkData()->GetCellData()->GetArray("Material")->SetComponent(i, 0, 32);
			break;

			// in the left atrial
		case 233:
			data.getVtkData()->GetCellData()->GetArray("Material")->SetComponent(i, 0, 33);
			break;
		case 188:
			data.getVtkData()->GetCellData()->GetArray("Material")->SetComponent(i, 0, 33);
			break;
		case 130:
			data.getVtkData()->GetCellData()->GetArray("Material")->SetComponent(i, 0, 33);
			break;
		case 99:
			data.getVtkData()->GetCellData()->GetArray("Material")->SetComponent(i, 0, 33);
			break;
		default:
			break;
		} // switch
	}

	if (data.getInputType() == DataFormat::vtp) {
		for (vtkIdType i = 0; i < data.getVtkData()->GetNumberOfCells(); i++) {
			switch ((int)data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->GetComponent(i, 0)) {
			case 159:
				data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->SetComponent(i, 0, 32);
				break;
			case 175:
				data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->SetComponent(i, 0, 32);
				break;
			case 131:
				data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->SetComponent(i, 0, 32);
				break;
			case 187:
				data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->SetComponent(i, 0, 32);
				break;
			case 176:
				data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->SetComponent(i, 0, 32);
				break;

			case 233:
				data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->SetComponent(i, 0, 33);
				break;
			case 188:
				data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->SetComponent(i, 0, 33);
				break;
			case 130:
				data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->SetComponent(i, 0, 33);
				break;
			case 99:
				data.getVtkData()->GetCellData()->GetArray("MaterialEndo")->SetComponent(i, 0, 33);
				break;
			default:
				break;
			} // switch
		}
	}
} // FiberOrientation::setMaterialtoOldMaterial


vtkSmartPointer<vtkPoints> FiberOrientation::testPointsandMove(DataFormat &data, vtkSmartPointer<vtkPoints> points, Material::Mat inMaterial) {
	for (int i = 0; i < points->GetNumberOfPoints(); i++) {
		if (!::isnan(points->GetPoint(i)[0]) && !::isnan(points->GetPoint(i)[1]) && !::isnan(points->GetPoint(i)[2])) {
			vector<double> point(3);
			point.at(0) = points->GetPoint(i)[0];
			point.at(1) = points->GetPoint(i)[1];
			point.at(2) = points->GetPoint(i)[2];

			if (!Methods::isPointinMaterial(data, point, inMaterial)) {
				vector<double> newPoint = Methods::findClosedPointinMaterialwithoutOrientation(data, point, inMaterial);
				double newPointArray[3] = { newPoint.at(0), newPoint.at(1), newPoint.at(2) };
				points->SetPoint(i, newPointArray);
			}
		}
	}
	return points;
}

/*! To update the cellpointer locator. It can not update otherwise as to write out the data and read in the data. But
 the cell-ids are different after that.
 \param data Pointer to the orginal mesh
 */
DataFormat FiberOrientation::update(DataFormat data) {
	UniversalWriter wirter;
	UniversalReader reader;

	if (data.getInputType() == DataFormat::vtu) {
		string fullstorePath = "temp.vtu";
		wirter.write(data, fullstorePath, false);

		DataFormat newdata = reader.read("temp.vtu", false);
		if (!ConfigGeneral::debug) {
			remove("temp.vtu");
			remove("temp_orginalOrientation.vtu");
		}
		data.setVtkData(newdata.getVtkData());
		data.setCentrePoints(newdata.getCentrePoints());

		return data;
	}
	else if (data.getInputType() == DataFormat::vtp) {
		string fullstorePath = "temp.vtp";
		wirter.write(data, fullstorePath, false);

		DataFormat newdata = reader.read("temp.vtp", false);
		if (!ConfigGeneral::debug) {
			remove("temp.vtp");
			remove("temp_orginalOrientation.vtu");
		}
		data.setVtkData(newdata.getVtkData());
		data.setCentrePoints(newdata.getCentrePoints());
		return data;
	}
	else { exit(-1); }
} // FiberOrientation::update

/*! The main calculation function for the fiber orientation.
 \param data Pointer to the orginal mesh with marked insite and outsite
 \param loadPointPath the path of the seed points txt-file
 \param filename the path for store the meshs with the fiber orientation
 \param freeBridgePath the path of the free bridge points txt-file
 */
void FiberOrientation::calculateFiberOriented(DataFormat &data, string loadPointPath, string filename, string freeBridgePath) {
	FiberOrientation fiberOrientation;

	fiberOrientation.readOutFiberPoints(loadPointPath);

	fiberOrientation.setResolution(Methods::getResolution(data));
	cout << fiberOrientation.getResolution() << endl;

	if (!ConfigFiberorientation::noInitialClean) {
		if (ConfigFiberorientation::rightAtrium) {
			if (data.getInputType() != DataFormat::vtp) {
				Methods::cleanCase(data, fiberOrientation.getRPoint(2));
			}
			else if (data.getInputType() == DataFormat::vtp) {
				Methods::cleanCaseVTP(data, fiberOrientation.getRPoint(3), fiberOrientation.getLPointEpi(0));
			}
		}
		else {
			Methods::cleanCase(data, fiberOrientation.getLPointEpi(2));
		}
	}

	// copy the orginal mesh to turn the fiber into the atrial wall
	DataFormat smoothData;
	smoothData = DataFormat::deepCopy(data);
	Methods::init(smoothData);

	data.setSurfaceForBridges(Methods::extractSurface(data));

	if (ConfigFiberorientation::rightAtrium) {
		if (data.getInputType() != DataFormat::vtp) {
			Methods::uncoupleRightFromLeft(data);
		}
	}

	// substitut the matrial in the surface mesh
	if (data.getInputType() == DataFormat::vtp) {
		for (vtkIdType i = 0; i < data.getVtkData()->GetNumberOfCells(); i++) {
			if (data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(i, 0) == Material::subendo_right_Atrial_Cardium) {
				data.getVtkData()->GetCellData()->GetArray("Material")->SetComponent(i, 0, Material::Vorhof_rechts);
			}
			if (data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(i, 0) == Material::Vorhof_links_Endo) {
				data.getVtkData()->GetCellData()->GetArray("Material")->SetComponent(i, 0, Material::Vorhof_links);
			}
		}
	}

	Methods::init(data);

	// substitut the matrial right endo in the volume mesh
	if (ConfigFiberorientation::rightAtrium) {
		if (data.getInputType() != DataFormat::vtp) {
			for (vtkIdType i = 0; i < data.getVtkData()->GetNumberOfCells(); i++) {
				if (data.getVtkData()->GetCellData()->GetArray("Material")->GetComponent(i, 0) == Material::subendo_right_Atrial_Cardium) {
					data.getVtkData()->GetCellData()->GetArray("Material")->SetComponent(i, 0, Material::Vorhof_rechts);
				}
			}
		}

		fiberOrientation.setRPoints(fiberOrientation.testPointsandMove(data, fiberOrientation.getRPoints(), Material::Vorhof_rechts));
	}


	vtkSmartPointer<vtkPoints> LPoints = vtkSmartPointer<vtkPoints>::New();
	LPoints->DeepCopy(fiberOrientation.getLPointsEpi());
	fiberOrientation.setLPointsEndo(LPoints);

	fiberOrientation.setLPointsEndo(fiberOrientation.testPointsandMove(data, fiberOrientation.getLPointsEndo(), Material::Vorhof_links_Endo));
	fiberOrientation.setLPointsEpi(fiberOrientation.testPointsandMove(data, fiberOrientation.getLPointsEpi(), Material::Vorhof_links));

	if ((ConfigFiberorientation::growPathTransmural || ConfigFiberorientation::growPathTransmuralPectEndo ||
		ConfigFiberorientation::growPathPectfromEndo || ConfigFiberorientation::markCoronary) &&
		(data.getInputType() != DataFormat::vtp)) {
		vtkSmartPointer<vtkPolyData> surface = Methods::triangulateSurface(Methods::extractSurface(data));
		if (data.getVtkData()->GetCellType(data.getVtkData()->GetNumberOfCells() / 2) != 10) {
			surface = Methods::smoothSurface(surface, 1000);
		}

		surface->GetCellData()->RemoveArray("Material");
		surface->GetCellData()->GetArray("RawMaterial")->SetName("Material");
		data.setSurfacewithNormals(surface);
		data.setCentrePointsSurface(Methods::calcCellCentroids(surface));
	}
	else if (data.getInputType() == DataFormat::vtp) {
		ConfigFiberorientation::growPathTransmural = false;
		ConfigFiberorientation::growPathTransmuralPectEndo = false;

		vtkSmartPointer<vtkPolyData> surface = Methods::triangulateSurface(Methods::extractSurface(data));

		surface->GetCellData()->RemoveArray("Material");
		surface->GetCellData()->GetArray("RawMaterial")->SetName("Material");
		data.setSurfacewithNormals(surface);
		data.setCentrePointsSurface(Methods::calcCellCentroids(surface));
	}
	else {
		ConfigFiberorientation::growPathTransmural = false;
		ConfigFiberorientation::growPathTransmuralPectEndo = false;
		ConfigFiberorientation::growPathPectfromEndo = false;
	}


	if (ConfigGeneral::intermediateData && !ConfigFiberorientation::noInitialClean) {
		Methods::writeIntermediateData(data, filename + "_CleanUncouple");
	}

	if (ConfigFiberorientation::rightAtrium) {
		fiberOrientation.setFiberOrientationRightAtrium(data);
	}

	if (ConfigGeneral::intermediateData && ConfigFiberorientation::rightAtrium) {
		Methods::writeIntermediateData(data, filename + "_RightAtriumFiberOrientation");
	}

	fiberOrientation.setFiberOrientationLeftEndoAtrium(data);

	if (ConfigGeneral::intermediateData) {
		Methods::writeIntermediateData(data, filename + "_LeftEndoAtriumFiberOrientation");
	}

	fiberOrientation.setFiberOrientationLeftEpiAtrium(data);

	if (ConfigGeneral::intermediateData) {
		Methods::writeIntermediateData(data, filename + "_LeftEpiAtriumFiberOrientation");
	}

	fiberOrientation.setFiberOrientationLeftVeins(data);

	if (ConfigGeneral::intermediateData) {
		Methods::writeIntermediateData(data, filename + "_LeftVeinsFiberOrientation");
	}

	fiberOrientation.setFiberOrientationLeftAppendes(data);

	if (ConfigGeneral::intermediateData || ConfigGeneral::debug) {
        std::string storename = filename + "_LeftAppendageAtriumFiberOrientation";
        if (ConfigGeneral::debug) {
            storename = "LeftAppendageAtriumFiberOrientation_debug";
        }
		Methods::writeIntermediateData(data, storename);
	}

	fiberOrientation.closeOrientation(data);
    
    if (ConfigGeneral::intermediateData || ConfigGeneral::debug) {
        std::string storename = filename + "_CloseFiberOrientation";
        if (ConfigGeneral::debug) {
            storename = "CloseFiberOrientation_debug";
        }
        Methods::writeIntermediateData(data, storename);
    }

	if (ConfigFiberorientation::rightAtrium) {
		fiberOrientation.setBridges(data, freeBridgePath);

        if (ConfigGeneral::intermediateData || ConfigGeneral::debug) {
            std::string storename = filename + "_Bridges";
            if (ConfigGeneral::debug) {
                storename = "Bridges_debug";
            }
            Methods::writeIntermediateData(data, storename);
        }
	}
	if (!ConfigFiberorientation::noEndCleanup) {
		if (ConfigFiberorientation::rightAtrium) {
			Methods::endCleanCase(data, fiberOrientation.getRPoint(2));
		}
		else {
			Methods::endCleanCase(data, fiberOrientation.getLPointEpi(2));
		}
	}
	fiberOrientation.closeOrientation(data);

    if (ConfigGeneral::intermediateData || ConfigGeneral::debug) {
        std::string storename = filename + "_CloseAndClean";
        if (ConfigGeneral::debug) {
            storename = "CloseAndClean_debug";
        }
        Methods::writeIntermediateData(data, storename);
    }

	if (ConfigFiberorientation::rightAtrium) {
		fiberOrientation.setSinusNode(data);

		if (ConfigGeneral::intermediateData) {
			Methods::writeIntermediateData(data, filename + "_Sinus");
		}
	}

	fiberOrientation.smoothOrientationToSurface(data, smoothData);

	if (ConfigGeneral::intermediateData) {
		Methods::writeIntermediateData(data, filename + "_SmoothFiberOrientationToSurface");
	}

	Methods::setThetaPhi(data);

	if (ConfigGeneral::intermediateData) {
		Methods::writeIntermediateData(data, filename + "_PhiThetaNormalFormat");
	}

	fiberOrientation.setMaterialtoOldMaterial(data);

	if (ConfigGeneral::intermediateData) {
		Methods::writeIntermediateData(data, filename + "_SubstitutionMatrial");
	}

	Methods::deleteArrays(data);

	if (ConfigGeneral::intermediateData) {
		Methods::writeIntermediateData(data, filename + "_Finish");
	}
} // FiberOrientation::calculateFiberOriented

/*!
 \page FiberOrientation

 \section DESCRIPTION_FiberOrientation DESCRIPTION
 Definition file for the annotation of the fiber orientation and insertion of the automatic bridges.

 \section SOURCE_FiberOrientation SOURCE

 FiberOrientation.cpp

 \section SEEALSO_FiberOrientation SEE ALSO
 \ref Methods \ref DataFormat \ref Reader \ref Writer \ref VcgTools \ref Material \ref Config

 \section CHANGELOG_FiberOrientation CHANGELOG
 V1.0.0 - 01.09.2015 (Andreas Wachter): Starting with changelog\n
 */

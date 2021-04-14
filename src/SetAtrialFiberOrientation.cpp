/*! \file SetAtrialFiberOrientation.cpp

   \brief Tool to set automatic the fiber orientation and interatrial bridges with the 22 seed points and free bridge points.

   \version 1.1.0

   \date Andreas Wachter 01.09.15

   \author Andreas Wachter\n
   Institute of Biomedical Engineering\n
   Kit\n
   http://www.ibt.kit.de\n

   \sa Synopsis \ref SetAtrialFiberOrientation
 */

#include "SetAtrialFiberOrientation.h"

int main(int argc, char*const argv[]) {

	for (int i = 0; i < argc; i++) {
		if (strcmp(argv[i], "-help") == 0 || argc == 0) {
			SetAtrialFiberOrientation::printOptions();
			exit(0);
		}
	}
	std::string freeBridgePath = "";
	std::string intermediateResultPrefix = "";

	bool prefix = false;

	try {
		if (argc > 3) {
			std::string fullLoadPath = argv[1];
			std::string fullStorePath = argv[2];
			std::string fullPointPath = argv[3];

			for (int i = 4; i < argc; i++) {
				if ((strcmp(argv[i], "-intermediateResultPrefix") == 0) && (i + 1 < argc)) {
					ConfigGeneral::intermediateData = true;
					intermediateResultPrefix = argv[++i];
					prefix = true;
				}
				else if (strcmp(argv[i], "-writeIntermediateData") == 0) {
					ConfigGeneral::intermediateData = true;
				}
				else if (strcmp(argv[i], "-noOpenVTPBridgeInAtrial") == 0) {
					ConfigFiberorientation::openVTPBridgeInAtrial = false;
				}
				else if (strcmp(argv[i], "-noBachmannBundle") == 0) {
					ConfigFiberorientation::bachmannBundle = false;
				}
				else if (strcmp(argv[i], "-noRightAtrium") == 0) {
					ConfigFiberorientation::rightAtrium = false;
				}
				else if (strcmp(argv[i], "-noWideBachmann") == 0) {
					ConfigFiberorientation::wideBachmann = false;
				}
				else if (strcmp(argv[i], "-noUpperPosteriorBridge") == 0) {
					ConfigFiberorientation::upperPosteriorBridge = false;
				}
				else if (strcmp(argv[i], "-noMiddlePosteriorBridge") == 0) {
					ConfigFiberorientation::middlePosteriorBridge = false;
				}
				else if (strcmp(argv[i], "-noCoronarySinusBridge") == 0) {
					ConfigFiberorientation::coronarySinusBridge = false;
				}
				else if (strcmp(argv[i], "-lowerAnteriorBridge") == 0) {
					ConfigFiberorientation::lowerAnteriorBridge = true;
				}
				else if (strcmp(argv[i], "-upperAnteriorBridge") == 0) {
					ConfigFiberorientation::upperAnteriorBridge = true;
				}
				else if (strcmp(argv[i], "-noIntercavalBundle") == 0) {
					ConfigFiberorientation::noIntercavalBundle = true;
				}
				else if (strcmp(argv[i], "-writeViewFiberOrientation") == 0) {
					ConfigFiberorientation::writeViewFiberOrientation = true;
				}
				else if (strcmp(argv[i], "-debug") == 0) {
					ConfigGeneral::debug = true;
				}
				else if ((strcmp(argv[i], "-numberOfPectinateMuscles") == 0) && (i + 1 < argc)) {
					ConfigFiberorientation::numbersofPectinates = std::atoi(argv[++i]);
				}
				else if ((strcmp(argv[i], "-growPathOffestLeft") == 0) && (i + 1 < argc)) {
					ConfigFiberorientation::growPathOffsetLeft = std::stod(argv[++i]);
				}
				else if (strcmp(argv[i], "-growPathTransmural") == 0) {
					ConfigFiberorientation::growPathTransmural = true;
				}
				else if (strcmp(argv[i], "-growPathPectfromEndo") == 0) {
					ConfigFiberorientation::growPathPectfromEndo = true;
				}
				else if (strcmp(argv[i], "-growPathTransmuralPectEndo") == 0) {
					ConfigFiberorientation::growPathTransmuralPectEndo = true;
				}
				else if (strcmp(argv[i], "-markCoronary") == 0) {
					ConfigFiberorientation::markCoronary = true;
				}
				else if ((strcmp(argv[i], "-freeBridge") == 0) && (i + 1 < argc)) {
					ConfigFiberorientation::freeBridge = true;
					freeBridgePath = argv[++i];
				}
				else if (strcmp(argv[i], "-noInitialClean") == 0) {
					ConfigFiberorientation::noInitialClean = true;
				}
				else if (strcmp(argv[i], "-noEndCleanup") == 0) {
					ConfigFiberorientation::noEndCleanup = true;
				}
				else if (strcmp(argv[i], "-growPectToRAppendage") == 0) {
					ConfigFiberorientation::growPectToRAppendage = true;
				}
				else if (strcmp(argv[i], "-growCristaDepOnWallthick") == 0) {
					ConfigFiberorientation::growCristaDepOnWallthick = true;
				}
				else {
					std::cerr << "argument: " << argv[i] << " not supported" << std::endl;
					exit(-1);
				}
			}


			long pathPos = fullStorePath.find_last_of("/");

			if (ConfigGeneral::intermediateData && !prefix) {
				long extensionPos = fullStorePath.find_last_of(".");
				intermediateResultPrefix = fullStorePath.substr(pathPos + 1, extensionPos - (pathPos + 1));
			}

			std::string interMediaFilename = fullStorePath.substr(0, pathPos + 1) + intermediateResultPrefix;

			UniversalReader reader;
			DataFormat data;

			data = reader.read(fullLoadPath);

			if (data.getInputType() == DataFormat::vtp) {
				data.setDoubleLayer(false);
			}
			else {
				data.setDoubleLayer(true);
			}

			std::cout << "Number of cells: " << data.getVtkData()->GetNumberOfCells() << std::endl;
			std::cout << "Number of points: " << data.getVtkData()->GetNumberOfPoints() << std::endl;


			FiberOrientation::calculateFiberOriented(data, fullPointPath, interMediaFilename, freeBridgePath);


			UniversalWriter wirter;
			wirter.write(data, fullStorePath);
		}
		else {
			SetAtrialFiberOrientation::printOptions();
		}
	}
	catch (int e) {
		std::cerr << argv[0] << " Error: " << e << "\n";
		exit(-1);
	}
	return 0;
}

void  SetAtrialFiberOrientation::printOptions() {
	std::cerr << "\t" << " SetAtrialFiberOrientation: This programm can used for setting semiautomatic atrial fieberorientation." << "\n";

	std::cerr << "\t" << "<" << "path for input datafile *.vtu | vtp | vtk" << ">" << "\n";
	std::cerr << "\t" << "<" << "path for output datafile *.vtu|vtp|vtk" << ">" << "\n";
	std::cerr << "\t" << "<" << "path for seedpoints file *.txt" << ">" << "\n";

	std::cerr << "\t" << "[" << "-intermediateResultPrefix" << "] " << "(default: inputdataname): prefix for intermediate results" << "\n";
	std::cerr << "\t" << "[" << "-writeIntermediateData" << "]: " << "for saving intermediate Results" << "\n";
	std::cerr << "\t" << "[" << "-noOpenVTPBridgeInAtrial" << "]: " << "open atrial Bridges in VTP" << "\n";
	std::cerr << "\t" << "[" << "-noRightAtrium" << "]" << "\n";
	std::cerr << "\t" << "[" << "-noBachmannBundle" << "]" << "\n";
	std::cerr << "\t" << "[" << "-noWideBachmann" << "]: " << "if you wouldn't have wide Bachmann Bundle from the right to the left atrial appendeage" << "\n";
	std::cerr << "\t" << "[" << "-noUpperPosteriorBridge" << "]" << "\n";
	std::cerr << "\t" << "[" << "-noMiddlePosteriorBridge" << "]" << "\n";
	std::cerr << "\t" << "[" << "-noCoronarySinusBridge" << "]" << "\n";
	std::cerr << "\t" << "[" << "-lowerAnteriorBridge" << "]" << "\n";
	std::cerr << "\t" << "[" << "-upperAnteriorBridge" << "]" << "\n";
	std::cerr << "\t" << "[" << "-writeViewFiberOrientation" << "]: " << "write out the orginal fiber orientation in the cell for illustration" << "\n";
	std::cerr << "\t" << "[" << "-numberOfPectinateMuscles N" << "] " << "(default: 15): number of pectinate muscles" << "\n";
	std::cerr << "\t" << "[" << "-growPathOffestLeft N" << "] " << "(default: 0): additiv Offset for the path dilatation in left Atrium" << "\n";
	std::cerr << "\t" << "[" << "-writeViewFiberOrientation" << "]: " << "growPath transmural in the atria" << "\n";
	std::cerr << "\t" << "[" << "-growPathTransmuralPectEndo" << "]: " << "grow endocardial Pectinate muscles in the atria" << "\n";
	std::cerr << "\t" << "[" << "-growPathPectFromEndo" << "]: " << "grow radial endocardial Pectinate muscles in the atria from the endocardial surface" << "\n";
	std::cerr << "\t" << "[" << "-growPectToRAppendage" << "]: " << "grow Pectinate muscles in the right atria appendeage" << "\n";
	std::cerr << "\t" << "[" << "-growCristaDepOnWallthick" << "]: " << "grow Crista Teriminals depending on avarage atrial wall thickness" << "\n";
	std::cerr << "\t" << "[" << "-markCoronary" << "]: " << "mark the coronary sinus" << "\n";
	std::cerr << "\t" << "[" << "-debug" << "]: " << "if you would have more information when program abort" << "\n";
	std::cerr << "\t" << "[" << "-freeBridge PATH" << "]: " << "if you would set free define bridges + definitionfilepath" << "\n";
	std::cerr << "\t" << "[" << "-noInitialClean" << "]: " << "no inital clean of the mesh e.g. when they are disconnected" << "\n";
	std::cerr << "\t" << "[" << "-noEndCleanup" << "]: " << "no end clean of the mesh e.g. when they are disconnected" << "\n";
	std::cerr << "\t" << "[" << "-debug" << "]: " << "if you would have more information when program abort" << "\n";
}

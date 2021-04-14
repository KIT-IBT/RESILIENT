/*! \file Config.cpp

   \brief Here are stored the default configuration settings for RESILIENT.
   \version 1.0.0

   \date Andreas Wachter 01.09.15

   \author Andreas Wachter\n
   Institute of Biomedical Engineering\n
   Kit\n
   http://www.ibt.kit.de\n

   \sa Synopsis \ref Config
 */

#include "Config.h"

bool ConfigFiberorientation::openVTPBridgeInAtrial     = true;
bool ConfigFiberorientation::wideBachmann              = true;
bool ConfigFiberorientation::upperPosteriorBridge      = true;
bool ConfigFiberorientation::coronarySinusBridge       = true;
bool ConfigFiberorientation::middlePosteriorBridge     = true;
bool ConfigFiberorientation::lowerAnteriorBridge       = false;
bool ConfigFiberorientation::upperAnteriorBridge       = false;
bool ConfigFiberorientation::noIntercavalBundle        = false;
bool ConfigFiberorientation::bachmannBundle            = true;
bool ConfigFiberorientation::rightAtrium               = true;
int  ConfigFiberorientation::numbersofPectinates       = 15;
bool ConfigFiberorientation::writeViewFiberOrientation = false;
bool ConfigFiberorientation::freeBridge                = false;
bool ConfigFiberorientation::noInitialClean            = false;
bool ConfigFiberorientation::noEndCleanup              = false;
bool ConfigFiberorientation::markCoronary              = false;
bool ConfigFiberorientation::growPathTransmural 	   = false;
bool ConfigFiberorientation::growCristaDepOnWallthick  = false;
bool ConfigFiberorientation::growPathTransmuralPectEndo = false;
bool ConfigFiberorientation::growPathPectfromEndo      = false;
bool ConfigFiberorientation::growPectToRAppendage      = false;

/*! The value are in mm.*/
double ConfigFiberorientation::growPathOffsetLeft = 0;

/*! The value are in mm. It is for marked the fiber orientation transmural.*/
double ConfigFiberorientation::maxWitdthAtrialWall = 15;

bool ConfigGeneral::intermediateData                             = false;
bool ConfigGeneral::debug                                        = false;

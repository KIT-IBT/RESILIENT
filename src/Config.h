/*!Here are stored the default configuration settings for RESILIENT*/

//
//  Config.h
//  RESILIENT
//
//  Created by Andreas Wachter on 01.09.15.
//  Copyright (c) 2015 IBT. All rights reserved.
//

#ifndef _RESILIENT_Config_
#define _RESILIENT_Config_

/*!Here are stored the default configuration settings for the fiber orientation tool and the ablation tool.*/

class ConfigFiberorientation {
 public:
  static bool openVTPBridgeInAtrial;
  static bool wideBachmann;
  static bool upperPosteriorBridge;
  static bool coronarySinusBridge;
  static bool lowerAnteriorBridge;
  static bool upperAnteriorBridge;
  static bool middlePosteriorBridge;
  static bool bachmannBundle;
  static bool noIntercavalBundle;
  static int numbersofPectinates;
  static bool writeViewFiberOrientation;
  static bool freeBridge;
  static double growPathOffsetLeft;
  static bool growPathTransmural;
  static bool growPathTransmuralPectEndo;
  static bool growPathPectfromEndo;
  static bool growCristaDepOnWallthick;
  static bool growPectToRAppendage;
  static bool markCoronary;
  static double maxWitdthAtrialWall;
  static bool rightAtrium;
  static bool noInitialClean;
  static bool noEndCleanup;

};

/*!Here are stored the default configuration settings for the fiber orientation tool and the ablation tool.*/

class ConfigGeneral {
 public:
  static bool intermediateData;
  static bool debug;
};

#endif /* defined(_RESILIENT_Config_) */

//

//  avergingOrientation.h
//  RESILIENT
//
//  Created by Andreas Wachter on 04.11.14.
//  Copyright (c) 2014 IBT. All rights reserved.
//

#ifndef _RESILIENT_avergingOrientation_
#define _RESILIENT_avergingOrientation_

#include <iostream>
#include <vector>
#include <list>
#include <tuple>
#include "Headers.h"

class AveragingOrientation 
{
public:
	static std::vector<double> averageOrientation(std::vector<double> vector1, std::vector<double> vector2, double distanceToVector1, double distanceToVector2);
	static std::vector<double> AveragingNVectors(std::list<std::tuple<int, double, std::vector<double>, std::vector<double>>> vectorList);
};

#endif /* defined(_RESILIENT_avergingOrientation_) */


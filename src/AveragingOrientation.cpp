/*! \file AveragingOrientation.cpp

 \brief For calculating the interpolated the fiber orientation.
 \version 1.0.0

 \date Andreas Wachter 01.09.15

 \author Andreas Wachter\n
 Institute of Biomedical Engineering\n
 Kit\n
 http://www.ibt.kit.de\n

 \sa Synopsis \ref AveragingOrientation
 */

#include "AveragingOrientation.h"

 /*! For calculation the weigthed interpolated orientation between two fiber orientations.
  \param vector1 Fiber orientation from cell 1 as normed vector.
  \param vector2 Fiber orientation from cell 2 as normed vector.
  \param distanceToVector1 The distance between the current calculated cell and the cell 1.
  \param distanceToVector1 The distance between the current calculated cell and the cell 2.
  \return The interpolated fiber orientation as normed vector with x,y,z component.
  */

std::vector<double> AveragingOrientation::averageOrientation(std::vector<double> vector1, std::vector<double> vector2, double distanceToVector1, double distanceToVector2) {
	double sumofDistance = distanceToVector1 + distanceToVector2;

	std::list<std::tuple<int, double, std::vector<double>, std::vector<double>>> differenceVectorList;
	double distance1 = 1 - distanceToVector1 / sumofDistance;
	std::tuple<int, double, std::vector<double>, std::vector<double>> insert1(0, distance1, vector1, vector1);

	differenceVectorList.push_back(insert1);

	double distance2 = 1 - distanceToVector2 / sumofDistance;
	std::tuple<int, double, std::vector<double>, std::vector<double>> insert2(0, distance2, vector2, vector2);
	differenceVectorList.push_back(insert2);
	return AveragingNVectors(differenceVectorList);
}


/*! Define an tuple that contains: material, disatnce, fiber orientation and main fiber orientation of the bundle*/
typedef std::tuple<int, double, std::vector<double>, std::vector<double>> mytuple;

/*! Define an bool equation. If the second element of the tupele is bigger than the second element of the second tupel
 -> than true otherwise false. */
bool mycompare(const mytuple &lhs, const mytuple &rhs) {
	return std::get<1>(lhs) > std::get<1>(rhs);
}

/*! For calculation the weigthed interpolated orientation between many fiber orientations.
 \param vectorList Is a list with the tupels. A tuple contains: material, disatnce, fiber orientation and main fiber
 orientation of the bundle
 \return The interpolated fiber orientation as normed vector with x,y,z component.
 */
std::vector<double> AveragingOrientation::AveragingNVectors(std::list<std::tuple<int, double, std::vector<double>, std::vector<double>>> vectorList) {
	if (vectorList.size() > 1) {
		vectorList.sort(mycompare);
		std::tuple<int, double, std::vector<double>, std::vector<double>> vector1 = vectorList.front();

		double vector1Array[3] = { std::get<2>(vector1).at(0), std::get<2>(vector1).at(1), std::get<2>(vector1).at(2) };
		double vector1MainorientationArray[3] = { std::get<3>(vector1).at(0), std::get<3>(vector1).at(1), std::get<3>(vector1).at(2) };
		vtkMath::Normalize(vector1Array);
		vtkMath::MultiplyScalar(vector1Array, std::get<1>(vector1));
		vectorList.pop_front();
		int k = 0;
		int j = 0;

		long size = vectorList.size();

		while (!vectorList.empty()) {
			std::tuple<int, double, std::vector<double>, std::vector<double>> vector2 = vectorList.front();
			vectorList.pop_front();
			size--;
			double vector2Array[3] = { std::get<2>(vector2).at(0), std::get<2>(vector2).at(1), std::get<2>(vector2).at(2) };
			double vector2MainorientationArray[3] = { std::get<3>(vector2).at(0), std::get<3>(vector2).at(1), std::get<3>(vector2).at(2) };

			vtkMath::Normalize(vector2Array);

			double dotProduct = vtkMath::Dot(vector1MainorientationArray, vector2MainorientationArray);
			double normVector1 = vtkMath::Norm(vector1MainorientationArray);
			double scalarProduct = dotProduct / normVector1;

			if (((scalarProduct < 0.1) && (scalarProduct > -0.1)) && (j < size) && !vectorList.empty()) {
				vectorList.push_back(vector2);
				size++;
				j++;
			}
			else {
				if (scalarProduct > 0) {
					vtkMath::MultiplyScalar(vector2Array, std::get<1>(vector2));
					vtkMath::Add(vector1Array, vector2Array, vector1Array);

					vector1MainorientationArray[0] = vector1Array[0];
					vector1MainorientationArray[1] = vector1Array[1];
					vector1MainorientationArray[2] = vector1Array[2];
				}
				else if (scalarProduct < 0) {
					vtkMath::MultiplyScalar(vector2Array, -1);
					vtkMath::MultiplyScalar(vector2Array, std::get<1>(vector2));
					vtkMath::Add(vector1Array, vector2Array, vector1Array);
					vector1MainorientationArray[0] = vector1Array[0];
					vector1MainorientationArray[1] = vector1Array[1];
					vector1MainorientationArray[2] = vector1Array[2];
				}
				else if (scalarProduct == 0) {
					if (k > size) {
						// cout << "Vector orthogonal? Vectors only would be add!" << endl;
						while (!vectorList.empty()) {
							std::tuple<int, double, std::vector<double>, std::vector<double>> vector2 = vectorList.front();
							vectorList.pop_front();
							double vector2Array[3] = { std::get<2>(vector2).at(0), std::get<2>(vector2).at(1), std::get<2>(vector2).at(2) };
							vtkMath::MultiplyScalar(vector2Array, std::get<1>(vector2));
							vtkMath::Add(vector1Array, vector2Array, vector1Array);

							vector1MainorientationArray[0] = vector1Array[0];
							vector1MainorientationArray[1] = vector1Array[1];
							vector1MainorientationArray[2] = vector1Array[2];
						}
					}
					else {
						vectorList.push_back(vector2);
						k++;
						size++;
					}
				}
			}
		}

		vtkMath::Normalize(vector1Array);
		std::vector<double> averageingVector = { vector1Array[0], vector1Array[1], vector1Array[2] };

		return averageingVector;
	}
	else {
		return std::get<2>(vectorList.front());
	}
} // AveragingOrientation::AveragingNVectors

/*!
 \page AveragingOrientation

 \section DESCRIPTION_AveragingOrientation DESCRIPTION
 For calculating the interpolated the fiber orientation.

 \section SOURCE_FiberOrientation SOURCE

 AveragingOrientation.cpp

 \section SEEALSO_AveragingOrientation SEE ALSO
 \ref Methods \ref Matrial

 \section CHANGELOG_AveragingOrientation CHANGELOG
 V1.0.0 - 01.09.2015 (Andreas Wachter): Starting with changelog\n
 */

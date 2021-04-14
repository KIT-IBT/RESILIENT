/*!All tissue classes and functions for tissues.*/

/*! \file Material.h

 \brief All tissue classes and functions for tissues.
 \version 1.0.0

 \date Andreas Wachter 01.09.15

 \author Andreas Wachter\n
 Institute of Biomedical Engineering\n
 Kit\n
 http://www.ibt.kit.de\n

 \sa Synopsis \ref Material
 */

#ifndef _RESILIENT_Material_
#define _RESILIENT_Material_

#include <iostream>
class Material {
public:
	/*!All tissue classes.*/
	enum Mat {
		undef = 0,
		Ventrikel_rechts = 30,
		Ventrikel_links = 31,
		Vorhof_rechts = 32,
		Vorhof_links = 33,
		Sinus_Knoten = 34,
		Crista_Terminalis_2 = 72,
		Pektinatmuskeln = 74,

		subepi_Atrial_Cardium = 99,
		Crista_Terminalis = 100,
		Pectinate_Muscle = 101,
		Bachmann_Bundle_right_atrium_or_all = 102,
		Bachmann_Bundle_left_atrium = 103,
		Inferior_Isthmus_right_atrium = 104,
		Coronary_Sinus_ostium_tissue = 105,
		Right_Atrial_Appendage = 111,
		Left_Atrial_Appendage = 112,

		subendo_right_Atrial_Cardium = 129,
		subendo_left_Atrial_Cardium = 130,
		not_right_Isthmus = 131,

		Blood_right_atrium = 150,
		Blood_left_atrium = 151,
		SVC = 159,
		IVC = 160,
		Tricuspid_Valve_Ring = 161,
		Mitral_Valve_Ring = 162,
		RIPV = 166,
		RSPV = 167,
		LIPV = 168,
		LSPV = 169,

		IntercavanalBundel = 175,
		Tricusline = 176,
		otherwise_right_or_all = 187,
		otherwise_left = 188,

		Interatrial_Bridge_left = 190,
		Interatrial_Bridge_right = 191,
		Middle_Posterior_Bridge_left = 192,
		Middle_Posterior_Bridge_right = 193,
		Upper_Posterior_Bridge_left = 194,
		Upper_Posterior_Bridge_right = 195,
		Lower_Anterior_Bridge_left = 196,
		Lower_Anterior_Bridge_right = 197,
		Coronary_Sinus_Bridge_left = 198,
		Coronary_Sinus_Bridge_right = 199,
		Upper_Anterior_Bridge_left = 200,
		Upper_Anterior_Bridge_right = 201,

		Vorhof_links_Endo = 233,

		insert_new_Material = 300,

		testMaterialRight = 500,
		testMaterialLeft = 501,
		testMaterialLeftEndo = 502
	};

	/*! Test if a tissue class is in left epicard.
	\param material Tissue class that would be tested.
	\return True if the tissue is in left epicard, else False.
	*/
	static bool isInLeftEpi(int material) {
		if ((material == Vorhof_links) || (material == subepi_Atrial_Cardium) ||
			(material == Bachmann_Bundle_left_atrium) || (material == Left_Atrial_Appendage) || (material == RIPV) ||
			(material == RSPV) || (material == LIPV) || (material == LSPV) || (material == Mitral_Valve_Ring) ||
			(material == otherwise_left) || (material == Upper_Anterior_Bridge_left) ||
			(material == Middle_Posterior_Bridge_left) || (material == Upper_Posterior_Bridge_left) ||
			(material == Lower_Anterior_Bridge_left) || (material == Coronary_Sinus_Bridge_left) ||
			(material == insert_new_Material)) {
			return true;
		}

		return false;
	}

	/*! Test if a tissue class is in left atrial.
	\param material Tissue class that would be tested.
	\return True if the tissue is in left atrial, else False.
	*/
	static bool isInLeft(int material) {
		if ((material == Vorhof_links) || (material == subepi_Atrial_Cardium) ||
			(material == Bachmann_Bundle_left_atrium) || (material == Left_Atrial_Appendage) || (material == RIPV) ||
			(material == RSPV) || (material == LIPV) || (material == LSPV) || (material == Mitral_Valve_Ring) ||
			(material == otherwise_left) || (material == Upper_Anterior_Bridge_left) ||
			(material == Middle_Posterior_Bridge_left) || (material == Upper_Posterior_Bridge_left) ||
			(material == Lower_Anterior_Bridge_left) || (material == Coronary_Sinus_Bridge_left) ||
			(material == Vorhof_links_Endo) || (material == insert_new_Material) || (material == subendo_left_Atrial_Cardium)) {
			return true;
		}

		return false;
	}

	/*! Test if a tissue class is in left endocard.
	\param material Tissue class that would be tested.
	\return True if the tissue is in left endocard, else False.
	*/
	static bool isInLeftEndo(int material) {
		if ((material == Vorhof_links_Endo) || (material == Left_Atrial_Appendage)) {
			return true;
		}

		return false;
	}

	/*! Test if a tissue class is in endocard.
	\param material Tissue class that would be tested.
	\return True if the tissue is in endocard, else False.
	*/
	static bool isInEndo(int material) {
		if ((material == Vorhof_links_Endo) || (material == subendo_right_Atrial_Cardium)) {
			return true;
		}

		return false;
	}

	/*! Test if a tissue class is in epicard.
	\param material Tissue class that would be tested.
	\return True if the tissue is in epidocard, else False.
	*/
	static bool isInEpi(int material) {
		if ((material == Vorhof_links) || (material == subepi_Atrial_Cardium) ||
			(material == Bachmann_Bundle_left_atrium) || (material == Left_Atrial_Appendage) || (material == RIPV) ||
			(material == RSPV) || (material == LIPV) || (material == LSPV) || (material == Mitral_Valve_Ring) ||
			(material == otherwise_left) || (material == Upper_Anterior_Bridge_left) ||
			(material == Middle_Posterior_Bridge_left) || (material == Upper_Posterior_Bridge_left) ||
			(material == Lower_Anterior_Bridge_left) || (material == Coronary_Sinus_Bridge_left) ||
			(material == insert_new_Material) || (material == Vorhof_rechts) || (material == Sinus_Knoten) || (material == IntercavanalBundel) ||
			(material == SVC) || (material == Crista_Terminalis) || (material == Crista_Terminalis_2) ||
			(material == Pectinate_Muscle) || (material == Bachmann_Bundle_right_atrium_or_all) ||
			(material == Inferior_Isthmus_right_atrium) || (material == Right_Atrial_Appendage) ||
			(material == not_right_Isthmus) || (material == Tricuspid_Valve_Ring) || (material == Tricusline) ||
			(material == otherwise_right_or_all) || (material == Upper_Anterior_Bridge_right) ||
			(material == Middle_Posterior_Bridge_right) || (material == Upper_Posterior_Bridge_right) ||
			(material == Lower_Anterior_Bridge_right) || (material == Coronary_Sinus_Bridge_right) ||
			(material == insert_new_Material) || (material == Pektinatmuskeln)) {
			return true;
		}

		return false;
	}


	/*! Test if a tissue class is in right atrial.
	\param material Tissue class that would be tested.
	\return True if the tissue is in right atrial, else False.
	*/
	static bool isInRight(int material) {
		if ((material == Vorhof_rechts) || (material == Sinus_Knoten) || (material == IntercavanalBundel) ||
			(material == SVC) || (material == Crista_Terminalis) || (material == Crista_Terminalis_2) ||
			(material == Pectinate_Muscle) || (material == Bachmann_Bundle_right_atrium_or_all) ||
			(material == Inferior_Isthmus_right_atrium) || (material == Right_Atrial_Appendage) ||
			(material == not_right_Isthmus) || (material == Tricuspid_Valve_Ring) || (material == Tricusline) ||
			(material == otherwise_right_or_all) || (material == Upper_Anterior_Bridge_right) ||
			(material == Middle_Posterior_Bridge_right) || (material == Upper_Posterior_Bridge_right) ||
			(material == Lower_Anterior_Bridge_right) || (material == Coronary_Sinus_Bridge_right) || (material == Coronary_Sinus_ostium_tissue) ||
			(material == insert_new_Material) || (material == subendo_right_Atrial_Cardium) || (material == Pektinatmuskeln)) {
			return true;
		}
		return false;
	}

	/*! Test if a tissue class is in right atrial without the bridges.
	\param material Tissue class that would be tested.
	\return True if the tissue is in right atrial without the bridges, else False.
	*/
	static bool isInRightWithoutBridges(int material) {
		if ((material == Vorhof_rechts) || (material == Sinus_Knoten) || (material == IntercavanalBundel) ||
			(material == SVC) || (material == Crista_Terminalis) || (material == Crista_Terminalis_2) ||
			(material == Pectinate_Muscle) || (material == Inferior_Isthmus_right_atrium) ||
			(material == Right_Atrial_Appendage) || (material == not_right_Isthmus) || (material == Tricuspid_Valve_Ring) ||
			(material == Tricusline) || (material == otherwise_right_or_all) ||
			(material == subendo_right_Atrial_Cardium) || (material == Coronary_Sinus_ostium_tissue) ||
			(material == Pektinatmuskeln)) {
			return true;
		}
		return false;
	}

	/*! Test if a tissue class is in right atrial without the right endocard.
	\param material Tissue class that would be tested.
	\return True if the tissue is in right atrial without the right endocard, else False.
	*/
	static bool isInRightWithoutEndo(int material) {
		if ((material == Vorhof_rechts) || (material == Sinus_Knoten) || (material == IntercavanalBundel) ||
			(material == SVC) || (material == Crista_Terminalis) || (material == Crista_Terminalis_2) ||
			(material == Pectinate_Muscle) || (material == Bachmann_Bundle_right_atrium_or_all) ||
			(material == Inferior_Isthmus_right_atrium) || (material == Right_Atrial_Appendage) ||
			(material == not_right_Isthmus) || (material == Tricuspid_Valve_Ring) || (material == Tricusline) ||
			(material == otherwise_right_or_all) || (material == Upper_Anterior_Bridge_right) ||
			(material == Middle_Posterior_Bridge_right) || (material == Upper_Posterior_Bridge_right) || (material == Coronary_Sinus_ostium_tissue) ||
			(material == Lower_Anterior_Bridge_right) || (material == Coronary_Sinus_Bridge_right) ||
			(material == insert_new_Material) || (material == Pektinatmuskeln)) {
			return true;
		}
		return false;
	}

	/*! Test if a tissue class is in left atrial without the bridges.
	\param material Tissue class that would be tested.
	\return True if the tissue is in left atrial without the bridges, else False.
	*/
	static bool isInLeftWithoutBridge(int material) {
		if ((material == Vorhof_links) || (material == subepi_Atrial_Cardium) || (material == Left_Atrial_Appendage) ||
			(material == RIPV) || (material == RSPV) || (material == LIPV) || (material == LSPV) ||
			(material == Mitral_Valve_Ring) || (material == otherwise_left) || (material == subendo_left_Atrial_Cardium) ||
			(material == Vorhof_links_Endo)) {
			return true;
		}

		return false;
	}

	/*! Test if a tissue class is in right atrial without the bridges and the blood.
	\param material Tissue class that would be tested.
	\return True if the tissue is in right atrial without the bridges and the blood, else False.
	*/
	static bool isInRightWithoutBridgesAndWithBlood(int material) {
		if ((material == Vorhof_rechts) || (material == Sinus_Knoten) || (material == IntercavanalBundel) ||
			(material == SVC) || (material == Crista_Terminalis) || (material == Crista_Terminalis_2) ||
			(material == Pectinate_Muscle) || (material == Inferior_Isthmus_right_atrium) ||
			(material == Right_Atrial_Appendage) || (material == not_right_Isthmus) || (material == Tricuspid_Valve_Ring) ||
			(material == Tricusline) || (material == otherwise_right_or_all) ||
			(material == subendo_right_Atrial_Cardium) || (material == Coronary_Sinus_ostium_tissue) ||
			(material == Blood_right_atrium) || (material == Pektinatmuskeln)) {
			return true;
		}
		return false;
	}

	/*! Test if a tissue class is in left atrial without the bridges and the blood.
	\param material Tissue class that would be tested.
	\return True if the tissue is in left atrial without the bridges and the blood, else False.
	*/
	static bool isInLeftWithoutBridgeAndWithBlood(int material) {
		if ((material == Vorhof_links) || (material == subepi_Atrial_Cardium) || (material == Left_Atrial_Appendage) ||
			(material == RIPV) || (material == RSPV) || (material == LIPV) || (material == LSPV) ||
			(material == Mitral_Valve_Ring) || (material == otherwise_left) || (material == subendo_left_Atrial_Cardium) ||
			(material == Blood_left_atrium)) {
			return true;
		}

		return false;
	}

	/*! Test if a tissue class is in right atrial with the blood.
	\param material Tissue class that would be tested.
	\return True if the tissue is in right atrial with the blood, else False.
	*/
	static bool isInRightWithBlood(int material) {
		if ((material == Vorhof_rechts) || (material == Sinus_Knoten) || (material == IntercavanalBundel) ||
			(material == SVC) || (material == Crista_Terminalis) || (material == Crista_Terminalis_2) ||
			(material == Pectinate_Muscle) || (material == Bachmann_Bundle_right_atrium_or_all) ||
			(material == Inferior_Isthmus_right_atrium) || (material == Right_Atrial_Appendage) ||
			(material == not_right_Isthmus) || (material == Tricuspid_Valve_Ring) || (material == Tricusline) ||
			(material == otherwise_right_or_all) || (material == Upper_Anterior_Bridge_right) ||
			(material == Middle_Posterior_Bridge_right) || (material == Upper_Posterior_Bridge_right) ||
			(material == Lower_Anterior_Bridge_right) || (material == Coronary_Sinus_Bridge_right) || (material == Coronary_Sinus_ostium_tissue) ||
			(material == insert_new_Material) || (material == subendo_right_Atrial_Cardium) || (material == Blood_right_atrium) ||
			(material == Pektinatmuskeln)) {
			return true;
		}
		return false;
	}

	/*! Test if a tissue class is in left atrial with the blood.
	\param material Tissue class that would be tested.
	\return True if the tissue is in left atrial with the blood, else False.
	*/
	static bool isInLeftWithBlood(int material) {
		if ((material == Vorhof_links) || (material == subepi_Atrial_Cardium) ||
			(material == Bachmann_Bundle_left_atrium) || (material == Left_Atrial_Appendage) || (material == RIPV) ||
			(material == RSPV) || (material == LIPV) || (material == LSPV) || (material == Mitral_Valve_Ring) ||
			(material == otherwise_left) || (material == Upper_Anterior_Bridge_left) ||
			(material == Middle_Posterior_Bridge_left) || (material == Upper_Posterior_Bridge_left) ||
			(material == Lower_Anterior_Bridge_left) || (material == Coronary_Sinus_Bridge_left) ||
			(material == Vorhof_links_Endo) || (material == insert_new_Material) || (material == subendo_left_Atrial_Cardium) ||
			(material == Blood_left_atrium)) {
			return true;
		}

		return false;
	}
}; // class Material

#endif /* defined(_RESILIENT_Material_) */

/*!
 \page Material

 \section DESCRIPTION_Material DESCRIPTION
 All tissue classes and functions for tissues.

 \section SOURCE_Material SOURCE

 Material.h

 \section CHANGELOG_Material CHANGELOG
 V1.0.0 - 01.09.2015 (Andreas Wachter): Starting with changelog\n
 */

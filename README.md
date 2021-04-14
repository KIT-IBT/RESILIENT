# RESILIENT
A RulE baSed atrIaL fIber gENeraTor.

## Overview
RESILIENT is a rule-based artial fiber generator for generating bi-atrial fibers and interatrial bridges. The algortihm was developed at IBT and published in Wachter et al. 2015 (see below). The algrothim works with triangle, terahedron and voxel meshes. The fibers are based on 22 seedpoints: 9 in the right atrium and 13 in the left atrium. They have to be selected by the user as described in Wachter et al.

![localization of the seedpoints](https://user-images.githubusercontent.com/70153727/114750043-d701d500-9d53-11eb-9d02-7608baddf1b3.jpg)

## Dependencies
To install and run RESILIENT, you need:
- VTK (https://gitlab.kitware.com/vtk/vtk)
- GMP (https://gmplib.org)
- Cork (https://github.com/gilbo/cork)
- Tetgen (https://www.tetgen.org)

Install these using your package manager (apt, homebrew, macports...) or download the source code directly. 
When compiling cork from source, make sure to give the correct paths to gmp in the "makeConstants" file.
The standard path where we look for cork is a subfolder "cork" here. For tetgen, it's "tetgen1.6.0".


## Install
We use CMake to configure the build process. If your dependencies are in non-standard locations, you might need to pass these paths to CMake.
```
git clone https://github.com/KIT-IBT/RESILIENT.git
cd RESILIENT
cmake .
make
sudo make install
```
## Run

### RESILIENT
```
RESILIENT: This programm setting semi-automatic the atrial fiber and interatrial bridge
<path for input datafile *.vtu | vtp | vtk>
<path for output datafile *.vtu | vtp | vtk>
<path for seedpoints file *.txt>

[-intermediateResultPrefix] (default: inputdataname): prefix for intermediate results
[-writeIntermediateData]: for saving intermediate Results
[-noOpenVTPBridgeInAtrial]: open atrial Bridges in VTP
[-noRightAtrium]
[-noBachmannBundle]
[-noWideBachmann]: if you wouldn't have wide Bachmann Bundle from the right to the left atrial appendeage
[-noUpperPosteriorBridge]
[-noMiddlePosteriorBridge]
[-noCoronarySinusBridge]
[-lowerAnteriorBridge]
[-upperAnteriorBridge]
[-writeViewFiberOrientation]: write out the orginal fiber orientation in the cell for illustration
[-numberOfPectinateMuscles N] (default: 15): number of pectinate muscles
[-growPathOffestLeft N] (default: 0): additiv Offset for the path dilatation in left Atrium
[-writeViewFiberOrientation]: growPath transmural in the atria
[-growPathTransmuralPectEndo]: grow endocardial Pectinate muscles in the atria
[-growPathPectFromEndo]: grow radial endocardial Pectinate muscles in the atria from the endocardial surface
[-growPectToRAppendage]: grow Pectinate muscles in the right atria appendeage
[-growCristaDepOnWallthick]: grow Crista Teriminals depending on avarage atrial wall thickness
[-markCoronary]: mark the coronary sinus
[-debug]: if you would have more information when program abort
[-freeBridge PATH]: if you would set free define bridges + definitionfilepath
[-noInitialClean]: no inital clean of the mesh e.g. when they are disconnected
[-noEndCleanup]: no end clean of the mesh e.g. when they are disconnected
[-debug]: if you would have more information when program abort


```
### TestTetrahedralize
```
TestTetrahedralize: This programm can used for testing if a bridge can tetralizied with tetgen.
<path for input datafile *.vtu | vtp | vtk>
```
### FindAndMarkSeedPoints
```
FindAndMarkSeedPoints: Mark the localization of the seedpoint in the mesh by using the material class.
<path for input datafile *.vtu | vtp | vtk>
<path for output datafile *.vtu|vtp|vtk>
<path for seedpoints file *.txt>
```

## Example



## Debug

## Citation
When using this program, please cite
Andreas Wachter, Axel Loewe, Martin W. Krueger, Olaf Dössel and Gunnar Seemann. "Mesh structure-independent modeling of patient-specific atrial fiber orientation". Current Directions in Biomedical Engineering, 2015, 1(1), 409-412.](https://doi.org/10.1515/cdbme-2015-0099)  

### Bibtex
```
@article{Wachter_2015
author = {Andreas Wachter and Axel Loewe and Martin W Krueger and Olaf Dössel and Gunnar Seemann},
doi = {doi:10.1515/cdbme-2015-0099},
url = {https://doi.org/10.1515/cdbme-2015-0099},
title = {Mesh structure-independent modeling of patient-specific atrial fiber orientation},
journal = {Current Directions in Biomedical Engineering},
number = {1},
volume = {1},
year = {2015},
pages = {409--412}
}
```

## Contact 
Andreas Wachter, Luca Azzolin, Axel Loewe
Institute of Biomedical Engineering
Karlsruhe Institute of Technology (KIT)
www.ibt.kit.edu
publications@ibt.kit.edu
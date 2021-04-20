# RESILIENT
RulE baSed atrIaL fIber gENeraTor.

## Overview
RESILIENT is a rule-based artial fiber generator for generating bi-atrial fibers and interatrial bridges. The algortihm was developed at IBT and published in Wachter et al. 2015 [[1](#1)]. The algrothim works with triangle, terahedron and voxel meshes. The fibers are based on 22 seedpoints: 9 in the right atrium and 13 in the left atrium. They have to be selected by the user as described in [[1](#1)].

## Dependencies
To install and run RESILIENT, you need:
- VTK (https://gitlab.kitware.com/vtk/vtk)
- GMP (https://gmplib.org)
- Cork (https://github.com/gilbo/cork)
- Tetgen (http://www.tetgen.org)

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


## Preprocessing steps

### Volume mesh (VTK/VTU) pre-processing
The atrial mesh requires special endo/epi materials and should be divided into left and right atrium for the algorithm.
The setpum can be part of the left atrium only or split into left/right. Both options should work.

```
Materials
 Material 32    for the right epicardium
 Material 33    for the left epicardium
 Material 129   for the right endocardium
 Material 233   for the left endocardium
```

### Surface mesh (VTK/VTP) pre-processing
The atrial mesh requires special epi materials and should be divided into left and right atrium for the algorithm. The artia should be closed e.g. the valves should be filled with blood (material: 150 right/ 151 left). 
The setpum can be part of the left atrium only or split into left/right. Both options should work.

```
Materials
 Material 32    for the right epicardium
 Material 33    for the left epicardium
```

### SeedPoints
The seed point are subdivides in 9 Points for the right atrium and 13 Points for the left atrium. 

![localization of the seedpoints](https://user-images.githubusercontent.com/70153727/114750043-d701d500-9d53-11eb-9d02-7608baddf1b3.jpg)

This point must manually annotated in every mesh by the user. It is important that the point do not lie in a extrema like a maximum or minimum on the surface.

The first three points on the right atrium describe the orifice of the superior cava vein. "R1" is arranged on the backward side of the right atrium at the borderline between the right and the left atrium. "R2" is arranged on the front side of the right atrium at the borderline between the right and the left atrium. "R3" is arranged on the back side of the right atrium at the beginning of the right appendage.  With these three points it should be possible to construct a ring. That means that the shortest way in the material between the points "R2" and "R3“ don’t have to cross the shortest way between "R1" and "R2" or between "R3" and "R2". 
The points "R4" and "R5" describe the orifice of the inferior cava vein. "R4" is arranged on the right side of the orifice and "R5" on the left side on the front side of the right atrium. The shortest way between these points have to be on the front side of the orifice. The point "R6" marks the peak of the right appendage over the opening of the right appendage. The last three points describe the ring of the tricuspidal valve. Point "R7" is arranged below of the right appendage. "R9" is the nearest point of the orifice to the left atrium. "R8" is arranged on the other side of the orifice. It is important the shortest way in material between the points "R8" and "R9" does not cross the shortest ways between "R7" and "R8" or between "R7" and "R9".
The first three points on the left atrium describe the ring of the orifice of the mitral valve. "L1" is arranged on the left lateral side below the left appendage. "L3" is on the right lateral side of the left atrium, ca. the nearest point of the orifice to the right atrium. It is important that the shortest way in material between these two points is in the front wall of the left atrium. "L2" is arranged ca. in the middle of the of the backward side of the orifice. Point "L4" is on the top of the front side of the left atrium, it is arranged at the beginning of the right superior pulmonal vein. "L5" is also arranged on the top of the front wall at the beginning of the left superior pulmonal vein. Normally it is between the left appendage and the left superior pulmonal vein. The next two points are on the top wall of the left atrium. "L6" marks a point between the right pulmonal veins and "L7" marks a point between both left pulmonal veins. The points "L8" and "L9" are located on the top of the posterior wall of the left atrium. "L8" marks the beginning of the right and "L9" the beginning of the left inferior pulmonal vein. The points "L10" and "L11" describe the left appendage. "L11" marks the peak of the left appendage over the opening of the right appendage and "L10" is arranged on the front wall of the left atrium and describe the excrescence of the left appendage on the front wall. The point "L12" is arranged in the left lateral wall over the left appendage between both left pulmonal veins, nearly direct below of the point "L7". It is important that the points "L5", "L7" and "L12" construct a ring of the orifice of the left superior pulmonal vein. The point "L13" is arranged in the right lateral wall of the left atrium below the right inferior pulmonal vein. It is important that the points "L6", "L8" and "L13" construct a ring of the orifice of the right inferior pulmonal vein.

For each three points on a common opening, make sure that the shortest path between two of those points is the path that does not include the third point. For example: L1L2L3: the path from L1 to L2 in clockwise direction (on the image above) is shorter than the path between L2 and L1 also in clockwise direction, since L3 is in between.   

The position of the seedpoints taken by the algorithm can be checked in the mesh with the tool: *FindAndMarkSeedPoints*.

#### Syntax for the seedpoints.txt
The points must be save in an TXT-file as the following syntax. 

```
 SCV1 G 12.2 -23.4 57.5
 SCV2 G 13 -25 79
 ...
 ...
 SCV9 G 48 27 -30
 LV1 G 59 -45 -98
 LV2 G 97  -45 -108
 ...
 ...
 LV13 G -105 23 -39
```
### Free atrial bridges
If you want an addititonal bridges between the atria you can define as follow in an TXT-file.

#### Syntax Freebridge
Point left Atria (XYZ) point right Atria (XYZ) material bridge radius
```
xl yl zl xr yr zr Material Bridgeradius
```

## Run
### RESILIENT
```
RESILIENT rParameter1 rParameter2 rParameter3 oParameter1 oParameter2 [...]
```
If you need help with parameter selection just run RESILENT in the command line. Then the options will pop up.

#### Requiered parameters:
```
<path for input datafile *.vtu | vtp | vtk>
<path for output datafile *.vtu | vtp | vtk>
<path for seedpoints file *.txt>
```
#### Optional parameters:
```
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

## Substituion of tissue class
Some tissue classes are substituted at the end of the calculation for the simulation.

|        Region        |  new created tissue class  |  for Simulation relevant tissue class  |
|:--------------------:|:--------------------------:|:--------------------------------------:|
|          SVC         |             159            |                   32                   |
|  IntercavanalBundel  |             175            |                   32                   |
|   not right Isthmus  |             131            |                   32                   |
|      other right     |             187            |                   32                   |
|      Tricusline      |             176            |                   32                   |
|   left atrium endo   |             233            |                   33                   |
|      other left      |             188            |                   33                   |
|        subendo       |             130            |                   33                   |
|        subepi        |             99             |                   33                   |

## Example
### Example for a vtp/vtu/vtk mesh:
If you don't know if the point are good for the mesh it makes sense to start RESILENT with the optional Parameter:
```
RESILENT inputMesh.vtp/vtu/vtk inputMesh_wafo.vtp/vtu/vtk seedPoint.txt -writeIntermediateData -debug 
```

### Surface mesh
```
RESILIENT ./example/surface_mesh/mesh_material.vtp ./example/surface_mesh/mesh_with_fiber.vtp ./example/surface_mesh/seedPoint.txt -noInitialClean
```

### Volume mesh
```
RESILIENT ./example/volume_mesh/mesh_material.vtu ./example/volume_mesh/mesh_with_fiber.vtu ./example/volume_mesh/seedPoint.txt -noInitialClean
```

## Debug and problems
### Fiber calculation
If the fiber could not be calculated, run RESILIENT again and use the "writeIntermediateData" flag and "debug" to get additional information and find the error. In most cases, one or two points are misplaced, so the base path has an incorrect path. Find them and replace them to fix the problem. For this, there is an array "basePath" in the network file where the base paths are stored. Compare them with the following figures.

#### Basic paths
The basic paths that would be calculated in the process of RESILIENT. This can help you if you to find the possible error in the fiber orientation calculation.

The basic paths in the right atrial, the left epicard and the automatic bridges:
![BasedPathEpicard](https://user-images.githubusercontent.com/70153727/115027351-e4d66800-9ec3-11eb-908f-52d2b552b1c1.jpg)

The basic paths in the right atrial and the left endocard:
![BasedPathEndo](https://user-images.githubusercontent.com/70153727/115027335-de47f080-9ec3-11eb-81cc-8f49b358f705.jpg)

Another problem can be that the dilated paths do not isolate a region. As a result, areas of the mesh may be annotated. In most cases you can fix the problem with replacement of the points or with the additional dilatation flag. 

### Bridge insertion
The bridges in the tetrahedron mesh sometimes can't be inserted. So the user has to edit the bridge manually. To do this, it is best to use the programs *Blender* and *Meshlab* in sequence. Mostly self-cutting and multiple edges are the problem, fix that. When you are done, you can test it with the tool *TestTetrahedralize*.

## Citation
When using this program, please cite

<a id="1">[1]</a> [Andreas Wachter, Axel Loewe, Martin W. Krueger, Olaf Dössel and Gunnar Seemann. "Mesh structure-independent modeling of patient-specific atrial fiber orientation". Current Directions in Biomedical Engineering, 2015, 1(1), 409-412.](https://doi.org/10.1515/cdbme-2015-0099)  

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

Institute of Biomedical Engineering, Karlsruhe Institute of Technology (KIT)

www.ibt.kit.edu 

<publications@ibt.kit.edu>

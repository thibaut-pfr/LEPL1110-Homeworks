# LEPL1110 - Finite Element Method

This repository contains the solutions to the homework problems of the 3rd year Finite Element Method course, as taught at UCLouvain by Vincent Legat in 2024-2025.

### Dependencies

* Cmake installed
* GMSH SDK installed for your machine (https://gmsh.info/)

### Installing

* You will need to add a folder called "build" for each homework that will contain the executables. Do this at the same level as the "src" folders.
* For each homework using gmsh, also create a folder called "gmsh", in which you put the download folder of the sdk.
* Verify that the version and minor of gmsh as described at line 7-8 of the $`\texttt{CMakeLists.txt}`$ is the same as the one you downloaded.

### Executing program

* To execute each homework, run the following commands from said "build" folder: 
```
$ cmake ..
$ make
$ ./myFem
```
## Authors

* Thibaut Peiffer
* Mathieu Mil-Homens Cavaco

# MTGoalAdaptiveFEM
## About
--------

MTGoalAdaptiveFEM is a parallel C++ program targeted at the forward modeling of 2D frequency-domain MT electromagnetic method using finite element method. The main features of MTGoalAdaptiveFEM including:
- Modular design
- Support Arbitrary anisotropic media
- Support the calculation using high-order shape functions

## Prerequisites
----------------

MTGoalAdaptiveFEM uses the open source library [dealii.II](https://www.dealii.org/), which is a C++ software library supporting the creation of finite element codes and an open community of users and developers.

The easiest way to install the dependencies of MTGoalAdaptiveFEM is using spack.

## Building
-----------

Once installed all the dependencies and unpacked the source code of MTGoalAdaptiveFEM into a directory `/path/to/mt2d`. Then configure, compile MTGoalAdaptiveFEM with
```shell
$ cd /path/to/mt2d
$ mkdir build
$ cd build
$ spack load deal.ii
$ cmake -DCMAKE_BUILD_TYPE=Release ..
$ make
```
Note that MTGoalAdaptiveFEM has only been tested on Linux and macOS systems. For Windows users, we recommend using WSL1/WSL2 or a virtual machine.

## Usage
--------

To run the forward modeling process, the user needs to provide three files: the model file, the receiver file and the configuration file. Then use the flowing command to run MTGoalAdaptiveFEM:
```shell
$ ./mt2d path/to/modelname
```

## Configuration File
---------------------

Before running MTGoalAdaptiveFEM, we need create some configuration files. The configuration files contain the model file, observation file and option file. Lines prefixed by the character # are considered as comments. 

### Model File

Creating a model is simple.You only need to give the approprate size of each element in the y and z directions and the coordinates of the two diagonal vertices. Set different materia_id for elements in different regions. The resistivity of each region is relative to the sequence of materia_id. We provide a script for the model creation. Modify it according to your requirements. Compiling the script, 
```shell
$ python ./scripts/mkdmdl.py path/to/modelname
```
you will get a model file with extension name `.vtk`. The mesh file can be visualized by [Paraview](https://www.paraview.org/) or [VisIt](https://wci.llnl.gov/simulation/computer-codes/visit).
### Part 2 - Transmitters
This block lists the number of transmitters and gives the parameters.
```
One line: <# of transmitters>
Following lines list # of transmitters:
  <x> <y> <z> <azimuth> <dip> <current> <dipole length> # <x> <y> <z> is the center of transmitter
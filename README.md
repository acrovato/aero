# Aero
Adrien Crovato  
ULiege, 2016-2017

## Features and limitations

### What is Aero?  
Aero was developed at the University of Liège during the academic year 2016-2017. 
The purpose of Aero was to compute the aerodynamic loads on an arbitrary wing for preliminary aircraft design computations.  

To achieve this goal, a Field Panel Method (FPM) solving the Full Potential Equation was developed.  
The FPM is an extension of the panel method. The body is discretized in panels onto which singularities are superimposed and a Cartesian grid is added around the body. The Cartesian grid also contain singularities (field sources) modelling the compressibility of the fluid.  

Additional theory as well as references can be found in the documentation under the [doc](doc/) folder.

### Is Aero finished?
Not really.  
The Field Panel Method works well for incompressible or subcritical flows. However as soon as shocks start to appear, the accuracy of the solution quickly drops out. Moreover, the FPM scales as N^2, where N is the number of field cells. The Fast Multipole Method could be implemented to reduce the computational cost.
  
## Compilation

Links to packages  

  - [CMake](https://cmake.org "https://cmake.org")  
  - [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page "http://eigen.tuxfamily.org/index.php?title=Main_Page ")  
  - [gmsh](http://geuz.org/gmsh/ "http://geuz.org/gmsh/")

### Linux (gcc)

Required packages
```bash
sudo apt-get update
sudo apt-get install build-essential
sudo apt-get install cmake
sudo apt-get install libeigen3-dev
sudo apt-get install gmsh
```
Compilation
```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j 4
```
Test
```bash
./aero ../IO/N12.cfg ../IO/N12.pts
```

### OSX (Clang)

Required packages  
  - CMake  
  - Eigen3  
  - gmsh  
Update Eigen path in CMakeLists.txt accordingly  
Compilation
```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j 4
```
Test
```bash
./aero ../IO/N12.cfg ../IO/N12.pts
```

### Windows (minGW, not tested)
Required packages  
  - CMake  
  - Eigen3  
  - gmsh  
Update Eigen path in CMakeLists.txt accordingly

## Usage

### Wing meshing

The mesh file (second input file, see below) is ASCII formatted and can be written by hand. If you have access to [MATLAB](https://www.mathworks.com "https://www.mathworks.com"), a meshing script is given to facilitate the generation of the mesh file.  
In order to define a wing, the user must change the following parameters:  

Files (ASCII Selig formatted) where the airfoils are defined:  

```MATLAB
% Folders and paths
nAirf = 3; % number of airfoil used to generate the wing
fname = cell(1,nAirf); % name of the airfoils data file
fpath = ''; % folder where airfoils data are
fname{1,1}(1,:) = 'nasaSC20714';
fname{1,2}(1,:) = 'nasaSC20712';
fname{1,3}(1,:) = 'nasaSC20712';
outFname = 'mesh.dat'; % output file name
```  

Geometry of the wing:  

```MATLAB
% Wing geometry (multi planform defintion for HALF wing)
meanC4AirfSweep = 0*pi/180; % mean quarter-chord sweep (0 if airfoil coordinates are provided STREAMWISE, mean sweep angle if CHORDWISE)
span = [2; 4]; % panform span
taper = [0.6; 0.4]; % planform taper
sweepLE = [25; 25]*pi/180; % leading edge planform sweep
% sweepLE = atan(tan(sweepC4)+(1-taper)./((AR)*(1+taper))); %if quarter chord sweep provided instead
dihedral = [5; 2]*pi/180; % planform dihedral
twist = [0.; 0.; -1.]*pi/180; % planform twist
rootChord = 2; % root chord length
```  

Mesh parameters:  

```MATLAB
% Mesh parameters
meshType = 1; % mesh distribution type (1=half-cosine, 2=linear, 3=geometric)
prog = 1.1; % if mesh type is geometric, defines progression
nSpanStat = [5; 11]; % number of spanwise stations (spanwise panels - 1) for EACH planform
nPoints = 101; % number of chordwise points (chordwise panels - 1), MUST BE ODD AND DIVISIBLE BY 4 when 1 is substracted
```

### Input files

Aero needs two input files: a configuration .cfg file and a (surface) grid .pts file.  

  - The cfg file is order in sections (marked by `$`). The order of the lines and sections should not be changed.  Below is a commented example file. `b` stands for boolean, `i` stands for int and `d` stands for double.  

    - `Config file` Title line
    - `$sym` Symmetry section
    - `XZ_SYM = b` Symmetry about xz plane (0=no, 1=yes) 
    - `$num` Numerical parameters section
    - `boxTol = d` Tolerance on the minimal size box enclosing the geometry (default: 1e-3)
    - `resRed = i` Tolerance on the relative source change (stop criterion) (default: 5)
    - `$geo` Geometry parameters
    - `S_REF = d` Reference area (projected on the z-plane)
    - `$flow` Flow parameters
    - `MACH = d` Mach number
    - `AoA = d` Angle of attack
    - `$grid` Cartesian grid parameters
    - `X_DIV = i` Number of cell in the x direction
    - `Y_DIV = i` Number of cell in the y direction
    - `Z_DIV = i` Number of cell in the z direction
    - `$box` Cartesian grid corner definition (reference looking in the +y direction)
    - `BOX_X0Y0Z0 = d d d` Coordinates of front lower left corner
    - `BOX_X1Y0Z0 = d d d` Coordinates of front lower right corner
    - `BOX_X0Y0Z1 = d d d` Coordinates of front upper left corner
    - `BOX_X1Y0Z1 = d d d` Coordinates of front upper right corner
    - `BOX_X0Y1Z0 = d d d` Coordinates of back lower left corner
    - `BOX_X1Y1Z0 = d d d` Coordinates of back lower right corner
    - `BOX_X0Y1Z1 = d d d` Coordinates of back upper left corner
    - `BOX_X1Y1Z1 = d d d` Coordinates of back upper right corner
    - `$sub` Sub-panelling parameters
    - `NS = i` Number of sub-panel per panel (16)
    - `NC = d` Cut-off percrentage of the distance cell-panel in the normal direction (default: 0.1)
    - `LC = d` Cut-off percrentage of the distance cell-panel in the longitudinal direction (default: 0.75)

	 
  - The pts file is ordered in two sections and contains the corner points of the panels on the wing.  
    The first line is the title. Then, the first section `$size` contains the number of chordwise and spanwise points. The second section `$points` contains the panel corner points. They are ordered in Selig format in the chordwise direction: from suction side TE to pressure side TE. When an airfoil (section of the wing) has been fully defined, the next section in the spanwise direction is described.  
    All the points are given in terms of x (flow), y (span), z (vertical) coordinates.

Two test cases are given as examples: a rectangular naca0012 wing and the weel-known Onera M6 wing.

### Post-processing and output files

Aero writes four files after each run.  

  - sp.dat: ASCII formatted file containing the pressure coefficient at the center of each panel  
  - fv.dat: ASCII formatted file containing the field variables at the center of each field cell  
  - Cp.pos: ASCII formatted gmsh file displaying the pressure coefficient on the wing  
  - M.pos: ASCII formatted gmsh file displaying the Mach number in the field

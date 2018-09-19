# Aero
Adrien Crovato  
ULiege, 2016-2017

## Features and limitations

### What is Aero?  
Aero was developed at the University of Liège during the academic year 2016-2017. 
The purpose of Aero was to compute the aerodynamic loads on an arbitrary wing for preliminary aircraft design computations.  

To achieve this goal, a Field Panel Method (FPM) solving the Full Potential Equation was developed.  
The FPM is an extension of the panel method. The body is discretized in panels onto which singularities are superimposed and a Cartesian grid is added around the body. The Cartesian grid also contain singularities (field sources) modeling the compressibility of the fluid.  

Additional theory as well as references can be found in the documentation under the [doc](doc/) folder.

### Should I use Aero?
Several cases should be considered:  

  - [x] Subcritical flow (M_local < 0.95)  
    The FPM will perform well, but the gain in accuracy will not be relevant compared to the increase in CPU cost when compared to the Panel Method (PM). A Fast Multipole Method (FMM) has to be implemented to decrease the computational time and make it close to a PM. Such methods have already been applied to the PM, so the extension should not be a problem.  
  - [x] Supercritical flow - really weak shock ( 0.95 < M_local < 1.15)  
    The FPM will give a better solution than the PM. The cost will not be that high. However, to remain competitive, an acceleration technique such as the Fast Multipole Method should be considered.  
  - [ ] Supercritical flow - weak shock (1.15 < M_local < 1.3)  
    The FPM results will be better than the PM, but the shock will be smeared and displaced compared to traditional field methods (Finite Elements/Volumes). Even if the FMM is implemented, the lack in shock resolution would still pose a serious problem for transonic aircraft design.  
  - [ ] Supercrtical flow - strong shock (M_local > 1.3)  
    The Full Potential Equation should not be used!  

Test cases illustrating thoses results can be found in the documentation under the [doc](doc/) folder.

### Cite us!
If you use this work, please acknowledge the authors:  
```text 
"Aero: a Field Panel Method for Aerodynamic Loads Computation in Preliminary Aircraft Design"
```  

The main paper about this work can be found in the 2018 [ICAS proceedings](https://www.icas.org/Papers_previous_congresses.html "https://www.icas.org/Papers_previous_congresses.html") or on [Orbi](http://hdl.handle.net/2268/227902 "http://hdl.handle.net/2268/227902"):
```text
A. Crovato, G. Dimitriadi, V.E. Terrapon, University of Liège. *Higher Fidelity Transonic Aerodynamic Modeling in Preliminary Aircraft Design*, 31st Congress of the International Council of the Aeronautical Sciences (ICAS), 9-14th September 2018, Belo Horionte, Brazil.
```
  
## Compilation
Aero is really simple to compile. It just needs to be linked to the pre-compiled (header based) Eigen library. All the configuration is handled by CMake, so you only need to adapt some path inside the CMakeLists.txt. Additionnaly, you may want to download gmsh to view some results during the post-processing.  

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
Path setting (in CMakeLists.txt)  
```cmake
IF (UNIX AND NOT APPLE)
include_directories(include; /usr/include/eigen3)
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
  - CMake, Eigen3 and gmsh  
  - Clang (through Command Line Tools)
```bash
xcode-select --install
```
Path setting (in CMakeLists.txt)  
```cmake
IF (APPLE)
include_directories(include; /path/to/eigen)
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

### Windows (MinGW)
Required packages  
  - CMake, Eigen3 and gmsh  
  - MinGW  
      - download from https://sourceforge.net/projects/mingw/files/ (mingw-get)  
      - install the manager  
      - install `mingw32-base` and `mingw32-gcc-g++`  
      - add `\path\to\MinGW\mingw32\bin` to the `PATH` Environment Variable in Windows  
Path setting (in CMakeLists.txt)  
```cmake
ELSEIF (WIN32)
include_directories(include; /path/to/eigen)
```
Compilation
```posh
md build
cd build
cmake -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release ..
mingw32-make -j 4
```
Test
```posh
.\aero.exe ..\IO\N12.cfg ..\IO\N12.pts
```

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

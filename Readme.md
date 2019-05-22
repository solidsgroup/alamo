# Getting Started #

**Note**: this README page is also the Doxygen main page, the Github readme page, and the Docs main page.
You can view it by running `make docs` in the root directory, then opening `docs/doxygen/html/index.html` or `docs/build/html/index.html` in a web browser. 

## Compiling Alamo ##

This section describes how to compile and install Alamo and its dependencies.

### Dependencies ###

Alamo requires `eigen3` and `AMReX` to compile. 
The following are instructions for how to install these dependencies on your platform.

* The website for `eigen3` is http://eigen.tuxfamily.org. 
  Download the source and store it in a directory (e.g. /home/myusername/eigen3/). 
  This is all you need to do.
  (The directory must be named eigen3)
* The website for `amrex` is https://github.com/AMReX-Codes/amrex
  To compile and install, clone the AMReX repository with 
  
      git clone https://github.com/AMReX-Codes/amrex.git
    
  Change into the AMReX directory with 
  
      cd amrex
  
  Configure AMReX by typing
  
      ./configure --dim=3 --debug=yes --prefix=/home/myusername/amrex/
  
  Note that you can compile in 3D with `--dim=2`, you can compile in non-debug mode with `--debug=no`, 
  and you can change the installation directory by changing the path in `--prefix`.
  Next, compile AMReX by typing
  
      make
      
  This may take some time. You can speed it up by compiling in parallel with `make -j2`(where 2 == two processors).
  Finally, install by typing
  
      make install
  

### Building - Using Makefile ###

To build alamo, in the alamo directory type

    ./configure

To see a full list of options type 

    ./configure --help

To explicitly specify AMReX and Eigen locations (which you probably should), type 

    ./configure --amrex /path/to/amrex --eigen /path/to/eigen

where `/path/to/amrex` and `/path/to/eigen` should match the `--prefix` and the location of the Eigen directory, respectively.
(The script will check this and error out if they do not.)
To specify the spatial dimension,

    ./configure --dim #OfDimensions

The default is to compile in 3 dimensions.
If you specify AMReX directly, it will check to make sure that AMReX is also compiled in the same number of dimensions.
Immediately after running the `./configure` script, type

    make

This again can be sped up using `make -j#OfProcessors`.
For additional help, type 

    make help

to get an extensive help message.
Finally, type 

    make clean

to clear out extra files that were generated.

### Building - Using CMAKE (Not up to Date)###

To build ALAMO using CMAKE, you must have CMAKE installed. 
1. Create a build directory.
   **Do not use the ALAMO root directory as your build directory**.
   The `.gitignore` file is configured to ignore the `./build` directory, but you can choose whatever directory you wish.
2. Within the build directory, run CMAKE with the following:

       cmake /path/to/alamo/root -DAMREX=/path/to/amrex
	   
   Note that `/path/to/amrex` must contain two directories: `lib` and `include`
3. Within the build directory, type

       make
		
   to generate executables.

## Testing ##

Upon successful compilation, run tests by 

    ./bin/test-3d-debug-g++

The output will indicate wether the tests pass or fail.
If you are committing changes, you should always make sure the tests pass in 2 and 3 dimensions before committing.

## Common Error Messages ##

The following are some common error messages and problems encountered.

* `MLLinOp: grids not coarsenable between AMR levels`:
  This is a conflict in the **multigrid solver** because the grid size is not a power of 2.
  Solve by changing the domain dimensions (`amr.n_cell`) so that they are powers of two.

* `static_cast<long>(i) < this->size() failed`:
  One common reason this happens is if Dirichlet/Neumann boundaries are specified but no boundary values are provided.

## Generating Documentation ##

Generating documentation requires the following packages:

* Doxygen (on Ubuntu: `sudo apt install doxygen`)
* Sphinx (on Ubuntu: `sudo apt install python-sphinx`)
* Breathe (on Ubuntu: `sudo apt install python-breathe`)
* M2R (on Ubuntu: `pip install m2r`)
* RTD theme (on Ubuntu: `pip install sphinx_rtd_theme`)
* GraphViz (on Ubuntu: `sudo apt install graphviz`)

To generate the documentation, type

    make docs

(You do not need to run `./configure` before generating documentation.)
Documentation will be generated in `docs/build/html` and can be viewed using a browser.

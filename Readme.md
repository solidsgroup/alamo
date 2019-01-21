## ALAMO - Getting Started ##

**Note**: this README page is also the Doxygen main page and the Github readme page.
You can view it by running `doxygen` in the root directory, then opening `doc/html/index.html` in a web browser.

## Compiling Alamo ##

This section describes how to compile and install Alamo and its dependencies.

### Dependencies ###

Alamo requires `eigen3` and `AMReX` to compile. 
The following are instructions for how to install these dependencies on your platform.

* The website for `eigen3` is http://eigen.tuxfamily.org. 
  Download the source and store it in a directory (e.g. /home/myusername/eigen3/). 
  This is all you need to do. 
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

    make AMREX=/path/to/amrex EIGEN=/path/to/eigen

where `/path/to/amrex` and `/path/to/eigen` should match the `--prefix` and the location of the Eigen directory, respectively.
To generate documentation, type

    make docs
    
(Note that this requires the Sphinx and Doxygen libraries.)
For additional help, type 

    make help

to get an extensive help message.
Finally, type 

    make clean

to clear out extra files that were generated.

### Building - Using CMAKE ###

To build ALAMO using CMAKE, you must have CMAKE installed. 
1. Create a build directory. **Do not use the ALAMO root directory as your build directory**. The `.gitignore` file is configured to ignore the `./build` directory, but you can choose whatever directory you wish.
2. Within the build directory, run CMAKE with the following:

       cmake /path/to/alamo/root -DAMREX=/path/to/amrex
	   
   Note that `/path/to/amrex` must contain two directories: `lib` and `include`
3. Within the build directory, type

       make
		
   to generate executables

## Testing ##

Upon successful compilation, run tests by 

    ./bin/test

## Common Error Messages ##

The following are some common error messages and problems encountered.

* `MLLinOp: grids not coarsenable between AMR levels`:
  This is a conflict in the **multigrid solver** because the grid size is not a power of 2.
  Solve by changing the domain dimensions (`amr.n_cell`) so that they are powers of two.

* `static_cast<long>(i) < this->size() failed'
  One common reason this happens is if Dirichlet/Neumann boundaries are specified but
  no boundary values are provided. 

## ALAMO - Getting Started ##

**Note**: this README page is also the Doxygen main page and the Github readme page.
You can view it by running `doxygen` in the root directory, then opening `doc/html/index.html` in a web browser.

### Building - Using CMAKE ###

To build ALAMO using CMAKE, you must have CMAKE installed. 
1. Create a build directory. **Do not use the ALAMO root directory as your build directory**. The `.gitignore` file is configured to ignore the `./build` directory, but you can choose whatever directory you wish.
2. Within the build directory, run CMAKE with the following:

       cmake /path/to/alamo/root -DAMREX=/path/to/amrex
	   
   Note that `/path/to/amrex` must contain two directories: `lib` and `include`
3. Within the build directory, type

       make
		
   to generate executables

### Building - Using Makefile ###

Unlike with CMAKE, you must have MPI installed already and AMReX already in your path. 
Once you have done this, you can build by typing

    make
	
in the root directory.

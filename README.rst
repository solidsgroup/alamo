.. getting-started:

Note: this README page is also the Doxygen main page, the Github readme page, 
and the Docs main page.
You can view it by running :code:`make docs` in the root directory, then opening 
:code:`docs/doxygen/html/index.html` or :code:`docs/build/html/index.html` in a web browser. 



Downloading Alamo
=================

Download alamo using git:

.. code-block::

    git clone git@github.com:solidsgroup/alamo.git
    
If you do not have a Github account and/or you have not uploaded your public SSH key, this will probably throw an error.
You can download alamo using HTTPS instead,

.. code-block::
    
    https://github.com/solidsuccs/alamo.git 

Note, however, that you will not be able to push anything using HTTPS authentication.
The :code:`master` branch is the most stable and is what is checked out by default.
The :code:`develompent` branch is generally stable, and includes the latest functionality.
To switch to :code:`development`, in the alamo directory,

.. code-block::
    
    git checkout development



.. literalinclude:: .git/workflows/dependencies.sh
   :caption: Dependencies
   :language: bash
    
Installing dependencies
=======================

Alamo depends primarily on AMReX, but has some of its own dependencies, such as Eigen.
The two officially supported Linux distributions are currently: the latest version of Ubuntu, and Ubuntu 20.04. (This is subject to change).
If you are using Ubuntu (or a debian-based distribution) you can install all necessary dependencies by running the dependencies script:

.. code-block::

    sudo bash .github/workflows/dependencies.sh

If you are using a non-debian based system such as RedHat, you can install the corresponding packages in :code:`dependencies.sh`. 
Windows and Mac OS are not supported.
To configure and run Alamo on a high-performance computing (HPC) cluster, please see :ref:`install_hpc`.

Setting the default MPI
=======================

Alamo can be compiled with either MPICH or OpenMPI.
(mvapich2 is semi-supported but not regularly tested.)
**MPICH currently does not work on Ubuntu 24.04**, so OpenMPI is now recommended.
On Ubuntu, you can check to see which version of mpi is being used by the :code:`update-alternatives` command.

.. code-block::

    $> sudo update-alternatives --config mpi

    There are 2 choices for the alternative mpi (providing /usr/bin/mpicc).
    
      Selection    Path                    Priority   Status
    ------------------------------------------------------------
      0            /usr/bin/mpicc.openmpi   50        auto mode
    * 1            /usr/bin/mpicc.mpich     40        manual mode
      2            /usr/bin/mpicc.openmpi   50        manual mode
    
    Press <enter> to keep the current choice[*], or type selection number:     

In this case, mpich is installed along with openmpi, but mpich is currently the default.
You can press 0 to switch to openmpi.
Do the same thing for :code:`mpirun`:

.. code-block::

    $> sudo update-alternatives --config mpirun
    
Now your system should be properly configured.

If you are using an HPC, it is easier to switch between versions using the :code:`module` command.
See the documentation for your specific platform to see how to load mpich or mvapich.    

Configuring
===========

To compile alamo, you must first run the configure script. 
This is done simply by running the following in the alamo directory 
(note that AMReX download is triggered by this command, so it may take a couple minutes to complete depending on your internet connection)

.. code-block::

    ./configure

By default, alamo will configure in 3D production mode. 
To compile in 2D debug mode, 

.. code-block::

    ./configure --dim=2 --debug

Additionally, the configuration step is when you specify which compiler will be used for compilation. By default, the configure script will use the GNU C++ Compiler (g++). To specify a different compiler, such as the Clang C++ Compiler, use the following command line argument, replacing clang++ with the supported compiler you'd like to use:

.. code-block::

   ./configure --comp clang++

There are many compilation options available for Alamo, and they must all be specified at configure time.
For a complete listing of the Alamo configuration options, type

.. code-block::

    ./configure --help

.. NOTE:: 
    The configure script produces output designed to assist in determining compile issues with Alamo.
    Whenever you request help with alamo, please always include the complete output of the configure script.

Compiling
=========

Once you have configured Alamo, compile it by

.. code-block::

    make

If you are on a platform with multiple cores, you can compile in parallel (for instance, with 4 cores) with

.. code-block::

    make -j4

The alamo exectuable will be stored in :code:`./bin/` and name according to the options specified at configure time.
For instance, if you are using GCC to make Alamo in 2D using debug mode, the alamo executable will be :code:`./bin/alamo-2d-debug-g++`.
You can work with multiple versions of Alamo at the same time without having to re-compile the entire code base.
All you need to do is re-run the configure script, and previous versions of Alamo and AMReX will be saved automatically.

.. WARNING::
    There is an issue with GNU Make that can cause I/O errors during parallel builds.
    You may get the following error:

    .. code-block::

        make[1]: write error: stdout

    To continue the build, just issue the :code:`make` command again and it should continue normally.
    You can also add the :code:`--output-sync=target` option which may help eliminate the issue.

Testing
=======

Upon successful compilation, run tests by

.. code-block::

    make test

This will run the unit tests and regression tests for all compiled production versions of Alamo.
If you have only run in 2D, only 2D tests will be generated.
If you are a developer and you are preparing to merge your branch into :code:`development`, you should perform a complete test via

.. code-block::

    ./configure --dim=2
    make
    ./configure --dim=3
    make
    make test

For a full description of the Alamo regression test system, please see 


Common Error Messages
=====================

The following are some common error messages and problems encountered:

* :code:`MLLinOp: grids not coarsenable between AMR levels`
  This is a conflict in the **multigrid solver** because the grid size is not a power of 2.
  Solve by changing the domain dimensions (`amr.n_cell`) so that they are powers of two.

* :code:`static_cast<long>(i) < this->size() failed`
  One common reason this happens is if Dirichlet/Neumann
  boundaries are specified but no boundary values are provided.

* :code:`error: lvalue required as left operand of assignment`
  This can happen when using the :code:`()` operator with a :code:`Isotropic` :code:`Matrix4`-type object.
  Because this data structure only stores two constants, it is not possible to define any of the values using
  indices. 
  (Similarly, you cannot set an :code:`Isotropic` 4-matrix to a :code:`Cubic` 4-matrix since the Cubic
  matrix has lower symmetry).
  If you get this error, you should use a lower-symmetry 4-matrix.

* :code:`Inconsistent box arrays`
  This is known to happen when using an :code:`Operator::Elastic` inside an :code:`Integrator`, e.g. in :code:`TimeStepBegin`.
  Typically this happens when the Elastic operator is not initialized within the routine in which it is used - i.e.e if it is declared as a member variable inside the :code:`Integrator` - derived class.
  (The reason is that there are AMReX-specific functions that only get called by the constructor.)
  The fix is to initialize the operator object inside the routine in which it is used - either by making the member variable a pointer and using the :code:`new` keyword, or by just creating the variable inside the function.
  
  

Generating this documentation
=============================

Generating documentation requires the following packages:

* Doxygen (on Ubuntu: :code:`sudo apt install doxygen`)
* Sphinx (on Ubuntu: :code:`sudo apt install python3-sphinx`)
* Breathe (on Ubuntu: :code:`sudo apt install python3-breathe`)
* M2R (on Ubuntu: :code:`python3 -m pip install m2r`)
* RTD theme (on Ubuntu: :code:`python3 -m pip install sphinx_rtd_theme`)
* GraphViz (on Ubuntu: :code:`sudo apt install graphviz`)

To generate the documentation, type

.. code-block::

    make docs

(You do not need to run :code:`./configure` before generating documentation.)
Documentation will be generated in `docs/build/html` and can be viewed using a browser.

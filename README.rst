.. raw:: html

   <p align="center">
       <img src="https://www.solids.group/wp-content/uploads/2025/03/alamo3-inv.png" alt="alamo" width="400">
   </p>

   <p align="center">
       <img src="https://github.com/solidsgroup/alamo/actions/workflows/linux.yml/badge.svg?branch=development">
       <img src="https://github.com/solidsgroup/alamo/actions/workflows/coverage.yml/badge.svg?branch=development">
       <img src="https://img.shields.io/github/last-commit/solidsuccs/alamo/development.svg?label=last%20commit%20%28development%29">
       <img src="https://img.shields.io/github/contributors/solidsuccs/alamo.svg">
       <img src="https://img.shields.io/github/issues-pr/solidsuccs/alamo.svg">
       <img src="https://img.shields.io/github/issues/solidsuccs/alamo.svg">
       <img src="https://zenodo.org/badge/DOI/10.5281/zenodo.10381768.svg">
   </p>

.. getting-started:

Alamo is a high-performance scientific code that uses block-structured adaptive mesh refinement
to solve such problems as: the ignition and burn of solid rocket propellant, plasticity, damage
and fracture in materials undergoing loading, and the interaction of compressible flow with
eroding solid materials. Alamo is powered by AMReX, and provides a set of unique methods,
models, and algorithms that enable it to solve solid-mechanics problems (coupled to other
physical behavior such as fluid flow or thermal diffusion) using the power of block-structured
adaptive mesh refinement.

`Alamo documentation <https://solidsgroup.github.io/alamo/docs/>`_

Downloading Alamo
-----------------

Download alamo using git:

.. code-block::

    git clone git@github.com:solidsgroup/alamo.git
    
If you do not have a Github account and/or you have not uploaded your public SSH key, this will probably throw an error.
You can download alamo using HTTPS instead,

.. code-block::
    
    https://github.com/solidsuccs/alamo.git 

Note, that you will not be able to push anything using HTTPS authentication.
The :code:`master` branch is the most stable and is what is checked out by default.
The :code:`develompent` branch is generally stable, and includes the latest functionality.
To switch to :code:`development`, in the alamo directory,

.. code-block::
    
    git checkout development
    
Installing dependencies
-----------------------

Alamo is routinely run and tested on Ubuntu and MacOS.
You can use the
`System Install Scripts <https://solidsgroup.github.io/alamo/docs/GettingStarted.html#system-install-scripts>`_
to install all necessary dependencies for your system.

Setting default MPI
-------------------

It may be necessary to use a specific MPI distribution.
On Ubuntu, you can change the distribution with the following:

.. code-block::

    $> sudo update-alternatives --config mpi

    There are 2 choices for the alternative mpi (providing /usr/bin/mpicc).
    
      Selection    Path                    Priority   Status
    ------------------------------------------------------------
    * 0            /usr/bin/mpicc.openmpi   50        auto mode
      1            /usr/bin/mpicc.mpich     40        manual mode
      2            /usr/bin/mpicc.openmpi   50        manual mode
    
    Press <enter> to keep the current choice[*], or type selection number:     

Do the same thing for mpirun.

.. code-block::

    $> sudo update-alternatives --config mpirun
    
Remember to run :code:`make realclean` every time you switch mpi versions. 

Configuring
-----------

To compile alamo, you must first run the configure script. 
This is done simply by running the following in the alamo directory 
(note that AMReX download is triggered by this command, so it may take a couple minutes to complete depending on your internet connection)

.. code-block::

    ./configure

By default, alamo will configure in 3D production mode. 
To compile in  2D debug mode, 

.. code-block::

    ./configure --dim=2 --debug

There are multiple compilation options available for Alamo, and they must all be specified at configure time.
For a complete listing of the Alamo configuration options, type

.. code-block::

    ./configure --help


.. NOTE:: 
    The configure script produces output designed to assist in determining compile issues with Alamo.
    Whenever you request help with alamo, please always include the complete output of the configure script.

Compiling
---------

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

Unit Testing
------------

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

Regression Testing
------------------

Alamo contains several `Regression Tests <https://solidsgroup.github.io/alamo/docs/Tests.html>`_ that are routinely tested
and checked with CI.
These are checked for `Performance <https://lookerstudio.google.com/s/id-e_zDzO8w>`_ 
and `Code Coverage <https://solidsgroup.github.io/alamo/cov/>`_


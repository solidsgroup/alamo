.. role:: cpp(code)
   :language: c++
.. role:: bash(code)
   :language: bash

400 Grain Benchmark Test
========================

(For background see Runnels *et al*: https://arxiv.org/abs/2001.04789)

Instructions to run
-------------------

1. Clone alamo

   .. code-block:: 
		   
        git clone git@github.com:solidsuccs/alamo.git
        cd alamo

2. Configure for 3D*

   .. code-block:: 
		   
        ./configure

   .. Note::
      python3 must be installed for the configure script to work.

3. Make

   .. code-block:: 
		   
        ./configure

   .. warning::
      If make fails in parallel, try making again.
      There are some known issues when building amrex and alamo together in
      parallel.
      Build succeeds if the last message is "DONE"

4. Run the test suite (assuming gcc is used)

   .. code-block:: 
		   
        ./bin/test-3d-g++
   
5. Run the benchmark:

   .. code-block:: 
		   
        mpirun -np $NPROCS ./bin/alamo-3d-g++ tests/400GrainStampede2/input
   
   By default output will go to :code:`tests/400GrainStampede2/input`.
   To change this, either modify the input file or run with
   
   .. code-block:: 
		   
        mpirun -np $NPROCS ./bin/alamo-3d-g++ tests/400GrainStampede2/input plot_file = /path/to/desired/output


.. note::
   The benchmark is passed if speedup meets or exceeds that in Figure 16 of the paper.

.. note::
   To run with 10 grains instead of 100, run with the argument

   .. code-block::

      ic.voronoi.number_of_grains = 10
   
   

Input file
----------

.. literalinclude:: ../../../tests/400GrainStampede2/input
   :caption: 400 Grain benchmark input file (tests/Eshelby/input)
   :language: makefile

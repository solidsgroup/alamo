.. _install_hpc:

:fas:`server` HPC 
=========================

This page provides reference information for compiling and running Alamo on high-performance computing (HPC) clusters.

.. WARNING::
    
    These instructions are for reference only and may not always work. The Alamo developers do not manage the software on these clusters, and configurations may change over time. If you encounter outdated instructions, please `open an issue on GitHub <https://github.com/solidsgroup/Alamo/issues>`_.

Managing dependencies on an HPC cluster
---------------------------------------

Compiling and running Alamo on an HPC cluster differs slightly from doing so on a local machine, primarily due to dependency management. HPC clusters often provide multiple versions of software dependencies for users. To manage these dependencies efficiently, HPC clusters commonly use `Environment Modules <https://modules.sourceforge.net/>`_ (or simply *modules*), which help users to load and unload software as needed. Most HPC clusters provide the required modules for compiling Alamo.

To load an environment module:

.. code-block:: bash

   module load module_1 module_2 ...

To unload an environment module:

.. code-block:: bash

   module unload module_1 module_2 ...

To unload all modules:

.. code-block:: bash

   module purge

While Environment Modules are widely used, another tool called `Spack <https://spack.readthedocs.io/en/latest/index.html>`_ was specifically developed for dependency management in shared computing environments. You may encounter either system while working on an HPC cluster. The instructions on this page assume the use of Environment Modules.

Configuring Alamo on an HPC cluster
-----------------------------------

Configuring Alamo on an HPC cluster is similar to configuring it on a local machine. However, you must ensure that a compiler, :code:`mpich`, and Python 3 are available via modules or other means. Additionally, Alamo relies on the `Eigen library <https://eigen.tuxfamily.org/index.php?title=Main_Page>`_, which can either be loaded as a module or installed during configuration (preferred) with the :code:`--get-eigen` flag:

.. code-block:: bash

   ./configure --get-eigen ...

Compiling Alamo on an HPC cluster
---------------------------------

Compiling code can be resource-intensive and time-consuming. Running large, multithreaded operations on the `login node <https://www.hpc.iastate.edu/guides/introduction-to-hpc-clusters/what-is-an-hpc-cluster>`_ of an HPC cluster is generally discouraged. To avoid this, compile within an interactive job or by submitting a batch job. If the cluster uses the `Slurm Workload Manager <https://slurm.schedmd.com/overview.html>`_, an interactive job can be started with `salloc <https://slurm.schedmd.com/salloc.html>`_:

.. code-block:: bash

   salloc --nodes=1 --cpus-per-task=16 --mem=16G --time=10:00

This command requests a single node with 16 cores and 16 GB of memory for 10 minutes. Once the resources are allocated, your shell will reload, and you can compile with:

.. code-block:: bash

   make -j16  # Replace 16 with the number of cores requested

Alternatively, if interactivity is not required, submit a non-interactive compilation job using `srun <https://slurm.schedmd.com/srun.html>`_:

.. code-block::

   srun --nodes=1 --cpus-per-task=16 --mem=16G --time=10:00 make -j16

Adjust resource requests as needed. To reduce wait times, specify a shorter duration using the :code:`--time` flag.

Running Alamo on an HPC cluster
-------------------------------

To verify that a simulation starts correctly, an interactive job may suffice. However, for full simulations, you should submit a batch job to the cluster's workload manager. For Slurm-based clusters, use the `sbatch <https://slurm.schedmd.com/sbatch.html>`_ command to submit a job script:

.. code-block::

   sbatch /path/to/job_script

The sections below have examples of job scripts that can be modified to suit your needs.

Reference scripts for Nova
--------------------------

.. NOTE::
    
    For more information on Nova or HPC in general, Iowa State University provides an `extensive online guide <https://www.hpc.iastate.edu/guides>`_.

By default, Python and GCC modules are loaded when you log in to Nova. However, additional modules are required for certain tasks:

+------------------------+--------------------+
| Task                   | Required Module(s) |
+========================+====================+
| configuring w/ g++     | mpich              |
+------------------------+--------------------+
| configuring w/ clang++ | mpich llvm         |
+------------------------+--------------------+
| compiling w/ g++       | mpich              |
+------------------------+--------------------+
| compiling w/ clang++   | mpich llvm         |
+------------------------+--------------------+
| running Alamo          | mpich              |
+------------------------+--------------------+

The scripts below automatically handle module management.

The configuration and compilation scripts below can be run line-by-line or in a Bash script. If you do the latter, remember to make the file executable (:code:`chmod +x /path/to/file`).

gcc Configure and Compile Script
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   # /bin/bash
   module purge
   module load mpich
   ./configure --get-eigen
   srun --nodes=1 --cpus-per-task=16 --mem=16G --time=10:00 make -j16

clang++ Configure and Compile Script
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   # /bin/bash
   module purge
   module load mpich llvm
   ./configure --get-eigen --comp clang++
   srun --nodes=1 --cpus-per-task=16 --mem=16G --time=10:00 make -j16

Alamo Simulation Slurm Job Script
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash
    
    #!/bin/bash
    #SBATCH --time=24:00:00
    #SBATCH --nodes=1
    #SBATCH --ntasks-per-node=36
    #SBATCH --job-name="alamo"
    #SBATCH --mail-user=your_email@iastate.edu
    #SBATCH --mail-type=BEGIN
    #SBATCH --mail-type=END
    #SBATCH --mail-type=FAIL

    module purge
    module load mpich
    srun --mpi=mpix ./path/to/alamo/executable /path/to/input/file

The script starts a parallel job on Nova. Modify these parameters as needed:


* :code:`--time`: The *wall clock time limit* (maximum job duration)
* :code:`--nodes`: Number of nodes requested
* :code:`--ntasks-per-node`: Number of tasks per node
* :code:`--mail-user`: Email for notifications (remove if not needed)
* **Executable path**: Example: :code:`./bin/alamo-2d-clang++`
* **Input file path**: Specify the input file for Alamo

.. NOTE::

   Iowa State University provides a `Slurm job script generator for Nova <https://www.hpc.iastate.edu/guides/nova/slurm-script-generator-for-nova>`_, which can help generate job scripts.

Nova uses a *fair-share scheduling system* to prioritize job execution based on requested resources and past usage. To reduce wait times, request only necessary resources and set reasonable time limits.

Slurm automatically determines the number of cores based on the :code:`--nodes` and :code:`--ntasks-per-node` values. Refer to the `Nova hardware guide <https://www.hpc.iastate.edu/guides/nova>`_ for appropriate values.

Once modifications are made, submit the job with `sbatch <https://slurm.schedmd.com/sbatch.html>`_:

.. code-block:: bash

   sbatch /path/to/job_script.sh


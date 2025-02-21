.. _install_hpc:

:fas:`server` HPC 
=========================

The information on this page is a reference for compiling and running Alamo on high-performance computing (HPC) clusters.

.. WARNING::
    
    These instructions are for reference only and may not always work. The Alamo developers do not manage the software on these clusters, which may change over time. If you encounter outdated instructions, please `open an issue on GitHub <https://github.com/solidsgroup/Alamo/issues>`_.

Managing dependencies on an HPC cluster
---------------------------------------

Compiling and running Alamo on a high-performance computing (HPC) cluster differs slightly from doing so on a machine locally, primarily due to dependency management. HPC clusters often require many software dependencies for their users, sometimes with several versions of the same dependency. To streamline this, HPC clusters commonly use `Environment Modules <https://modules.sourceforge.net/>`_ (or simply *modules*), which help manage and load the necessary tools. Most HPC clusters likely provide all the required modules for compiling Alamo.

To load an environment module:

.. code-block::

   module load module_1 module_2 ...

To unload an environment module:

.. code-block::

   module unload module_1 module_2 ...

To completely unload all modules:

.. code-block:

   module purge

While Environment Modules have many uses outside HPC, a new tool called `Spack <https://spack.readthedocs.io/en/latest/index.html>`_ was developed specifically to tackle the challenges related to dependency management in shared computing environments. While working on HPC clusters, you may encounter either system. The resources available on this page are tailored to Environment Modules.

Before configuring or compiling Alamo, see the sections below for information on which modules must be loaded for each HPC cluster.

Configuring Alamo on an HPC cluster
-----------------------------------

Configuring Alamo on an HPC cluster is nearly the same as configuring Alamo on a local machine. However, you must be certain a compiler, mpich, and Python3 is available via modules or otherwise. Additionally, Alamo relies in the `package eigen <https://eigen.tuxfamily.org/index.php?title=Main_Page>`_. This dependency can either be loaded as a module or installed during configuration (preferred) with the :code:`--get-eigen` command-line argument:

.. code-block:

   ./configure --get-eigen ...

Compiling Alamo on an HPC cluster
---------------------------------

Compiling code can be computationally expensive and time-consuming. Running large, multithreaded operations on the `login node <https://www.hpc.iastate.edu/guides/introduction-to-hpc-clusters/what-is-an-hpc-cluster>`_ of an HPC cluster is generally discouraged. To avoid this, it's recommended to compile within an interactive job or by submitting a job. If the cluster uses the `Slurm Workload Manager <https://slurm.schedmd.com/overview.html>`_, this can be accomplished with the `salloc <https://slurm.schedmd.com/salloc.html>`_ command:

.. code-block::

   salloc --nodes=1 --cpus-per-task=16 --mem=16G --time=10:00

This command will block until the requested resources are allocated. The example above requests a single node with 16 cores and 16 GB of memory. Once the resources have been allocated, your shell will reload, and you can then run :code:`make -j16`, replacing 16 with the number cores you requested. Alternatively, if you don't need interactivity, the compile command can be submitted with a `srun <https://slurm.schedmd.com/srun.html>`_ command:

.. code-block::

   srun --nodes=1 --cpus-per-takes=16 --mem=16G --time=10:00 make -j16

Once the resources are available, your compilation job will run. For either compilation method, adjust the resources as needed. To speed up resource allocation, request a shorter time with the :code:`--time` flag.

Running Alamo on an HPC cluster
-------------------------------

To check whether a simulation will start correctly, an interactive job may be sufficient; however, to run a simulation to completion, you should submit a batch job to the cluster's workload manager. For the `Slurm Workload Manager <https://slurm.schedmd.com/overview.html>`_, jobs are defined by a bash script known as a job script and submitted with the `sbatch <https://slurm.schedmd.com/sbatch.html>`_ command:

.. code-block::

   sbatch /path/to/job_script

The sections below have examples of job scripts that can be modified to suit your needs.

Nova Scripts
------------

.. NOTE::
    
    For more information on Nova or HPC in general, Iowa State University offers an `extensive guide online <https://www.hpc.iastate.edu/guides>`_.

The Python and gcc modules are loaded by default when you log in to Nova. However, the following modules are required to perform certain tasks:

+-----------------------+--------------------+
|         Task          | Required Module(s) |
+=======================+====================+
|   configuring (g++)   |       mpich        |
+-----------------------+--------------------+
| configuring (clang++) |     mpich llvm     |
+-----------------------+--------------------+
|    compiling (g++)    |       mpich        |
+-----------------------+--------------------+
|  compiling (clang++)  |     mpich llvm     |
+-----------------------+--------------------+
|     running Alamo     |       mpich        |
+-----------------------+--------------------+



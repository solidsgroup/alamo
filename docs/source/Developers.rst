.. role:: cpp(code)
   :language: c++

.. _developers:

=======================================
:icon:`developer_guide` Developer Guide
=======================================

.. NOTE::

    This page is a work in progress and is not yet comprehensive.


Conventions
===========

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Namespaces
----------

Classes and structures are organized using namespaces.
Namespaces should be brief but descriptive, and two nested namespaces should correspond to two nested subsets.
For instance :code:`Shape::Polygon::Triangle` is acceptable whereas :code:`Shape::Circle::Round` is not (since "round" is not a subset of "Circle").

The directory structure must reflect the namespace structure.
For instance :code:`Shape::Polygon::Triangle` must be declared in :code:`src/Shape/Polygon/Triangle.H`.

Regular and Pure Virtual classes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Pure virtual classes must have the same name as the namespace in which they are contained.
For instance :code:`Model::Elastic::Elastic` is pure virtual, whereas :code:`Model::Elastic::Isotropic` must inherit from :code:`Model::Elastic::Elastic` and is not pure virtual.

Include Guards
--------------

Include guards should be all uppercase, and should be the same as the path (relative to `src/`) with `/` replaced by underscores.
For instance: a file called `src/Integrator/Integrator.H` should have the include guard 

.. code-block:: cpp

   #ifndef INTEGRATOR_INTEGRATOR_H
   #define INTEGRATOR_INTEGRATOR_H
   ...
   #endif

Member and Argument Variable Names
----------------------------------

Members of a class should have names beginning with :code:`m_`; argument variables should begin with :code:`a_`.
For instance:

.. code-block:: cpp

   class MyClass
   {
   public:
      MyClass(a_myvariable) : m_myvariable(a_myvariable)
      {}
   private:
      m_myvariable;
   }



Python (In development)
=======================

Alamo supports a (currently limited) Python interface.
It can be run as simply as

.. code-block:: bash

      python alamo.py

where :code:`alamo.py` is a python script that you write.
An example (that currently works) is

.. code-block:: cpp

      import alamo

      alamo.Util.Initialize()

      mytest = alamo.Test.Operator.Elastic()

      mytest.setBounds([1.0,1.0,1.0])

      for lev in [1,2,3]:
          print("Levels of refinement: " + str(lev))
      mytest.Define([32,32,32],lev,1,mytest.Grid.XYZ)
      failed = mytest.TrigTest(0,0,1,"")
      if failed: print("Test failed")
      else: print ("Test passed")

      alamo.Util.Finalize()

Note that the C++ namespace structure is mirrored in the Python bindings.

.. WARNING::
   Note that object constructors require :code:`()`, e.g. :code:`test = alamo.Test.Operator.Elastic()`.
   If you forget the :code:`()` the python interpreter **will not complain!**
   Instead, it will give you an extremely cryptic message about member methods not having the right number of arguments.
   A typical error message might look like this:

   .. code-block:: bash

           Traceback (most recent call last):
           File "alamo.py", line 7, in <module>
           mytest.setBounds([1.0,1.0,1.0])
           TypeError: unbound method Boost.Python.function object must be called with Elastic instance as first argument (got list instance instead)

Compiling Alamo Python interface
--------------------------------

The python interface requires all code to be compiled using the :code:`-fPIC` flag.
To compile AMReX,

.. code-block:: bash

        ./compile [whatever flags you would normally use
        make USE_COMPILE_PIC=TRUE
        make install

To compile Alamo, you need the Boost python library.
On ubuntu, install with

.. code-block:: bash

        sudo apt install libboost-python-dev

Then to compile Alamo, 

.. code-block:: bash

        ./configure [whatever flags you would normally use] --python
        make python

This will create a file called :code:`alamo.so` in the Alamo root directory.
(Note that this file must be in your python path (or your current directory) to use.)

Extending the Alamo Python interface
------------------------------------

Python bindings are currently limited and you will likely need to add them if you want them for a specific component.
Alamo uses the Boost Python library to create the Python bindings.
All of this code is located in the :code:`py/` directory, and mimics the directory structure of :code:`src/`.
These are all C++ files, but use the :code:`.cpy` extension.
To add bindings (e.g. to an object) you must

1. Create an appropriately named file (e.g. :code:`py/Test/Operator/Elastic.cpy` to bind :code:`Test::Operator::Elastic`).
2. Define a function that encapsulates the bindings (e.g. :code:`void exportTestOperatorElastic`).
3. Write the Boost Python bindings in the function.
   Keep in mind that you need to specify namespace explicitly.
4. Include the file in :code:`py/alamo.cpy`
5. Call the function in :code:`py/alamo.cpy`
   

Restart files
=============

You can restart a simulation from any point at which you have dumped an output file.
For instance, suppose you have run alamo with this command 

.. code ::

    ./bin/alamo-2d-g++ input

and it produced the following output.

.. code ::

    ./output
        00000cell 00000node
        00010cell 00010node
        00020cell 00020node
        00030cell 00030node
        00040cell 00040node
        00050cell 00050node
        00060cell 00060node
        00070cell 00070node
        celloutput.visit
        nodeoutput.visit
        metadata
    
You may wish to *restart* the simulation without starting from the beginning.
Perhaps the simulation was fine at :code:`00060` but became unstable at :code:`00070`, and you want to continue from that point with a different parameter value.
To do that, you can simply run:

.. code ::

    ./bin/alamo-2d-g++ input restart=./output/00060cell

for the simulation to pick up where it left off, but with updated values (or updated code) that may have changed since you first ran it.

It is important to distinguish what the restart function does and does not do.

The restart function *DOES*:
    - Import all of the data from fields named "XXX" "YYY" etc into corresponding fields with the same names in the new simulation.
    - Set the clock equal to the time corresponding to the output file.
    - Set the current timestep equal to that corresponding to the output file.

The restart function DOES NOT:
    - Set any other simulation parameters. The input file is responsible for setting all simulation parameters. It is up to you to use them (or change them) as needed.
    - Initialize non-matching data. If you try to restart using data containing fields "AAA" and "BBB" but your simulation expects fields "BBB" and "CCC", only "BBB" will be initialized - the rest is undefined.
    - Continue to output to the same directory. A different directory will be created. Only the metadata file will record where the restard data came from.
   

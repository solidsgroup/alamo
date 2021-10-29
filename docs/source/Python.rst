.. _building-python:
.. role:: cpp(code)
   :language: c++

Python Interface
================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Running Alamo with Python
-------------------------

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

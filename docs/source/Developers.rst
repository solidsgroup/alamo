.. role:: cpp(code)
   :language: c++

.. _developers:

=================================
:fas:`code;fa-fw` Developer Guide
=================================

:fas:`rocket;fa-fw` Step-by-step development guide
====================================================

This is a quick tutorial intended for new Alamo developers to get up to speed with contributing to the Alamo code base.

#. **Create Github credentials.**
    #. If you do not have one already, create a new GitHub account, and email your username to brunnels@iastate.edu to request push access
    #. Follow `these Github instructions <https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account>`__ 
       to add an SSH key to your account.
    #. :bdg-danger:`important` Make sure your repository points to the SSH version of Alamo. 
       You can check this with the command

      .. code-block:: bash
         
          git remote show origin
      
      If the output begins with

      .. code-block:: bash

         Fetch URL: git@github.com:solidsgroup/alamo.git
         Push  URL: git@github.com:solidsgroup/alamo.git

      then you have configured correctly :icon:`check`.
      However, if the output begins with

      .. code-block:: bash

         Fetch URL: https://github.com/solidsuccs/alamo.git
         Push  URL: https://github.com/solidsuccs/alamo.git

      this means that you are using HTTPS authentication, which will not allow you to push.
      To fix, run

      .. code-block:: bash

         git remote remove origin
         git remote add origin git@github.com:solidsgroup/alamo.git

      Run a quick :code:`git pull` to make sure that you still have access.
      If you get an authentication error, this likely means that your SSH key is not configured correctly.
      
#. **Create a new branch.**
   You can create a branch `from the terminal <https://git-scm.com/book/en/v2/Git-Branching-Basic-Branching-and-Merging>`_, 
   or you can `use the Github online interface <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-and-deleting-branches-within-your-repository>`_.
   Follow these guidlines when naming your branch

   - Use combinations of lowercase letters and hyphens, avoid mixed case, numbers, and underscores (:code:`gas-combustion` good, :code:`GasCombustion2` bad).
   - Be as descriptive as possible and avoid overly general names (:code:`gb-damage-model` good, :code:`model1` / :code:`mymodel` /  :code:`foo` bad). Long is ok.
   - Keep it professional. No profanity or excessive whimsy in branch names.

   If you have created a branch (e.g. :code:`mytestbranch`) locally, make sure to run

   .. code-block:: bash

      git push --set-upstream origin mytestbranch

   to make pushing easier.

#. **Create a capability input file**.
   This is an input file that is designed to demonstrate the functionality of your new code.
   The capability input file should always be located in the alamo root directory.
   The :code:`plot_file` should always be :code:`output`.
   Once you are ready to merge your code into development, convert your input file to a regression test and remove it from the root. 
   (See the Autodoc and Autotest section).

#. **Use an EditorConfig-compliant text editor.**
   Alamo uses `EditorConfig <https://editorconfig.org/>`_ to maintain code formatting standards.
   `VS Code <https://code.visualstudio.com/>`_ is recommended, as it supports multiple keybinding schemes and automatically
   supports editorconfig.

#. **Write your code.**
   As you are writing, keep the following in mind:

   ======================== ===================================================================================================
   :bdg-success:`always`    Commit your code - at least once per day you are coding, even if the code is not working.
                            There are `multiple great tutorials <https://www.atlassian.com/git/tutorials/saving-changes>`_
                            if you are not sure how to commit your code.
   ------------------------ ---------------------------------------------------------------------------------------------------
   :bdg-success:`always`    Write meaningful commit messages. This will help you more than you may realize when you are 
                            trying to debug your code later.
   ------------------------ ---------------------------------------------------------------------------------------------------
   :bdg-success:`always`    Follow the development guide below and stick to the specified convention.
   ------------------------ ---------------------------------------------------------------------------------------------------
   :bdg-primary:`regularly` Merge the latest development into your code (see below for more details).
                            The more frequently you do this, the less painful your life will be later.
   ------------------------ ---------------------------------------------------------------------------------------------------
   :bdg-primary:`regularly` Run the regression test suite with :code:`make test`.
                            This helps to ensure that you are not breaking your (or someone else's) code.
   ------------------------ ---------------------------------------------------------------------------------------------------
   :bdg-primary:`regularly` Ensure that your commits pass the `commit tests <https://github.com/solidsgroup/alamo/actions>`_.
                            Sometimes this is not possible, especially if your code is in active development, but the more you
                            keep your code up to standard, the easier your life will be later.
   ------------------------ ---------------------------------------------------------------------------------------------------
   :bdg-primary:`regularly` Document your code.
                            Use the autodoc system to ensure a base level of documentation.
   ------------------------ ---------------------------------------------------------------------------------------------------
   :bdg-warning:`rarely`    Break backwards compatibility.
                            If you are implementing an improved version of a model that is better than the legacy model,
                            leave the option to run using the existing model.
                            This is important because publications may have used the legacy model, and may need the legacy
                            model to reproduce important results.
                            Backwards compatibility should only be broken when implementing fundamental new features.
                            Even then, changes required to reproduce legacy results should be kept minimal.
   ------------------------ ---------------------------------------------------------------------------------------------------
   :bdg-danger:`never`      Modify the GitHub actions to get your branch to pass the tests.
                            Email failure notifications may be annoying, but they are there for a reason.
   ------------------------ ---------------------------------------------------------------------------------------------------
   :bdg-danger:`never`      Commit large or unnecessary files to the repository.
                            Use :code:`git status` liberally, before and after :code:`git commit`, before you push.
                            Edit your `gitignore <https://www.atlassian.com/git/tutorials/saving-changes/gitignore>`_ file
                            as needed.
   ------------------------ ---------------------------------------------------------------------------------------------------
   :bdg-danger:`never`      Commit sensitive content to the repository.
                            This is usually an issue only if your project is export controlled, which is rare.
                            But when in doubt, ask.
   ======================== ===================================================================================================

#. **Merge your branch into development**
   To include your changes into the development branch, submit a
   `pull request <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request>`_
   from your branch into development.
   This will trigger a review before the merge can be approved.
   (Note that all tests must pass before the merge will be reviewed.)
   The more frequently you merge in the latest from development, the easier your merge will be.   
   
   :bdg-primary-line:`How do I know my code is ready to merge?` 
   You should always merge your code once you know that it works, but sometimes this is a gray area.
   You should always merge your code if you are submitting a publication, or if you are a PhD student about to graduate.
   

:fas:`user-graduate;fa-fw` Tutorial: A new Integrator
===========================================================

Integrators are the basic workhorse in Alamo.
They are called Integrators because they usually (although not always) are intended to integrate a PDE in time.
In this tutorial, we will create a new integrator based on the existing HeatConduction integrator.

#. Create a copy of :code:`./src/Integrator/HeatConduction.H` inside the same directory, called `MyPDE.H`
   (We will use the working name "MyPDE" here - you can replace with your own name.)

#. Rename all instances of HeatConduction with MyPDE inside :code:`./src/Integrator/MyPDE.H`.
   This includes the names throughout the code as well as the `include guard <https://en.wikipedia.org/wiki/Include_guard>`__
   in the first two lines.

#. Include the file in :code:`./src/alamo.cc`

   .. code-block:: cpp

      #include "Integrator/MyPDE.H"

   and add a caller inside the main function:

   .. code-block:: cpp

      //existing
      else if (program == "thermoelastic")integrator = new Integrator::ThermoElastic(pp);
      //new
      else if (program == "mypde")        integrator = new Integrator::MyPDE(pp);
   
#. Finally copy the Heat Conduction example file to the root directory

   .. code-block:: bash

      cp ./tests/HeatConduction/input ./input
   
   In the input file change the :code:`plot_file` to `output` and :code:`alamo.program` to :code:`mypde`.

#. You should now be able to compile and run the code.
   The result will be the same as for the heat conduction example, but will be based on the newly copied integrator.

#. You can now begin alternating the code to achieve different outputs.
   The HeatConduction integrator has extensive, line-by-line documentation.
      

:fas:`rectangle-list;fa-fw` Conventions
=======================================

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



:fab:`python;fa-fw` Python (In development)
===========================================

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
   

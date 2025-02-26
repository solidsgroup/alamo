.. _autodoc:

============================================
:fas:`swatchbook;fa-fw` Autodoc and Autotest
============================================

Alamo includes a self-contained automatic documentation (autodoc) and automatic test (autotest) 
system.
The purpose of the custom autodoc and autotest systems is to:

1. Make it easy for developers to add documentation and regression tests
2. Keep documentation and tests accessible and up-to-date
3. Prevent old or out-of-date documentation.

Autodoc system
==============

The autodoc system parses the source files stored in :code:`./src/` to generate the content of the
Inputs portion of this documentation.
There are two ways to create compliant documentation: via parser comments, and via header file 
comments.

Parser comments
---------------

The :code:`IO::ParmParse` class (inherited from the AMReX :code:`ParmParse` class) is used to read
all input paramters from the input file, usually though lines like this:

.. code:: cpp

    pp.query("lame",lambda);
    pp.query("shear",mu);

To generate the autodocumentation, simply add a standard C++ comment following the query line:

.. code:: cpp

    pp.query("lame",lambda); // Lame parameter
    pp.query("shear",mu);    // Shear modulus

or immediately preceeding the query line:

.. code:: cpp

    // Lame parameter
    pp.query("lame",lambda);
    // Shear modulus
    pp.query("shear",mu);

Note that the parser object must be called :code:`pp` for the autodoc system to locate the query.
As long as the convention is followed, the content of the comments will be automatically scraped
and included in the Inputs section of the documentation. 
(For instance: :ref:`Model::Solid::Linear::Isotropic`)


Header comments
---------------

General comments to document a class, or collection of methods, goes at the top of the associated 
header file.
Use standard C++ comments (:code:`//`).
These comments will get scraped by the autodoc system, and will be formatted along with the rest
of the inputs.

For example, consider the documentation for the linear isotropic material model class.
The following leading comment is formatted at :ref:`Model::Solid::Linear::Isotropic`

.. literalinclude:: ../../src/Model/Solid/Linear/Isotropic.H
    :language: cpp
    :lines: 1-50


ReStructuredText
----------------

All comments are formatted using 
`restructuredtext <https://docutils.sourceforge.io/docs/ref/rst/restructuredtext.html>`_
markup.
You can make use of this is if you like, but it is not required for good documentation.

Documentation generation
------------------------

To generate this documentation, complete with the scraped markup, run

.. code:: bash

    make docs

in the alamo root directory.
The file :code:`alamo/docs/requirements.txt` contains a list of the necessary packages.
You can install them using pip.
Once the documenation generation is complete, you can view it by

.. code:: bash

    google-chrome docs/build/html/index.html


Autotest system
===============

Alamo contains a complete regression testing system to ensure that established capabilities do not get lost during subsequent development.
Regression tests are run automatically by github, and must all pass prior to merge into the development or master branches.
For information on how to run the regression tests, see :ref:`Testing`.

Test requirements
-----------------

Remember that your test will be run automatically on GitHub and by other users, so avoid creating long, memory-intensive tests.
Instead, design your tests so that they run long enough to identify errors, but short enough so that they can be performed tractable.
Make sure that your tests are not redundant with existing tests as well.

Creating a test
---------------

#. Create a subdirectory in the :code:`./tests` with your test name, e.g. "MyTest" [:code:`./tests/MyTest`]
#. Place your test input file in this directory and name it `input` [:code:`./tests/input`].
#. At the top of the test input file, add

   .. code-block:: cpp

      #@ [sectionname]

   where :code:`sectionname` can be any unique handle name.
   The :code:`#@` directive is ignored by alamo, but all comments beginning with it will be parsed by the autotest system.

#. Test your test by running

   .. code-block:: bash

      ./scripts/runtests ./tests/MyTest

   If alamo has been compiled in 3D, this should cause your test to run.
   The output of your test will be stored in

   .. code-block:: bash

      ./tests/MyTest/output_YYYY-MM-DD_HH.MM.SS_hostname_sectionname

   
Once you have done this, your test will be executed along with all other regression tests by GitHub.
If your test fails to complete - for instance, it segfaults or aborts abnormally - this will trigger a "failed test"

Input file directives
---------------------

There are several directives that you can use for your inputs.
These are written using the #@ comments at the top of the input file.

.. code-block:: make

    #@  [2D-serial-5levels]        | each [...] specifies a run of the test 
    #@  dim    = 2                 | specify two dimensions
    #@  nprocs = 1                 | specify running in serial (the default)
    #@  check  = false             | for this test we do not want to run the verification
    #@   
    #@  [3D-parallel-4levels]      | another test
    #@  dim    = 3                 | specify running in 3D (the default)
    #@  nprocs = 4                 | specify run in parallel with 4 processors
    #@  args   = amr.max_level=4   | specify an additional argument to be added to input file
    #@  benchmark-beaker = 16.10   | timing test: benchmark-id record current average time to
    #@  benchmark-statler = 11.36  |     run on platform "id". The autotest system will let you 
    #@  benchmark-github = 22.75   |     know if future changes slow the test down.
    #@  ignore = myargument        | tell alamo to ignore certain arguments

For additional examples, see the Tests section.

Test verification
-----------------

Successfully running a test does not tell you much unless you know that the correct answer was produced.
To automatically verify your test, add an executable file called :code:`test` to the test directory:

.. code-block:: bash

   ./tests/MyTest/test

The test script must be executable, and can be written in a language of your choice, as long as you have the appropriate shebang at the top.
There are two requirements on the test script:

#. It must take a single argument, the test directory name.
#. It produces a zero error code if the test passes, and a nonzero error code if the test fails.

It is up to you to make sure that the test accurately assesses the output.
You can store reference data inside the test directory; e.g.

.. code-block:: bash

   ./tests/MyTest/reference/stress_xx.dat
   ./tests/MyTest/reference/stress_xy.dat
   ....

as long as the data files are reasonably small in size and, of course, are text-based.
For example tests, see the existing tests in the repository.



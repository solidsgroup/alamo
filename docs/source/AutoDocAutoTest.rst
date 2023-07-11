.. _autodoc:

========================================
:icon:`play_lesson` Autodoc and Autotest
========================================

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

Alamo contains a complete regression testing system to ensure that established capabilities
do not get lost during further development.
Regression tests are run automatically by github, and must all pass prior to merge into the
development or master branches.
For information on how to run the regression tests, see :ref:`Testing`.


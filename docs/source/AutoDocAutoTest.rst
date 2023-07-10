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
:ref:`Inputs` portion of this documentation.
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







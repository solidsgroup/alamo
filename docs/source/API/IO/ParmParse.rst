.. role:: cpp(code)
   :language: c++

IO::ParmParse
-------------

:code:`IO::ParmParse` inherits from :code:`amrex:ParmParse` which is documented at https://amrex-codes.github.io/amrex/docs_html/Basics.html#parmparse -- see that website for general usage.

General features of the :code:`IO::ParmParse` (instead of the AMReX version) are:

    *  Ability to read :code:`Set::Vector` and :code:`Set::Matrix` with :code:`pp.queryarr("name",obj)`.
    *  Ability to read arbitrary class with :code:`pp.queryclass` method - more on this below.

Writing class parse files
~~~~~~~~~~~~~~~~~~~~~~~~~

Suppose you have the following class:

.. code-block:: cpp

    class Circle
    {
        Set::Scalar radius;
        Set::Vector x0;
    }

and you would like to parse the following input file:

.. code-block:: bash

    parser.circle1.radius = 2.0
    parser.circle1.x0 = 1.0 2.0 3.0
    parser.circle2.radius = 4.0
    parser.circle2.x0 = 5.0 6.0 7.0

where there may be other :code:`parser` entries that you are reading in.

.. code-block:: cpp

    class Circle
    {
        Set::Scalar radius;
        Set::Vector x0;
    public:
        static void Parse(Circle & value, IO::ParmParse &pp)
        {
            pp.query("radius",value.radius);
            pp.queryarr("x0",value.x0);
        }
    }

then, in the application code, you would do the following:

.. code-block:: cpp

    IO::ParmParse("parser"); 
    // other parse code
    Circle circle1, circle2; // Note that we want two circle objects
    pp.queryclass("circle1",circle1);
    pp.queryclass("circle2",circle2);

Note that the :code:`.circle1.` and :code:`.circle2.` prefixes are appended automatically.
If you want to use no prefix, e.g. you want to read the following:

.. code-block:: bash

    parser.radius = 2.0
    parser.x0 = 1.0 2.0 3.0

then you can simply call

.. code-block:: cpp

    pp.queryclass(circle1);


Parsing abstract pointers to derived classes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Suppose the :code:`Circle` class inherits from an abstract :code:`Shape` class:

.. code-block:: cpp

    class Circle : public Shape
    ...

and in your code, you have the following pointer

.. code-block:: cpp

    Shape * shape;

that you would like to initialize as a :code:`Circle`. 
You would do the following:

.. code-block:: cpp

    Shape * shape;
    shape = new Circle();
    pp.queryclass("circle1",static_cast<Circle*>(shape));

The :code:`static_cast` is necessary for the compiler to know to call the :code:`Circle` parsing function. 
(Otherwise there is nothing that tells the compiler that :code:`shape` points to a :code:`Circle` object.)





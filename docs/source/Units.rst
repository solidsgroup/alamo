============================================
:fas:`ruler;fa-fw` Unit system
============================================


Introduction
------------

Many scientific and engineering codes accept numerical inputs for physical parameters but ignore their dimensions entirely. This can lead to subtle bugs, difficult-to-interpret outputs, and significant user error, sometimes resulting in
`massive system level failure <https://en.wikipedia.org/wiki/Mars_Climate_Orbiter>`_
Alamo takes a unit-aware approach by **explicitly tracking the physical dimensions of inputs**, which

- enforces dimensional consistency throughout the simulation.
- allows inputs to be expressed in arbitrary physical units.
- simplifies post-processing by producing outputs in user-specified units.

This design helps prevent unit mismatch errors and promotes reproducible, well-documented simulations.
At the time of this writing, enforcement is not yet universal, and the system remains fully backwards compatible with legacy input files that omit units.
However, unit-checking and enforcement will be increasingly integrated and eventually required of all new code development.

Specifying Units in Input Files
-------------------------------

System-level units are specified using the following commands.
Any valid length and time scale is acceptable.
(The default units are meters and seconds.)

.. code-block:: none

   system.length = mm                 # ✅ System output will be in millimeters
   system.time = ms                   # ✅ System output will be in milliseconds

Units are specified in the input file by appending them to the numerical value using a single underscore (:code:`_`).
Whitespace is not permitted.

.. code-block:: none

   timestep = 0.5_s                   # ✅ Valid
   timestep = 0.5_us                  # ✅ Valid: prefixes are acceptable
   timestep = 0.5                     # ⚠️ Valid, but will change based on system.time

   viscosity = 0.001_kg/m/s           # ✅  valid - compount units are acceptable
   density   = 1.0_g/cm^3             # ✅  valid - can use exponents
   velocity = 50_mph                  # ✅  valid - reduced units like mph, N, Pa are accepted
   velocity = 50_furlong/fortnight    # ✅  valid - any units can work :)

   timestep = 0.5s                    # ❌ missing underscore
   timestep = 0.5 s                   # ❌ missing underscore
   timestep = 0.5_m                   # ❌ will error if units are incorrect

This syntax ensures that the parser treats the entire token as a single string and avoids ambiguity.
There must be exactly one underscore separating the number and the unit expression.

Units can be **arbitrary**, as long as they are dimensionally equivalent to the expected quantity. For example, the following are all equivalent inputs for a velocity parameter:

.. code-block:: none

   velocity = 5.0_m/s                 # ✅  valid
   velocity = 5000.0_mm/s             # ✅  si-prefixes are applied at read time
   velocity = 196.85_in/s             # ✅  imperial units 
   velocity = 1.666667_ly/year        # ✅  less conventional but acceptable
   velocity = 1.0_furlong/fortnight   # ✅  any compatible unit system will work

Even uncancelled units like :code:`gal/ft^2` will be accepted—provided they are dimensionally consistent with the expected quantity (e.g., length in this case).

To determine what unit is expected for a given input, consult the **Inputs documentation**. 
Most common compound units like `Pa`, `N`, and `J` are supported out of the box.
Additional relationships (e.g., `1 cal = 4.184 J`) can be easily added by the user (see below).

Units API
--------------

Units are encapsulated in the :code:`Unit` class (:code:`Unit/Unit.H`).
This class is a standalone header that follows the Alamo convention but can be used outside of Alamo.

(The unit class **is not intended to be a high performance data type**: that is, it is used to read in and convert inputs - and possibly to recombine input parameters at parse time - but not to be used in power loops.
Unit resolution involves the use of :code:`std::map` structures, which are not intended for performance.)

The `Unit` class is tightly integrated with the input parser. A typical pattern for querying a parameter with a unit is:

.. code-block:: c++

   pp.query_default("density", rho, "0.0_kg/m^3", Unit::Density());

This call will:

- Parse the string and extract the numeric value and unit.
- Check that the unit is dimensionally consistent with density.
- Convert the value to system units (as defined by `system.length`, `system.time`, etc.).
- Store the result in :code:`rho`.

The `Unit` class supports full dimensional algebra, including:

- Multiplication and division: :code:`Unit::Density() / Unit::Acceleration("m/s^2")` yields `Unit::Force()`.
- Parsing unit expressions from strings: :code:`Unit::Parse("kPa")`, :code:`Unit::Parse("3_m^2/s")`, :code:`Unit::Parse(3.0, "m^2/s")` etc.
- Comparison of dimensional consistency: :code:`u1.isType(u2)`.

This makes it easy to perform internal unit manipulations or validate user-specified expressions.

New relationships between units can be easily added by adding additional entries to the :code:`Unit::compound` map.
For example, the following defines force relationships:

.. code-block:: c++

   {"N",         {1.0,     "kg*m/s^2" }}, // Define the Newton: N = 1 kg*m/s^2
   {"kN",        {1e3,     "N"        }}, // Define kilonewton  kN = 1e3 N
   {"lbf",       {4.44822162, "N"     }}, // Pound-force in terms of Newtons
   {"dyn",       {1e-5,    "N"        }}, // dyn in terms of Newtons

Relationships are resolved recursively; as long as the new unit is defined in terms of existing units, any form is acceptable.


=============
Running Alamo
=============

This section describes how to run :code:`Integrator`-based applications.


.. NOTE::

    This section is a work in progress and is not yet comprehensive.

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

.. role:: cpp(code)
   :language: c++

Simulation Browser Analysis (SimBA)
===================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Reproducibility is critical to ensuring the integrity of simulation results.
SimBA is a set of python scripts designed to manage, organize, and view Alamo simulations through a web browser.
It integrates with Alamo's metadata recording functionality and uses a database to manage simulation information.

Running SimBA
-------------

Suppose you have a directory containing a collection of simulation results.

.. code-block:: bash

   output
   output.old.1111111
   output.old.2222222
   output.old.3333333

Begin by running :code:`./simba/database.py output*`.
This should produce the following output:

.. code-block:: bash
   
    ADDED TABLE: simulation_data
    ├╴Inserting: output
    ├╴Inserting: output.old.1111111
    ├╴Inserting: output.old.2222222
    ├╴Inserting: output.old.3333333
    └╴Done

and should create a database called :code:`results.db` in your directory.
Note that the database name and table name are all customizable.

Now, you can start the web interface with :code:`./simba/web.py`.
Open a browser and go to :code:`localhost:5000`.
You should see each of your simulation entries listed.
You can click on each entry to see greater detail, and you can edit the "Description" and "Tags" fields to record information about the simulation.
You can also delete simulations (and the corresponding directory).

The SimBA web interface will automatically recognize images stored within the simulation directory, and will display them under the simulation entry.

Updating SimBA Records
----------------------

Simulation directories are never overwritten - following the AMReX implementation, old directories are renamed if there is a naming conflict.
SimBA recognizes this and accounts for it.
In the above example, suppose that another simulation was run so that :code:`output` was renamed to be :code:`output.old.4444444`.
Running SimBA again produces the following output

.. code-block:: bash
   
    ADDED TABLE: simulation_data
    ├╴Inserting: output
    ├╴Updating:  output.old.1111111 ( record already exists )
    ├╴Updating:  output.old.2222222 ( record already exists )
    ├╴Updating:  output.old.3333333 ( record already exists )
    ├╴Moving:    output --> output.old.4444444
    └╴Done

Alamo records a unique identifier for each simulation that is stored in the :code:`metadata` file - this way simulations can be moved around and renamed without losing their designation in the SimBA database.


.. NOTE::
   SimBA does need to be updated using the :code:`database.py` script - the web interface does not
   do this automatically (yet).


Sharing results on the network
------------------------------

An advantage of SimBA is that it allows the posting of results to be viewed over the internet.
For instance, suppose you want someone else on your network to be able to view your SimBA page and the outputs that you have generated.
To do this, start the web interface this way:

.. code-block:: bash

        ./simba/web -i 123.456.789.10

where :code:`123.456.789.10` is your computer's IP address, which you can find (on Linux) using :code:`ifconfig`.

.. DANGER::
   Putting your results online allows other people to have interactive access with your database and your
   results!
   Unless you are on a private network, you should use the :code:`-s` option to run in "Safe Mode" - this
   disables the deleting/editing capabilities.

Once SimBA is online, other people on your network can view the page by going to :code:`123.456.789.10:5000` in their browser.


   

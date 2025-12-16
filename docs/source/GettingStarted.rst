.. toctree::
   :maxdepth: 2
   :caption: Contents:

.. _getting-started:

====================================
:fas:`forward;fa-fw` Getting Started
====================================



.. include:: ../../README.rst
 

System Install Scripts
----------------------

The Alamo system install scripts automate the setup of all necessary dependencies
for running Alamo on supported platforms, including Ubuntu, MacOS, and Windows. The
scripts are actively maintained in the repository under :file:`.github/workflows` and
leverage **GitHub Actions**, a continuous integration (CI) service that automatically
builds, tests, and deploys code whenever changes are pushed. For more information about
GitHub Actions, see `GitHub Actions <https://docs.github.com/en/actions>`_.

Every commit to the :code:`development` and :code:`master` branches triggers these
install scripts via GitHub Actions, ensuring that installation works correctly across
all supported platforms and that any issues introduced by recent changes are
immediately detected. Users can check the current status of these workflows on the
`Alamo GitHub Actions page <https://github.com/solidsgroup/alamo/actions>`_. Each page
in this documentation also includes **status badges** that indicate whether the latest
workflow for that platform and branch completed successfully, allowing users to
quickly verify the health of the scripts.

The scripts themselves include platform-specific workflows as well as shared workflows
for common setup tasks, and can be inspected or modified if custom behavior is needed.


.. toctree::
   :maxdepth: 1

   Install/Ubuntu
   Install/MacOS
   Install/Windows
           

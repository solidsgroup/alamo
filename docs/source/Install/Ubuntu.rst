.. _install_ubuntu:

:fab:`ubuntu;fa-b` Ubuntu
=========================

.. |ubuntu-badge| image:: https://github.com/solidsgroup/alamo/actions/workflows/install.yml/badge.svg
.. |ubuntu-badge-devel| image:: https://github.com/solidsgroup/alamo/actions/workflows/install.yml/badge.svg?branch=development
.. |ubuntu-badge-master| image:: https://github.com/solidsgroup/alamo/actions/workflows/install.yml/badge.svg?branch=master

+---------------+-----------------------+
| recent        | |ubuntu-badge|        |
+---------------+-----------------------+
| development   | |ubuntu-badge-devel|  |
+---------------+-----------------------+
| master        | |ubuntu-badge-master| |
+---------------+-----------------------+

Install dependencies
~~~~~~~~~~~~~~~~~~~~

.. tab-set::

    .. tab-item:: 24.04

        .. literalinclude:: ../../../.github/workflows/dependencies-ubuntu-24.04.sh
           :language: bash
           :lines: 2-

    .. tab-item:: 22.04

        .. literalinclude:: ../../../.github/workflows/dependencies-ubuntu-22.04.sh
           :language: bash
           :lines: 2-


Build code and run unit tests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: ../../../.github/workflows/build-and-test.sh
   :language: bash
   :lines: 2-

Build this documentation
~~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: ../../../.github/workflows/build-docs.sh
   :language: bash
   :lines: 2-
      

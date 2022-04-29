.. role:: cpp(code)
   :language: c++

===========
Conventions
===========

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Namespaces
==========

Classes and structures are organized using namespaces.
Namespaces should be brief but descriptive, and two nested namespaces should correspond to two nested subsets.
For instance :code:`Shape::Polygon::Triangle` is acceptable whereas :code:`Shape::Circle::Round` is not (since "round" is not a subset of "Circle").

The directory structure must reflect the namespace structure.
For instance :code:`Shape::Polygon::Triangle` must be declared in :code:`src/Shape/Polygon/Triangle.H`.

Regular and Pure Virtual classes
================================

Pure virtual classes must have the same name as the namespace in which they are contained.
For instance :code:`Model::Elastic::Elastic` is pure virtual, whereas :code:`Model::Elastic::Isotropic` must inherit from :code:`Model::Elastic::Elastic` and is not pure virtual.

Include Guards
==============

Include guards should be all uppercase, and should be the same as the path (relative to `src/`) with `/` replaced by underscores.
For instance: a file called `src/Integrator/Integrator.H` should have the include guard 

.. code-block:: cpp

   #ifndef INTEGRATOR_INTEGRATOR_H
   #define INTEGRATOR_INTEGRATOR_H
   ...
   #endif

Member and Argument Variable Names
==================================

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

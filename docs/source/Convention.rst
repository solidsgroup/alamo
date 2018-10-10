Conventions
===========

Namespaces
----------

Classes and structures are organized using namespaces.
Namespaces should be brief but descriptive, and two nested namespaces should correspond to two nested subsets.
For instance `Shape::Polygon::Triangle` is acceptable whereas `Shape::Circle::Round` is not (since "round" is not a subset of "Circle").

The directory structure must reflect the namespace structure.
For instance `Shape::Polygon::Triangle` must be declared in `src/Shape/Polygon/Triangle.H`.

Regular and Pure Virtual classes
--------------------------------

Pure virtual classes must have the same name as the namespace in which they are contained.
For instance `Model::Elastic::Elastic` is pure virtual, whereas `Model::Elastic::Isotropic` must inherit from `Model::Elastic::Elastic` and is not pure virtual.



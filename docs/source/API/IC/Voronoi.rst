Voronoi
-------
Create a series of patches using the Voronoi Tesselation

For a thorough discussion of the Voronoi tesselation algorithm, see
[this wikipedia article](https://en.wikipedia.org/wiki/Voronoi_diagram)

Let `N` be the number of patches in the tesselation, 
and let :math:`\alpha_{1},\ldots,\alpha_{n}` be a set of scalar values.
There are two methods for applying the tesselation IC:

  * Partition: This requires an `N`-component multifab. 
    (There is currently no check for this.)
    Each component of the fab corresponds to one of the patches zero elsewhere.
    (This is mostly useful for multiphase field microstructure)

  * Values: This uses a 1-component fab. It sets the value of the ith fab 
    equal to :math:`\alpha_{i}`.

 Note on periodicity: This IC uses the `amrex::Geometry` object to determine
 if it is periodic, and if so, it will generate a periodic tesselation.



.. doxygenclass:: IC::Voronoi
   :project: alamo
   :members:
   :protected-members:
   :private-members:
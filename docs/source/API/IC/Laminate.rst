
Laminate
--------

 Parameter List
 ~~~~~~~~~~~~~~

.. list-table::
    :header-rows: 1

    * - Parameter name
      - Type
      - Description
    * - :code:`number_of_inclusions`
      - Single value
      -  How many laminates (MUST be greater than or equal to 1). Default = 1
    * - :code:`center`
      - Array
      -  (x,y,[z]) values for the center point of the laminate
    * - :code:`thickness`
      - Array
      -  thickness of the laminate
    * - :code:`orientation`
      - Array
      -  Vector normal to the interface of the laminate
    * - :code:`eps`
      - Array
      -  Diffuse thickness
    * - :code:`mollifier`
      - Single value
      - None

.. doxygenclass:: IC::Laminate
   :project: alamo
   :members:
   :protected-members:
   :private-members:

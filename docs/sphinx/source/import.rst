Importing Wind Models
#####################

In Python it is possible to import a wind model from an arbitrary grid. Currently
supported by Python are 1D spherical grids, or 2D cylindrical or polar grids.

Using an Imported Grid
======================

In order for Python to be able to interpret and read in the provided grid, Python
expects the input files to be formatted in a certain standard. The input format
for the different coordinate systems is mostly standardised, with only small
differences relating to the coordinates for the different systems.

.. admonition :: Quantity Locations

    The grid coordinates and velocities of a cell are defined at the inner vertex
    of cells. But, the electron and radiation temperature, as well as the mass
    density are defined at the center of cells.

Required Parameters
-------------------

To use an imported model, when choosing a wind type, you must select `imported`,
e.g.,

`Wind.type(SV,star,hydro,corona,kwd,homologous,shell,imported) imported`

The coordinate system must then be specified using the parameter, e.g.,

`Wind.coord_system(spherical,cylindrical,polar,cyl_var) spherical`

Finally, the path to the model grid to read in must be specified using, e.g.,

`Wind.model2import 1d_model.import.txt`

Spherical Grids
---------------

Using a spherical coordinate system, a 1D spherically symmetric model can be
read into Python.

To read in a grid of this type, the following columns are required for each cell:

* i                        :  the element number for each cell
* :math:`r`                :  the radial coordinate in CGS
* :math:`v_{r}`            :  the radial velocity in CGS
* :math:`\rho`             :  the mass density in CGS
* :math:`T_{e}` (optional) :  the electron temperature in Kelvin
* :math:`T_{r}` (optional) :  the radiation temperature in Kelvin

.. admonition :: Grid Coordinates

    The radial coordinates of the cells must be constantly increasing in size.

Cylindrical Grids
-----------------

Using cylindrical coordinates, a 2.5D model can be read into Python.

.. admonition :: Grid Coordinates

    Note that the grid coordinates and the velocity is specified in Cartesian
    coordinates.

To read in a grid of this type, the following columns are required for each cell:

* i                        :  the i element number (row)
* j                        :  the j element number (column)
* inwind                   :  a flag indicating whether the cell is in the wind or not
* :math:`x`                :  the x coordinate in CGS
* :math:`z`                :  the z coordinate in CGS
* :math:`v_x`              :  the velocity in the x direction in CGS
* :math:`v_y`              :  the velocity in the y direction in CGS
* :math:`v_z`              :  the velocity in the z direction in CGS
* :math:`\rho`             :  the mass density in CGS
* :math:`T_{e}` (optional) :  the electron temperature in Kelvin
* :math:`T_{r}` (optional) :  the radiation temperature in Kelvin

.. admonition :: Unstructed/non-linear Grids

    In principle, it is possible to read in an unstructed or non-linear 
    cylindrical grid, i.e. where the cells are not regularly spaced, however,
    Python has been designed for structured grids with regular grid spacing, and
    as such there may be undefined behaviour for unstructed grids.    

Polar Grids
-----------

Using polar coordinates, a 2.5D model can be read into Python.

.. admonition :: Cartesian Velocity

    The velocity in for the polar grid is required to be in Cartesian
    coordinates due to conventions within the Python programming style. As such,
    any polar velocity components must first be projected into their Cartesian
    equivalent.


* i                        :  the i element number (row)
* j                        :  the j element number (column)
* inwind                   :  a flag indicating whether the cell is in the wind or not
* :math:`r`                :  the radial coordinate in CGS
* :math:`\theta`           :  the :math:`\theta` coordinate in degrees
* :math:`v_x`              :  the velocity in the x direction in CGS
* :math:`v_y`              :  the velocity in the y direction in CGS
* :math:`v_z`              :  the velocity in the z direction in CGS
* :math:`\rho`             :  the mass density in CGS
* :math:`T_{e}` (optional) :  the electron temperature in Kelvin
* :math:`T_{r}` (optional) :  the radiation temperature in Kelvin

.. admonition :: :math:`\theta`-cells

    The :math:`\theta` range should extend from at least 0 to 90°. It is possible
    to extend beyond 90°, but these cells should not be inwind.

.. admonition :: Unstructed/non-linear Grids

    In principle, it is possible to read in an unstructed or non-linear 
    polar grid, i.e. where the cells are not regularly spaced, however,
    Python has been designed for structured grids with regular grid spacing, and
    as such there may be undefined behaviour for unstructed grids.    

Guard Cells and Setting Values for `inwind` 
-------------------------------------------

The `inwind` flag is used to mark if a grid cell is either in the wind or not
in the wind. The following enumerator flags are used,

.. code :: c

    W_IGNORE      = -2   // ignore this grid cell
    W_NOT_INWIND  = -1   // this cell is not in the wind
    W_ALL_INWIND  =  0   // this cell is in the wind

Whilst it is possible to set in `inwind = 1` for a grid cell, that is that the
cell is partially in the wind, Python will instead set these cells with
`inwind = -2` and ignore these grid cells.

Spherical
^^^^^^^^^

Three guard cells are expected. One guard cell is expected at the inner edge of
wind and two are expected at the outer edge of the wind. Guard cells should still
have a velocity, but the mass density and temperatures should be zero. 

An example of a correctly formatted spherical grid is below.

+---+-------------------+-------------------+---------+------+
| i |                 r |                 v |     rho |   t_e|
+---+-------------------+-------------------+---------+------+
| 0 | 1208000000000000.0|  1100258151.526268|      0.0|   0.0|
+---+-------------------+-------------------+---------+------+
| 1 | 1236000000000000.0|  1100258151.526268| 7.41e-14| 40000|
+---+-------------------+-------------------+---------+------+
| 2 | 1263000000000000.0| 1124299782.0866106| 6.34e-14| 40000|
+---+-------------------+-------------------+---------+------+
| 3 | 1291000000000000.0|  1951614716.074871| 1.32e-15| 40000|
+---+-------------------+-------------------+---------+------+
| 4 | 1320000000000000.0|  1994041122.946064|      0.0|   0.0|
+---+-------------------+-------------------+---------+------+
| 5 | 1350000000000000.0| 2050609665.4409878|      0.0|   0.0|
+---+-------------------+-------------------+---------+------+

Cylindrical
^^^^^^^^^^^

For cylindrical grids, the outer boundaries of the wind should have two layers 
of  guard cells in the same way as the a spherical grid, as above. For these
cells, and all cells which do not make up the wind, an inwind value of -1 or -2 
should be set. 

.. figure:: images/import_cylindrical_inwind.png
    :width: 700px
    :align: center

    A colour plot of the inwind variable for the cv_standard.pf example. Here, a
    SV model is being imposed on a cylindrical coordinate grid.

Polar
^^^^^

For polar grids, the outer boundaries of the wind should have two layers of 
guard cells in the same way as the a spherical grid, as above. For these cells, 
and all cells which do not make up the wind, an inwind value of -1 or -2 should be set. 

In this example, the theta cells extend beyond 90°. But, as they are not inwind,
Python is happy to include these cells. For a stellar wind in polar coordinates,
these extra :math:`\theta` cells extending beyond 90° are required. 

.. figure:: images/import_polar_inwind.png
    :width: 700px
    :align: center

    A colour plot of the inwind variable for the rtheta.pf example. Here, a SV
    model is being imposed on an polar coordinate grid.

.. figure:: images/import_stellar_polar_inwind.png
    :width: 700px
    :align: center

    A colour plot of the inwind variable for a stellar wind imposed on a polar
    coordinate grid. Important to note is the "halo" of inwind = -1 cells 
    surrounding the inwind cells. The cells with inwind = 1 will be set to
    inwind = -2 when imported into Python and ignored.

Setting Wind Temperatures
-------------------------

Reading in a temperature is optional when importing a model. However, if one
temperature value for a cell is provided, then Python assumes that this is
the electron temperature and the radiation temperature will be initialised as,

.. math ::
    T_{r} = 1.1 T_{e}.

However, if two temperature values are provided for the cells, then the first
temperature will be assumed as being the electron temperature and the second
will be the radiation temperature.

If no temperature is provided with the imported model, then the radiation 
temperature will be initialised using the parameter, e.g.,

`Wind.t.init 40000`

The electron temperature is then initialised using the Lucy approximation,

.. math ::
    T_{e} = 0.9 T_{r}


Maximum and Minimum Wind Radius
--------------------------------

The maximum and minimum spherical extent of the wind is calculated automatically
by Python, and does not take into account guard cells when it is doing this.

Tools for creating Imported Grids
=================================

Some tools to convert Python `root.wind_save` files into models which can be
imported exist in `$PYTHON/py_progs` are are named,

* import_spherical.py
* import_cyl.py
* import_rtheta.py

Using these on the example parameter files can be a good way to figure out the
expected standard for imported model grids.

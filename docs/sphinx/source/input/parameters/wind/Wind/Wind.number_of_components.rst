Wind.number_of_components
=========================
While most simple description of a wind consist of a single region of space, Python can calculate
radiative transfer through more complicated structres, where one region of space is described with one
prescription and another region of space with a second prescription. For example, one might want to place
a disk atmoosphere between the disk and a wind.  This parameter describes the number of components (aka domains)
of the wind.

Type
  Integer

Values
  Greater than 0

File
  `python.c <https://github.com/agnwinds/python/blob/master/source/python.c>`_


Parent(s)
  * :ref:`System_type`: ``star``, ``binary``, ``agn``


Child(ren)
  * :ref:`Wind.t.init`

  * :ref:`Wind.coord_system`

  * :ref:`Diag.adjust_grid`

  * :ref:`Wind.radmax`

  * :ref:`Wind.filling_factor`

  * :ref:`Wind.dim.in.z_or_theta.direction`

  * :ref:`Wind.type`

  * :ref:`Wind.dim.in.x_or_r.direction`


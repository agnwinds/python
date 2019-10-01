Wind.filling_factor
===================
The volume filling factor of the outflow. The implementation
of clumping (microclumping) is described in
Matthews et al. (2016), 2016MNRAS.458..293M. Asked once per domain.

Type
  Double

Values
  0 < f <= 1, where 1 is a fully smooth wind.

File
  `setup_domains.c <https://github.com/agnwinds/python/blob/master/source/setup_domains.c>`_


Parent(s)
  * :ref:`Wind.number_of_components`: Greater than 0. Once per domain.



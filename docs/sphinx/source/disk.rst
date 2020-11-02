The Disk
########

The disk is normally treated as infinitely thin and defined by an inner boundary and an outer boundary.  It assumed to be in  Keplerian  rotation about 
the central object in the system.   The temperature distribution of the disk
is normally assumed to be that of a standard Shakura-Sunyaev disk, with a hard
boundary at its inner edge.   Options are provided for reading in a non-standard
temperature distribution.

An option is provide for a vertically extended disk, whose thickness increases
as with distance from the central object object.   

Vertically Extended disk (Details)
##################################

.. figure:: images/vertical_disk.png
    :width: 300px
    :align: center

In defining a vertically extended disk in the context of parameterized 
models, such as  KWD of SV, one needs to decide how to tranlated values from
a parameterized wind on a flat disk to a parameterized wind on verticallye extended
disk.   The choices we have made are (intended to be) as follows:

* The temperature of a vertically extended disk are given by the distance from the central object. 
* The magnitude of the poloidal velocity at the footpoint is the same as the flat 
  disk that underines it.
* The density at the base of the wind is defined as the same as the flat disk that underlies it.
* For the KWD disk, the streamlines are the streamlines that reflect the focus position.  
  For the SV model, the streamline direction is also determined by the footpoint on the flat disk.

(Note that the fact that we take most parameters from the flat disk at the same cylindrical distance from the central object, but the streamline directions from the footpoints introduces a small error
in the mass loss rate of the wind that is generated compared to specified mass loss rate.)

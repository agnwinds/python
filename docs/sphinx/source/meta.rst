Meta-documentation
##################

How to document Python
======================

This documentation is written in ReStructured Text, and parsed by Sphinx.
We're trying to maintain a roughly consistent format for the documentation.

Parameter documentation
=======================

Parameters are documented in a consistent way. They have a set of properties.
Not every parameter will have all properties but you should fill them all in where possible.
A full example outline is:

.. literalinclude:: reference_rst.txt
   :language: rst

The sections we expect are entered as a definition list.
A definition list consists of titles followed by a definition block indented by 2 characters.
The headings, in the order we expect, are:

Name
  The parameter name, as used by Python input files.

Description
  A description of the parameter and its function.
  This can include links to other pages and parameters, using the format

  .. literalinclude:: reference_rst.txt
     :language: rst
     :lines: 4-4

Type
  This is whether the parameter is an integer, float, or enumerator (a list of choices).

Unit
  This is the unit. It can be something like `cm`, `m` or even derived from other parameters
  (e.g. `Central_object.radius`).

Values
  If the parameter is an integer or float, this should describe the range of values it can take.
  For example, `Greater than 0` or `0-1`.

  If the variable type is `Enumerator`, then instead it should include a nested definition list of
  the possible choices. Where each choice implies a different set of possible children
  (e.g. :ref:`Wind.type``) then each choice should have its own Children definition list.

File
  The file the parameter is found in. This is a link to the file on the `master` branch.

Child(ren)
  If the parameter implies any others.
  For example, :ref:`Spectrum.no_observers` has child parameters :ref:`Spectrum.angle`.

Parent(s)
  If the parameter depends on another.
  For example, :ref:`KWD.rmax` is only required for a specific choice of :ref:`Wind.type`.




Old bit to update
=================

It will read the files in the *docs/parameters* directory and build a *HTML* page,
then automatically open it in your browser. You may see some errors during the creation::

   while scanning a block scalar
    in "/Users/amsys/python/docs/parameters/sv.diskmin.yaml", line 5, column 7
   expected chomping or indentation indicators, but found 'c'
    in "/Users/amsys/python/docs/parameters/sv.diskmin.yaml", line 5, column 8

This occurs during the *YAML* step when a file has an invalid entry, and that parameter file will be
skipped as a consequence. Typically it is cause when a text line is malformed. For inline text
(e.g. ``values: condition``) it occurs when the string starts with a non-alphanumeric character like
``>`` (so ```values: >0`` is invalid).

You may also see it when using block text (e.g. ``description:|``) if the following lines are not
correctly indented e.g.::

   description:|
   The text for values has not been indented correctly
   and so will throw an error.

   description:|
     This text has been indented
     and so will not throw an error

A common warning is::

   WARNING: Unknown target name:

This occurs during the *Sphinx* step when it cannot link to the parent for a variable. If a
parameter has a parent provided as a ``key:value`` pair e.g. ``reverb.type: Anything above 3``.

*Sphinx* will automatically try to link to ``reverb.type`` as its parent. If the name is misspelled
(e.g. ``rverb.type``) or the parameter simply doesn't exist any more, it will throw a warning and
generate a dead link.

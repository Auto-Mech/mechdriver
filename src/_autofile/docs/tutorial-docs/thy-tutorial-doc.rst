.. _thy-tutorial-doc:

Theory Tutorial
==========================

We can create a theory directory manager inside a given species
directory.

.. code-block:: python

    >>> pfx = fs[-1].path(['InChI=1S/H', 0, 2])
    >>> tfs = autofile.fs.theory(pfx)

The theory filesystem has only one layer, which can be accessed using either
``0`` or ``-1`` for the index, and takes method, basis, and orbital type as its
locator values.

.. code-block:: python

    >>> tfs[-1].create(['b3lyp', '6-31g*', 'U'])
    >>> tfs[-1].create(['b3lyp', '6-31g*', 'R'])
    >>> tfs[-1].path(['b3lyp', '6-31g*', 'U'])
    '/home/avcopan/SPC/H/YZCKVEUIGOORGS/0/2/UHFFFAOYSA-N/ezvlpJU'
    >>> tfs[-1].path(['b3lyp', '6-31g*', 'R'])
    '/home/avcopan/SPC/H/YZCKVEUIGOORGS/0/2/UHFFFAOYSA-N/ezvlpJR'

The theory directory manager allows for the reading and writing of various
files within a given directory. One does this through the file attribute.

.. code-block:: python

    >>> tfs[-1].file
    namespace(energy=<...>, geometry=<...>, hessian=<...>, zmatrix=<...>)

The file attribute is a namespace of several file I/O managers. I have cut out
the object identifiers above to make the printed value more readable, but they
are all ``autofile.system.model.DataSeriesFile objects``.

Tip: If you want a readable print-out of what the files are in a given layer,
you can use the following.

.. code-block:: python

    >>> tfs[-1].file.__dict__.keys()
    dict_keys(['energy', 'geometry', 'hessian', 'zmatrix'])

Otherwise, the files for each layer are also listed in the function docstrings
for this module.

As an example, let us do some I/O with an energy file.

First, we'll check that the file doesn't exist yet.

.. code-block:: python

    >>> tfs[-1].file.energy.exists(['b3lyp', '6-31g*', 'U'])
    False

Notice that we need the same three locators! The argument doesn't change.

Let's write a made-up energy value to the file.

.. code-block:: python

    >>> tfs[-1].file.energy.write(5.7, ['b3lyp', '6-31g*', 'U'])

Now the file exists.

.. code-block:: python

    >>> tfs[-1].file.energy.exists(['b3lyp', '6-31g*', 'U'])
    True

The path to this file is as follows.

.. code-block:: python

    >>> tfs[-1].file.energy.path(['b3lyp', '6-31g*', 'U'])
    '/home/avcopan/SPC/H/YZCKVEUIGOORGS/0/2/UHFFFAOYSA-N/ezvlpJU/geom.ene'

We can confirm that our made-up value was correctly stored by reading it back
out.

.. code-block:: python

    >>> tfs[-1].file.energy.read(['b3lyp', '6-31g*', 'U'])
    5.7


|
|
|

.. note::
    Move on to the next tutorial :ref:`cnf-tutorial-doc` to learn the conformer system

    Or return to the tutorial hub :ref:`tutorial-hub` to check out more tutorials

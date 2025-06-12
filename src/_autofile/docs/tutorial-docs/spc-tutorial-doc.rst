.. _spc-tutorial-doc:

Species Tutorial
==========================

To begin, let us create a species filesystem in the current directory, ``'.'``.

.. code-block:: python

    >>> import autofile
    >>> fs = autofile.fs.species('.')

``fs`` is simply a tuple of two directory managers (``DataSeries`` objects),
one for the trunk directory and one for the leaf directory.

.. code-block:: python

    >>> fs
    (DataSeries('/home/avcopan', species_trunk), DataSeries('/home/avcopan', species_leaf))

Since this is a tuple, the two data series can be accessed using indices.

.. code-block:: python

    >>> fs[0]
    DataSeries('/home/avcopan', species_trunk)
    >>> fs[1]
    DataSeries('/home/avcopan', species_leaf)
    >>> fs[-1]
    DataSeries('/home/avcopan', species_leaf)

In this case there are only two elements, so the leaf data series can be
accessed as either ``fs[1]`` or ``fs[-1]``.

These ``DataSeries`` objects has several methods, which we will explore in
turn: ``exists()``, ``create()``, ``path()``, and ``existing()``.
Let us first check whether the trunk directory exists.

.. code-block:: python

    >>> fs[0].exists([])
    False

Since it doesn't exist, create it.

.. code-block:: python

    >>> fs[0].create([])

Now, we will find that it does exist.

.. code-block:: python

    >>> fs[0].exists([])
    True

The empty list argument to each of these functions is the sequence of "locator
values" for accessing this directory.
The species trunk directory doesn't take any locator values, so the list is
empty.

We can print the path the this trunk directory as follows.

.. code-block:: python

    >>> fs[0].path([])
    '/home/avcopan/SPC'

Obviously, the path on your system will be different.

Now, we can create the species leaf directories, which go inside the trunk
directory. 
The manager for the leaf directories is the second (and final) element of
``fs``.

Let's create some directories for atoms.

.. code-block:: python

    >>> fs[-1].create(['InChI=1S/H', 0, 2])
    >>> fs[-1].create(['InChI=1S/He', 0, 1])
    >>> fs[-1].create(['InChI=1S/O', 0, 3])
    >>> fs[-1].create(['InChI=1S/O', 0, 1])

We can see that the species leaf directory takes three locator values: 1. the
inchi, 2. the charge, and 3. the multiplicity.
We need these three values every time we want to access the file for a
particular species.

If you wish to see which directories have already been created, you can can use
the ``DataSeries.existing()`` method to retrieve a full list.

.. code-block:: python

    >>> fs[-1].existing()
    (['InChI=1S/H', 0, 2], ['InChI=1S/He', 0, 1], ['InChI=1S/O', 0, 1], ['InChI=1S/O', 0, 3])

This method is useful for traversing a file system after it has been created.

Let's take a look at the paths for each leaf directory:

.. code-block:: python

    >>> fs[-1].path(['InChI=1S/H', 0, 2])
    '/home/avcopan/SPC/H/YZCKVEUIGOORGS/0/2/UHFFFAOYSA-N'
    >>> fs[-1].path(['InChI=1S/He', 0, 1])
    '/home/avcopan/SPC/He/SWQJXJOGLNCZEY/0/1/UHFFFAOYSA-N'
    >>> fs[-1].path(['InChI=1S/O', 0, 3])
    '/home/avcopan/SPC/O/QVGXLLKOCUKJST/0/3/UHFFFAOYSA-N'
    >>> fs[-1].path(['InChI=1S/O', 0, 1])
    '/home/avcopan/SPC/O/QVGXLLKOCUKJST/0/1/UHFFFAOYSA-N'

Note that there is no correspondence between the number of locators and the
number of directories.

|
|
|

.. note::
    Move on to the next tutorial :ref:`thy-tutorial-doc` to learn the theory system and begin
    writing and reading molecular data.

    Or return to the tutorial hub :ref:`tutorial-hub` to check out more tutorials

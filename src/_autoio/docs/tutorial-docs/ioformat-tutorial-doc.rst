.. _ioformat-tutorial-doc:

IOformat Tutorial
==========================

IOformat helps us keep consistant formatting throughout AutoMech. We can practice working with strings:

.. code-block:: python

    >>> import ioformat
    >>> mystring = 'I love AutoMech'
    >>> ioformat.indent(mystring, 4)
    '    I love AutoMech'
    >>> ioformat.addchar(mystring, ':) ', side='pre')
    ':) I love AutoMech'        
    >>> ioformat.addchar(mystring, '!', side='post')
    'I love AutoMech !'

ioformat can also be used to remove unwanted characters from strings:

.. code-block:: python

    >>> uglystring = 'AutoMech    '
    >>> ioformat.remove_whitespace_from_string(uglystring)
    'AutoMech'

Paths must be formatted very particularly and ioformat gives us that control:

.. code-block:: python

    >>> path = ioformat.pathtools.current_path() 
    >>> tutorial_path = 'automech-tutorial'
    >>> new_path = ioformat.pathtools.prepare_path((path, tutorial_path), make=True)
    >>> new_path
    '/gpfs/fs1/home/elliott/automech-tutorial'

ioformat optionally can create this path and allows us to go to it

.. code-block:: python

    >>> ioformat.pathtools.go_to(new_path)
    >>> ioformat.pathtools.write_file(mystring, new_path, 'myfile.txt') 
    >>> ioformat.pathtools.read_file(new_path, 'myfile.txt')
    'I love AutoMech'

|
|
|

.. note::
    Move on to the next tutorial :ref:`autoparse-tutorial-doc` to ...

    Or return to the tutorial hub :ref:`base-tutorial-hub` to check out more tutorials

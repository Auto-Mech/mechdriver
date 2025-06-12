.. _autoparse-tutorial-doc:

Autoparse Tutorial
==========================

Autoparse allows for user-friendly regex usage to parse for information.  Consider if you
had the following string:

.. code-block:: python

    >>> mystring = ' ___A_a_ * & b ___c_ C 1d 2 __e_ D__'

We can use autoparse to work with the string

.. code-block:: python

    >>> import autoparse
    >>> autoparse.find.split_words(mystring)
    ('___A_a_', '*', '&', 'b', '___c_', 'C', '1d', '2', '__e_', 'D__')

Or we can capture specific patterns

.. code-block:: python

    >>> pattern = autoparse.pattern.capturing(autoparse.pattern.UPPERCASE_LETTER)
    >>> autoparse.find.all_captures(pattern, mystring)
    ('A', 'C', 'D')
    >>> pattern = autoparse.pattern.capturing(autoparse.pattern.DIGIT)
    >>> autoparse.find.first_capture(pattern, mystring)
    '1'

We can also use autoparse to make sure a string is formatted as expected

.. code-block:: python

    >>> message = 'Greetings, user'
    >>> autoparse.find.starts_with('Greet', message) 
    True
    >>> autoparse.find.ends_with('goodbye', message) 
    False
    >>> mynumber = '  400 '
    >>> autoparse.find.is_number(mynumber)
    True

Now that you've got the basics, lets see how autoparse can help us parse real
data we'll encounter in AutoMech.  We are defining a search pattern that will 
parse the cartesian coordinates out of a file.

.. code-block:: python

    >>> atom_symbols = (                                         
    ...     autoparse.pattern.LETTER +                                  
    ...     autoparse.pattern.maybe(autoparse.pattern.LETTER)           
    ... )
    >>> number = autoparse.pattern.FLOAT                   
    >>> xyz_lines = autoparse.pattern.LINESPACES.join([     
    ...     autoparse.pattern.capturing(atom_symbols),      
    ...     autoparse.pattern.capturing(number),           
    ...     autoparse.pattern.capturing(number),           
    ...     autoparse.pattern.capturing(number),           
    ... ])


Now we can put in an example xyz string   

.. code-block:: python


    >>> xyz_string = """
    ... dummy line 1
    ... another dummy line
    ... 
    ... that dummy line was blank
    ... 6                    
    ... charge: 0, mult: 1                   
    ... F    1.584823  -0.748487  -0.427122  
    ... C    0.619220   0.190166  -0.271639  
    ... C   -0.635731  -0.183914  -0.180364  
    ... Cl  -1.602333   0.736678  -0.026051  
    ... H    0.916321   1.229946  -0.227127  
    ... H   -0.882300  -1.224388  -0.229636  
    ... let's end on a dummy line
    ... """                                  

    >>> autoparse.cast(autoparse.find.all_captures(xyz_lines, xyz_string))
    (('F', 1.584823, -0.748487, -0.427122), ('C', 0.61922, 0.190166, -0.271639), ('C', -0.635731, -0.183914, -0.180364), ('Cl', -1.602333, 0.736678, -0.026051), ('H', 0.916321, 1.229946, -0.227127), ('H', -0.8823, -1.224388, -0.229636))

|
|
|

.. note::
    Move on to the next tutorial :ref:`autoread-tutorial-doc` to ...

    Or return to the tutorial hub :ref:`base-tutorial-hub` to check out more tutorials

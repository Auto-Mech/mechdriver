
theory.dat
----------

Details various electronic structure methods.

Theory levels are given as a series of blocks like so:

General format::

    level <level_name>
            <keyword1 = value>
            <keyword2 = value>
              …
    end level

Above, <level_name> serves as a descriptor of the method that is used in other parts of the code. This string can be whatever the user wishes. It must contain no whitespace.

Any number of these level blocks can be given in theory.dat. The only blocks used by the code are those specified in run.dat and models.dat.


.. csv-table:: keywords for theory blocks
    :file: tables/thy_keys.csv
    :header-rows: 1
    :widths: 10, 10, 10, 10


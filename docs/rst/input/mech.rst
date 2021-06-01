
mechanism.dat
=============

List of all the reactions that AutoMech may run.

Currently, the input file takes the form of a CHEMKIN input file. The full details of the format can be found in the CHEMKIN manual.
However, a simplified version will suffice for running AutoMech. The file can take the form::

    REACTIONS
      Reac1+Reac2=Prod1+Prod2         1.00 0.00 0.00
      Reac3=Prod3                     1.00 0.00 0.00
      Reac4=Prod3+Prod5               1.00 0.00 0.00
      …
    END

Above, are names of chemical species undergoing the reaction. These names must correspond to the names given in the species.csv file.

The numbers included in the above example are included for the input parser to work (but are also ignored). Note that rate-constant fit parameters normally included in the CHEMKIN file are ignored.

The reactions do not need to be in any particular order, as AutoMech can sort the reactions into PESs itself.

The only reactions that will be run by the code, will be those specified in run.dat.

Commenting
----------

Reaction lines preceded by a ‘!’ are ignored by the parser. (maybe include #?)

In accordance with chemkin rules


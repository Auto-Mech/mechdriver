
Overview
========

MechAnalyzer
------------

MechAnalyzer is a collection of tools to process, analyze, and construct chemical reaction
mechanisms.

Many of the tasks can be considered pre-process and post-process tasks.

Herein is defined various Python objects for representing data structures useful. Also contains
useful functionality for handling these mechanism objects which are used by the workflow code.

Mechanalyzer is packaged with, and utilizes the ratefit and thermfit packages.

RateFit
-------

This package handles conversion between numerical rate constants calculated from master
equation simulations and functional form rate expressions used to represent the reactions in
mechanisms.

ThermFit
--------

This package handles the calculation of calculation of heats-of-formation which are used
alongside partition function to calculate the coefficients of NASA polynomials used to
represent the thermochem of species in mechanisms.

Most importantly, it also contains code to determine bases for thermochem via CBH methods.



Workflow
========

Hiearchy of Input
-----------------

For complex automated processes, a central challenge involves the code being able to decide whether
it should use information from a specified input or the database, this allows for grabbing user having the ability t

To simplify the logic of the process, we have instituted a strict hierarchy of input information

(1) User input
    If the user has specified a value or process in one of the MechDriver inputs, we will use that value
(2) Permament storage
    IF no input has been provided, then we search the save filesystem for the values
(3) Estimate a value
    If no value is available via input or fileystem we will use internal procedures to estimate what it could be
(4) Assume a default
    If there is no estimation procedure we assume the value is some default.

keep in mind this is also for certain situations where a value must be present in the save fileystem.
A bit different when certian drivers DEMAND a value to be there (ktp vs. es)


One minor exception is specifying models. 

If user specifies 1DHR model for a PES, we will build the MESS with all input for the hindered rotor potentials,
however if there are no HRs we just do harmonic freqs for species as this is not inconsistent treatment



Object Description
------------------


Kinetics
~~~~~~~~

rxn_ktp_dct
Description: rate constants for each reaction as a function of T and P

.. code-block:: python

    # Description
    {rxn1: ktp_dct1, rxn2: ...}

    # Type
    dict[tuple(tuple(str)): dict[]

    # Example
    {(('CH4', 'H'), ('CH3', 'H2'), ()): {}}

ktp_dct
Description: rate constants as a function of T and P

.. code-block:: python

    {pressure1: (temps1, kts1), pressure2: ...}

rxn_param_dct
description: reactions and their accompanying rate expressions

.. code-block:: python


    # Description
    {rxn1: (param_tup1, param_tup2, ...), rxn2: ...}

    # Example::
    {((‘H’, ‘O2’), (‘OH’, ‘O’), (None,)):
        ([1e15, 0, 15000], None, None, None, None, None),
        ([1e10, 0, 5000], None, None, None, None, None)}

.. code-block:: python

    # Example
    {(('H', 'O2'), ('OH', 'O'), (None,)):
        (([1E+15, 0.00, 25000], None, None, None,
            {0.1: [1E+15, 0.00, 25000],
            1.0: [1E+16, 0.00, 25000],
            10.0: [1E+17, 0.00, 25000],
            100.0: [1E+18, 0.00, 25000]}, None),
        ([1E+15, 0.00, 25000], None, None, None,
            {0.1: [1E+15, 0.00, 25000],
             1.0: [1E+16, 0.00, 25000],
           100.0: [1E+18, 0.00, 25000]}, None),)}


Thermo
~~~~~~

aligned_rxn_ktp_dct

.. code-block:: python

    {rxn1: [ktp_dct1, ktp_dct2, ...], rxn2: ...}

aligned_rxn_ratio_dct

.. code-block:: python

    {rxn1: [ratio_dct1, ratio_dct2, ...], rxn2: ...}

aligned_spc_thermo_dct

.. code-block:: python

    {spc1: [thermo_array1, thermo_array2, ...] , spc2: ...}

see thermo_array (clickable)

aligned_spc_diff_dct

.. code-block:: python

    {spc1: [diff_array1, diff_array2, ...] , spc2: ...}

nasa7_params

.. code-block:: python

    [NEED]

ratio_dct:
Description: similar structure to a ktp_dct, except give a ratio of k(T,P) values relative to another ktp_dct

.. code-block:: python

    {pressure1: (temps1, ratios1), pressure2: ...}

spc_nasa7_dct

.. code-block:: python

    {spc1: nasa7_params1, spc2: ...}

spc_thermo_dct

.. code-block:: python

    {spc1: thermo_array1, spc2: ...}

thermo_array 
each item is a 1xN numpy array 

.. code-block:: python

    [temps, h, cp, s, g]


Parameters
~~~~~~~~~~

param_tup 
Description: rate expression for a reaction

.. code-block:: python

    (highp _params, lowp_params, troe_params, cheb_dct, plog_dct, collider_dct)

highp_params
Description: Arrhenius parameters for the high-pressure limit
Note: highp_params should only ever contain a single Arrhenius expression

.. code-block:: python

    [A, n, Ea]

lowp_params
Description: Arrhenius parameters for the low-pressure limit. Only used for Lindemann and Troe expressions.
Note: lowp_params should only ever contain a single Arrhenius expression

.. code-block:: python

    [A, n, Ea]

troe_params
Description: Troe parameters

.. code-block:: python

    [alpha, T***, T*, T**]

cheb_dct
Description: Chebyshev parameters

.. code-block:: python

    {'t_limits': [tmin, tmax],
     'p_limits': [pmin, pmax],
     'alpha_elm': cheb_coeffs,
     'a_units': units of the output rate coefficient}

cheb_coeffs 
Description: Chebyshev polynomial coefficients
type: Numpy array 
shape (N, M), where N is the number of basis functions along the temperature axis and M is the number of basis functions along the pressure axis
Note: N, M is the same order that these parameters are defined in the Chemkin CHEB command
For N=2, M=3, this Numpy array would look like [[a,b,c], [d,e,f]]
format: Numpy_array[[coeff1, coeff2, ...], [...], ...]
units: the units of the output rate constant are given by the ‘a_units’ value, which is a str that can be either ‘moles’ or ‘molecules’ (see the cheb_dct entry)

plog_dct
Description: PLOG parameters
Note: for pressures with more than one Arrhenius expression, duplicates are described by multiple param_tuples (see the rxn_param_dct entry)

.. code-block:: python

    {pressure1: highp _params1, pressure2: ...}

units:
A is on a molar basis
n is relative to a reference temp of 1 Kelvin
Ea is in cal/mol
alpha is dimensionless
T***, T*, and T** are in Kelvin
T** is optional; it can either be omitted from the array or specified as None
a_units: the units of the output rate constant are given by the ‘a_units’ value, which is a str that can be either ‘moles’ or ‘molecules’


Physical Values
---------------

| Activation Energy
|     units = kcal/mol
|     float
| 
| Pressures
|     units = atmospheres
|     numpy array of shape (N,)
|
| Rate Constants (T-dependent)
|     units = mol, cm, s; values determined by molecularity of the reaction
|     numpy array of shape (N,)
|     
| Temperatures
|     units = Kelvin
|     numpy array of shape (N,)


Basic
~~~~~

spc
description: spc_name
type: str
format: spc
rxn
description: reaction name
type: tuple
format: ((rct1, rct2, ...), (prd1, prd2, ...), (third_bod1, third_bod2, ...))
each entry (e.g., rct1) is a species (see the spc entry)
the third bodies are confusing in that they have a ‘+’ in front and may also have ‘()’ enclosing them
can be a generic third body instead of a specific species: ‘+M’ or ‘(+M)’
spc_ident_dct
description: species and their accompanying chemically unique descriptions 
type: dct
format: {spc1: ident_dct1, spc2: ...}
ident_dct
description: chemically unique description of a spc
type: dct
format: {‘smiles’: SMILES, ‘inchi’: InChI, ‘inchikey’: InChI_key, ‘mult’: multiplicity, ‘charge’: charge, ‘sens’: sensitivity, ‘fml’: fml_dct}
SMILES is a str
InChI is a str
InChI_key is a str
multiplicity is an int
charge is an int
sensitivity is a float
fml_dct is a dct describing the chemical formula of a species
for example, for formaldehyde, the fml_dct would be {‘C’: 1, ‘H’: 2, ‘O’: 1}


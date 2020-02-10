""" lib of objects
"""

class PES():
    """ Potential Energy Surface
    """
    def __init__():
        self.formula
        self.pes_idx 
        self.sub_pes_idx_lst
        self.channels

    # Rates and Thermo Data
    self.temperatures
    self.pressures
    self.energy_transfer_params


class Reaction():
    """ Reaction objects
    """

    def __init__():
        # Molecular Species
        self.reacs
        self.prods
        self.nwells
        self.wells
        # Class
        self.reaction_class
        self.sadpt
        # Treatment
        self.model
        # Code Variables
        self.kickoff


class Species():
    """ Species objects
    """
    # Species ID
    self.name
    self.smiles
    self.inchi
    self.hind_inc
    # Molecular Info
    self.multiplicity
    self.charge
    self.geometry
    self.zmatrix
    self.tors_names
    self.elec_levels
    self.symmetry
    self.electronic_energy
    self.zpve
    # Mechanism Info
    self.sensitivity
    
    # Model info
    geo_lvl
    har_lvl
    tors_lvl

    # Code Stuff
    mc_nsamp

    # Filesystem

    class TS():
        """ 
        """
        self.dist_info
        self.frm_bnd_key
        self.brk_bnd_key
        self.grid
        self.saddle = False,True
        self.low_mult
        self.high_mult

    class Wells():



"""
  Read the save filesystem for all of the required information specified by
    (1) the models specified for partition function and
    (2) the electronic structure levels
  in order to write portions of MESS strings for species and reaction paths
  and calculate electronic and zero-point vibrational energies.
"""

import automol
import autofile
import autorun
from phydat import phycon
from mechanalyzer.inf import spc as sinfo
from mechlib import filesys
from mechlib.amech_io import printer as ioprinter
from mechroutines.models import ene
from mechroutines.models import typ
from mechroutines.models import etrans
from mechroutines.models import _rot as rot
from mechroutines.models import _tors as tors
from mechroutines.models import _symm as symm
from mechroutines.models import _vib as vib
from mechroutines.models import _flux as flux
from mechroutines.models import _util as util
from mechroutines.thermo import basis
import thermfit

# General readers
def read_spc_data(spc_dct, spc_name,
                  pes_mod_dct_i, spc_mod_dct_i,
                  run_prefix, save_prefix, chn_basis_ene_dct,
                  calc_chn_ene=True, spc_locs=None):
    """ Reads all required data from the SAVE filesystem for a given species.
        Also sets the writer for appropriately formatting the data into
        an MESS input file string.

        All of the data that is read is determined by the models that
        are described in the pes and spc model dictionaries.

        Info and basis species stored in dicts.

        :param spc_dct:
        :type spc_dct:
        :param spc_name: mechanism name of species
        :type spc_name: str
        :param pes_mod_dct_i: keyword dict of specific PES model
        :type pes_mod_dct_i: dict[]
        :param spc_mod_dct_i: keyword dict of specific species model
        :type spc_mod_dct_i: dict[]
        :param run_prefix: root-path to the run-filesystem
        :type run_prefix: str
        :param save_prefix: root-path to the save-filesystem
        :type save_prefix: str
        :param chn_basis_ene_dct: basis species <names> for mechanism species
        :type chn_basis_ene_dct: dict[]
        :rtype: (dict[], dict[])
    """

    ioprinter.obj('line_plus')
    ioprinter.reading(
        'filesystem info for {}'.format(spc_name), newline=1)

    vib_model = spc_mod_dct_i['vib']['mod']
    tors_model = spc_mod_dct_i['tors']['mod']
    spc_dct_i = spc_dct[spc_name]
    if typ.is_atom(spc_dct_i):
        inf_dct = atm_data(
            spc_dct, spc_name, pes_mod_dct_i, spc_mod_dct_i,
            run_prefix, save_prefix)
        writer = 'atom_block'
    else:
        if vib_model == 'tau' or tors_model == 'tau':
            inf_dct = tau_data(
                spc_dct_i, spc_mod_dct_i,
                run_prefix, save_prefix, saddle=False)
            writer = 'tau_block'
        else:
            inf_dct, chn_basis_ene_dct = mol_data(
                spc_name, spc_dct,
                pes_mod_dct_i, spc_mod_dct_i, chn_basis_ene_dct,
                run_prefix, save_prefix, calc_chn_ene=calc_chn_ene,
                spc_locs=spc_locs, zrxn=None)
            writer = 'species_block'

    # Add writer to inf dct
    inf_dct['writer'] = writer

    return inf_dct, chn_basis_ene_dct


def read_ts_data(spc_dct, tsname, rcts, prds,
                 pes_mod_dct_i, spc_mod_dct_i,
                 run_prefix, save_prefix, chn_basis_ene_dct):
    """ Reads all required data from the SAVE filesystem for a transition state.
        Also sets the writer for appropriately formatting the data into
        an MESS input file string.

        All of the data that is read is determined by the models that
        are described in the pes and spc model dictionaries.

        :param spc_dct:
        :type spc_dct:
        :param tsname: mechanism name of transition state
        :type tsname: str
        :param rcts: mechanism names of reactants connected to transition state
        :type rcts: tuple(str)
        :param prds: mechanism names of products connected to transition state
        :type prds: tuple(str)
        :param pes_mod_dct_i: keyword dict of specific PES model
        :type pes_mod_dct_i: dict[]
        :param spc_mod_dct_i: keyword dict of specific species model
        :type spc_mod_dct_i: dict[]
        :param run_prefix: root-path to the run-filesystem
        :type run_prefix: str
        :param save_prefix: root-path to the save-filesystem
        :type save_prefix: str
        :param chn_basis_ene_dct: basis species <names> for mechanism species
        :type chn_basis_ene_dct: dict[]
        :rtype: (dict[], dict[])
    """

    ioprinter.obj('line_plus')
    ioprinter.reading(
        'Reading filesystem info for {}'.format(tsname), newline=1)

    ts_dct = spc_dct[tsname]
    reac_dcts = [spc_dct[name] for name in rcts]
    prod_dcts = [spc_dct[name] for name in prds]

    ts_mod = spc_mod_dct_i['ts']
    ts_sadpt, ts_nobar = ts_mod['sadpt'], ts_mod['nobar']

    # Get all of the information for the filesystem
    if not automol.par.is_radrad(ts_dct['class']):

        # Set up the saddle point keyword
        sadpt = True
        search = ts_dct.get('ts_search')
        if search is not None:
            if 'vtst' in search:
                sadpt = False

        # Build MESS string for TS at a saddle point
        if ts_sadpt == 'pst':
            inf_dct = pst_data(
                ts_dct, reac_dcts,
                spc_mod_dct_i,
                run_prefix, save_prefix)
            writer = 'pst_block'
        elif ts_sadpt == 'rpvtst':
            inf_dct = rpvtst_data(
                ts_dct, reac_dcts, spc_mod_dct_i,
                run_prefix, save_prefix, sadpt=sadpt)
            writer = 'rpvtst_block'
        else:
            pf_filesystems = filesys.models.pf_filesys(
                spc_dct[tsname], spc_mod_dct_i,
                run_prefix, save_prefix, True, name=tsname)
            [cnf_fs, _, min_cnf_locs, _, _] = pf_filesystems['harm']
            cnf_path = cnf_fs[-1].path(min_cnf_locs)
            zma_fs = autofile.fs.zmatrix(cnf_path)
            zrxn = zma_fs[-1].file.reaction.read((0,))

            inf_dct, chn_basis_ene_dct = mol_data(
                tsname, spc_dct,
                pes_mod_dct_i, spc_mod_dct_i,
                chn_basis_ene_dct,
                run_prefix, save_prefix, zrxn=zrxn)
            writer = 'species_block'
    else:

        # Build MESS string for TS with no saddle point
        if ts_nobar == 'pst':
            if len(rcts) == 2:
                inf_dct = pst_data(
                    ts_dct, reac_dcts,
                    spc_mod_dct_i,
                    run_prefix, save_prefix)
            else:
                inf_dct = pst_data(
                    ts_dct, prod_dcts,
                    spc_mod_dct_i,
                    run_prefix, save_prefix)
            writer = 'pst_block'
        elif ts_nobar == 'rpvtst':
            inf_dct = rpvtst_data(
                ts_dct, reac_dcts,
                spc_mod_dct_i,
                run_prefix, save_prefix, sadpt=False)
            writer = 'rpvtst_block'
        elif ts_nobar == 'vrctst':
            inf_dct = flux_data(
                ts_dct, spc_mod_dct_i)
            writer = 'vrctst_block'

    # Add writer to inf dct
    inf_dct['writer'] = writer

    return inf_dct, chn_basis_ene_dct


# Data Readers
def atm_data(spc_dct, spc_name, pes_mod_dct_i, spc_mod_dct_i,
             run_prefix, save_prefix):
    """ Reads all required data from the SAVE filesystem for an atom.
        Stores data into an info dictionary.

        All of the data that is read is determined by the models that
        are described in the pes and spc model dictionaries.

        :param spc_dct:
        :type spc_dct:
        :param pes_mod_dct_i: keyword dict of specific PES model
        :type pes_mod_dct_i: dict[]
        :param spc_mod_dct_i: keyword dict of specific species model
        :type spc_mod_dct_i: dict[]
        :param run_prefix: root-path to the run-filesystem
        :type run_prefix: str
        :param save_prefix: root-path to the save-filesystem
        :type save_prefix: str
        :rtype: dict[]
    """

    spc_dct_i = spc_dct[spc_name]
    # Set up all the filesystem objects using models and levels
    pf_filesystems = filesys.models.pf_filesys(
        spc_dct_i, spc_mod_dct_i, run_prefix, save_prefix, False)

    ioprinter.info_message(
        'Obtaining the geometry...', newline=1)
    geom = rot.read_geom(pf_filesystems)

    ioprinter.info_message(
        'Obtaining the electronic energy...', newline=1)
    ene_chnlvl = ene.read_energy(
        spc_dct_i, pf_filesystems, spc_mod_dct_i,
        run_prefix, read_ene=True, read_zpe=False)

    hf0k, hf0k_trs, _, _ = basis.enthalpy_calculation(
        spc_dct, spc_name, ene_chnlvl,
        {}, pes_mod_dct_i, spc_mod_dct_i,
        run_prefix, save_prefix,
        pforktp='ktp', zrxn=None)

    # Create info dictionary
    inf_dct = {
        'geom': geom,
        'sym_factor': 1.0,
        'freqs': tuple(),
        'mess_hr_str': '',
        'mass': util.atom_mass(spc_dct_i),
        'elec_levels': spc_dct_i['elec_levels'],
        'ene_chnlvl': hf0k,
        'ene_reflvl': None,
        'ene_tsref': hf0k_trs,
        'zpe_chnlvl': None
    }

    return inf_dct


def mol_data(spc_name, spc_dct,
             pes_mod_dct_i, spc_mod_dct_i,
             chn_basis_ene_dct,
             run_prefix, save_prefix, calc_chn_ene=True, zrxn=None,
             spc_locs=None):
    """ Reads all required data from the SAVE filesystem for a molecule.
        Stores data into an info dictionary.

        All of the data that is read is determined by the models that
        are described in the pes and spc model dictionaries.

        :param spc_dct:
        :type spc_dct:
        :param pes_mod_dct_i: keyword dict of specific PES model
        :type pes_mod_dct_i: dict[]
        :param spc_mod_dct_i: keyword dict of specific species model
        :type spc_mod_dct_i: dict[]
        :param run_prefix: root-path to the run-filesystem
        :type run_prefix: str
        :param save_prefix: root-path to the save-filesystem
        :type save_prefix: str
        :rtype: dict[]
    """

    spc_dct_i = spc_dct[spc_name]
    print('spc_name test:', spc_name, spc_dct_i)
    ene_chnlvl = None
    ene_reflvl = None
    zpe = None
    hf0k = None
    hf0k_trs = None
    hf0k = None

    # Initialize all of the elements of the inf dct
    geom, sym_factor, freqs, imag, elec_levels = None, None, None, None, None
    allr_str, mdhr_dat = '', ''
    xmat, rovib_coups, rot_dists = None, None, None

    # Set up all the filesystem objects using models and levels
    pf_filesystems = filesys.models.pf_filesys(
        spc_dct_i, spc_mod_dct_i, run_prefix, save_prefix,
        zrxn is not None, name=spc_name, spc_locs=spc_locs)

    # Obtain rotation partition function information
    ioprinter.info_message(
        'Obtaining info for rotation partition function...', newline=1)
    geom = rot.read_geom(pf_filesystems)

    if typ.nonrigid_rotations(spc_mod_dct_i):
        rovib_coups, rot_dists = rot.read_rotational_values(pf_filesystems)

    # Obtain vibration partition function information
    ioprinter.info_message(
        'Preparing internal rotor info building partition functions...',
        newline=1)
    rotors = tors.build_rotors(
        spc_dct_i, pf_filesystems, spc_mod_dct_i)
    ioprinter.info_message(
        'Obtaining the vibrational frequencies and zpves...', newline=1)
    freqs, imag, zpe, _, tors_strs, _, _, _ = vib.full_vib_analysis(
        spc_dct_i, pf_filesystems, spc_mod_dct_i,
        run_prefix, zrxn=zrxn)
    allr_str = tors_strs[0]

    # ioprinter.info_message('zpe in mol_data test:', zpe)
    if typ.anharm_vib(spc_mod_dct_i):
        xmat = vib.read_anharmon_matrix(pf_filesystems)

    # Obtain symmetry factor
    ioprinter.info_message(
        'Determining the symmetry factor...', newline=1)

    zma = None
    if zrxn:
        [_, cnf_save_path, _, _, _] = pf_filesystems['harm']
        # Build the rotors
        if cnf_save_path:
            zma_fs = autofile.fs.zmatrix(cnf_save_path)
            zma = zma_fs[-1].file.zmatrix.read([0])

    sym_factor = symm.symmetry_factor(
        pf_filesystems, spc_mod_dct_i, spc_dct_i, rotors, grxn=zrxn, zma=zma)

    # Obtain electronic energy levels
    elec_levels = spc_dct_i['elec_levels']

    # Obtain energy levels
    ioprinter.info_message(
        'Obtaining the electronic energy + zpve...', newline=1)
    if calc_chn_ene:
        chn_ene = ene.read_energy(
            spc_dct_i, pf_filesystems, spc_mod_dct_i,
            run_prefix, read_ene=True, read_zpe=False, saddle=zrxn is not None)
        ene_chnlvl = chn_ene + zpe

        zma = None
        # Determine info about the basis species used in thermochem calcs
        hf0k, hf0k_trs, chn_basis_ene_dct, _ = basis.enthalpy_calculation(
            spc_dct, spc_name, ene_chnlvl,
            chn_basis_ene_dct, pes_mod_dct_i, spc_mod_dct_i,
            run_prefix, save_prefix, zrxn=zrxn)

    ene_reflvl = None

    #  Build the energy transfer section strings
    if zrxn is None:
        ioprinter.info_message(
            'Determining energy transfer parameters...', newline=1)
        well_info = sinfo.from_dct(spc_dct_i)
        # ioprinter.debug_message('well_inf', well_info)
        # bath_info = ['InChI=1S/N2/c1-2', 0, 1]  # how to do...
        bath_info = ['InChI=1S/Ar', 0, 1]  # how to do...
        etrans_dct = etrans.build_etrans_dct(spc_dct_i)

        edown_str, collid_freq_str = etrans.make_energy_transfer_strs(
            well_info, bath_info, etrans_dct)
    else:
        edown_str, collid_freq_str = None, None

    # Create info dictionary
    keys = ['geom', 'sym_factor', 'freqs', 'imag', 'elec_levels',
            'mess_hr_str', 'mdhr_dat',
            'xmat', 'rovib_coups', 'rot_dists',
            'ene_chnlvl', 'ene_reflvl', 'zpe_chnlvl', 'ene_tsref',
            'edown_str', 'collid_freq_str']
    vals = [geom, sym_factor, freqs, imag, elec_levels,
            allr_str, mdhr_dat,
            xmat, rovib_coups, rot_dists,
            hf0k, ene_reflvl, zpe, hf0k_trs,
            edown_str, collid_freq_str]
    inf_dct = dict(zip(keys, vals))

    return inf_dct, chn_basis_ene_dct


# VRCTST
def flux_data(ts_dct, spc_mod_dct_i):
    """ Read a VRC-TST flux file from the SAVE filesystem.

        :param ts_dct: species dict entry for a transition state
        :type ts_dct: dict[]
        :param spc_mod_dct_i: keyword dict of specific species model
        :type spc_mod_dct_i: dict[]
        :rtype: dict[str: str]
    """

    # Read the flux file from the filesystem
    _, ts_save_path, _, _ = filesys.models.set_rpath_filesys(
        ts_dct, spc_mod_dct_i['rpath'][1])

    flux_str = flux.read_flux(ts_save_path)

    # Create info dictionary
    inf_dct = {'flux_str': flux_str}

    return inf_dct


# VTST
def rpvtst_data(ts_dct, reac_dcts, spc_mod_dct_i,
                run_prefix, save_prefix, sadpt=False):
    """ Pull all of the neccessary information from the
        filesystem for a species
    """

    zrxn = ts_dct['zrxn']

    # Set up all the filesystem objects using models and levels
    if sadpt:
        # Set up filesystems and coordinates for saddle point
        # Scan along RxnCoord is under THY/TS/CONFS/cid/Z
        pf_filesystems = filesys.models.pf_filesys(
            ts_dct, spc_mod_dct_i, run_prefix, save_prefix, True)
        tspaths = pf_filesystems['harm']
        [_, cnf_save_path, min_locs, _, cnf_run_fs] = tspaths
        ts_run_path = cnf_run_fs[-1].path(min_locs)

        # Set TS reaction coordinate
        frm_name = 'IRC'
        scn_vals = filesys.models.get_rxn_scn_coords(cnf_save_path, frm_name)
        scn_vals.sort()
        scn_ene_info = spc_mod_dct_i['ene'][1][0][1]  # fix to be ene lvl
        scn_prefix = cnf_save_path
    else:
        # Set up filesystems and coordinates for reaction path
        # Scan along RxnCoord is under THY/TS/Z
        tspaths = filesys.models.set_rpath_filesys(
            ts_dct, spc_mod_dct_i['rpath'][1])
        ts_run_path, _, _, thy_save_path = tspaths

        # Set TS reaction coordinate
        scn_vals = filesys.models.get_rxn_scn_coords(thy_save_path, frm_name)
        scn_vals.sort()
        scn_ene_info = spc_mod_dct_i['rpath'][1][0]
        scn_prefix = thy_save_path

    # Modify the scn thy info
    ioprinter.debug_message('scn thy info', scn_ene_info)
    ioprinter.info_message('scn vals', scn_vals)
    mod_scn_ene_info = filesys.inf.modify_orb_restrict(
        filesys.inf.get_spc_info(ts_dct), scn_ene_info)
    # scn thy info [[1.0, ['molpro2015', 'ccsd(t)', 'cc-pvdz', 'RR']]]

    # Need to read the sp vals along the scan. add to read
    ref_ene = 0.0
    enes, geoms, grads, hessians, _, _ = filesys.read.potential(
        [frm_name], [scn_vals],
        scn_prefix,
        mod_scn_ene_info, ref_ene,
        constraint_dct=None,   # No extra frozen treatments
        read_geom=True,
        read_grad=True,
        read_hess=True)
    script_str = autorun.SCRIPT_DCT['projrot']
    freqs = autorun.projrot.pot_frequencies(
        script_str, geoms, grads, hessians, ts_run_path)

    # Get the energies and zpes at R_ref
    if not sadpt:
        _, ene_hs_sr_ref, ene_hs_mr_ref = ene.rpath_ref_idx(
            ts_dct, scn_vals, frm_name, scn_prefix,
            spc_mod_dct_i['ene'],
            spc_mod_dct_i['rpath'][1])
    fr_idx = len(scn_vals) - 1
    zpe_ref = (sum(freqs[(fr_idx,)]) / 2.0) * phycon.WAVEN2KCAL

    # Get the reactants and infinite seperation energy
    reac_ene = 0.0
    ene_hs_sr_inf = 0.0
    for dct in reac_dcts:
        pf_filesystems = filesys.models.pf_filesys(
            dct, spc_mod_dct_i, run_prefix, save_prefix, False)
        new_spc_dct_i = {
            'ene': spc_mod_dct_i['ene'],
            'harm': spc_mod_dct_i['harm'],
            'tors': spc_mod_dct_i['tors']
        }
        reac_ene += ene.read_energy(
            dct, pf_filesystems, new_spc_dct_i, run_prefix,
            read_ene=True, read_zpe=True, saddle=sadpt)

        ioprinter.debug_message('rpath', spc_mod_dct_i['rpath'][1])
        new_spc_dct_i = {
            'ene': ['mlvl', [[1.0, spc_mod_dct_i['rpath'][1][2]]]],
            'harm': spc_mod_dct_i['harm'],
            'tors': spc_mod_dct_i['tors']
        }
        ene_hs_sr_inf += ene.read_energy(
            dct, pf_filesystems, new_spc_dct_i, run_prefix,
            read_ene=True, read_zpe=False)

    # Scale the scn values
    if sadpt:
        scn_vals = [val / 100.0 for val in scn_vals]
    # scn_vals = [val * phycon.BOHR2ANG for val in scn_vals]

    # Grab the values from the read
    inf_dct = {}
    inf_dct['rpath'] = []
    pot_info = zip(scn_vals, enes.values(), geoms.values(), freqs.values())
    for rval, pot, geo, frq in pot_info:

        # Scale the r-values

        # Get the relative energy (edit for radrad scans)
        zpe = (sum(frq) / 2.0) * phycon.WAVEN2KCAL
        if sadpt:
            zero_ene = (pot + zpe) * phycon.KCAL2EH
        else:
            ioprinter.debug_message('enes')
            ioprinter.debug_message('reac ene', reac_ene)
            ioprinter.debug_message('hs sr', ene_hs_sr_ref)
            ioprinter.debug_message('inf', ene_hs_sr_inf)
            ioprinter.debug_message('hs mr', ene_hs_mr_ref)
            ioprinter.debug_message('pot R', pot * phycon.KCAL2EH)
            ioprinter.debug_message('zpe', zpe)
            ioprinter.debug_message('zpe ref', zpe_ref)

            elec_ene = (
                ene_hs_sr_ref - ene_hs_sr_inf -
                ene_hs_mr_ref + pot * phycon.KCAL2EH
            )
            zpe_pt = zpe - zpe_ref
            zero_ene = reac_ene + (elec_ene + zpe_pt * phycon.KCAL2EH)
            ioprinter.info_message('elec ene', elec_ene)
            ioprinter.info_message('zero ene', zero_ene)

        # ENE
        # ene = (reac_ene +
        #        ene_hs_sr(R_ref) - ene_hs_sr(inf) +
        #        ene_ls_mr(R_ref) - ene_hs_mr(R_ref) +
        #        ene_ls_mr(R) - ene_ls_mr(R_ref))
        # ene = (reac_ene +
        #        ene_hs_sr(R_ref) - ene_hs_sr(inf) -
        #        ene_hs_mr(R_ref) + ene_ls_mr(R))
        # inf_sep_ene = reac_ene + hs_sr_ene - hs_mr_ene
        # inf_sep_ene_p = (reac_ene +
        #                  hs_sr_ene(R_ref) - ene_hs_sr(inf) +
        #                  ls_mr_ene(R_ref) - hs_mr_ene(R_ref))
        # ene = inf_sep_ene_p + ene_ls_mr(R) - ene_ls_mr(R_ref)
        # ZPE
        # zpe = zpe(R) - zpe(inf)
        # or
        # zpe = zpe_ls_mr(R) - zpe_ls_mr(R_ref)

        # Set values constant across the scan
        elec_levels = ts_dct['elec_levels']

        # Create info dictionary and append to lst
        keys = ['rval', 'geom', 'freqs', 'elec_levels', 'ene_chnlvl']
        vals = [rval, geo, frq, elec_levels, zero_ene]
        inf_dct['rpath'].append(dict(zip(keys, vals)))

    # Calculate and store the imaginary mode
    if sadpt:
        _, imag, _ = vib.read_harmonic_freqs(
            pf_filesystems, run_prefix, zrxn=zrxn)
        ts_idx = scn_vals.index(0.00)
    else:
        imag = None
        ts_idx = 0
    inf_dct.update({'imag': imag})
    inf_dct.update({'ts_idx': ts_idx})

    return inf_dct


# PST
def pst_data(ts_dct, reac_dcts,
             spc_mod_dct_i, run_prefix, save_prefix):
    """ Determine the parameters used to set a phase space theory (PST)
        description of transition states.

        Simply reads a geometries for reactants and uses it along
        with various PST conditions from the user to set the n and Cn
        paramater for the Vn = Cn/R^n potential model for PST in MESS.

        :param ts_dct: species dict entry for a transition state
        :type ts_dct: dict[]
        :param reac_dcts: species dict entries for connected reactants
        :type reac_dcts:
        :param spc_mod_dct_i: keyword dict of specific species model
        :type spc_mod_dct_i: dict[]
        :param run_prefix: root-path to the run-filesystem
        :type run_prefix: str
        :param save_prefix: root-path to the save-filesystem
        :type save_prefix: str

    """

    # Get the k(T), T, and n values to get a Cn
    kt_pst = ts_dct.get('kt_pst', 4.0e-10)  # cm3/s
    temp_pst = ts_dct.get('temp_pst', 300.0)  # K
    n_pst = ts_dct.get('n_pst', 6.0)  # unitless

    ioprinter.info_message(
        'Determining parameters for Phase Space Theory (PST)',
        'treatment that yields k({} K) = {}'.format(temp_pst, kt_pst),
        newline=1)
    ioprinter.info_message(
        'Assuming PST model potential V = C0 / R^{}'.format(n_pst),
        indent=1)

    # Obtain the reduced mass of the reactants
    ioprinter.info_message(
        'Reading reactant geometries to obtain reduced mass...', newline=1)
    geoms = []
    for dct in reac_dcts:
        pf_filesystems = filesys.models.pf_filesys(
            dct, spc_mod_dct_i, run_prefix, save_prefix, False)
        geoms.append(rot.read_geom(pf_filesystems))
    mred = automol.geom.reduced_mass(geoms[0], geoms[1])

    cn_pst = automol.reac.pst_cn(kt_pst, n_pst, mred, temp_pst)

    # Create info dictionary
    keys = ['n_pst', 'cn_pst']
    vals = [n_pst, cn_pst]
    inf_dct = dict(zip(keys, vals))

    return inf_dct


# TAU
def tau_data(spc_dct_i,
             spc_mod_dct_i,
             run_prefix, save_prefix, saddle=False):
    """ Read the filesystem to get information for TAU
    """

    # Set up all the filesystem objects using models and levels
    pf_filesystems = filesys.models.pf_filesys(
        spc_dct_i, spc_mod_dct_i, run_prefix, save_prefix, saddle)
    [harm_cnf_fs, _,
     harm_min_locs, harm_save, _] = pf_filesystems['harm']
    # [tors_cnf_fs, _, tors_min_locs, _, _] = pf_filesystems['tors']

    # Get the conformer filesys for the reference geom and energy
    if harm_min_locs:
        geom = harm_cnf_fs[-1].file.geometry.read(harm_min_locs)
        min_ene = harm_cnf_fs[-1].file.energy.read(harm_min_locs)

    # Set the filesystem
    tau_save_fs = autofile.fs.tau(harm_save)

    # Get the rotor info
    rotors = tors.build_rotors(spc_dct_i, pf_filesystems, spc_mod_dct_i)

    run_path = filesys.models.make_run_path(pf_filesystems, 'tors')
    tors_strs = tors.make_hr_strings(rotors)
    [_, hr_str, flux_str, prot_str, _] = tors_strs

    # Use model to determine whether to read grads and hessians
    vib_model = spc_mod_dct_i['vib']['mod']
    freqs = ()
    harm_zpve = 'fix broken'
    proj_zpve = 'fix broken'
    # _, _, proj_zpve, harm_zpve = vib.tors_projected_freqs_zpe(
    #     pf_filesystems, hr_str, prot_str, run_prefix, zrxn=None)
    zpe_chnlvl = proj_zpve * phycon.EH2KCAL

    # Set reference energy to harmonic zpve
    db_style = 'directory'
    reference_energy = harm_zpve * phycon.EH2KCAL
    if vib_model == 'tau':
        if db_style == 'directory':
            tau_locs = [locs for locs in tau_save_fs[-1].existing()
                        if tau_save_fs[-1].file.hessian.exists(locs)]
        elif db_style == 'jsondb':
            tau_locs = [locs for locs in tau_save_fs[-1].json_existing()
                        if tau_save_fs[-1].json.hessian.exists(locs)]
    else:
        if db_style == 'directory':
            tau_locs = tau_save_fs[-1].existing()
        elif db_style == 'jsondb':
            tau_locs = tau_save_fs[-1].json_existing()

    # Read the geom, ene, grad, and hessian for each sample
    samp_geoms, samp_enes, samp_grads, samp_hessians = [], [], [], []
    for locs in tau_locs:

        # ioprinter.info_message('Reading tau info at path {}'.format(
        #     tau_save_fs[-1].path(locs)))

        if db_style == 'directory':
            geo = tau_save_fs[-1].file.geometry.read(locs)
        elif db_style == 'jsondb':
            geo = tau_save_fs[-1].json.geometry.read(locs)

        geo_str = autofile.data_types.swrite.geometry(geo)
        samp_geoms.append(geo_str)

        if db_style == 'directory':
            tau_ene = tau_save_fs[-1].file.energy.read(locs)
        elif db_style == 'jsondb':
            tau_ene = tau_save_fs[-1].json.energy.read(locs)
        rel_ene = (tau_ene - min_ene) * phycon.EH2KCAL
        ene_str = autofile.data_types.swrite.energy(rel_ene)
        samp_enes.append(ene_str)

        if vib_model == 'tau':
            if db_style == 'directory':
                grad = tau_save_fs[-1].file.gradient.read(locs)
            elif db_style == 'jsondb':
                grad = tau_save_fs[-1].json.gradient.read(locs)
            grad_str = autofile.data_types.swrite.gradient(grad)
            samp_grads.append(grad_str)

            if db_style == 'directory':
                hess = tau_save_fs[-1].file.hessian.read(locs)
            elif db_style == 'jsondb':
                hess = tau_save_fs[-1].json.hessian.read(locs)
            hess_str = autofile.data_types.swrite.hessian(hess)
            samp_hessians.append(hess_str)

    # Read a geometry, grad, and hessian for a reference geom if needed
    ref_geom, ref_grad, ref_hessian = [], [], []
    if vib_model != 'tau':

        # Get harmonic filesystem information
        [harm_save_fs, _, harm_min_locs, _, _] = pf_filesystems['harm']

        # Read the geometr, gradient, and Hessian
        geo = harm_save_fs[-1].file.geometry.read(harm_min_locs)
        geo_str = autofile.data_types.swrite.geometry(geo)
        ref_geom.append(geo_str)

        grad = harm_save_fs[-1].file.gradient.read(harm_min_locs)
        grad_str = autofile.data_types.swrite.gradient(grad)
        ref_grad.append(grad_str)

        hess = harm_save_fs[-1].file.hessian.read(harm_min_locs)
        hess_str = autofile.data_types.swrite.hessian(hess)
        ref_hessian.append(hess_str)

    # Obtain symmetry factor
    ioprinter.info_message('Determining the symmetry factor...', newline=1)
    sym_factor = symm.symmetry_factor(
        pf_filesystems, spc_mod_dct_i, spc_dct_i, rotors,
    )

    # Create info dictionary
    keys = ['geom', 'sym_factor', 'elec_levels', 'freqs', 'flux_mode_str',
            'samp_geoms', 'samp_enes', 'samp_grads', 'samp_hessians',
            'ref_geom', 'ref_grad', 'ref_hessian',
            'zpe_chnlvl', 'reference_energy']
    vals = [geom, sym_factor, spc_dct_i['elec_levels'], freqs, flux_str,
            samp_geoms, samp_enes, samp_grads, samp_hessians,
            ref_geom, ref_grad, ref_hessian,
            zpe_chnlvl, reference_energy]
    inf_dct = dict(zip(keys, vals))

    return inf_dct

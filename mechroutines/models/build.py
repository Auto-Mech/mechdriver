"""
  Read the save filesystem for all of the required information specified by
    (1) the models specified for partition function and
    (2) the electronic structure levels
  in order to write portions of MESS strings for species and reaction paths
  and calculate electronic and zero-point vibrational energies.
"""

import os
import automol
import elstruct
import autofile
import autorun
from phydat import phycon
from mechanalyzer.inf import spc as sinfo
from mechanalyzer.inf import thy as tinfo
from mechlib import filesys
from mechlib.amech_io import printer as ioprinter
from mechlib.amech_io._path import job_path
from mechroutines.models import ene
from mechroutines.models import typ
from mechroutines.models import etrans
from mechroutines.models import _rot as rot
from mechroutines.models import _tors as tors
from mechroutines.models import _symm as symm
from mechroutines.models import _vib as vib
from mechroutines.models import _util as util
from mechroutines.thermo import basis
from mechroutines.es.ts import ts_zma_locs
# import thermfit


# General readers
def read_spc_data(spc_dct, spc_name,
                  pes_mod_dct_i, spc_mod_dct_i,
                  run_prefix, save_prefix, chn_basis_ene_dct,
                  calc_chn_ene=True,
                  calc_ene_trans=True,
                  spc_locs=None):
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
    ioprinter.reading(f'filesystem info for {spc_name}', newline=1)
    vib_model = spc_mod_dct_i['vib']['mod']
    tors_model = spc_mod_dct_i['tors']['mod']
    spc_dct_i = spc_dct[spc_name]
    if typ.is_atom(spc_dct_i):
        inf_dct = atm_data(
            spc_dct, spc_name,
            pes_mod_dct_i, spc_mod_dct_i,
            run_prefix, save_prefix)
        writer = 'atom_block'
    else:
        if vib_model == 'tau' or 'tau' in tors_model:
            inf_dct = tau_data(
                spc_dct_i,
                spc_mod_dct_i,
                run_prefix, save_prefix, saddle=False)
            writer = 'tau_block'
        else:
            inf_dct, chn_basis_ene_dct = mol_data(
                spc_name, spc_dct,
                pes_mod_dct_i, spc_mod_dct_i, chn_basis_ene_dct,
                run_prefix, save_prefix,
                calc_chn_ene=calc_chn_ene,
                calc_ene_trans=calc_ene_trans,
                spc_locs=spc_locs, zrxn=None)
            writer = 'species_block'

    # Add writer to inf dct
    inf_dct['writer'] = writer

    return inf_dct, chn_basis_ene_dct


def read_ts_data(spc_dct, tsname, rcts, prds,
                 pes_mod_dct_i, spc_mod_dct_i,
                 run_prefix, save_prefix, chn_basis_ene_dct,
                 spc_locs=None):
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
    ioprinter.reading(f'Reading filesystem info for {tsname}', newline=1)
    ts_dct = spc_dct[tsname]
    reac_dcts = [spc_dct[name] for name in rcts]
    prod_dcts = [spc_dct[name] for name in prds]

    ts_mod = spc_mod_dct_i['ts']
    ts_sadpt, ts_nobar = ts_mod['sadpt'], ts_mod['nobar']
    # Override the models specification with species.dat, if avail
    search = ts_dct.get('ts_search')
    if search is not None:
        if search == 'rpvtst':
            ts_sadpt, ts_nobar = 'rpvtst', 'rpvtst'
        elif search == 'vrctst':
            ts_sadpt, ts_nobar = 'vrctst', 'vrctst'
        elif search == 'pst':
            ts_sadpt, ts_nobar = 'pst', 'pst'
        elif search == 'sadpt':
            ts_sadpt, ts_nobar = 'sadpt', 'sadpt'

    # Get all of the information for the filesystem
    rxn_class = ts_dct['class']
    if automol.ReactionInfo.is_radical_radical(rxn_class) and automol.ReactionInfo.is_low_spin(rxn_class):

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
            inf_dct, chn_basis_ene_dct = rpvtst_nobar(
                tsname, spc_dct, reac_dcts,
                pes_mod_dct_i, spc_mod_dct_i,
                chn_basis_ene_dct,
                run_prefix, save_prefix)
            writer = 'rpvtst_block'
        elif ts_nobar == 'vrctst':
            inf_dct = flux_data(
                ts_dct, spc_mod_dct_i)
            writer = 'vrctst_block'
    else:  # write the TS for a well-defined saddle point

        # Set up the saddle point keyword
        if search is not None:
            if 'vtst' in search:
                ts_sadpt = False

        # Build MESS string for TS at a saddle point
        if ts_sadpt == 'pst':
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
        elif ts_sadpt == 'rpvtst':
            inf_dct, chn_basis_ene_dct = rpvtst_sadpt(
                tsname, spc_dct,
                pes_mod_dct_i, spc_mod_dct_i,
                chn_basis_ene_dct,
                run_prefix, save_prefix)
            writer = 'rpvtst_block'
        else:
            print("Obtaining a ZRXN object from conformer any TS, "
                  "shouldn't matter")
            print('-----')
            pf_filesystems = filesys.models.pf_filesys(
                spc_dct[tsname], spc_mod_dct_i,
                run_prefix, save_prefix, True, name=tsname)
            [cnf_fs, _, min_cnf_locs, _, _] = pf_filesystems['harm']
            cnf_path = cnf_fs[-1].path(min_cnf_locs)
            zma_fs = autofile.fs.zmatrix(cnf_path)
            zma_locs = ts_zma_locs(spc_dct, tsname, zma_fs)
            zrxn = zma_fs[-1].file.reaction.read(zma_locs)
            print('-----')

            inf_dct, chn_basis_ene_dct = mol_data(
                tsname, spc_dct,
                pes_mod_dct_i, spc_mod_dct_i,
                chn_basis_ene_dct,
                run_prefix, save_prefix,
                calc_ene_trans=False,
                zrxn=zrxn, spc_locs=spc_locs)
            writer = 'species_block'
    # Add writer to inf dct
    inf_dct['writer'] = writer

    return inf_dct, chn_basis_ene_dct


# Data Readers
def atm_data(spc_dct, spc_name,
             pes_mod_dct_i, spc_mod_dct_i,
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
        spc_dct, run_prefix, read_ene=True, read_zpe=False)

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
             run_prefix, save_prefix,
             calc_chn_ene=True,
             calc_ene_trans=True,
             zrxn=None,
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
    ene_chnlvl = None
    ene_reflvl = None
    zpe = None
    hf0k = None
    hf0k_trs = None

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

    # Obtain vibration partition function information
    ioprinter.info_message(
        'Preparing internal rotor info building partition functions...',
        newline=1)
    ioprinter.info_message(
        'Obtaining the vibrational frequencies and zpves...', newline=1)
    vib_anal_dct = vib.full_vib_analysis(
        spc_dct_i, pf_filesystems, spc_mod_dct_i,
        spc_dct, run_prefix, zrxn=zrxn)

    freqs = vib_anal_dct['fund_proj_RTimagTors']
    imag = vib_anal_dct['harm_imag']
    zpe = vib_anal_dct['anharm_zpe']
    tors_strs = vib_anal_dct['mess_tors_strs']
    rotors = vib_anal_dct['rotors']
    if typ.anharm_core(spc_mod_dct_i):
        xmat = vib_anal_dct['x_mat']
    if typ.nonrigid_rotations(spc_mod_dct_i):
        rovib_coups = vib_anal_dct['rovib_mat']
        rot_dists = vib_anal_dct['rot_dists']

    # Get the torsion strings
    allr_str = tors_strs[0]
    mdhr_dat = tors_strs[4]

    # Obtain symmetry factor
    ioprinter.info_message(
        'Determining the symmetry factor...', newline=1)

    zma = None
    zma_locs = (0,)
    if zrxn:
        [_, cnf_save_path, _, _, _] = pf_filesystems['harm']
        # Build the rotors
        if cnf_save_path:
            zma_fs = autofile.fs.zmatrix(cnf_save_path)
            zma_locs = ts_zma_locs(spc_dct, spc_name, zma_fs)
            zma = zma_fs[-1].file.zmatrix.read(zma_locs)

    racemic = True
    ioprinter.info_message('Setting symmetry factors as racemic=', racemic)
    sym_factor = symm.symmetry_factor(
        pf_filesystems, spc_mod_dct_i, spc_dct_i, rotors, grxn=zrxn, zma=zma,
        racemic=racemic)

    # Obtain electronic energy levels
    elec_levels = spc_dct_i['elec_levels']

    # Obtain energy levels
    ioprinter.info_message(
        'Obtaining the electronic energy + zpve...', newline=1)
    if calc_chn_ene:
        chn_ene = ene.read_energy(
            spc_dct_i, pf_filesystems, spc_mod_dct_i,
            spc_dct, run_prefix, read_ene=True, read_zpe=False, saddle=zrxn is not None)
        ene_chnlvl = chn_ene + zpe

        zma = None
        # Determine info about the basis species used in thermochem calcs
        hf0k, hf0k_trs, chn_basis_ene_dct, _ = basis.enthalpy_calculation(
            spc_dct, spc_name, ene_chnlvl,
            chn_basis_ene_dct, pes_mod_dct_i, spc_mod_dct_i,
            run_prefix, save_prefix, zrxn=zrxn)

    ene_reflvl = None

    #  Build the energy transfer section strings
    if calc_ene_trans and zrxn is None:
        ioprinter.info_message(
            'Determining energy transfer parameters...', newline=1)
        etrans_dct = etrans.etrans_dct_for_species(
            spc_dct_i, pes_mod_dct_i)
        spc_info = sinfo.from_dct(spc_dct_i, canonical=True)
        bath_info = etrans.set_bath(spc_dct, etrans_dct)

        edown_str, collid_freq_str = etrans.make_energy_transfer_strs(
            spc_info, bath_info, etrans_dct)
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
        ts_dct, spc_mod_dct_i['rpath']['geolvl'][1])

    # Set the prefix assuming locs of 00 for TS and VRC
    ts_save_prefix = os.path.join(ts_save_path, '00')
    vrc_locs = (0,)

    vrc_fs = autofile.fs.vrctst(ts_save_prefix)
    if vrc_fs[-1].file.vrctst_flux.exists(vrc_locs):

        flux_str = vrc_fs[-1].file.vrctst_flux.read(vrc_locs)

        vrc_path = vrc_fs[-1].file.vrctst_flux.path(vrc_locs)
        ioprinter.info_message(f'Reading flux file from {vrc_path}')
        ioprinter.warning_message('We have assumed VRC locs of 00')
    else:
        flux_str = None
        ioprinter.warning_message(f'No flux file at {ts_save_prefix}')

    # Create info dictionary
    inf_dct = {'flux_str': flux_str,
               'elec_levels': ts_dct['elec_levels']}

    return inf_dct


# VTST
def rpvtst_sadpt(tsname, spc_dct,
                 pes_mod_dct_i, spc_mod_dct_i,
                 chn_basis_ene_dct,
                 run_prefix, save_prefix):
    """ Pull all of the neccessary information from the
        filesystem for a species
    """

    ts_dct = spc_dct[tsname]

    # Set up filesystems and coordinates for saddle point
    # Scan along RxnCoord is under THY/TS/CONFS/cid/Z
    pf_filesystems = filesys.models.pf_filesys(
        ts_dct, spc_mod_dct_i, run_prefix, save_prefix,
        saddle=True, name=tsname)
    cnf_save_path = pf_filesystems['harm'][1]

    print(f'Using conformer from path {cnf_save_path}')

    # Set TS reaction coordinate and the scan values
    scn_vals, frm_name = filesys.models.get_rxn_scn_coords(
        cnf_save_path, coord_name='IRC', zma_locs=(0,))
    scn_vals.sort()
    scn_prefix = cnf_save_path

    # Modify the scn thy info
    # Assumes no composite ene being used for variational treats
    scn_ene_info = spc_mod_dct_i['ene']['lvl1'][1][1]
    mod_scn_ene_info = tinfo.modify_orb_label(
        scn_ene_info, sinfo.from_dct(ts_dct))
    ioprinter.warning_message(
        'using energy {scn_ene_info}, '
        'assuming no composite energy for variational for now.')

    # Need to read the sp vals along the scan. add to read
    ioprinter.info_message('Reading molecular data across the potential')
    ioprinter.info_message(f'Scan potential name = {frm_name}')
    enes, geoms, freqs, _ = _rpath_ene_data(
        ts_dct, frm_name, scn_vals,
        None, spc_mod_dct_i,
        mod_scn_ene_info,
        scn_prefix, run_prefix, save_prefix)

    # Grab the values from the read
    inf_dct = {}
    inf_dct['rpath'] = []
    pot_info = zip(scn_vals, enes.values(), geoms.values(), freqs.values())
    for rval, pot, geo, frqs in pot_info:

        ioprinter.info_message(
            f'\nDetermining variational values for R = {rval} Bohr')

        # Scale the R value if needed
        if frm_name == 'IRC':
            rval /= 100.0

        # Determine energy values
        zero_ene = _rpath_rel_ene(pot, frqs)
        hf0k, hf0k_trs, chn_basis_ene_dct, _ = basis.enthalpy_calculation(
            spc_dct, tsname, zero_ene,
            chn_basis_ene_dct, pes_mod_dct_i, spc_mod_dct_i,
            run_prefix, save_prefix, zrxn=ts_dct['zrxn'])

        # Set electronic energy values for all points
        elec_levels = ts_dct['elec_levels']

        # Create info dictionary and append to lst
        keys = ['rval', 'geom', 'freqs', 'elec_levels',
                'ene_chnlvl', 'ene_tsref']
        vals = [rval, geo, frqs, elec_levels,
                hf0k, hf0k_trs]
        inf_dct['rpath'].append(dict(zip(keys, vals)))

    # Calculate and store the imaginary mode and the index of the saddle point
    _, imag, _, _ = vib.read_harmonic_freqs(
        pf_filesystems, run_prefix, zrxn=ts_dct['zrxn'])
    inf_dct.update({'imag': imag})
    inf_dct.update({'ts_idx': scn_vals.index(0.00)})

    return inf_dct, chn_basis_ene_dct


def rpvtst_nobar(tsname, spc_dct, reac_dcts,
                 pes_mod_dct_i, spc_mod_dct_i,
                 chn_basis_ene_dct,
                 run_prefix, save_prefix):
    """ get info for barrierless transition state
    """

    ts_dct = spc_dct[tsname]

    # Set up filesystems and coordinates for reaction path
    # Scan along RxnCoord is under THY/TS/Z
    _, _, _, thy_save_path = filesys.models.set_rpath_filesys(
        ts_dct, spc_mod_dct_i['rpath'][1])

    # Set TS reaction coordinate
    scn_vals, frm_name = filesys.models.get_rxn_scn_coords(
        thy_save_path, coord_name=None, zma_locs=(0,))
    scn_vals.sort()
    scn_prefix = thy_save_path

    # Modify the scn thy info
    scn_ene_info = spc_mod_dct_i['rpath'][1][0]
    mod_scn_ene_info = tinfo.modify_orb_label(
        scn_ene_info, sinfo.from_dct(ts_dct))
    ioprinter.warning_message(
        'using energy {scn_ene_info}, '
        'assuming no composite energy for variational for now.')

    # Need to read the sp vals along the scan. add to read
    ioprinter.info_message('Reading molecular data across the potential')
    ioprinter.info_message(f'Scan potential name = {frm_name}')
    enes, geoms, freqs, mref = _rpath_ene_data(
        ts_dct, frm_name, scn_vals,
        reac_dcts, spc_mod_dct_i,
        mod_scn_ene_info,
        scn_prefix, run_prefix, save_prefix)
    ene_hs_sr_inf, ene_hs_sr_ref, ene_hs_mr_ref, reac_ene, zpe_ref = mref

    # Grab the values from the read
    inf_dct = {}
    inf_dct['rpath'] = []
    pot_info = zip(scn_vals, enes.values(), geoms.values(), freqs.values())
    for rval, pot, geo, frq in pot_info:

        ioprinter.info_message(
            f'\nDetermining variational values for R = {rval} Bohr')

        # Determine energy values
        zero_ene = _rpath_rel_ene(
            pot, frq,
            ene_hs_sr_ref=ene_hs_sr_ref,
            ene_hs_sr_inf=ene_hs_sr_inf,
            ene_hs_mr_ref=ene_hs_mr_ref,
            reac_ene=reac_ene,
            zpe_ref=zpe_ref)

        hf0k, hf0k_trs, chn_basis_ene_dct, _ = basis.enthalpy_calculation(
            spc_dct, tsname, zero_ene,
            chn_basis_ene_dct, pes_mod_dct_i, spc_mod_dct_i,
            run_prefix, save_prefix, zrxn=ts_dct['zrxn'])

        # Set electronic energy values for all points
        elec_levels = ts_dct['elec_levels']

        # Create info dictionary and append to lst
        keys = ['rval', 'geom', 'freqs', 'elec_levels',
                'ene_chnlvl', 'ene_tsref']
        vals = [rval, geo, frq, elec_levels,
                hf0k, hf0k_trs]
        inf_dct['rpath'].append(dict(zip(keys, vals)))

    # Calculate and store the imaginary mode
    inf_dct.update({'imag': None})
    inf_dct.update({'ts_idx': 0})

    return inf_dct, chn_basis_ene_dct


def _rpath_ene_data(ts_dct, frm_name, scn_vals,
                    reac_dcts, spc_mod_dct_i,
                    mod_scn_ene_info,
                    scn_prefix, run_prefix, save_prefix):
    """ Read molecular data along a scan including what is needed
        to get the multireference infinite separation energy value
    """

    fml_str = 'RPATH'
    vib_path = job_path(run_prefix, 'PROJROT', 'FREQ', fml_str)

    # Read all of the data along the scan
    ref_ene = 0.0
    enes, geoms, grads, hessians, _, _ = filesys.read.potential(
        [frm_name], [scn_vals], scn_prefix,
        mod_scn_ene_info, ref_ene,
        None,   # No extra frozen treatments
        read_geom=True,
        read_grad=True,
        read_hess=True,
        read_energy_backstep=False,
        remove_bad_points=False)

    script_str = autorun.SCRIPT_DCT['projrot']
    freqs = autorun.projrot.pot_frequencies(
        script_str, geoms, grads, hessians, vib_path)

    # Read all data needed to get multiref inf sep ene values
    # Based on if scan is a multireference method
    if elstruct.par.Method.is_multiref(mod_scn_ene_info[1]):
        ioprinter.info_message(
            'Scan ene method ({mod_scn_ene_info}) is multireference method. '
            '\n Reading all values to get infinite separation energy')

        # Get the energies and zpes at R_ref
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
                spc_dct, read_ene=True, read_zpe=True, saddle=False)

            ioprinter.debug_message('rpath', spc_mod_dct_i['rpath'][1])
            new_spc_dct_i = {
                'ene': ['mlvl', [[1.0, spc_mod_dct_i['rpath'][1][2]]]],
                'harm': spc_mod_dct_i['harm'],
                'tors': spc_mod_dct_i['tors']
            }
            ene_hs_sr_inf += ene.read_energy(
                dct, pf_filesystems, new_spc_dct_i, run_prefix,
                spc_dct, read_ene=True, read_zpe=False, saddle=True)

        mref_data = (ene_hs_sr_inf, ene_hs_sr_ref, ene_hs_mr_ref,
                     reac_ene, zpe_ref)
    else:
        mref_data = None

    return enes, geoms, freqs, mref_data


def _rpath_rel_ene(pot, frqs,
                   ene_hs_sr_ref=None,
                   ene_hs_sr_inf=None,
                   ene_hs_mr_ref=None,
                   reac_ene=None,
                   zpe_ref=None):
    """ Calculate the relative energy
    """

    # Get the relative energy (edit for radrad scans)
    zpe = (sum(frqs) / 2.0) * phycon.WAVEN2KCAL

    if ene_hs_sr_ref is None:
        # Single reference energy
        zero_ene = (pot + zpe) * phycon.KCAL2EH
    else:
        # Use energy components for a multireference energy
        ioprinter.debug_message('Energies to get total ene:')
        ioprinter.debug_message("""
        elec_ene = (ene_hs_sr_ref - ene_hs_sr_inf -
                    ene_hs_mr_ref + pot)""")
        ioprinter.debug_message("""
            zero_ene = (reac_ene +
                        (elec_ene + zpe_pt * phycon.KCAL2EH))""")
        ioprinter.debug_message('reac ene', reac_ene)
        ioprinter.debug_message('ene_hs_sr_ref', ene_hs_sr_ref)
        ioprinter.debug_message('ene_hs_sr_inf', ene_hs_sr_inf)
        ioprinter.debug_message('ene_hs_mr_ref', ene_hs_mr_ref)
        ioprinter.debug_message('pot R', pot * phycon.KCAL2EH)
        ioprinter.debug_message('zpe', zpe)
        ioprinter.debug_message('zpe ref', zpe_ref)

        elec_ene = (
            ene_hs_sr_ref - ene_hs_sr_inf -
            ene_hs_mr_ref + pot * phycon.KCAL2EH
        )
        zpe_pt = zpe - zpe_ref
        zero_ene = reac_ene + (elec_ene + zpe_pt * phycon.KCAL2EH)

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

    return zero_ene


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
        f'treatment that yields k({temp_pst} K) = {kt_pst}',
        newline=1)
    ioprinter.info_message(
        f'Assuming PST model potential V = C0 / R^{n_pst}',
        indent=1)

    # Obtain the reduced mass of the reactants
    ioprinter.info_message(
        'Reading reactant geometries to obtain reduced mass...', newline=1)
    geoms = []

    for dct in reac_dcts:
        pst_mod_dct_i = {'vib': spc_mod_dct_i['vib']}
        pf_filesystems = filesys.models.pf_filesys(
            dct, pst_mod_dct_i, run_prefix, save_prefix, False)
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

    # Set up model and basic thy objects
    spc_info = sinfo.from_dct(spc_dct_i, canonical=True)
    thy_info = spc_mod_dct_i['vib']['geolvl'][1][1]
    mod_thy_info = tinfo.modify_orb_label(
        thy_info, spc_info)

    vib_model = spc_mod_dct_i['vib']['mod']

    # Set up reference conformer filesys
    pf_filesystems = filesys.models.pf_filesys(
        spc_dct_i, spc_mod_dct_i, run_prefix, save_prefix, saddle)
    [harm_save_fs, _, harm_min_locs, _, _] = pf_filesystems['harm']

    # Obtain all values from initial reference conformer
    vib_anal_dct = vib.full_vib_analysis(
        spc_dct_i, pf_filesystems, spc_mod_dct_i,
        run_prefix, zrxn=None)
    freqs = vib_anal_dct['fund_proj_RTimagTors']
    # imag = vib_anal_dct['harm_imag']
    zpe = vib_anal_dct['anharm_zpe']
    tors_strs = vib_anal_dct['mess_tors_strs']
    rotors = vib_anal_dct['rotors']
    harm_freqs = vib_anal_dct['harm_proj_RTimag']

    harm_zpve = 0.5 * sum(harm_freqs) * phycon.WAVEN2EH

    ioprinter.info_message('Determining the symmetry factor...', newline=1)
    sym_factor = symm.symmetry_factor(
        pf_filesystems, spc_mod_dct_i, spc_dct_i, rotors,
    )

    zpe_chnlvl = zpe * phycon.EH2KCAL
    ref_ene = harm_zpve * phycon.EH2KCAL

    ref_geom = [harm_save_fs[-1].file.geometry.read(harm_min_locs)]
    ref_grad = [harm_save_fs[-1].file.gradient.read(harm_min_locs)]
    ref_hessian = [harm_save_fs[-1].file.hessian.read(harm_min_locs)]

    min_cnf_ene = filesys.read.energy(
        harm_save_fs, harm_min_locs, mod_thy_info)

    # Set up the TAU filesystem objects, get locs, and read info
    _, tau_save_fs = filesys.build_fs(
        run_prefix, save_prefix, 'TAU',
        spc_locs=spc_info, thy_locs=mod_thy_info[1:])

    db_style = 'jsondb'
    vib_model = spc_mod_dct_i['vib']['mod']
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

    ioprinter.info_message(
        'Reading data for the Monte Carlo samples from db.json'
        f'at path {tau_save_fs[0].path()}')
    samp_geoms, samp_enes, samp_grads, samp_hessians = [], [], [], []
    tot_locs = len(tau_locs)
    for idx, locs in enumerate(tau_locs):

        if db_style == 'directory':
            geo = tau_save_fs[-1].file.geometry.read(locs)
        elif db_style == 'jsondb':
            geo = tau_save_fs[-1].json.geometry.read(locs)

        # geo_str = autofile.data_types.swrite.geometry(geo)
        samp_geoms.append(geo)

        if db_style == 'directory':
            tau_ene = tau_save_fs[-1].file.energy.read(locs)
        elif db_style == 'jsondb':
            tau_ene = tau_save_fs[-1].json.energy.read(locs)
        rel_ene = (tau_ene - min_cnf_ene) * phycon.EH2KCAL
        # ene_str = autofile.data_types.swrite.energy(rel_ene)
        samp_enes.append(rel_ene)

        if vib_model == 'tau':
            if db_style == 'directory':
                grad = tau_save_fs[-1].file.gradient.read(locs)
            elif db_style == 'jsondb':
                grad = tau_save_fs[-1].json.gradient.read(locs)
            # grad_str = autofile.data_types.swrite.gradient(grad)
            samp_grads.append(grad)

            if db_style == 'directory':
                hess = tau_save_fs[-1].file.hessian.read(locs)
            elif db_style == 'jsondb':
                hess = tau_save_fs[-1].json.hessian.read(locs)
            # hess_str = autofile.data_types.swrite.hessian(hess)
            samp_hessians.append(hess)

        # Print progress message (every 150 geoms read)
        if idx % 149 == 0:
            print(f'Read {idx+1}/{tot_locs} samples...')

    # Determine the successful conformer ratio
    inf_obj = tau_save_fs[0].file.info.read()
    excluded_volume_factor = len(samp_geoms) / inf_obj.nsamp
    print('excluded volume factor test:',
          excluded_volume_factor, len(samp_geoms), inf_obj.nsamp)

    # Create info dictionary
    keys = ['geom', 'sym_factor', 'elec_levels',
            'freqs', 'flux_mode_str',
            'samp_geoms', 'samp_enes', 'samp_grads', 'samp_hessians',
            'ref_geom', 'ref_grad', 'ref_hessian',
            'zpe_chnlvl', 'ref_ene', 'excluded_volume_factor']
    vals = [ref_geom[0], sym_factor, spc_dct_i['elec_levels'],
            freqs, tors_strs[2],
            samp_geoms, samp_enes, samp_grads, samp_hessians,
            ref_geom, ref_grad, ref_hessian,
            zpe_chnlvl, ref_ene, excluded_volume_factor]
    inf_dct = dict(zip(keys, vals))

    return inf_dct

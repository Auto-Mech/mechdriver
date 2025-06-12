"""
Writes the full sections of a MESS input file
"""

import os
from ioformat import build_mako_str
from ioformat import indent
from ioformat import remove_trail_whitespace
from phydat import phycon
from mess_io.writer._mol_inf import core_rigidrotor
from mess_io.writer._spc import molecule
from mess_io.writer._rxnchan import species
from mess_io.writer import _format as messformat


# OBTAIN THE PATH TO THE DIRECTORY CONTAINING THE TEMPLATES #
SRC_PATH = os.path.dirname(os.path.realpath(__file__))
TEMPLATE_PATH = os.path.join(SRC_PATH, 'templates')
SECTION_PATH = os.path.join(TEMPLATE_PATH, 'sections')


# Write the full input file strings
def messrates_inp_str(globkey_str, rxn_chan_str,
                      energy_trans_str=None, well_lump_str=None,
                      use_short_names=False):
    """ Combine various MESS strings together to combined MESS rates
    """

    # Create dictionary to fill template
    full_rxn_inp_keys = {
        'globkey_str': globkey_str,
        'energy_trans_str': energy_trans_str,
        'well_lump_str': well_lump_str,
        'rxn_chan_str': rxn_chan_str,
        'use_short_names': use_short_names
    }

    mess_inp_str = build_mako_str(
        template_file_name='full_rates_inp.mako',
        template_src_path=SECTION_PATH,
        template_keys=full_rxn_inp_keys)

    return remove_trail_whitespace(mess_inp_str)


def messpf_inp_str(globkey_str, spc_str):
    """ Combine various MESS strings together to combined MESSPF
    """
    return '\n'.join([globkey_str, spc_str]) + '\n'


def messhr_inp_str(geo, hind_rot_str,
                   float_type='double'):
    """ Special MESS input string to calculate frequencies and ZPVEs
        for hindered rotors
    """

    global_pf_str = global_pf_input(
        temperatures=(100.0, 200.0, 300.0, 400.0, 500),
        rel_temp_inc=0.001,
        atom_dist_min=0.6,
        float_type=float_type)
    dat_str = molecule(
        core=core_rigidrotor(geo, 1.0),
        freqs=(1000.0,),
        elec_levels=((0.0, 1.0),),
        hind_rot=hind_rot_str,
    )
    spc_str = species(
        spc_label='Tmp',
        spc_data=dat_str,
        zero_ene=0.0
    )

    return messpf_inp_str(global_pf_str, spc_str)


# Write individual sections of the input file
def global_rates_input_v1(
        temperatures, pressures,
        calculation_method='well-reduction',
        model_ene_limit=800.0,
        ene_stepover_temp=0.2, excess_ene_temp=None,
        well_extension=0.001,
        chem_eig_max=None,
        well_reduction_thresh=10.0,
        ground_ene_shift_max=None,
        ped_spc_lst=None,
        hot_enes_dct=None,
        micro_out_params=None,
        float_type='double',
        ktp_outname='rate.out',
        ke_outname='ke.out',
        ped_outname='ped.out'):
    """ Writes the global keywords section of the MESS input file by
        formatting input information into strings a filling Mako template.

        This is for the original MESS binary code.

        :param temperatures: List of temperatures (in K)
        :type temperatures: float
        :param pressures: List of pressures (in atm)
        :type pressures: float
        :param reduction_method:
        :type reduction_method: str
        :param well_reduction_thresh:
        :type well_reduction_thresh: float
        :param excess_ene_temp:
        :type excess_ene_temp: str/float
        :param float_type: precision of floats used in MESS calculation
        :type float_type: str
        :return global_str: String for section
        :rtype: string
    """

    assert calculation_method in ('direct', 'well-reduction'), (
        f'calculation_method is {calculation_method}, '
        'not direct or well_reduction')
    assert float_type in ('double', 'quadruple'), (
        f'float_type is {float_type}, not double or quadruple')

    # Format temperature and pressure lists
    temperature_list = '  '.join(str(val) for val in temperatures)
    pressure_list = '  '.join(str(val) for val in pressures)

    # Format the other keywords as needed
    if excess_ene_temp is not None:
        assert isinstance(excess_ene_temp, float), (
            'ExcessEnergyOverTemperature value must be a float'
        )
        excess_ene_temp_str = f'{excess_ene_temp:.2f}'
    else:
        excess_ene_temp_str = None

    if well_extension is not None:
        if well_extension != 'auto':
            assert isinstance(well_extension, float), (
                'WellExtension value must be a float'
            )
            well_extension_str = f'{well_extension:.4f}'
        else:
            well_extension_str = ''
    else:
        well_extension_str = None

    if chem_eig_max is not None:
        chem_eig_max = f'{chem_eig_max:.2f}'

    if ped_spc_lst is not None:
        ped_spc_str = messformat.format_ped_species(ped_spc_lst)
    else:
        ped_spc_str = None

    if hot_enes_dct is not None:
        nhot, hot_ene_str = messformat.format_hot_enes(hot_enes_dct)
        if ped_spc_lst is not None:
            calculation_method = 'well-reduction'
    else:
        nhot, hot_ene_str = 0, None

    well_reduction_thresh_str = f'{well_reduction_thresh:.2f}'

    if ground_ene_shift_max is not None:
        ground_ene_shift_max = f'{ground_ene_shift_max*phycon.EH2KCAL:.2f}'

    if micro_out_params is not None:
        assert (len(micro_out_params) == 3 and
                all(isinstance(x, float) for x in micro_out_params)), (
            f'{micro_out_params} is not a tuple/list of three floats')

    # Create dictionary to fill template
    globrxn_keys = {
        'temperatures': temperature_list,
        'pressures': pressure_list,
        'model_ene_limit': f'{model_ene_limit:.2f}',
        'ene_stepover_temp': f'{ene_stepover_temp:.2f}',
        'excess_ene_temp': excess_ene_temp_str,
        'calculation_method': calculation_method,
        'well_reduction_thresh': well_reduction_thresh_str,
        'well_extension': well_extension_str,
        'ground_ene_shift_max': ground_ene_shift_max,
        'chem_eig_max': chem_eig_max,
        'hot_ene_str': hot_ene_str,
        'nhot': nhot,
        'ped_spc_str': ped_spc_str,
        'micro_out_params': micro_out_params,
        'float_type': float_type,
        'ktp_outname': ktp_outname,
        'ke_outname': ke_outname,
        'ped_outname': ped_outname,
    }

    return build_mako_str(
        template_file_name='global_rates.mako',
        template_src_path=SECTION_PATH,
        template_keys=globrxn_keys)


def global_rates_input_v2(
        temperatures, pressures,
        ref_temperature=1000.0, ref_pressure=1.0,
        model_ene_limit=800.0,
        ene_stepover_temp=0.2, ene_cutoff_temp=20.0, excess_ene_temp=10.0,
        chem_tol=1.0e-10, chem_thresh=0.1,
        well_pojection_thresh=0.1, well_reduction_thresh=10.0,
        time_propagation_limit=50.0, time_propagation_step=0.02,
        well_extension=0.001,
        ground_ene_shift_max=None,
        ped_spc_lst=None, hot_enes_dct=None,
        micro_out_params=None,
        float_type='double',
        ktp_outname='rate.out',
        ke_outname='ke.out',
        ped_outname='ped.out'):
    """ Writes the global keywords section of the MESS input file by
        formatting input information into strings a filling Mako template.

        This is for the new MESS binary code.

        :param temperatures: List of temperatures (in K)
        :type temperatures: float
        :param pressures: List of pressures (in atm)
        :type pressures: float
        :return global_str: String for section
        :rtype: string
    """

    assert float_type in ('double', 'quadruple'), (
        f'float_type is {float_type}, not double or quadruple')

    # Format temperature and pressure lists
    temperature_list = '  '.join(str(val) for val in temperatures)
    pressure_list = '  '.join(str(val) for val in pressures)

    if ped_spc_lst is not None:
        ped_spc_str = messformat.format_ped_species(ped_spc_lst)
    else:
        ped_spc_str = None

    if hot_enes_dct is not None:
        nhot, hot_ene_str = messformat.format_hot_enes(hot_enes_dct)
    else:
        nhot, hot_ene_str = 0, None

    if ground_ene_shift_max is not None:
        ground_ene_shift_max = f'{ground_ene_shift_max*phycon.EH2KCAL:.2f}'

    if micro_out_params is not None:
        assert (len(micro_out_params) == 3 and
                all(isinstance(x, float) for x in micro_out_params)), (
            f'{micro_out_params} is not a tuple/list of three floats')

    # Create dictionary to fill template
    globrxn_keys = {
        'temperatures': temperature_list,
        'pressures': pressure_list,
        'ref_temperature': f'{ref_temperature:.2f}',
        'ref_pressure': f'{ref_pressure:.2f}',
        'model_ene_limit': f'{model_ene_limit:.2f}',
        'ene_stepover_temp': f'{ene_stepover_temp:.2f}',
        'ene_cutoff_temp': f'{ene_cutoff_temp:.2f}',
        'excess_ene_temp': f'{excess_ene_temp:.2f}',
        'chem_tol': f'{chem_tol:.2e}',
        'chem_thresh': f'{chem_thresh:.2f}',
        'well_projection_thresh': f'{well_pojection_thresh:.2f}',
        'well_reduction_thresh': f'{well_reduction_thresh:.2f}',
        'time_propagation_limit': f'{time_propagation_limit:.2f}',
        'time_propagation_step': f'{time_propagation_step:.2f}',
        'well_extension': f'{well_extension:.4f}',
        'ground_ene_shift_max': ground_ene_shift_max,
        'hot_ene_str': hot_ene_str,
        'nhot': nhot,
        'ped_spc_str': ped_spc_str,
        'micro_out_params': micro_out_params,
        'float_type': float_type,
        'ktp_outname': ktp_outname,
        'ke_outname': ke_outname,
        'ped_outname': ped_outname
    }

    return build_mako_str(
        template_file_name='global_rates_v2.mako',
        template_src_path=SECTION_PATH,
        template_keys=globrxn_keys)


def global_pf_input(temperatures=(),
                    temp_step=100, ntemps=30,
                    rel_temp_inc=0.001, atom_dist_min=1.13384,
                    float_type='double'):
    """ Writes the global keywords section of the MESS input file by
        formatting input information into strings a filling Mako template.

        :param temperatures: List of temperatures (in K)
        :type temperatures: list(float)
        :param temp_step: temperature step (in K)
        :type temp_step: float
        :param ntemps: number of temperature values on grid
        :type ntemps: int
        :param rel_temp_inc: increment for temps
        :type rel_temp_inc: float
        :param atom_dist_min: cutoff for atom distances (Bohr)
        :type atom_dist_min: float
        :param float_type: precision of floats used in MESS calculation
        :type float_type: str
        :return global_pf_str: string for section
        :rtype: string
    """

    assert float_type in ('double', 'quadruple'), (
        f'float_type is {float_type}, not double or quadruple')

    if temperatures:
        temperature_list = '  '.join(str(val) for val in temperatures)
        temp_step = None
        ntemps = None
    else:
        temperature_list = ''

    # Convert the atom distance minimum
    atom_dist_min = f'{atom_dist_min * phycon.BOHR2ANG:.2f}'

    # Create dictionary to fill template
    globpf_keys = {
        'temperatures': temperature_list,
        'temp_step': temp_step,
        'ntemps': ntemps,
        'rel_temp_inc': rel_temp_inc,
        'atom_dist_min': atom_dist_min,
        'float_type': float_type
    }

    return build_mako_str(
        template_file_name='global_pf.mako',
        template_src_path=SECTION_PATH,
        template_keys=globpf_keys)


def global_energy_transfer_input(edown_str, collid_freq_str):
    """ Writes the global energy transfer section of the MESS input file by
        formatting input information into strings a filling Mako template.

        :param edown_str: String for the energy down parameters
        :type edown_str: str
        :param collid_freq_str: String for the collisional freq parameters
        :type collid_freq_str: str
        :rtype: str
    """

    edown_str = indent(edown_str, 2)
    collid_freq_str = indent(collid_freq_str, 2)

    # Create dictionary to fill template
    glob_etrans_keys = {
        'edown_str': edown_str,
        'collid_freq_str': collid_freq_str
    }

    return build_mako_str(
        template_file_name='global_etrans.mako',
        template_src_path=SECTION_PATH,
        template_keys=glob_etrans_keys)


# Write data to output file formats if you want to make a formatted output file
def pf_output(fml_str, temps, logq, dq_dt, d2q_dt2, svals=None, cpvals=None):
    """ Writes partition function data into a string that is formatted like the
        output file
    """

    mess_out_str = (
        "Natural log of the partition function, "
        "its derivatives, entropy, and thermal capacity:\n"
        f"T, K          {fml_str}\n"
        "               Z_0          Z_1          Z_2 "
        "S, cal/mol/K C, cal/mol/K"
    )
    for idx, _ in enumerate(temps):
        mess_out_str += '\n'
        mess_out_str += f'{temps[idx]:>8.3f}    '
        mess_out_str += f'{logq[idx]:>8.6f}    '
        mess_out_str += f'{dq_dt[idx]:>8.8f}    '
        mess_out_str += f'{d2q_dt2[idx]:>8.6e}    '
        if svals is not None:
            mess_out_str += f'{svals[idx]:>8.6f}    '
        if cpvals is not None:
            mess_out_str += f'{cpvals[idx]:>8.6f}    '

    return mess_out_str

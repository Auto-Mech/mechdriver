""" Calculate Total Energy Distribution and related quantities
"""

import automol
import intder_io
from autorun._run import from_input_string


INPUT_NAME = 'intder.inp'
OUTPUT_NAMES = ('intder.out',)
SCRIPT_NAME = 'run_intder.sh'


# Specialized Runners
def reaction_coordinate_check_idxs(script_str, run_dir,
                                   geo, zma, hess, zrxn,
                                   input_name=INPUT_NAME,
                                   output_names=OUTPUT_NAMES):
    """ Calculate internal coordinate rep of mode from TED
    """

    if not automol.zmat.dummy_keys(zma):

        ted_mode_dct, ted_mode_freq = ted_zmatrix_coordinates(
            script_str, run_dir, geo, zma, hess, 0,
            input_name=input_name,
            output_names=output_names)

        if ted_mode_dct is not None:
            # Get the zmat names corresponding to frm/brk keys; remove Nones
            ted_keys = frozenset(ted_mode_dct.keys())
            tsg = automol.reac.ts_graph(zrxn)
            frm_keys = automol.graph.ts.forming_bond_keys(tsg)
            brk_keys = automol.graph.ts.breaking_bond_keys(tsg)
            zrxn_keys = frozenset().union(*(frm_keys, brk_keys))
            print('all keys test', zrxn_keys)

            print(f'TED Comparisons for {automol.reac.class_(zrxn)}:')

            print(f'- Freq: {ted_mode_freq} cm-1')

            _strs = []
            for idxs, pct in ted_mode_dct.items():
                _idx_str = "-".join((str(idx+1) for idx in idxs))
                _strs.append(f'{_idx_str} ({pct})')
            tedname_str = ', '.join(_strs)
            print(f'- TED: {tedname_str}'.format(tedname_str))

            _strs = []
            for keys in zrxn_keys:
                _key_str = "-".join((str(key+1) for key in keys))
                _strs.append(_key_str)
            rkey_str = ', '.join(_strs)
            print(f'- RXN: {rkey_str}  (only R)')

            if ted_keys & zrxn_keys:
                unmatched_keys = zrxn_keys - ted_keys
                _strs = []
                for keys in unmatched_keys:
                    _key_str = "-".join((str(key+1) for key in keys))
                    _strs.append(_key_str)
                _unmatched_str = ', '.join(_strs)
                print(f'- Unmatched ({len(unmatched_keys)}): {_unmatched_str}')
                print('- Overlap of coordinates found, possible success')
                success = True
            else:
                print('- No similarity of coords, likely something is wrong')
                success = False
        else:
            print('INTDER had some error, skipping TED check')
            success = True
    else:
        print('Z-Matrix has dummy atoms, cannot do TED check')
        success = True

    return success


def reaction_coordinate_check(script_str, run_dir,
                              geo, zma, hess, zrxn,
                              input_name=INPUT_NAME,
                              output_names=OUTPUT_NAMES):
    """ Calculate internal coordinate rep of mode from TED
    """

    if not automol.zmat.dummy_keys(zma):

        ted_mode_names_dct, ted_mode_freq = ted_zmatrix_coordinates(
            script_str, run_dir, geo, zma, hess, 0,
            input_name=input_name,
            output_names=output_names)

        if ted_mode_names_dct is not None:
            # Get the zmat names corresponding to frm/brk keys; remove Nones
            ted_names = tuple(ted_mode_names_dct.keys())

            rxn_names = automol.reac.zmatrix_coordinate_names(zrxn, zma)
            rxn_names = (tuple(x for x in rxn_names[0] if x is not None) +
                         tuple(x for x in rxn_names[1] if x is not None))

            print(f'TED Comparisons for {automol.reac.class_(zrxn)}:')
            print(f'- Freq: {ted_mode_freq} cm-1')
            tedname_str = ' '.join(
                f'{name} ({pct})' for name, pct in ted_mode_names_dct.items()
            )
            print(f'- TED: {tedname_str}'.format(tedname_str))
            rname_str = ' '.join(rxn_names)
            print(f'- RXN: {rname_str}  (only R)')

            _set_ted_names, _set_rxn_names = set(ted_names), set(rxn_names)
            if _set_ted_names & _set_rxn_names:
                _unmatched = _set_rxn_names - _set_ted_names
                _unmatched_str = ' '.join(_unmatched)
                print(f'- Unmatched ({len(_unmatched)}): {_unmatched_str}')
                print('- Overlap of coordinates found, possible success')
                success = True
            else:
                print('- No similarity of coords, likely something is wrong')
                success = False
        else:
            print('INTDER had some error, skipping TED check')
            success = True
    else:
        print('Z-Matrix has dummy atoms, cannot do TED check')
        success = True

    return success


def ted_zmatrix_coordinates(script_str, run_dir,
                            geo, zma, hess, mode_idx,
                            input_name=INPUT_NAME,
                            output_names=OUTPUT_NAMES):
    """ Calculate internal coordinate rep of mode from TED

        Note that geo and zma atom ordering must match.
        Currently, calculation does not support dummy atoms/linear segments
    """

    # Run INTDER
    output_strs = direct(script_str, run_dir, zma, geo, hess,
                         script_name=SCRIPT_NAME,
                         input_name=input_name,
                         output_names=output_names)
    output_str = output_strs[0]

    # Get modes from TED output
    if output_strs is not None:
        # Read the internal coordinates and TED assignments from output
        intl_coords = intder_io.reader.internal_coordinates(output_str)
        ted_assign = intder_io.reader.ted_assignments(output_str)

        # Get the indices
        # ted_mode_dct_zma = intder_io.ted_zmatrix_coordinates(
        #     zma, mode_idx, intl_coords, ted_assign)
        # Get the Z-Matrix names
        ted_mode_dct_zma = intder_io.ted_coordinate_indices(
            mode_idx, intl_coords, ted_assign)

        # Get the frequency of the mode
        ted_mode_freq = ted_assign[mode_idx][0]
    else:
        ted_mode_dct_zma = None
        ted_mode_freq = None

    return ted_mode_dct_zma, ted_mode_freq


# Generalized Runner
def direct(script_str, run_dir, zma, geo, hess,
           script_name=SCRIPT_NAME,
           input_name=INPUT_NAME,
           output_names=OUTPUT_NAMES):
    """ Generates an input file for a ProjRot job, runs it directly, and
        obtains all of the possible output file strings
    """

    input_str = intder_io.writer.input_file(geo, zma)
    aux_dct = {'file15': intder_io.writer.cartesian_hessian_file(hess)}

    return from_input_string(
        script_str, run_dir, input_str,
        aux_dct=aux_dct,
        script_name=script_name,
        input_name=input_name,
        output_names=output_names)

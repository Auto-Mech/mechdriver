""" Script to read in a PES from a MESS input file and then plot it.

    PLot generated is ReactionCoordinate vs. Relative Energy.
"""

import sys
import os
import argparse

import ioformat
import automol.util
import mechanalyzer.plotter
import mechanalyzer.parser


def plot_single_sccs(
        spc_file, mech_file, cwd, pes_idx, ccs_idx,
        i=None, j=None):
    """ plot the graph of a single sccs
    """
    # Parse the input based on the initial type
    inp_spc_str = ioformat.pathtools.read_file(cwd, spc_file)
    inp_mech_str = ioformat.pathtools.read_file(cwd, mech_file)
    if inp_spc_str is None and inp_mech_str is None:
        print('ERROR: Input species.csv or mechanism.dat not found')
        sys.exit()
    mech_spc_dct = mechanalyzer.parser.spc.build_spc_dct(
        inp_spc_str, 'csv')
    pes_dct = mechanalyzer.parser.pes.pes_dictionary(
        inp_mech_str, 'chemkin', mech_spc_dct)

    # Grab the channels of the PES dct that was requested
    print(f'Plotting PES-CCS: {pes_idx}-{ccs_idx}')
    chnls = ()
    for (form, pidx, sidx), sub_chnls in pes_dct.items():
        if pes_idx == -10:
            chnls += sub_chnls
        elif pidx == pes_idx-1 and sidx == ccs_idx-1:
            chnls += sub_chnls
            break
    spc_ichs, spc_names = (), ()
    conn_lst = ()
    for chnl in chnls:
        # Obtain the InChIs and SMILES for reagents from the spc dct
        rcts, prds = chnl[1][0], chnl[1][1]
        rct_ichs = tuple(mech_spc_dct[rct]['inchi'] for rct in rcts)
        prd_ichs = tuple(mech_spc_dct[prd]['inchi'] for prd in prds)
        rct_smis = tuple(mech_spc_dct[rct]['smiles'] for rct in rcts)
        prd_smis = tuple(mech_spc_dct[prd]['smiles'] for prd in prds)

        # Build connection list with smiles or NAME
        conn_lst += (('+'.join(rcts), '+'.join(prds)),)
        # conn_lst += (('+'.join(rct_smis), '+'.join(prd_smis)),)

        # Build master species list with inchis
        spc_ichs += rct_ichs + prd_ichs
        spc_names += rct_smis + prd_smis

    # Remove redundant smiles and names
    spc_ichs = automol.util.remove_duplicates_with_order(spc_ichs)
    spc_names = automol.util.remove_duplicates_with_order(spc_names)

    print('Connections:')
    for conn in conn_lst:
        print(f'{conn[0]} - {conn[1]}')

    # Produce the 2D image plot
    # img_path = os.path.join(cwd, f'structs_{options["pes"]}.pdf')
    # automol.chi.draw_grid(spc_ichs, names=spc_names, save_path=img_path)

    G = mechanalyzer.plotter.pes.sccs_graph(conn_lst, (i, j))
    # if len(G.nodes) > 0:
    #     mechanalyzer.plotter.pes.show_sccs(G, '{:g}_{:g}'.format(i, j))
    return G, conn_lst


def plot_all_sccs(cwd, options):
    g_lst = ()
    g_conns = ()
    label_lst = ()
    pes = options['pes'].split('_')
    if len(pes) > 2:
        pes_idx, ccs_idx, sccs_idx = pes
        sccs_idx_lst = ()
        if len(ccs_idx.split('-')) > 1:
            start_idx, end_idx = [int(x) for x in ccs_idx.split('-')]
            ccs_idx_lst = list(range(start_idx, end_idx + 1))
        else:
            ccs_idx_lst = [int(ccs_idx)]
        if len(sccs_idx.split('-')) > 1:
            start_idx, end_idx = [int(x) for x in sccs_idx.split('-')]
            sccs_idx_lst = list(range(start_idx, end_idx + 1))
        else:
            sccs_idx_lst = [int(sccs_idx)]
        for i in ccs_idx_lst:
            for j in sccs_idx_lst:
                if i in [9, 10, 11, 12, 13, 14]:
                    continue
                mech_file = options['mechanism'] + '_{:g}_{:g}'.format(i, j)
                spc_file = options['species'] + '_{:g}_{:g}'.format(i, j)
                if not os.path.exists(mech_file):
                    print('{} not found, skipping'.format(mech_file))
                    break
                G, conn_lst = plot_single_sccs(
                    spc_file, mech_file, cwd, -10, -10, i=i, j=j)
                g_lst += (G,)
                g_conns += (conn_lst,)
                label_lst += ((i, j),)
    else:
        pes_idx, ccs_idx = pes
        mech_file = options['mechanism']
        spc_file = options['species']
        if len(pes_idx.split('-')) > 1:
            start_idx, end_idx = [int(x) for x in pes_idx.split('-')]
            pes_idx_lst = list(range(start_idx, end_idx + 1))
        else:
            pes_idx_lst = [int(pes_idx)]
        if len(ccs_idx.split('-')) > 1:
            start_idx, end_idx = [int(x) for x in ccs_idx.split('-')]
            ccs_idx_lst = list(range(start_idx, end_idx + 1))
        else:
            ccs_idx_lst = [int(ccs_idx)]
        for pes_idx in pes_idx_lst:
            for ccs_idx in ccs_idx_lst:
                G, conn_lst = plot_single_sccs(
                        spc_file, mech_file, cwd,
                        pes_idx, ccs_idx, i=pes_idx, j=ccs_idx)
                g_lst += (G,)
                g_conns += (conn_lst,)
                label_lst += ((pes_idx, ccs_idx),)
    return g_lst, g_conns, label_lst


if __name__ == '__main__':
    # Set useful global variables
    oscwd = os.getcwd()

    # Parse the command line
    DSTR = 'Generates a graph plot of a PES from Chemkin mechanism'
    par = argparse.ArgumentParser(description=DSTR)
    par.add_argument('-s', '--species', default='species.csv',
                     help='species file (species.csv)')
    par.add_argument('-m', '--mechanism', default='mechanism.dat',
                     help='mechanism file (mechanism.dat)')
    par.add_argument('-p', '--pes', default='1_1',
                     help='PES_SUBPES to plot (1_1)')
    par.add_argument('-o', '--output', default='surface.pdf',
                     help='name of output plot file (surface.pdf)')
    opts = vars(par.parse_args())

    sccs_ret = plot_all_sccs(oscwd, opts)
    all_g, all_conns, labels_lst = sccs_ret
    if len(all_g) > 1:
        mechanalyzer.plotter.pes.show_pes(all_g, all_conns, labels_lst)

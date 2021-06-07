"""
"""


def test__build_fs():
    """ build fs
    """
    RUN_PREFIX = os.path.join(os.getcwd(), 'run')
    SAVE_PREFIX = os.path.join(os.getcwd(), 'save')

    RXN_LOCS = []
    SPC_LOCS = ['InChI=1S/H2O/h1H2', 0, 1]
    THY_LOCS = ['gaussian09', 'wb97xd', 'cc-pvtz', 'R']
    CNF_LOCS = [autofile.schema.generate_new_conformer_id()]
    ZMA_LOCS = [0]
    SCN_LOCS = [['D2'], [2.15]]

    CNF_FS = build_fs(
        RUN_PREFIX, SAVE_PREFIX, 'CONFORMER',
        spc_locs=SPC_LOCS,
        thy_locs=THY_LOCS[1:])

    CNF_RUN_PATH = CNF_FS[0][-1].path(CNF_LOCS)
    CNF_SAVE_PATH = CNF_FS[1][-1].path(CNF_LOCS)
    CNF_FS[0][-1].create(CNF_LOCS)
    CNF_FS[1][-1].create(CNF_LOCS)
    print('cnf fs', CNF_FS)
    print('cnf run', CNF_RUN_PATH)
    print('cnf save', CNF_SAVE_PATH)

    ZMA_FS = build_fs(
        CNF_RUN_PATH, CNF_SAVE_PATH, 'ZMATRIX')

    ZMA_RUN_PATH = ZMA_FS[0][-1].path(ZMA_LOCS)
    ZMA_SAVE_PATH = ZMA_FS[1][-1].path(ZMA_LOCS)
    ZMA_FS[0][-1].create(ZMA_LOCS)
    ZMA_FS[1][-1].create(ZMA_LOCS)
    print('zma run', ZMA_RUN_PATH)
    print('zma save', ZMA_SAVE_PATH)

    SCN_FS = build_fs(
        CNF_RUN_PATH, CNF_SAVE_PATH, 'SCAN',
        zma_locs=ZMA_LOCS)

    SCN_RUN_PATH = SCN_FS[0][-1].path(SCN_LOCS)
    SCN_SAVE_PATH = SCN_FS[1][-1].path(SCN_LOCS)
    SCN_FS[0][-1].create(SCN_LOCS)
    SCN_FS[1][-1].create(SCN_LOCS)
    print('scn run', SCN_RUN_PATH)
    print('scn save', SCN_SAVE_PATH)


if __name__ == '__main__':
    test__build_fs()

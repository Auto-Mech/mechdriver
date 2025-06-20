""" elstruct parameters
"""

from automol.util import sort_by_list
from elstruct import pclass
from elstruct import option


def standard_case(name):
    """ Reformat name in standard lowercase capitalization.

        :param name: name to reformat
        :type name: str
        :rtype: str
    """
    return name.lower() if isinstance(name, str) else name


class Module():
    """ elstruct module names """
    WRITER = 'writer'
    READER = 'reader'


class Program():
    """ Programs supported in elstruct """
    CFOUR2 = 'cfour2'
    GAUSSIAN09 = 'gaussian09'
    GAUSSIAN03 = 'gaussian03'
    GAUSSIAN16 = 'gaussian16'
    MOLPRO2021 = 'molpro2021'
    MOLPRO2015 = 'molpro2015'
    MRCC2018 = 'mrcc2018'
    NWCHEM6 = 'nwchem6'
    ORCA4 = 'orca4'
    PSI4 = 'psi4'
    QCHEM5 = 'qchem5'


def programs():
    """ List all electronic structure backend programs that are supported.
    """
    return pclass.all_values(Program)


def is_program(prog):
    """ Asssess if a program is supported as a backend program.

        :param prog: name of program to test
        :type prog: str
        :rtype: bool
    """
    return standard_case(prog) in programs()


class Reference():
    """ References for density functional and wavefunction methods. """
    RHF = 'rhf'
    UHF = 'uhf'
    ROHF = 'rohf'
    RKS = 'rks'
    UKS = 'uks'


class Method():
    """ Program specific names for various electronic structure methods
        as well as their references.

        (name, {program: singlet name,
                         multiplet name,
                         singlet orb modes,
                         multiplet orb modes})
    """

    HF = ('hf',
          {Program.CFOUR2: (
              'rhf', 'hf',
              ('R',), ('U', 'R')),
           Program.GAUSSIAN09: (
               'hf', 'hf',
               ('R',), ('U', 'R')),
           Program.GAUSSIAN03: (
               'hf', 'hf',
               ('R',), ('U', 'R')),
           Program.GAUSSIAN16: (
               'hf', 'hf',
               ('R',), ('U', 'R')),
           Program.MOLPRO2015: (
               'hf', 'hf',
               ('R',), ('U', 'R')),
           Program.MOLPRO2021: (
               'hf', 'hf',
               ('R',), ('U', 'R')),
           Program.MRCC2018: (
               'hf', 'hf',
               ('R',), ('U', 'R')),
           Program.ORCA4: (
               'hf', 'uhf',
               ('R',), ('U', 'R')),
           Program.PSI4: (
               'hf', 'hf',
               ('R',), ('U', 'R'))})

    class SemiEmpirical():
        PM3 = ('pm3',
            {Program.GAUSSIAN16: (
                'pm3', 'pm3',
                ('R',), ('U', 'R')),
            Program.GAUSSIAN09: (
                'pm3', 'pm3',
                ('R',), ('U', 'R'))})
        PM6 = ('pm6',
            {Program.GAUSSIAN16: (
                'pm6', 'pm6',
                ('R',), ('U', 'R')),
            Program.GAUSSIAN09: (
                'pm6', 'pm6',
                ('R',), ('U', 'R'))})

    class Corr():
        """ Correlated method names """
        MP2 = ('mp2',
               {Program.CFOUR2: (
                   'mp2', 'mp2',
                   ('R',), ('U', 'R')),
                Program.GAUSSIAN09: (
                    'mp2', 'mp2',
                    ('R',), ('U', 'R')),
                Program.GAUSSIAN03: (
                    'mp2', 'mp2',
                    ('R',), ('U', 'R')),
                Program.GAUSSIAN16: (
                    'mp2', 'mp2',
                    ('R',), ('U', 'R')),
                Program.MOLPRO2015: (
                    'mp2', 'ump2',
                    ('R',), ('U', 'R')),
                Program.MOLPRO2021: (
                    'mp2', 'ump2',
                    ('R',), ('U', 'R')),
                Program.MRCC2018: (
                    'mp2', 'mp2',
                    ('R',), ('U', 'R')),
                Program.ORCA4: (
                    'mp2', 'mp2',
                    ('R',), ('U', 'R')),
                Program.PSI4: (
                    'mp2', 'mp2',
                    ('R',), ('U', 'R'))})
        CCSD = ('ccsd',
                {Program.CFOUR2: (
                    'ccsd', 'ccsd',
                    ('R',), ('U', 'R',)),
                 Program.MOLPRO2015: (
                     'ccsd', 'uccsd',
                     ('R',), ('R', 'R')),
                 Program.MOLPRO2021: (
                     'ccsd', 'uccsd',
                     ('R',), ('R', 'R')),
                 Program.PSI4: (
                    'ccsd', 'ccsd',
                    ('R',), ('U', 'R',))})
        CCSD_T = ('ccsd(t)',
                  {Program.CFOUR2: (
                      'ccsd(t)', 'ccsd(t)',
                      ('R',), ('R',)),
                   Program.MOLPRO2015: (
                       'ccsd(t)', 'uccsd(t)',
                       ('R',), ('R',)),
                   Program.MOLPRO2021: (
                       'ccsd(t)', 'uccsd(t)',
                       ('R',), ('R',)),
                   Program.MRCC2018: (
                       'ccsd(t)', 'ccsd(t)',
                       ('R',), ('R', 'R')),
                   Program.PSI4: (
                        'ccsd(t)', 'ccsd(t)',
                        ('R',), ('U', 'R'))})
        CCSDT = ('ccsdt',
                 {Program.MOLPRO2015: (
                     'mrcc,method=ccsdt', 'mrcc,method=ccsdt',
                     ('R',), ('U', 'R')),
                  Program.MOLPRO2021: (
                     'mrcc,method=ccsdt', 'mrcc,method=ccsdt',
                     ('R',), ('U', 'R')),
                  Program.MRCC2018: (
                      'ccsdt', 'ccsdt',
                      ('R',), ('R', 'R'))})
        CCSDT_Q = ('ccsdt(q)',
                   {Program.MOLPRO2015: (
                       'mrcc,method=ccsdt(q)', 'mrcc,method=ccsdt(q)',
                       ('R',), ('U', 'R')),
                    Program.MOLPRO2021: (
                       'mrcc,method=ccsdt(q)', 'mrcc,method=ccsdt(q)',
                       ('R',), ('U', 'R')),
                    Program.MRCC2018: (
                        'ccsdt(q)', 'ccsdt(q)',
                        ('R',), ('R', 'R'))})
        MP2_F12 = ('mp2-f12',
                   {Program.MOLPRO2015: (
                       'mp2-f12', 'ump2-f12',
                       ('R',), ('R',)),
                    Program.MOLPRO2021: (
                       'mp2-f12', 'ump2-f12',
                       ('R',), ('R',))})
        MP2_F12_RI = ('mp2-f12-ri',
                      {Program.ORCA4: (
                        'mp2-f12-ri', 'mp2-f12-ri',
                        ('R',), ('U',))})
        MP2_F12D_RI = ('mp2-f12d-ri',
                       {Program.ORCA4: (
                        'mp2-f12d-ri', 'mp2-f12d-ri',
                        ('R',), ('U',))})
        CCSD_F12 = ('ccsd-f12',
                    {Program.MOLPRO2015: (
                        'ccsd-f12', 'uccsd-f12',
                        ('R',), ('R',)),
                     Program.MOLPRO2021: (
                        'ccsd-f12', 'uccsd-f12',
                        ('R',), ('R',))})
        CCSD_T_F12 = ('ccsd(t)-f12',
                      {Program.MOLPRO2015: (
                          'ccsd(t)-f12', 'uccsd(t)-f12',
                          ('R',), ('R',)),
                       Program.MOLPRO2021: (
                          'ccsd(t)-f12', 'uccsd(t)-f12',
                          ('R',), ('R',))})
        CCSD_T_F12_RI = ('ccsd(t)-f12-ri',
                         {Program.ORCA4: (
                            'ccsd(t)-f12/ri', 'ccsd(t)-f12/ri',
                            ('R',), ('U',))})

    class MultiRef():
        """ Multireference electronic structure methods
        """
        CASSCF = ('casscf',
                  {Program.MOLPRO2015: (
                      'casscf', 'casscf',
                      ('R',), ('R', 'R')),
                   Program.MOLPRO2021: (
                      'casscf', 'casscf',
                      ('R',), ('R', 'R'))})
        CASPT2 = ('caspt2',
                  {Program.MOLPRO2015: (
                      'rs2', 'rs2',
                      ('R',), ('R', 'R')),
                   Program.MOLPRO2021: (
                      'rs2', 'rs2',
                      ('R',), ('R', 'R'))})
        CASPT2I = ('caspt2i',
                   {Program.MOLPRO2015: (
                       'rs2', 'rs2',
                       ('R',), ('R', 'R')),
                    Program.MOLPRO2021: (
                       'rs2', 'rs2',
                       ('R',), ('R', 'R'))})
        CASPT2C = ('caspt2c',
                   {Program.MOLPRO2015: (
                       'rs2c', 'rs2c',
                       ('R',), ('R', 'R')),
                    Program.MOLPRO2021: (
                       'rs2c', 'rs2c',
                       ('R',), ('R', 'R'))})
        MRCISDQ = ('mrcisd_q',
                   {Program.MOLPRO2015: (
                       'mrci', 'mrci',
                       ('R',), ('R', 'R')),
                    Program.MOLPRO2021: (
                       'mrci', 'mrci',
                       ('R',), ('R', 'R'))})
        MRCIF12 = ('mrcisd_q-f12',
                   {Program.MOLPRO2015: (
                       'mrci-f12', 'mrci-f12',
                       ('R',), ('R', 'R')),
                    Program.MOLPRO2021: (
                       'mrci-f12', 'mrci-f12',
                       ('R',), ('R', 'R'))})

    class Dft():
        """ Density functional theory method names """
        BP86 = ('bp86',
                 {Program.PSI4: (
                     'BP86', 'BP86',
                     ('R',), ('U',)),
                  Program.GAUSSIAN09: (
                      'bp86', 'bp86',
                      ('R',), ('U',)),
                  Program.GAUSSIAN03: (
                      'bp86', 'bp86',
                      ('R',), ('U',)),
                  Program.GAUSSIAN16: (
                      'bp86', 'bp86',
                      ('R',), ('U',))})
        B3LYP = ('b3lyp',
                 {Program.PSI4: (
                     'B3LYP', 'B3LYP',
                     ('R',), ('U',)),
                  Program.GAUSSIAN09: (
                      'b3lyp', 'b3lyp',
                      ('R',), ('U',)),
                  Program.GAUSSIAN03: (
                      'b3lyp', 'b3lyp',
                      ('R',), ('U',)),
                  Program.GAUSSIAN16: (
                      'b3lyp', 'b3lyp',
                      ('R',), ('U',))})
        WB97XD = ('wb97xd',
                  {Program.PSI4: (
                      'WB97X-D', 'WB97X-D',
                      ('R',), ('U',)),
                   Program.QCHEM5: (
                       'wb97X-D', 'wb97X-D',
                       ('R',), ('U',)),
                   Program.GAUSSIAN09: (
                       'wb97xd', 'wb97xd',
                       ('R',), ('U',)),
                   Program.GAUSSIAN03: (
                       'wb97xd', 'wb97xd',
                       ('R',), ('U',)),
                   Program.GAUSSIAN16: (
                       'wb97xd', 'wb97xd',
                       ('R',), ('U',))})
        M062X = ('m062x',
                 {Program.PSI4: (
                     'M06-2X', 'M06-2X',
                     ('R',), ('U',)),
                  Program.GAUSSIAN09: (
                      'm062x', 'm062x',
                      ('R',), ('U',)),
                  Program.GAUSSIAN03: (
                      'm062x', 'm062x',
                      ('R',), ('U',)),
                  Program.GAUSSIAN16: (
                      'm062x', 'm062x',
                      ('R',), ('U',))})
        REVDSD = ('revdsd',
                 {Program.GAUSSIAN16: (
                      'dsdpbep86', 'udsdpbep86',
                      ('R',), ('U',))})        
        B2PLYPD3 = ('b2plypd3',
                    {Program.GAUSSIAN09: (
                        'b2plypd3', 'b2plypd3',
                        ('R',), ('U',)),
                     Program.GAUSSIAN03: (
                        'b2plypd3', 'b2plypd3',
                        ('R',), ('U',)),
                     Program.GAUSSIAN16: (
                         'b2plypd3', 'b2plypd3',
                         ('R',), ('U',))})

    class ModPrefix():
        """ Allowed Prefixes for methods
            (full prefix name, prefix used in method names)
        """
        ALL_ELEC = ('all-electron', 'ae-')
        DBOC = ('diagonal-Born-Oppenheimer-correction', 'dboc-')
        DF = ('density-fitting', 'df-')
        L_PNO = ('pair-natural-orbital-local', 'pno-l')
        REL_DKH = ('relativistic-dkh', 'dkh-')

    @classmethod
    def contains(cls, name):
        """ Assess if provided method is a part of this class.

            :param cls: class object
            :type cls: obj
            :param name: name of method
            :type name: str
        """

        name = standard_case(name)
        names = [row[0] for row in pclass.all_values(cls)]

        return name in names

    @classmethod
    def is_semi_empirical(cls, name):
        """ Assess if a method is a semi-empirical method.

            :param cls: class object
            :type cls: obj
            :param name: name of method
            :type name: str
        """

        name = standard_case(name)
        semi_emp_names = [row[0] for row in pclass.all_values(cls.SemiEmpirical)]

        return name in semi_emp_names

    @classmethod
    def is_correlated(cls, name):
        """ Assess if a method is a single-reference correlated method.

            :param cls: class object
            :type cls: obj
            :param name: name of method
            :type name: str
        """

        name = standard_case(name)
        corr_names = [row[0] for row in pclass.all_values(cls.Corr)]

        return name in corr_names

    @classmethod
    def is_multiref(cls, name):
        """ Assess if a method is a multi-reference correlated method.

            :param cls: class object
            :type cls: obj
            :param name: name of method
            :type name: str
        """

        name = standard_case(name)
        multiref_names = [row[0] for row in pclass.all_values(cls.MultiRef)]

        return name in multiref_names

    @staticmethod
    def is_casscf(name):
        """ Assess if a method is CASSCF.

            :param cls: class object
            :type cls: obj
            :param name: name of method
            :type name: str
        """
        return standard_case(name) == 'casscf'

    @classmethod
    def is_standard_dft(cls, name):
        """ Assess if a method corresponds to a density functional
            currently defined in elstruct.

            :param cls: class object
            :type cls: obj
            :param name: name of method
            :type name: str
        """

        name = standard_case(name)
        dft_names = [row[0] for row in pclass.all_values(cls.Dft)]

        return name in dft_names

    @staticmethod
    def is_nonstandard_dft(name):
        """ Assess if a method corresponds to a user-defined
            density functional (indicated by 'dft:<name>').

            :param name: name of method
            :type name: str
        """
        return name.lower().startswith('dft:')

    @classmethod
    def is_dft(cls, name):
        """ Assess if a method corresponds to a density functional
            (either standard or non-standard).
        """
        return cls.is_standard_dft(name) or cls.is_nonstandard_dft(name)

    @classmethod
    def nonstandard_dft_name(cls, name):
        """ Extract the name of a non-standard density functional
            (indicated by 'dft:<name>').

            :param cls: class object
            :type cls: obj
            :param name: name of method
            :type name: str
        """
        assert cls.is_nonstandard_dft(name)
        return name[4:]

    # Handle prefixes to methods
    @classmethod
    def ordered_prefix_lst(cls, pfx_lst):
        """ Returns a list of prefix corresponding to a certain order
            to be used by other functions
        """
        ord_lst = (cls.ModPrefix.ALL_ELEC[0],
                   cls.ModPrefix.DBOC[0],
                   cls.ModPrefix.DF[0],
                   cls.ModPrefix.L_PNO[0],
                   cls.ModPrefix.REL_DKH[0])
        return sort_by_list(pfx_lst, ord_lst)

    @classmethod
    def evaluate_method_type(cls, name):
        """ Analyze a method name and alter to return the core method
            and flags signaling if it is modified.

            All modifications presented as 'mod-' as a prefix to the core name

            Ex: ae-ccsd => ccsd, all_electron=True
            Ex: dkh-ae-ccsd => ccsd, all_electron=True, douglas_kroll=True

            :param name: name of method
            :type name: str
        """

        def _demodify_name(_name, mod):
            """ Assess if the electronic structure method modifier exists in
                the method name and remove it if so.
            """
            if mod in _name:
                _mod_name = _name.replace(mod, '')
                _has_mod = True
            else:
                _mod_name = _name
                _has_mod = False
            return _mod_name, _has_mod

        # Assess what method modifiers exist and remove prefix from name
        _core = name
        _core, has_ae = _demodify_name(_core, cls.ModPrefix.ALL_ELEC[1])
        _core, has_dkh = _demodify_name(_core, cls.ModPrefix.REL_DKH[1])
        _core, has_pnol = _demodify_name(_core, cls.ModPrefix.L_PNO[1])
        _core, has_dboc = _demodify_name(_core, cls.ModPrefix.DBOC[1])
        _core, has_df = _demodify_name(_core, cls.ModPrefix.DF[1])

        # Build list of prefixes in standard order
        pfxs = ()
        if has_ae:
            pfxs += (cls.ModPrefix.ALL_ELEC[0],)
        if has_dkh:
            pfxs += (cls.ModPrefix.REL_DKH[0],)
        if has_pnol:
            pfxs += (cls.ModPrefix.L_PNO[0],)
        if has_dboc:
            pfxs += (cls.ModPrefix.DBOC[0],)
        if has_df:
            pfxs += (cls.ModPrefix.DF[0],)

        ord_pfxs = cls.ordered_prefix_lst(pfxs)

        return _core, ord_pfxs


def program_methods_info(prog):
    """ List methods available for a given program, with their information.

        :param prog: electronic structure program name
        :type prog: str
        :rtype: dict[str: str]
    """
    prog = standard_case(prog)
    return {row[0]: row[1][prog] for row in pclass.all_values(Method)
            if prog in row[1]}


def program_methods(prog):
    """ List methods available for a given program.

        :param prog: electronic structure program name
        :type prog: str
        :rtype: tuple(str)
    """
    return tuple(sorted(program_methods_info(prog)))


def program_dft_methods(prog):
    """ List density functional theory methods available for a given program.

        :param prog: electronic structure program name
        :type prog: str
        :rtype: tuple(str)
    """
    prog_methods = program_methods(prog)
    return tuple(method for method in prog_methods
                 if Method.is_standard_dft(method))


def program_nondft_methods(prog):
    """ List Hartree-Fock wavefunction methods available for a given program.

        :param prog: electronic structure program name
        :type prog: str
        :rtype: tuple(str)
    """
    prog_methods = program_methods(prog)
    return tuple(method for method in prog_methods
                 if not Method.is_standard_dft(method))


def program_method_name(prog, method, singlet=True):
    """ Obtain the name of a given method specific to the program provided.

        :param prog: electronic structure program name
        :type prog: str
        :param method: electronic structure method name
        :type method: str
        :param singlet: Parameter specifying name for singlet species
        :type singlet: bool
        :rtype: tuple(str)
    """

    prog = standard_case(prog)

    if Method.is_nonstandard_dft(method):
        method = Method.nonstandard_dft_name(method)
    else:
        method = standard_case(method)
        prog_method_dct = program_methods_info(prog)
        assert method in prog_method_dct
        name = (prog_method_dct[method][0] if singlet else
                prog_method_dct[method][1])
        method = method if name is None else name

    return method


def program_method_orbital_types(prog, method, singlet):
    """ Obtain the available orbital modes of a given method
        specific to the program provided.

        :param prog: electronic structure program name
        :type prog: str
        :param method: electronic structure method name
        :type method: str
        :param singlet: Parameter specifying name for singlet species
        :type singlet: bool
        :rtype: tuple(str)
    """

    prog = standard_case(prog)
    method = standard_case(method)
    prog_method_dct = program_methods_info(prog)
    assert method in prog_method_dct
    orb_types = (prog_method_dct[method][2] if singlet else
                 prog_method_dct[method][3])

    return orb_types


def is_program_method(prog, method):
    """ Assess if the given method is valid for the specific program.

        :param prog: electronic structure program name
        :type prog: str
        :param method: electronic structure method name
        :type method: str
    """

    prog = standard_case(prog)
    method = standard_case(method)

    return method in program_methods(prog)


def is_program_method_orbital_type(prog, method, singlet, orb_type):
    """ is this a valid method for this program?
    """
    prog = standard_case(prog)
    method = standard_case(method)
    assert isinstance(singlet, bool)
    return orb_type in program_method_orbital_types(prog, method, singlet)


class Basis():
    """ Electronic structure basis sets, defined internally in elstruct,
        as well as dictionary that maps the name of the basis set into the
        name for all of the implemented electronic structure programs.

        (name, {program: name})
    """
    NONE = (None, {Program.CFOUR2: None,
                     Program.GAUSSIAN09: None,
                     Program.GAUSSIAN03: None,
                     Program.GAUSSIAN16: None,
                     Program.MOLPRO2015: None,
                     Program.MOLPRO2021: None,
                     Program.MRCC2018: None,
                     Program.NWCHEM6: None,
                     Program.ORCA4: None,
                     Program.PSI4: None})

    STO3G = ('sto-3g', {Program.CFOUR2: None,
                        Program.GAUSSIAN09: None,
                        Program.GAUSSIAN03: None,
                        Program.GAUSSIAN16: None,
                        Program.MOLPRO2015: None,
                        Program.MOLPRO2021: None,
                        Program.MRCC2018: None,
                        Program.NWCHEM6: None,
                        Program.ORCA4: None,
                        Program.PSI4: None})

    DEF2SV_P = ('def2-sv(p)', {Program.CFOUR2: None,
                               Program.GAUSSIAN09: 'def2svpp',
                               Program.GAUSSIAN03: 'def2svpp',
                               Program.GAUSSIAN16: 'def2svpp',
                               Program.MOLPRO2015: None,
                               Program.MOLPRO2021: None,
                               Program.MRCC2018: None,
                               Program.NWCHEM6: None,
                               Program.ORCA4: None,
                               Program.PSI4: None})

    DEF2SVP = ('def2-svp', {Program.CFOUR2: None,
                            Program.GAUSSIAN09: 'def2svp',
                            Program.GAUSSIAN03: 'def2svp',
                            Program.GAUSSIAN16: 'def2svp',
                            Program.MOLPRO2015: None,
                            Program.MOLPRO2021: None,
                            Program.MRCC2018: None,
                            Program.NWCHEM6: None,
                            Program.ORCA4: None,
                            Program.PSI4: None})

    DEF2TZVP = ('def2-tzvp', {Program.CFOUR2: None,
                              Program.GAUSSIAN09: 'def2tzvp',
                              Program.GAUSSIAN03: 'def2tzvp',
                              Program.GAUSSIAN16: 'def2tzvp',
                              Program.MOLPRO2015: None,
                              Program.MOLPRO2021: None,
                              Program.MRCC2018: None,
                              Program.NWCHEM6: None,
                              Program.ORCA4: None,
                              Program.PSI4: None})

    class Pople:
        """ Pople basis sets """
        P321 = ('3-21g', {Program.CFOUR2: None,
                          Program.GAUSSIAN09: None,
                          Program.GAUSSIAN03: None,
                          Program.GAUSSIAN16: None,
                          Program.MOLPRO2015: None,
                          Program.MOLPRO2021: None,
                          Program.MRCC2018: None,
                          Program.NWCHEM6: None,
                          Program.ORCA4: None,
                          Program.PSI4: None})
        P321S = ('3-21g*', {Program.PSI4: None})
        P631 = ('6-31g', {Program.CFOUR2: None,
                          Program.GAUSSIAN09: None,
                          Program.GAUSSIAN03: None,
                          Program.GAUSSIAN16: None,
                          Program.MOLPRO2015: None,
                          Program.MOLPRO2021: None,
                          Program.MRCC2018: None,
                          Program.NWCHEM6: None,
                          Program.ORCA4: None,
                          Program.PSI4: None})
        P631S = ('6-31g*', {Program.CFOUR2: None,
                            Program.GAUSSIAN09: None,
                            Program.GAUSSIAN03: None,
                            Program.GAUSSIAN16: None,
                            Program.MOLPRO2015: None,
                            Program.MOLPRO2021: None,
                            Program.MRCC2018: None,
                            Program.NWCHEM6: None,
                            Program.ORCA4: None,
                            Program.PSI4: None,
                            Program.QCHEM5: None})
        P631PS = ('6-31+g*', {Program.CFOUR2: None,
                              Program.GAUSSIAN09: None,
                              Program.GAUSSIAN03: None,
                              Program.GAUSSIAN16: None,
                              Program.MOLPRO2015: None,
                              Program.MOLPRO2021: None,
                              Program.MRCC2018: None,
                              Program.NWCHEM6: None,
                              Program.ORCA4: None,
                              Program.PSI4: None})
        P6311SS = ('6-311g**', {Program.CFOUR2: None,
                                Program.GAUSSIAN09: None,
                                Program.GAUSSIAN03: None,
                                Program.GAUSSIAN16: None,
                                Program.MOLPRO2015: None,
                                Program.MOLPRO2021: None,
                                Program.MRCC2018: None,
                                Program.NWCHEM6: None,
                                Program.ORCA4: None,
                                Program.PSI4: None})
        P6311PSS = ('6-311+g**', {Program.CFOUR2: None,
                                  Program.GAUSSIAN09: None,
                                  Program.GAUSSIAN03: None,
                                  Program.GAUSSIAN16: None,
                                  Program.MOLPRO2015: None,
                                  Program.MOLPRO2021: None,
                                  Program.MRCC2018: None,
                                  Program.NWCHEM6: None,
                                  Program.ORCA4: None,
                                  Program.PSI4: None})
        P6311PPSS = ('6-311++g**', {Program.CFOUR2: None,
                                    Program.GAUSSIAN09: None,
                                    Program.GAUSSIAN03: None,
                                    Program.GAUSSIAN16: None,
                                    Program.MOLPRO2015: None,
                                    Program.MOLPRO2021: None,
                                    Program.MRCC2018: None,
                                    Program.NWCHEM6: None,
                                    Program.ORCA4: None,
                                    Program.PSI4: None})

    class Dunning():
        """ Dunning basis sets """
        D = ('cc-pvdz', {Program.CFOUR2: None,
                         Program.GAUSSIAN09: None,
                         Program.GAUSSIAN03: None,
                         Program.GAUSSIAN16: None,
                         Program.MOLPRO2015: None,
                         Program.MOLPRO2021: None,
                         Program.MRCC2018: None,
                         Program.NWCHEM6: None,
                         Program.ORCA4: None,
                         Program.PSI4: None})
        T = ('cc-pvtz', {Program.CFOUR2: None,
                         Program.GAUSSIAN09: None,
                         Program.GAUSSIAN03: None,
                         Program.GAUSSIAN16: None,
                         Program.MOLPRO2015: None,
                         Program.MOLPRO2021: None,
                         Program.MRCC2018: None,
                         Program.NWCHEM6: None,
                         Program.ORCA4: None,
                         Program.PSI4: None})
        Q = ('cc-pvqz', {Program.CFOUR2: None,
                         Program.GAUSSIAN09: None,
                         Program.GAUSSIAN03: None,
                         Program.GAUSSIAN16: None,
                         Program.MOLPRO2015: None,
                         Program.MOLPRO2021: None,
                         Program.MRCC2018: None,
                         Program.NWCHEM6: None,
                         Program.ORCA4: None,
                         Program.PSI4: None})
        P = ('cc-pv5z', {Program.CFOUR2: None,
                         Program.GAUSSIAN09: None,
                         Program.GAUSSIAN03: None,
                         Program.GAUSSIAN16: None,
                         Program.MOLPRO2015: None,
                         Program.MOLPRO2021: None,
                         Program.MRCC2018: None,
                         Program.NWCHEM6: None,
                         Program.ORCA4: None,
                         Program.PSI4: None})
        CD = ('cc-pcvdz', {Program.CFOUR2: None,
                           Program.GAUSSIAN09: None,
                           Program.GAUSSIAN03: None,
                           Program.GAUSSIAN16: None,
                           Program.MOLPRO2015: None,
                           Program.MOLPRO2021: None,
                           Program.MRCC2018: None,
                           Program.NWCHEM6: None,
                           Program.ORCA4: None,
                           Program.PSI4: None})
        CT = ('cc-pcvtz', {Program.CFOUR2: None,
                           Program.GAUSSIAN09: None,
                           Program.GAUSSIAN03: None,
                           Program.GAUSSIAN16: None,
                           Program.MOLPRO2015: None,
                           Program.MOLPRO2021: None,
                           Program.MRCC2018: None,
                           Program.NWCHEM6: None,
                           Program.ORCA4: None,
                           Program.PSI4: None})
        CQ = ('cc-pcvqz', {Program.CFOUR2: None,
                           Program.GAUSSIAN09: None,
                           Program.GAUSSIAN03: None,
                           Program.GAUSSIAN16: None,
                           Program.MOLPRO2015: None,
                           Program.MOLPRO2021: None,
                           Program.MRCC2018: None,
                           Program.NWCHEM6: None,
                           Program.ORCA4: None,
                           Program.PSI4: None})
        CP = ('cc-pcv5z', {Program.CFOUR2: None,
                           Program.GAUSSIAN09: None,
                           Program.GAUSSIAN03: None,
                           Program.GAUSSIAN16: None,
                           Program.MOLPRO2015: None,
                           Program.MOLPRO2021: None,
                           Program.MRCC2018: None,
                           Program.NWCHEM6: None,
                           Program.ORCA4: None,
                           Program.PSI4: None})

        class Aug():
            """ Augmented Dunning basis sets """
            AD = ('aug-cc-pvdz', {Program.CFOUR2: None,
                                  Program.GAUSSIAN09: None,
                                  Program.GAUSSIAN03: None,
                                  Program.GAUSSIAN16: None,
                                  Program.MOLPRO2015: None,
                                  Program.MOLPRO2021: None,
                                  Program.MRCC2018: None,
                                  Program.NWCHEM6: None,
                                  Program.ORCA4: None,
                                  Program.PSI4: None})
            AT = ('aug-cc-pvtz', {Program.CFOUR2: None,
                                  Program.GAUSSIAN09: None,
                                  Program.GAUSSIAN03: None,
                                  Program.GAUSSIAN16: None,
                                  Program.MOLPRO2015: None,
                                  Program.MOLPRO2021: None,
                                  Program.MRCC2018: None,
                                  Program.NWCHEM6: None,
                                  Program.ORCA4: None,
                                  Program.PSI4: None})
            AQ = ('aug-cc-pvqz', {Program.CFOUR2: None,
                                  Program.GAUSSIAN09: None,
                                  Program.GAUSSIAN03: None,
                                  Program.GAUSSIAN16: None,
                                  Program.MOLPRO2015: None,
                                  Program.MOLPRO2021: None,
                                  Program.MRCC2018: None,
                                  Program.NWCHEM6: None,
                                  Program.ORCA4: None,
                                  Program.PSI4: None})
            A5 = ('aug-cc-pv5z', {Program.CFOUR2: None,
                                  Program.GAUSSIAN09: None,
                                  Program.GAUSSIAN03: None,
                                  Program.GAUSSIAN16: None,
                                  Program.MOLPRO2015: None,
                                  Program.MOLPRO2021: None,
                                  Program.MRCC2018: None,
                                  Program.NWCHEM6: None,
                                  Program.ORCA4: None,
                                  Program.PSI4: None})
            CAD = ('aug-cc-pcvdz', {Program.CFOUR2: None,
                                    Program.GAUSSIAN09: None,
                                    Program.GAUSSIAN03: None,
                                    Program.GAUSSIAN16: None,
                                    Program.MOLPRO2015: None,
                                    Program.MOLPRO2021: None,
                                    Program.MRCC2018: None,
                                    Program.NWCHEM6: None,
                                    Program.ORCA4: None,
                                    Program.PSI4: None})
            CAT = ('aug-cc-pcvtz', {Program.CFOUR2: None,
                                    Program.GAUSSIAN09: None,
                                    Program.GAUSSIAN03: None,
                                    Program.GAUSSIAN16: None,
                                    Program.MOLPRO2015: None,
                                    Program.MOLPRO2021: None,
                                    Program.MRCC2018: None,
                                    Program.NWCHEM6: None,
                                    Program.ORCA4: None,
                                    Program.PSI4: None})
            CAQ = ('aug-cc-pcvqz', {Program.CFOUR2: None,
                                    Program.GAUSSIAN09: None,
                                    Program.GAUSSIAN03: None,
                                    Program.GAUSSIAN16: None,
                                    Program.MOLPRO2015: None,
                                    Program.MOLPRO2021: None,
                                    Program.MRCC2018: None,
                                    Program.NWCHEM6: None,
                                    Program.ORCA4: None,
                                    Program.PSI4: None})
            CA5 = ('aug-cc-pcv5z', {Program.CFOUR2: None,
                                    Program.GAUSSIAN09: None,
                                    Program.GAUSSIAN03: None,
                                    Program.GAUSSIAN16: None,
                                    Program.MOLPRO2015: None,
                                    Program.MOLPRO2021: None,
                                    Program.MRCC2018: None,
                                    Program.NWCHEM6: None,
                                    Program.ORCA4: None,
                                    Program.PSI4: None})
            JT = ('jun-cc-pvtz', {Program.CFOUR2: None,
                                    Program.GAUSSIAN09: None,
                                    Program.GAUSSIAN03: None,
                                    Program.GAUSSIAN16: None,
                                    Program.MOLPRO2015: None,
                                    Program.MOLPRO2021: None,
                                    Program.MRCC2018: None,
                                    Program.NWCHEM6: None,
                                    Program.ORCA4: None,
                                    Program.PSI4: None})
            JQ = ('jun-cc-pvqz', {Program.CFOUR2: None,
                                    Program.GAUSSIAN09: None,
                                    Program.GAUSSIAN03: None,
                                    Program.GAUSSIAN16: None,
                                    Program.MOLPRO2015: None,
                                    Program.MOLPRO2021: None,
                                    Program.MRCC2018: None,
                                    Program.NWCHEM6: None,
                                    Program.ORCA4: None,
                                    Program.PSI4: None})
            ATd = ('aug-cc-pv(t+d)z', {Program.MOLPRO2015: None})

        class F12():
            """ Dunning F12 basis sets """
            DF = ('cc-pvdz-f12', {Program.MOLPRO2015: None,
                                  Program.MOLPRO2021: None,
                                  Program.ORCA4: None})
            TF = ('cc-pvtz-f12', {Program.MOLPRO2015: None,
                                  Program.MOLPRO2021: None,
                                  Program.ORCA4: None})
            QF = ('cc-pvqz-f12', {Program.MOLPRO2015: None,
                                  Program.MOLPRO2021: None,
                                  Program.ORCA4: None})
            CDF = ('cc-pcvdz-f12', {Program.MOLPRO2015: None,
                                    Program.MOLPRO2021: None})
            CTF = ('cc-pcvtz-f12', {Program.MOLPRO2015: None,
                                    Program.MOLPRO2021: None})
            CQF = ('cc-pcvqz-f12', {Program.MOLPRO2015: None,
                                    Program.MOLPRO2021: None})

    @classmethod
    def contains(cls, name):
        """ Assess if provided basis set is a part of this class.

            :param cls: class object
            :type cls: obj
            :param name: name of basis set
            :type name: str
        """
        name = standard_case(name)
        names = [row[0] for row in pclass.all_values(cls)]
        return name in names

    is_standard_basis = contains

    @staticmethod
    def is_nonstandard_basis(name):
        """ Assess if a basis set corresponds to user-defined
            basis set (indicated by 'basis:<name>').

            :param name: name of basis set
            :type name: str
        """
        return isinstance(name, str) and name.lower().startswith('basis:')

    @classmethod
    def nonstandard_basis_name(cls, name):
        """ Extract the aname of non-standard basis set
            (indicated by 'basis:<name>').
        """
        assert cls.is_nonstandard_basis(name)
        return name[6:]


def program_bases(prog):
    """ List basis sets available for a given program.

        :param prog: electronic structure program name
        :type prog: str
        :rtype: dict[str: str]
    """
    prog = standard_case(prog)
    return {row[0]: row[1][prog] for row in pclass.all_values(Basis)
            if prog in row[1]}


def program_basis_name(prog, basis):
    """ Obtain the name of a given basis set specific to the program provided.

        :param prog: electronic structure program name
        :type prog: str
        :param basis: electronic structure basis set
        :type basis: str
    """

    prog = standard_case(prog)

    if Basis.is_nonstandard_basis(basis):
        basis = Basis.nonstandard_basis_name(basis)
    else:
        basis = standard_case(basis)
        prog_bases = program_bases(prog)
        assert basis in prog_bases
        name = prog_bases[basis]
        basis = basis if name is None else name

    return basis


def is_program_basis(prog, basis):
    """ Assess if the given basis set is valid for the specific program.

        :param prog: electronic structure program name
        :type prog: str
        :param basis: electronic structure basis set
        :type basis: str
    """

    prog = standard_case(prog)
    basis = standard_case(basis)

    return basis in program_bases(prog)


class Job():
    """ Names of electronic structure jobs whose input (output)
        can be written (read).
    """

    ENERGY = 'energy'
    GRADIENT = 'gradient'
    HESSIAN = 'hessian'
    OPTIMIZATION = 'optimization'
    VPT2 = 'vpt2'
    IRCF = 'ircf'
    IRCR = 'ircr'
    MOLPROP = 'molecular_properties'

    @classmethod
    def contains(cls, name):
        """ Assess if provided method is a part of this class.

            :param cls: class object
            :type cls: obj
            :param name: name of job
            :type name: str
        """

        name = standard_case(name)
        names = pclass.all_values(cls)
        # names = [row for row in pclass.all_values(cls)]

        return name in names


class Error():
    """ Types of error messages that can be found in electronic structure
        program output files.
    """
    SCF_NOCONV = 'scf_noconv'
    MCSCF_NOCONV = 'mcscf_noconv'
    CC_NOCONV = 'cc_noconv'
    OPT_NOCONV = 'opt_noconv'
    IRC_NOCONV = 'irc_noconv'
    SYMM_NOFIND = 'symm_nofind'
    LIN_DEP_BASIS = 'linear_dependent_basis'


class Success():
    """ Types of sucess errors that can be found in electronic structure
        program output files.
    """
    SCF_CONV = 'scf_conv'
    CC_CONV = 'cc_conv'
    OPT_CONV = 'opt_conv'
    IRC_CONV = 'irc_conv'


class Option():
    """ Option values for the electronic structure program writer module
    """

    class Mol():
        """ Options for molecules """
        NOSYMM_ = option.create('no_symmetry')

    class Scf():
        """ SCF options (passed to `scf_options`) """
        MAXITER_ = option.create('scf_maxiter', ['num'])
        DIIS_ = option.create('scf_diis', ['bool'])

        class Guess():
            """ Initial SCF guess generation methods """
            CORE = option.create('scf_guess_core')
            HUCKEL = option.create('scf_guess_huckel')
            MIX = option.create('scf_guess_mix')

    class Casscf():
        """ CASSCF options to set active space (passed to `casscf_options`) """
        OCC_ = option.create('casscf_occ', ['num'])
        CLOSED_ = option.create('casscf_closed', ['num'])
        WFN_ = option.create(
            'casscf_wavefunction',
            ['nelec', 'sym', 'spin', 'charge', 'nstates'])

    class Corr():
        """ Options for correlated methods"""
        ALL_ELEC_ = option.create('all_electron')

    class MRCorr():
        """ Correlated multiref method options (passed to `corr_options`) """
        SHIFT_ = option.create('level_shift', ['num'])
        IPEA_ = option.create('ipea', ['num'])

    class Opt():
        """ Optimization options (passed to `job_options`) """
        MAXITER_ = option.create('opt_maxiter', ['num'])

        class Coord():
            """ Corrdinate system to perform optimization in """
            CARTESIAN = option.create('opt_coord_cartesian')
            ZMATRIX = option.create('opt_coord_zmatrix')
            REDUNDANT = option.create('opt_coord_redundant')

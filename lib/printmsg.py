"""
  Library of Runtime Messages
"""


def tsk_msg(tsk, spc, thy_info, ini_thy_info):
    """ print a task message
    """
    print('\nTask {} \t {}//{} \t Species {}'.format(
        tsk, '/'.join(thy_info), '/'.join(ini_thy_info), spc))


def sadpt_tsk_msg(tsk, sadpt, spc_dct, thy_info, ini_thy_info):
    """ print a task message for a TS
    """
    print('Task {} \t for {} \t {}//{} \t {} = {}'.format(
        tsk, sadpt,
        '/'.join(thy_info),
        '/'.join(ini_thy_info),
        '+'.join(spc_dct[sadpt]['reacs']),
        '+'.join(spc_dct[sadpt]['prods'])))


def ini_info_noavail_msg(tsk):
    """ print a message saying initial input info not available for task
    """
    print(
        'Initial level of theory for conformers must be ',
        'run before {} '.format(tsk))


def run_tsk_msg(tsk):
    """ print a message saying a task is running
    """
    print('running task {}'.format(tsk))


def program_header(driver):
    """ print the header for a program
    """
    header_dct = {
        'amech': AMECH_MSG,
        'ktp': KTP_MSG,
        'thermo': THM_MSG,
        'es': ES_MSG
    }
    print(header_dct[driver]+'\n\n')


AMECH_MSG = """
          ================================================================
          ==                        AUTOMECHANIC                        ==
          ===         Andreas Copan, Sarah Elliott, Kevin Moore,       ===
          ===     Daniel Moberg, Carlo Cavallotti, Yuri Georgievski,   ===
          ==       Ahren Jasper, Murat Keceli, Stephen Klippenstein     ==
          ================================================================"""

KTP_MSG = """
          ================================================================
          ==                          KTPDRIVER                         ==
          ===         Sarah Elliott, Kevin Moore, Andreas Copan,       ===
          ===      Daniel Moberg, Carlo Cavallotti, Yuri Georgievski,  ===
          ==            Ahren Jasper, Stephen Klippenstein              ==
          ================================================================"""

THM_MSG = """
          ================================================================
          ==                        THERMODRIVER                        ==
          ===         Sarah Elliott, Kevin Moore, Andreas Copan,       ===
          ===    Murat Keceli, Yuri Georgievski, Stephen Klippenstein   ==
          ================================================================"""

ES_MSG = """
          ================================================================
          ==                          ESDRIVER                          ==
          ====        Sarah Elliott, Andreas Copan, Kevin Moore,      ====
          ==            Carlo Cavolotti, Stephen Klippenstein           ==
          ================================================================"""


def random_cute_animal():
    import random
    msg = random.choice(["""\n\t\t
                                  _,--._
                                ,'      `.
                        |\     /          \     /|
                        )o),/ ( ,--,  ,--, ) \.(o(
                       /o/// /|            |\ \\ \\o\\
                      / / |\ \(   .----,   )/ /| \ \\
                      | | \o`-/    `--'    \-'o/ | |
                      \ \  `,'              `.'  / /
                   \.  \ `-'  ,'|   /\   |`.  `-' /  ,/
                    \`. `.__,' /   /  \   \ `.__,' ,'/
                     \o\     ,'  ,'    `.  `.     /o/
                      \o`---'  ,'        `.  `---'o/
                       `.____,'           `.____,'  """,

     """\n    
                          ,,,         ,,,
                        ;"   ^;     ;'   ",
                       ;    s$$$$$$$s      ;
                        ,  ss$$$$$$$$$$s  ,'
                        ;s$$$$$$$$$$$$$$$
                        $$$$$$$$$$$$$$$$$$
                       $$$$P""Y$$$Y""W$$$$$
                       $$$$  0"$$$"0  $$$$$
                       $$$$  .$$$$$.  $$$$
                        $$$$$$$$$$$$$$$$$
                         "Y$$$"'*'"$$$Y"
                            "$$b.d$$"        """,
    """\n
                                       _.---~-~-~~-..
                   ..       __.    .-~               ~-.
                   ((\     /   `}.~                     `.
                    \\\\\   {     }              /     \   \\
                (\   \\\\~~       }             |       }   \\
                 \`.-~-^~     }  ,-,.         |       )    \\
                 (___,    ) _}  (    :        |    / /      `.
                  `----._-~.     _\ \ |_       \   / /- _     -.
                         ~~----~~  \ \| ~~--~~~(  + /     ~-.   '--~.
                                   /  /         \  \         `~-..__ `~__
                                __/  /          _\  )               `~~---'
                              .<___.'         .<___/
    """])
    print(msg)
    print('\n\n\n')

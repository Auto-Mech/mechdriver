""" functions to switch safemode on/off
"""


SAFEMODE = True


def safemode_is_on():
    """ indicates whether or not the safemode global variable is on
    """
    return SAFEMODE


def turn_off_safemode():
    """ turn off the safemode global variable
    """
    global SAFEMODE
    SAFEMODE = False


def turn_on_safemode():
    """ turn off the safemode global variable
    """
    global SAFEMODE
    SAFEMODE = True

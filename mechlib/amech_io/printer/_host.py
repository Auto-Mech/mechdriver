""" Handle process and host info
"""

import autorun
from mechlib.amech_io.printer._print import message


HOST_MSG = """
=========================================================
HOST: {}
PID: {}
========================================================="""


def host_name():
    """ print the host the calculation is running on
    """
    host_node = autorun.host_node()
    pid = autorun.process_id()
    message(HOST_MSG.format(host_node, pid))

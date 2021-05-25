""" Handle process and host info
"""

import os
import subprocess
from mechlib.amech_io.printer._print import message


HOST_MSG = """
=========================================================
HOST: {}
PID: {}
========================================================="""


def host_name():
    """ print the host the calculation is running on
    """
    host_node = _get_host_node()
    pid = _get_pid()
    message(HOST_MSG.format(host_node, pid))


def _get_host_node():
    """ get the nodes
    """
    proc = subprocess.Popen(['hostname'], stdout=subprocess.PIPE)
    host_node = proc.stdout.read()
    host_node = host_node.decode('ascii')
    host_node = host_node.strip()

    return host_node


def _get_pid():
    """ get pid
    """
    return os.getpid()

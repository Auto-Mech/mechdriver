""" Handle process and host info
"""

import os
import subprocess


HOST_MSG = """
=========================================================
HOST: {}
PID: {}
========================================================="""


def print_host_name():
    """ print the host the calculation is running on
    """
    host_node = get_host_node()
    pid = get_pid()
    print(HOST_MSG.format(host_node, pid))


def get_host_node():
    """ get the nodes
    """
    proc = subprocess.Popen(['hostname'], stdout=subprocess.PIPE)
    host_node = proc.stdout.read()
    host_node = host_node.decode('ascii')
    host_node = host_node.strip()

    return host_node


def get_pid():
    """ get pid
    """
    return os.getpid()

""" Library used to obtain various information about the
    shell process and run node that is associated with the
    MechDriver calculations the user launched.
"""

import os
import subprocess


def get_host_node():
    """ Calls the BASH `hostname` command to obtain the name of the
        node server that MechDriver is running on.

        :rtype: str
    """
    proc = subprocess.Popen(['hostname'], stdout=subprocess.PIPE)
    host_node = proc.stdout.read()
    host_node = host_node.decode('ascii')
    host_node = host_node.strip()

    return host_node


def get_pid():
    """ Gets the shell process ID for the MechDriver process running.

        :rtype: int
    """
    return os.getpid()

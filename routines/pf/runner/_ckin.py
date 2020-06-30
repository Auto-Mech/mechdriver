"""
  Handle ckin paths
"""

def ckin_path():
    """ Set path and make directories
    """
    starting_path = pfrunner.get_starting_path()
    path = pfrunner.prepare_path(starting_path, 'ckin')
    if not os.path.exists(path):
        os.makedirs(path)

    return path

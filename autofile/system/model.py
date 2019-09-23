""" defines the filesystem model
"""
import os
import glob
import types
import shutil
import autofile.file


class DataFile():
    """ file manager for a given datatype """

    def __init__(self, name, writer_=(lambda _: _), reader_=(lambda _: _)):
        """
        :param name: the file name
        :type name: str
        :param writer_: writes data to a string
        :type writer_: callable[object->str]
        :param reader_: reads data from a string
        :type reader_: callable[str->object]
        """
        self.name = name
        self.writer_ = writer_
        self.reader_ = reader_

    def path(self, dir_pth):
        """ file path
        """
        return os.path.join(dir_pth, self.name)

    def exists(self, dir_pth):
        """ does this file exist?
        """
        pth = self.path(dir_pth)
        return os.path.isfile(pth)

    def write(self, val, dir_pth):
        """ write data to this file
        """
        assert os.path.exists(dir_pth)
        pth = self.path(dir_pth)
        val_str = self.writer_(val)
        autofile.file.write_file(pth, val_str)

    def read(self, dir_pth):
        """ read data from this file
        """
        assert self.exists(dir_pth)
        pth = self.path(dir_pth)
        val_str = autofile.file.read_file(pth)
        val = self.reader_(val_str)
        return val


class DataSeries():
    """ directory manager mapping locator values to a directory series
    """

    def __init__(self, prefix, map_, nlocs, depth, loc_dfile=None,
                 root_ds=None, removable=False):
        """
        :param map_: maps `nlocs` locators to a segment path consisting of
            `depth` directories
        :param info_map_: maps `nlocs` locators to an information object, to
            be written in the data directory
        """
        assert os.path.isdir(prefix)
        self.prefix = os.path.abspath(prefix)
        self.map_ = map_
        self.nlocs = nlocs
        self.depth = depth
        self.loc_dfile = loc_dfile
        self.root = root_ds
        self.removable = removable
        self.file = types.SimpleNamespace()

    def add_data_files(self, dfile_dct):
        """ add DataFiles to the DataSeries
        """
        dfile_dct = {} if dfile_dct is None else dfile_dct

        for name, dfile in dfile_dct.items():
            assert isinstance(name, str)
            assert isinstance(dfile, DataFile)
            dsfile = _DataSeriesFile(ds=self, dfile=dfile)
            setattr(self.file, name, dsfile)

    def path(self, locs=()):
        """ absolute directory path
        """
        if self.root is None:
            prefix = self.prefix
        else:
            root_locs = self._root_locators(locs)
            locs = self._self_locators(locs)
            prefix = self.root.path(root_locs)
        assert len(locs) == self.nlocs

        pth = self.map_(locs)
        assert _path_is_relative(pth)
        assert _path_has_depth(pth, self.depth)
        return os.path.join(prefix, pth)

    def exists(self, locs=()):
        """ does this directory exist?
        """
        pth = self.path(locs)
        return os.path.isdir(pth)

    def remove(self, locs=()):
        """ does this directory exist?
        """
        if self.removable:
            pth = self.path(locs)
            if self.exists(locs):
                shutil.rmtree(pth)
        else:
            raise ValueError("This data series is not removable")

    def create(self, locs=()):
        """ create a directory at this prefix
        """
        # recursively create starting from the first root directory
        if self.root is not None:
            root_locs = self._root_locators(locs)
            self.root.create(root_locs)

        # create this directory in the chain, if it doesn't already exist
        if not self.exists(locs):
            pth = self.path(locs)
            os.makedirs(pth)

            if self.loc_dfile is not None:
                locs = self._self_locators(locs)
                self.loc_dfile.write(locs, pth)

    def existing(self, root_locs=(), relative=False):
        """ return the list of locators for existing paths
        """
        if self.loc_dfile is None:
            raise ValueError("This function does not work "
                             "without a locator DataFile")

        pths = self.existing_paths(root_locs)
        locs_lst = tuple(self.loc_dfile.read(pth) for pth in pths)
        if not relative:
            locs_lst = tuple(map(list(root_locs).__add__, locs_lst))

        return locs_lst

    def existing_paths(self, root_locs=()):
        """ existing paths at this prefix/root directory
        """
        if self.root is None:
            prefix = self.prefix
        else:
            prefix = self.root.path(root_locs)

        pth_pattern = os.path.join(prefix, *('*' * self.depth))
        pths = filter(os.path.isdir, glob.glob(pth_pattern))
        pths = tuple(sorted(os.path.join(prefix, pth) for pth in pths))
        return pths

    # helpers
    def _self_locators(self, locs):
        """ locators for this DataSeriesDir
        """
        nlocs = len(locs)
        assert nlocs >= self.nlocs
        root_nlocs = nlocs - self.nlocs
        return locs[root_nlocs:]

    def _root_locators(self, locs):
        """ locators for the root DataSeriesDir, if there is one
        """
        nlocs = len(locs)
        assert nlocs >= self.nlocs
        root_nlocs = nlocs - self.nlocs
        return locs[:root_nlocs]


class FileSystem(types.SimpleNamespace):
    """ a collection of DataSeries
    """

    def __init__(self, dseries_dct):
        self.update(dseries_dct)

    def __iter__(self):
        for key, val in vars(self).items():
            yield key, val

    def update(self, dseries_dct):
        """ update the filesystem dataseries
        """
        for name, obj in dict(dseries_dct).items():
            assert isinstance(name, str)
            assert isinstance(obj, DataSeries)
            setattr(self, name, obj)


# helpers:
class _DataSeriesFile():
    """ file manager mapping locator values to files in a directory series
    """

    def __init__(self, ds, dfile):
        self.dir = ds
        self.file = dfile

    def path(self, locs=()):
        """ absolute file path
        """
        return self.file.path(self.dir.path(locs))

    def exists(self, locs=()):
        """ does this file exist?
        """
        return self.file.exists(self.dir.path(locs))

    def write(self, val, locs=()):
        """ write data to this file
        """
        self.file.write(val, self.dir.path(locs))

    def read(self, locs=()):
        """ read data from this file
        """
        return self.file.read(self.dir.path(locs))


def _path_is_relative(pth):
    """ is this a relative path?
    """
    return os.path.relpath(pth) == pth


def _path_has_depth(pth, depth):
    """ does this path have the given depth?
    """
    return len(_os_path_split_all(pth)) == depth


def _os_path_split_all(pth):
    """ grabbed this from the internet """
    allparts = []
    while 1:
        parts = os.path.split(pth)
        if parts[0] == pth:    # sentinel for absolute paths
            allparts.insert(0, parts[0])
            break
        elif parts[1] == pth:  # sentinel for relative paths
            allparts.insert(0, parts[1])
            break
        else:
            pth = parts[0]
            allparts.insert(0, parts[1])
    return allparts

""" defines the filesystem model
"""
import os
import glob
import types
import shutil
import itertools
import autofile.io_


class DataFile():
    """ file manager for a given datatype

        :param name: the file name
        :type name: str
        :param writer_: writes data to a string
        :type writer_: callable[object->str]
        :param reader_: reads data from a string
        :type reader_: callable[str->object]
        :param removable: Is this file removable?
        :type removable: bool

    """
    def __init__(self, name, writer_=(lambda _: _), reader_=(lambda _: _)):
        self.name = name
        self.writer_ = writer_
        self.reader_ = reader_
        self.removable = False

    def path(self, dir_pth):
        """ file path

        :param dir_pth: directory path
        :type dir_pth: str
        :returns: datafile path
        :return type: str
        """
        return os.path.join(dir_pth, self.name)

    def exists(self, dir_pth):
        """ does this file exist?

        :param dir_pth: directory path
        :type dir_pth: str
        :returns: existance of datafile
        :return type: bool
        """
        pth = self.path(dir_pth)
        return os.path.isfile(pth)

    def write(self, val, dir_pth):
        """ write data to this file

        :param val: value to be written
        :type val: int/float/str/tuple
        :param dir_pth: directory path
        :type dir_pth: str
        """
        assert os.path.exists(dir_pth), (
            f'No path exists: {dir_pth}'
        )
        pth = self.path(dir_pth)
        val_str = self.writer_(val)
        autofile.io_.write_file(pth, val_str)

    def read(self, dir_pth):
        """ read data from this file

        :param dir_pth: directory path
        :type dir_pth: str
        :returns: datafile contents
        :return type: int/float/str/tuple
        """
        assert self.exists(dir_pth), (
            f'Either requested file {self}',
            f'or requested path does not exist {dir_pth}'
        )

        pth = self.path(dir_pth)
        val_str = autofile.io_.read_file(pth)
        val = self.reader_(val_str)
        return val

    def remove(self, dir_pth):
        """ remove this file

        (only possible if `removable` attribute is set to `True`)

        :param dir_pth: directory path
        :type dir_pth: str
        """
        if self.removable:
            pth = self.path(dir_pth)
            os.remove(pth)
        else:
            raise ValueError("This data series is not removable")

    def __repr__(self):
        """ represent this object as a string

        """
        return f"DataFile('{self.name}')"


class DataSeries():
    """ directory manager mapping locator values to a directory series


        :param map_: maps `nlocs` locators to a segment path consisting of
            `depth` directories
        :param info_map_: maps `nlocs` locators to an information object, to
            be written in the data directory
    """

    def __init__(self, prefix, map_, nlocs, depth, loc_dfile=None,
                 root_ds=None, removable=False):
        self.prefix = os.path.abspath(prefix)
        self.map_ = map_
        self.nlocs = nlocs
        self.depth = depth
        self.loc_dfile = loc_dfile
        self.root = root_ds
        self.removable = removable
        self.file = types.SimpleNamespace()
        self.json_file = 'db.json'
        self.json = types.SimpleNamespace()

    def add_data_files(self, dfile_dct):
        """ add DataFiles to the DataSeries

        """
        dfile_dct = {} if dfile_dct is None else dfile_dct

        for name, dfile in dfile_dct.items():
            assert isinstance(name, str)
            assert isinstance(dfile, DataFile)
            dsfile = DataSeriesFile(dseries=self, dfile=dfile)
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
        """ remove this directory


        (only possible if `removable` attribute is set to `True`)
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
            os.makedirs(pth, exist_ok=True)

            if self.loc_dfile is not None:
                locs = self._self_locators(locs)
                self.loc_dfile.write(locs, pth)
        if self.loc_dfile is not None:
            try:
                pth = self.path(locs)
                locs = self._self_locators(locs)
                self.loc_dfile.write(locs, pth)
            except AssertionError:
                pass

    def existing(self, root_locs=(), relative=False, ignore_bad_formats=True):
        """ return the list of locators for existing paths

        """
        if self.nlocs == 0:
            # If there are no locators, this DataSeries only produces one
            # directory and we can just check if it exists or not
            if self.exists():
                locs_lst = ([],)
            else:
                locs_lst = ()
        else:
            # If there are one or more locators, we need to read in from the
            # locator data files
            assert self.nlocs > 0

            if self.loc_dfile is None:
                raise ValueError("This function does not work "
                                 "without a locator DataFile")

            root_nlocs = self.root_locator_count()

            # Recursion for when we have a root DataSeries
            if len(root_locs) < root_nlocs:
                locs_lst = tuple(itertools.chain(*(
                    self.existing(root_locs_)
                    for root_locs_ in self.root.existing(root_locs))))
            else:
                assert root_nlocs == len(root_locs), (
                    f'{root_nlocs} != {len(root_locs)}'
                )
                pths = self._existing_paths(root_locs)
                if ignore_bad_formats:
                    locs_lst = []
                    for pth in pths:
                        if self.loc_dfile.exists(pth):
                            try:
                                pth_loc = self.loc_dfile.read(pth)
                                locs_lst.append(pth_loc)
                            except ValueError as exception:
                                print(
                                    'currently allowing ' +
                                    f'exception {exception}' +
                                    ' in existing to avoid crashes from' +
                                    '  CONF/cid in RUN')
                            except KeyError as exception:
                                print(
                                    'currently allowing ' +
                                    f'exception {exception}' +
                                    ' in existing to avoid crashes from' +
                                    '  CONF/cid in RUN')
                else:
                    locs_lst = tuple(self.loc_dfile.read(pth) for pth in pths
                                     if self.loc_dfile.exists(pth))

                if not relative:
                    locs_lst = tuple(map(list(root_locs).__add__, locs_lst))

        return locs_lst

    def _existing_paths(self, root_locs=()):
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

    def json_path(self, json_layer=None):
        """ json file path

        """
        if self.root is None:
            prefix = self.prefix
        else:
            prefix = self.root.path()
        if json_layer:
            prefix = _remove_layer_from_path(prefix, json_layer)
        return os.path.join(prefix, self.json_file)

    def json_exists(self, json_layer=None):
        """ does this file exist?

        """
        pth = self.json_path(json_layer=json_layer)
        return os.path.isfile(pth)

    def json_existing(self, locs=(), json_layer=None):
        """ returns a list of locations (aka keys) in the json file

        """
        locs = []
        if self.json_exists(json_layer=json_layer):
            json_data = autofile.json_.read_json(
                self.json_path(json_layer=json_layer))
            keys = json_data.keys()
            dct = json_data
            for nested_key in locs:
                if nested_key in keys:
                    dct = dct[nested_key]
                    keys = dct.keys()
                else:
                    dct = {}
                    keys = []
            locs = []
            for key in keys:
                if isinstance(dct[key], dict):
                    locs.append(key)
            locs = tuple([key] for key in keys)

        return locs

    def map(self, locs):
        """ returns a list of mapped locations

        """
        if locs:
            locs = [self.map_(locs)]
        return locs

    def add_json_entries(self, entry_dct):
        """ add DataFiles to the DataSeries

        """
        entry_dct = {} if entry_dct is None else entry_dct

        for name, obj in entry_dct.items():
            assert isinstance(name, str)
            assert isinstance(obj, JSONObject)
            jsentry = JSONEntry(jseries=self, jobject=obj)
            setattr(self.json, name, jsentry)

    def json_create(self, json_layer=None):
        """ json file creation

        """
        if not self.json_exists():
            autofile.json_.write_json(
                {}, self.json_path(json_layer=json_layer))

    def root_locator_count(self):
        """ count the number of root locator values recursively

        """
        if self.root is None:
            root_nlocs = 0
        else:
            root_nlocs = self.root.nlocs + self.root.root_locator_count()
        return root_nlocs

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
        assert nlocs >= self.nlocs, (
            f'{nlocs} != {self.nlocs}'
        )
        root_nlocs = nlocs - self.nlocs
        return locs[:root_nlocs]

    def __repr__(self):
        """ represent this object as a string

        """
        return f"DataSeries('{self.prefix}', {self.map_.__name__})"


class DataSeriesFile():
    """ file manager mapping locator values to files in a directory series

    """

    def __init__(self, dseries, dfile):
        self.dir = dseries
        self.file = dfile
        self.removable = False

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

    def remove(self, locs=()):
        """ remove this file


        (only possible if `removable` attribute is set to `True`)
        """
        if self.removable:
            pth = self.path(locs)
            os.remove(pth)
        else:
            raise ValueError("This data series is not removable")

    def __repr__(self):
        """ represent this object as a string

        """
        return (
            "DataSeriesFile("
            f"{self.dir.prefix}, {self.dir.map_.__name__}, {self.file.name}"
            ")"
        )


class JSONObject():
    """ json manager """

    def __init__(
            self, name, json_prefix=(None, None), removable=False,
            writer_=(lambda _: _), reader_=(lambda _: _)):
        """

        :param name: the file name
        :type name: str
        :param json_prefix: static, top level keys
        :type json_prefix: tuple
        :param writer_: writes data to a string
        :type writer_: callable[object->str]
        :param reader_: reads data from a string
        :type reader_: callable[str->object]
        :param removable: Is this file removable?
        :type removable: bool
        """
        self.name = name
        self.json_prefix, self.json_layer = json_prefix
        self.removable = removable
        self.writer_ = writer_
        self.reader_ = reader_

    def add_layer(self, key):
        """ add a key to the json loc list

        """
        layered_key = key
        if self.json_layer:
            layered_key = self.json_prefix.copy()
            layered_key.append(self.json_layer)
            layered_key.extend(key)
        return layered_key

    def exists(self, key, path):
        """ check existance of a json

        """
        exists = True
        key = self.add_layer(key)
        json_data = read_json(path)
        keys = json_data.keys()
        dct = json_data
        for nested_key in key:
            if nested_key in keys:
                dct = dct[nested_key]
                keys = dct.keys()
            else:
                exists = False
        if exists:
            if self.name not in dct:
                exists = False
        return exists

    def read(self, key, path):
        """ read a key out of a json file

        """
        json_data = read_json(path)
        return self._read(key, json_data)

    def _read(self, key, json_data):
        """ read a key out of a json file

        """
        exists = True
        key = self.add_layer(key)
        keys = json_data.keys()
        dct = json_data
        for nested_key in key:
            if nested_key in keys:
                dct = dct[nested_key]
                keys = dct.keys()
            else:
                exists = False
        if exists:
            if self.name in dct:
                return self.reader_(dct[self.name])
        return None

    def read_all(self, keys, path):
        """ read a key out of a json file for

            many keys
        """
        json_data = read_json(path)
        ret = []
        for key in keys:
            ret.append(self._read(key, json_data))
        return ret

    def write(self, val, key, path):
        """ write a value for a key in a json

        """
        key = self.add_layer(key)
        current_json = read_json(path)
        keys = current_json.keys()
        dct = current_json
        for nested_key in key:
            if nested_key not in keys:
                dct[nested_key] = {}
            dct = dct[nested_key]
            keys = dct.keys()
        val = self.writer_(val)
        dct[self.name] = val
        write_json(current_json, path)

    def write_all(self, vals, all_keys, path):
        """ write values for multiple keys in a json

        """
        current_json = read_json(path)
        for key, val in zip(all_keys, vals):
            key = self.add_layer(key)
            keys = current_json.keys()
            dct = current_json
            for nested_key in key:
                if nested_key not in keys:
                    dct[nested_key] = {}
                dct = dct[nested_key]
                keys = dct.keys()
            val = self.writer_(val)
            dct[self.name] = val
        write_json(current_json, path)


class JSONEntry():
    """ json manager for a given datatype

    """
    def __init__(self, jseries, jobject):
        """ json manager mapping locator values to objects in a json file

        """
        self.jseries = jseries
        self.json = jobject
        self.removable = False

    def exists(self, key=('database_entry'), mapping=True):
        """ does this entry exist?

        """
        ret = False
        if mapping:
            ret = self.json.exists(
                self.jseries.map(key), self.jseries.json_path(
                    json_layer=self.json.json_layer))
        else:
            ret = self.json.exists(key, self.jseries.json_path(
                json_layer=self.json.json_layer))
        return ret

    def existing(self, key=()):
        """ returns the keys nested under this key

        """
        return self.jseries.json_existing(
            locs=self.json.add_layer(key),
            json_layer=self.json.json_layer)

    def write(self, val, key=('database_entry'), mapping=True):
        """ write data to this file

        """
        if not self.jseries.json_exists(json_layer=self.json.json_layer):
            self.jseries.json_create(json_layer=self.json.json_layer)
        if mapping:
            key = self.jseries.map(key)
        self.json.write(val, key, self.jseries.json_path(
            json_layer=self.json.json_layer))

    def write_all(self, vals, keys=(('database_entry')), mapping=True):
        """ write data to this file

        """
        if not self.jseries.json_exists(json_layer=self.json.json_layer):
            self.jseries.json_create(json_layer=self.json.json_layer)
        if mapping:
            new_keys = []
            for key in keys:
                new_keys.append(self.jseries.map(key))
        else:
            new_keys = keys
        self.json.write_all(vals, new_keys, self.jseries.json_path(
            json_layer=self.json.json_layer))

    def read(self, key=('database_entry'), mapping=True):
        """ read data from this file

        """
        ret = None
        if self.exists(key, mapping=mapping):
            if mapping:
                ret = self.json.read(
                    self.jseries.map(key),
                    self.jseries.json_path(json_layer=self.json.json_layer))
            else:
                ret = self.json.read(key, self.jseries.json_path(
                    json_layer=self.json.json_layer))
        return ret

    def read_all(self, keys=(('database_entry')), mapping=True):
        """ read data from this file

        """
        ret = []
        new_keys = []
        for key in keys:
            if self.exists(key, mapping=mapping):
                if mapping:
                    new_keys.append(self.jseries.map(key))
                else:
                    new_keys.append(key)
        if new_keys:
            ret = self.json.read_all(new_keys, self.jseries.json_path(
                json_layer=self.json.json_layer))
        return ret


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
        if parts[1] == pth:  # sentinel for relative paths
            allparts.insert(0, parts[1])
            break
        pth = parts[0]
        allparts.insert(0, parts[1])
    return allparts


def _remove_layer_from_path(path, json_layer):
    head, tail = os.path.split(path)
    if tail == json_layer:
        path = head
    return path


def read_json(path):
    """ read a json

    """

    return autofile.json_.read_json(path)


def write_json(json_data, path):
    """ write a json

    """
    autofile.json_.write_json(json_data, path)


def items(path):
    """ what are the parent keys in this json file

    """
    keys = []
    entries = []
    if os.path.isfile(path):
        keys, entries = zip(*read_json(path).items())
    return keys, entries


def _keys(path):
    """ return the keys of a json

    """
    keys_, _ = items(path)
    return keys_


def _entries(path):
    """ return the entries of a json

    """
    _, entries_ = items(path)
    return entries_

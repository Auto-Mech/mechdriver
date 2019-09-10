from autoparse.find    import first_capture
from autoparse.find    import all_captures
from autoparse.pattern import capturing
from autoparse.pattern import zero_or_more
from autoparse.pattern import one_or_more
from autoparse.pattern import one_of_these
from autoparse.pattern import maybe
from autoparse.pattern import UNSIGNED_INTEGER
from autoparse.pattern import STRING_START
from autoparse.pattern import VARIABLE_NAME
from autoparse.pattern import NONNEWLINE
from autoparse.pattern import NONSPACE
from autoparse.pattern import SPACE
from autoparse.pattern import WILDCARD
from autoparse.pattern import INTEGER
from autoparse.pattern import FLOAT
import autoparse.find as apf
import automol.geom
import automol.zmatrix
import automol.convert.zmatrix
import autoread.zmatrix
import autofile.file

from esdriver.load import _species
from esdriver.load import _theory
from esdriver.load import _run

import numbers
from collections import OrderedDict 
import logging
log   = logging.getLogger(__name__)



def load_logger(outfile=''):
    loglevel = logging.DEBUG
    logging.addLevelName(logging.ERROR, 'ERROR: ')
    logging.addLevelName(logging.WARNING, 'WARNING: ')
    logging.addLevelName(logging.DEBUG, 'Status: ')
    logging.addLevelName(logging.INFO, '')
    if outfile:
        logging.basicConfig(format='%(levelname)s%(message)s', level=loglevel, filename=outfile, filemode='w')
    else:
        logging.basicConfig(format='%(levelname)s%(message)s', level=loglevel)
    return
 
def load_params():
    tsks = _load_tsks() 
    mods = _load_mods()
    lvls = _load_lvls()
    spcs = _load_spcs()
    return tsks, mods, lvls, spcs


def _load_spcs(fname='species.dat'):
    return _spcdic(_load_file(fname))

def _load_lvls(fname='theory.dat'):
    return _lvldic(_load_file(fname))

def _load_mods(fname='run.dat'):
    return _moddic(_load_file(fname))

def _load_tsks(fname='run.dat'):
    return _tskdic(_load_file(fname))

def _spcdic(s):
    spcdic   = OrderedDict()
    spckeys = _species.get_defined_species(s)
    for spec in spckeys:
        spcdic[spec] = {}
        keys = _species.get_defined_keywords(s, spec)    
        prevkey = '' 
        for key in keys:
            if prevkey == 'geom':
                call  = _species.get_attr_call(key, True)
            else:
                call  = _species.get_attr_call(key)
            if call != 'error_call':
                prevkey = key
            elif prevkey == 'geom': 
                continue
            call  = getattr(_species, call) 
            result = call(s, spec)
            spcdic[spec][key] = result
        if not 'inchi' in spcdic[spec]:
            if 'geom' in spcdic[spec]:
                geo = spcdic[spec]['geom']
                if not '=' in geo:   #geo was given in cartesian
                    syms, xyzs = autoread.geom.read(geo)
                    geo = automol.geom.from_data(syms, xyzs, angstrom=True)
                    spcdic[spec]['geoobj'] = geo
                    spcdic[spec][ 'inchi'] = automol.geom.inchi(geo)
                else:                #geo was given in zmat
                    zmat = spcdic[spec]['geom']
                    #syms, key_mat, name_mat = 
                    syms, key_mat, name_mat, val_dct = autoread.zmatrix.read(zmat, sym_ptt=autoread.par.Pattern.ATOM_SYMBOL + maybe(UNSIGNED_INTEGER),  key_ptt=one_of_these([UNSIGNED_INTEGER, VARIABLE_NAME]))
                    key_dct = dict(map(reversed, enumerate(syms)))
                    key_dct[None] = 0
                    key_mat = [[key_dct[val]+1 if not isinstance(val, numbers.Real) else val
                                for val in row] for row in key_mat]
                    sym_ptt = STRING_START + capturing(autoread.par.Pattern.ATOM_SYMBOL)
                    syms = [apf.first_capture(sym_ptt, sym) for sym in syms]
                    zma = automol.zmatrix.from_data(
                         syms, key_mat, name_mat, val_dct,
                          one_indexed=True, angstrom=True, degree=True)
                    spcdic[spec]['geoobj'] = automol.convert.zmatrix.geometry(zma)
                    spcdic[spec][ 'inchi'] = automol.geom.inchi(automol.convert.zmatrix.geometry(zma))
            elif 'smiles' in spcdic[spec]:
                spcdic[spec][ 'inchi'] = automol.smiles.inchi(spcdic[spec]['smiles'])
                spcdic[spec]['geeobj'] = automol.inchi.geometry(spcdic[spec][ 'inchi'])
        if not 'geoobj' in spcdic[spec]:
            if 'inchi' in spcdic[spec]:
                spcdic[spec]['geoobj'] = automol.inchi.geometry(spcdic[spec]['inchi'])
            elif 'smiles' in spcdic[spec]:
                spcdic[spec]['geoobj'] = automol.inchi.geometry(automol.smiles.inchi(spcdic[spec]['smiles']))
            else:
                log.error('No geom, inchi, or smiles provided for species {}'.format(spec))
    return spcdic 

def _lvldic(s):
    lvldic = OrderedDict()
    lvlkeys = _theory.get_defined_lvls(s)
    for lvl in lvlkeys:
        lvldic[lvl] = {}
        lvldic[lvl]['orb_res'] = 'RU'
        keys = _theory.get_defined_keywords(s, lvl)     
        for key in keys:
            call_ = _theory.get_attr_call(key)
            call  = getattr(_theory, call_) 
            if 'key' in call_:
                result = call(s, lvl, key)
            else:
                result = call(s, lvl)
            lvldic[lvl][key] = result
    return lvldic

def _moddic(s):
    moddic = {}
    tskkeys   = _run.get_defined_tasks(s)
    for tsk in tskkeys:
        moddic[tsk] = {}
        keys = _run.get_defined_keywords(s, tsk)     
        availkeys = ['job', 'reactants', 'products', 'references', 'paths']
        for key in availkeys:
            moddic[tsk][key] = []
            if key == 'paths':
                moddic[tsk][key] = ['rundir', 'savedir']
        keys  = list(set(availkeys) & set(keys))
        for key in keys:
            call_ = _run.get_attr_call(key)
            call  = getattr(_run, call_) 
            if 'key' in call_:
                result = call(s, tsk, key)
            else:
                result = call(s, tsk)
            moddic[tsk][key] = result
    return moddic

def _tskdic(s):
    tskdic = OrderedDict()
    tskkeys   = _run.get_defined_tasks(s)
    for tsk in tskkeys:
        tskdic[tsk] = {}
        keys = _run.get_defined_keywords(s, tsk)     
        unavailkeys = ['job', 'reactants', 'products', 'references', 'paths']
        keys  = sorted(set(keys) - set(unavailkeys), key = keys.index)
        for key in keys:
            call_ = _run.get_attr_call(key)
            call  = getattr(_run, call_) 
            if 'key' in call_:
                result = call(s, tsk, key)
            else:
                result = call(s, tsk)
            tskdic[tsk][key] = result
    return tskdic

def _load_file(fname):
    return autofile.file.read_file(fname)   


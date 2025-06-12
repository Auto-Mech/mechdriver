""" test the autofile.info module
"""
import autofile.info


def test_():
    """ tests
    """
    inf_obj = autofile.info.Info(
        a=['b', 'c', 'd', 'e'], x=autofile.info.Info(y=1, z=2))
    assert dict(inf_obj) == {'a': ['b', 'c', 'd', 'e'], 'x': {'y': 1, 'z': 2}}
    assert autofile.info.object_(dict(inf_obj)) == inf_obj
    print(autofile.info.object_(
        {'a': ['b', 'c', 'd', 'e'], 'x': {'y': 1, 'z': 2}}))
    print(dict(autofile.info.object_(
        {'a': ['b', 'c', 'd', 'e'], 'x': {'y': 1, 'z': 2}})))

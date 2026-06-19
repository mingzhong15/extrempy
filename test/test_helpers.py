"""Shared mock setup for tests that import through extrempy/__init__.py.

Usage: in any test file, do:
    from test_helpers import setup_mocks
    setup_mocks()
    import sys; sys.path.insert(0, ...)
"""
import sys, types
from unittest.mock import MagicMock


def setup_mocks():
    """Set up sys.modules mocks for all heavy dependencies."""
    # numpy: needs arbitrary attributes (pi, ones, ...)
    sys.modules['numpy'] = MagicMock()
    sys.modules['numpy'].__path__ = []

    # matplotlib
    sys.modules['matplotlib'] = types.ModuleType('matplotlib')
    sys.modules['matplotlib'].__path__ = []
    sys.modules['matplotlib'].pyplot = MagicMock()
    sys.modules['matplotlib.pyplot'] = sys.modules['matplotlib'].pyplot

    # dpdata
    sys.modules['dpdata'] = types.ModuleType('dpdata')
    sys.modules['dpdata'].__path__ = []
    sys.modules['dpdata'].System = MagicMock()
    sys.modules['dpdata'].LabeledSystem = MagicMock()

    # dscribe
    sys.modules['dscribe'] = types.ModuleType('dscribe')
    sys.modules['dscribe'].__path__ = []
    sys.modules['dscribe.descriptors'] = types.ModuleType('dscribe.descriptors')
    sys.modules['dscribe'].descriptors = sys.modules['dscribe.descriptors']

    # scipy
    sys.modules['scipy'] = types.ModuleType('scipy')
    sys.modules['scipy'].__path__ = []
    sys.modules['scipy'].optimize = types.ModuleType('scipy.optimize')
    sys.modules['scipy'].optimize.curve_fit = MagicMock()
    sys.modules['scipy.optimize'] = sys.modules['scipy'].optimize
    sys.modules['scipy'].interpolate = types.ModuleType('scipy.interpolate')
    sys.modules['scipy'].interpolate.interp1d = MagicMock()
    sys.modules['scipy.interpolate'] = sys.modules['scipy'].interpolate

    # ase: package with build submodule
    sys.modules['ase'] = types.ModuleType('ase')
    sys.modules['ase'].__path__ = []
    sys.modules['ase'].build = types.ModuleType('ase.build')
    sys.modules['ase'].build.bulk = MagicMock()
    sys.modules['ase'].build.make_supercell = MagicMock()
    sys.modules['ase'].io = types.ModuleType('ase.io')
    sys.modules['ase'].io.write = MagicMock()
    sys.modules['ase'].Atoms = MagicMock()
    sys.modules['ase.build'] = sys.modules['ase'].build
    sys.modules['ase.io'] = sys.modules['ase'].io

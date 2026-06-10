import numpy as np
import os

from ase.build import bulk, make_supercell
from ase.io import write
from ase.data import chemical_symbols, atomic_numbers
from ase import Atoms


# ================================================================
#  Element Phase Data
#  Key: 'phases' ordered from low to high temperature
#  'a'/'c' in Angstrom; 'T_min'/'T_max' in Kelvin
#  References: CRC Handbook, ASM, NIST, Pearson's Crystal Data
# ================================================================

ELEMENT_PHASE_DATA = {
    # ================================================================
    #  Period 1
    # ================================================================
    'H':  {'symbol': 'H',  'z': 1,  'phases': [{'structure': 'hcp', 'label': 'H',
            'a': 3.75, 'c': 6.12, 'T_min': 0, 'T_max': 14, 'note': 'solid H at low T'}], 'Tm': 14},
    'He': {'symbol': 'He', 'z': 2,  'phases': [{'structure': 'hcp', 'label': 'He',
            'a': 3.57, 'c': 5.83, 'T_min': 0, 'T_max': 1, 'note': 'under pressure'}], 'Tm': None},

    # ================================================================
    #  Period 2
    # ================================================================
    'Li': {'symbol': 'Li', 'z': 3,  'phases': [
        {'structure': 'bcc', 'label': 'α-Li', 'a': 3.509, 'T_min': 0, 'T_max': 80, 'note': 'martensitic bcc-hcp near 80K, often ignored'},
        {'structure': 'bcc', 'label': 'β-Li', 'a': 3.509, 'T_min': 80, 'T_max': 454},
    ], 'Tm': 454, 'rt_structure': 'bcc'},

    'Be': {'symbol': 'Be', 'z': 4,  'phases': [
        {'structure': 'hcp', 'label': 'α-Be', 'a': 2.286, 'c': 3.584, 'T_min': 0, 'T_max': 1527},
        {'structure': 'bcc', 'label': 'β-Be', 'a': 2.55,  'T_min': 1527, 'T_max': 1560},
    ], 'Tm': 1560, 'rt_structure': 'hcp'},

    'B':  {'symbol': 'B',  'z': 5,  'phases': [{'structure': 'rhombohedral', 'label': 'β-B',
            'a': 10.93, 'c': 23.79, 'T_min': 0, 'T_max': 2349}], 'Tm': 2349},

    'C':  {'symbol': 'C',  'z': 6,  'phases': [{'structure': 'graphite', 'label': 'C-graphite',
            'a': 2.464, 'c': 6.711, 'T_min': 0, 'T_max': None, 'note': 'sublimes at ~3915 K'}], 'Tm': None},

    'N':  {'symbol': 'N',  'z': 7,  'phases': [], 'Tm': 63},

    'O':  {'symbol': 'O',  'z': 8,  'phases': [], 'Tm': 55},

    'F':  {'symbol': 'F',  'z': 9,  'phases': [], 'Tm': 54},

    'Ne': {'symbol': 'Ne', 'z': 10, 'phases': [], 'Tm': 25},

    # ================================================================
    #  Period 3
    # ================================================================
    'Na': {'symbol': 'Na', 'z': 11, 'phases': [
        {'structure': 'bcc', 'label': 'α-Na', 'a': 4.291, 'T_min': 0, 'T_max': 371},
    ], 'Tm': 371, 'rt_structure': 'bcc'},

    'Mg': {'symbol': 'Mg', 'z': 12, 'phases': [
        {'structure': 'hcp', 'label': 'Mg', 'a': 3.209, 'c': 5.211, 'T_min': 0, 'T_max': 923},
    ], 'Tm': 923, 'rt_structure': 'hcp'},

    'Al': {'symbol': 'Al', 'z': 13, 'phases': [
        {'structure': 'fcc', 'label': 'Al', 'a': 4.050, 'T_min': 0, 'T_max': 933},
    ], 'Tm': 933, 'rt_structure': 'fcc'},

    'Si': {'symbol': 'Si', 'z': 14, 'phases': [
        {'structure': 'diamond', 'label': 'Si', 'a': 5.431, 'T_min': 0, 'T_max': 1687},
    ], 'Tm': 1687, 'rt_structure': 'diamond'},

    'P':  {'symbol': 'P',  'z': 15, 'phases': [], 'Tm': 317},

    'S':  {'symbol': 'S',  'z': 16, 'phases': [], 'Tm': 388},

    'Cl': {'symbol': 'Cl', 'z': 17, 'phases': [], 'Tm': 172},

    'Ar': {'symbol': 'Ar', 'z': 18, 'phases': [], 'Tm': 84},

    # ================================================================
    #  Period 4
    # ================================================================
    'K': {'symbol': 'K', 'z': 19, 'phases': [
        {'structure': 'bcc', 'label': 'K', 'a': 5.328, 'T_min': 0, 'T_max': 337},
    ], 'Tm': 337, 'rt_structure': 'bcc'},

    'Ca': {'symbol': 'Ca', 'z': 20, 'phases': [
        {'structure': 'fcc', 'label': 'α-Ca', 'a': 5.588, 'T_min': 0, 'T_max': 716},
        {'structure': 'bcc', 'label': 'β-Ca', 'a': 4.485, 'T_min': 716, 'T_max': 1115},
    ], 'Tm': 1115, 'rt_structure': 'fcc'},

    'Sc': {'symbol': 'Sc', 'z': 21, 'phases': [
        {'structure': 'hcp', 'label': 'α-Sc', 'a': 3.309, 'c': 5.268, 'T_min': 0, 'T_max': 1610},
        {'structure': 'bcc', 'label': 'β-Sc', 'a': 3.73,  'T_min': 1610, 'T_max': 1814},
    ], 'Tm': 1814, 'rt_structure': 'hcp'},

    'Ti': {'symbol': 'Ti', 'z': 22, 'phases': [
        {'structure': 'hcp', 'label': 'α-Ti', 'a': 2.951, 'c': 4.686, 'T_min': 0, 'T_max': 1155},
        {'structure': 'bcc', 'label': 'β-Ti', 'a': 3.307, 'T_min': 1155, 'T_max': 1941},
    ], 'Tm': 1941, 'rt_structure': 'hcp'},

    'V': {'symbol': 'V', 'z': 23, 'phases': [
        {'structure': 'bcc', 'label': 'V', 'a': 3.027, 'T_min': 0, 'T_max': 2183},
    ], 'Tm': 2183, 'rt_structure': 'bcc'},

    'Cr': {'symbol': 'Cr', 'z': 24, 'phases': [
        {'structure': 'bcc', 'label': 'Cr', 'a': 2.885, 'T_min': 0, 'T_max': 2180},
    ], 'Tm': 2180, 'rt_structure': 'bcc', 'magnetic': {'type': 'antiferromagnetic', 'TN': 311}},

    'Mn': {'symbol': 'Mn', 'z': 25, 'phases': [
        {'structure': 'cbcc', 'label': 'α-Mn', 'a': 8.915, 'T_min': 0, 'T_max': 1000, 'note': 'complex bcc, 58 atoms/cell'},
        {'structure': 'cpcc', 'label': 'β-Mn', 'a': 6.315, 'T_min': 1000, 'T_max': 1368, 'note': 'primitive cubic, 20 atoms/cell'},
        {'structure': 'fcc',  'label': 'γ-Mn', 'a': 3.86,  'T_min': 1368, 'T_max': 1406},
        {'structure': 'bcc',  'label': 'δ-Mn', 'a': 3.08,  'T_min': 1406, 'T_max': 1519},
    ], 'Tm': 1519, 'rt_structure': 'cbcc'},

    'Fe': {'symbol': 'Fe', 'z': 26, 'phases': [
        {'structure': 'bcc', 'label': 'α-Fe', 'a': 2.866, 'T_min': 0, 'T_max': 1185},
        {'structure': 'fcc', 'label': 'γ-Fe', 'a': 3.647, 'T_min': 1185, 'T_max': 1667},
        {'structure': 'bcc', 'label': 'δ-Fe', 'a': 2.93,  'T_min': 1667, 'T_max': 1811},
    ], 'Tm': 1811, 'rt_structure': 'bcc', 'magnetic': {'type': 'ferromagnetic', 'Tc': 1043}},

    'Co': {'symbol': 'Co', 'z': 27, 'phases': [
        {'structure': 'hcp', 'label': 'ε-Co', 'a': 2.507, 'c': 4.070, 'T_min': 0, 'T_max': 695},
        {'structure': 'fcc', 'label': 'α-Co', 'a': 3.544, 'T_min': 695, 'T_max': 1768},
    ], 'Tm': 1768, 'rt_structure': 'hcp', 'magnetic': {'type': 'ferromagnetic', 'Tc': 1388}},

    'Ni': {'symbol': 'Ni', 'z': 28, 'phases': [
        {'structure': 'fcc', 'label': 'Ni', 'a': 3.524, 'T_min': 0, 'T_max': 1728},
    ], 'Tm': 1728, 'rt_structure': 'fcc', 'magnetic': {'type': 'ferromagnetic', 'Tc': 627}},

    'Cu': {'symbol': 'Cu', 'z': 29, 'phases': [
        {'structure': 'fcc', 'label': 'Cu', 'a': 3.615, 'T_min': 0, 'T_max': 1358},
    ], 'Tm': 1358, 'rt_structure': 'fcc'},

    'Zn': {'symbol': 'Zn', 'z': 30, 'phases': [
        {'structure': 'hcp', 'label': 'Zn', 'a': 2.665, 'c': 4.947, 'T_min': 0, 'T_max': 693},
    ], 'Tm': 693, 'rt_structure': 'hcp'},

    'Ga': {'symbol': 'Ga', 'z': 31, 'phases': [{'structure': 'orthorhombic', 'label': 'α-Ga',
            'a': 4.526, 'b': 4.520, 'c': 7.659, 'T_min': 0, 'T_max': 303}], 'Tm': 303},

    'Ge': {'symbol': 'Ge', 'z': 32, 'phases': [
        {'structure': 'diamond', 'label': 'Ge', 'a': 5.658, 'T_min': 0, 'T_max': 1211},
    ], 'Tm': 1211, 'rt_structure': 'diamond'},

    'As': {'symbol': 'As', 'z': 33, 'phases': [{'structure': 'rhombohedral', 'label': 'As',
            'a': 4.132, 'c': 10.939, 'T_min': 0, 'T_max': 1090, 'note': 'sublimes'}], 'Tm': None},

    'Se': {'symbol': 'Se', 'z': 34, 'phases': [], 'Tm': 494},

    'Br': {'symbol': 'Br', 'z': 35, 'phases': [], 'Tm': 266},

    'Kr': {'symbol': 'Kr', 'z': 36, 'phases': [], 'Tm': 116},

    # ================================================================
    #  Period 5
    # ================================================================
    'Rb': {'symbol': 'Rb', 'z': 37, 'phases': [
        {'structure': 'bcc', 'label': 'Rb', 'a': 5.710, 'T_min': 0, 'T_max': 312},
    ], 'Tm': 312, 'rt_structure': 'bcc'},

    'Sr': {'symbol': 'Sr', 'z': 38, 'phases': [
        {'structure': 'fcc', 'label': 'α-Sr', 'a': 6.085, 'T_min': 0, 'T_max': 830},
        {'structure': 'bcc', 'label': 'β-Sr', 'a': 4.87,  'T_min': 830, 'T_max': 1050},
    ], 'Tm': 1050, 'rt_structure': 'fcc'},

    'Y': {'symbol': 'Y', 'z': 39, 'phases': [
        {'structure': 'hcp', 'label': 'α-Y', 'a': 3.647, 'c': 5.731, 'T_min': 0, 'T_max': 1755},
        {'structure': 'bcc', 'label': 'β-Y', 'a': 4.10,  'T_min': 1755, 'T_max': 1799},
    ], 'Tm': 1799, 'rt_structure': 'hcp'},

    'Zr': {'symbol': 'Zr', 'z': 40, 'phases': [
        {'structure': 'hcp', 'label': 'α-Zr', 'a': 3.232, 'c': 5.148, 'T_min': 0, 'T_max': 1136},
        {'structure': 'bcc', 'label': 'β-Zr', 'a': 3.609, 'T_min': 1136, 'T_max': 2128},
    ], 'Tm': 2128, 'rt_structure': 'hcp'},

    'Nb': {'symbol': 'Nb', 'z': 41, 'phases': [
        {'structure': 'bcc', 'label': 'Nb', 'a': 3.300, 'T_min': 0, 'T_max': 2750},
    ], 'Tm': 2750, 'rt_structure': 'bcc'},

    'Mo': {'symbol': 'Mo', 'z': 42, 'phases': [
        {'structure': 'bcc', 'label': 'Mo', 'a': 3.147, 'T_min': 0, 'T_max': 2896},
    ], 'Tm': 2896, 'rt_structure': 'bcc'},

    'Tc': {'symbol': 'Tc', 'z': 43, 'phases': [
        {'structure': 'hcp', 'label': 'Tc', 'a': 2.737, 'c': 4.391, 'T_min': 0, 'T_max': 2430},
    ], 'Tm': 2430, 'rt_structure': 'hcp'},

    'Ru': {'symbol': 'Ru', 'z': 44, 'phases': [
        {'structure': 'hcp', 'label': 'Ru', 'a': 2.706, 'c': 4.282, 'T_min': 0, 'T_max': 2607},
    ], 'Tm': 2607, 'rt_structure': 'hcp'},

    'Rh': {'symbol': 'Rh', 'z': 45, 'phases': [
        {'structure': 'fcc', 'label': 'Rh', 'a': 3.803, 'T_min': 0, 'T_max': 2237},
    ], 'Tm': 2237, 'rt_structure': 'fcc'},

    'Pd': {'symbol': 'Pd', 'z': 46, 'phases': [
        {'structure': 'fcc', 'label': 'Pd', 'a': 3.891, 'T_min': 0, 'T_max': 1828},
    ], 'Tm': 1828, 'rt_structure': 'fcc'},

    'Ag': {'symbol': 'Ag', 'z': 47, 'phases': [
        {'structure': 'fcc', 'label': 'Ag', 'a': 4.085, 'T_min': 0, 'T_max': 1235},
    ], 'Tm': 1235, 'rt_structure': 'fcc'},

    'Cd': {'symbol': 'Cd', 'z': 48, 'phases': [
        {'structure': 'hcp', 'label': 'Cd', 'a': 2.979, 'c': 5.619, 'T_min': 0, 'T_max': 594},
    ], 'Tm': 594, 'rt_structure': 'hcp'},

    'In': {'symbol': 'In', 'z': 49, 'phases': [
        {'structure': 'bct',  'label': 'In', 'a': 3.252, 'c': 4.946, 'T_min': 0, 'T_max': 430},
    ], 'Tm': 430, 'rt_structure': 'bct'},

    'Sn': {'symbol': 'Sn', 'z': 50, 'phases': [
        {'structure': 'diamond', 'label': 'α-Sn', 'a': 6.489, 'T_min': 0, 'T_max': 286, 'note': 'gray tin, semiconductor'},
        {'structure': 'bct',     'label': 'β-Sn', 'a': 5.832, 'c': 3.181, 'T_min': 286, 'T_max': 505, 'note': 'white tin, metallic'},
    ], 'Tm': 505, 'rt_structure': 'bct'},

    'Sb': {'symbol': 'Sb', 'z': 51, 'phases': [{'structure': 'rhombohedral', 'label': 'Sb',
            'a': 4.308, 'c': 11.274, 'T_min': 0, 'T_max': 904}], 'Tm': 904},

    'Te': {'symbol': 'Te', 'z': 52, 'phases': [], 'Tm': 723},

    'I':  {'symbol': 'I',  'z': 53, 'phases': [], 'Tm': 387},

    'Xe': {'symbol': 'Xe', 'z': 54, 'phases': [], 'Tm': 161},

    # ================================================================
    #  Period 6
    # ================================================================
    'Cs': {'symbol': 'Cs', 'z': 55, 'phases': [
        {'structure': 'bcc', 'label': 'Cs', 'a': 6.141, 'T_min': 0, 'T_max': 302},
    ], 'Tm': 302, 'rt_structure': 'bcc'},

    'Ba': {'symbol': 'Ba', 'z': 56, 'phases': [
        {'structure': 'bcc', 'label': 'Ba', 'a': 5.023, 'T_min': 0, 'T_max': 1000},
    ], 'Tm': 1000, 'rt_structure': 'bcc'},

    'La': {'symbol': 'La', 'z': 57, 'phases': [
        {'structure': 'dhcp', 'label': 'α-La', 'a': 3.774, 'c': 12.171, 'T_min': 0, 'T_max': 613},
        {'structure': 'fcc',  'label': 'β-La', 'a': 5.303, 'T_min': 613, 'T_max': 1138},
        {'structure': 'bcc',  'label': 'γ-La', 'a': 4.26,  'T_min': 1138, 'T_max': 1193},
    ], 'Tm': 1193, 'rt_structure': 'dhcp'},

    'Ce': {'symbol': 'Ce', 'z': 58, 'phases': [
        {'structure': 'dhcp', 'label': 'β-Ce', 'a': 3.681, 'c': 11.857, 'T_min': 0, 'T_max': 280},
        {'structure': 'fcc',  'label': 'γ-Ce', 'a': 5.161, 'T_min': 280, 'T_max': 999},
        {'structure': 'bcc',  'label': 'δ-Ce', 'a': 4.12,  'T_min': 999, 'T_max': 1068},
    ], 'Tm': 1068, 'rt_structure': 'dhcp', 'note': 'γ→β at ~280K; α-Ce (collapsed fcc) appears below ~116K'},

    'Pr': {'symbol': 'Pr', 'z': 59, 'phases': [
        {'structure': 'dhcp', 'label': 'α-Pr', 'a': 3.673, 'c': 11.835, 'T_min': 0, 'T_max': 1068},
        {'structure': 'bcc',  'label': 'β-Pr', 'a': 4.13,  'T_min': 1068, 'T_max': 1208},
    ], 'Tm': 1208, 'rt_structure': 'dhcp'},

    'Nd': {'symbol': 'Nd', 'z': 60, 'phases': [
        {'structure': 'dhcp', 'label': 'α-Nd', 'a': 3.658, 'c': 11.799, 'T_min': 0, 'T_max': 1136},
        {'structure': 'bcc',  'label': 'β-Nd', 'a': 4.13,  'T_min': 1136, 'T_max': 1297},
    ], 'Tm': 1297, 'rt_structure': 'dhcp'},

    'Pm': {'symbol': 'Pm', 'z': 61, 'phases': [
        {'structure': 'dhcp', 'label': 'α-Pm', 'a': 3.65,  'c': 11.65,  'T_min': 0, 'T_max': 1163},
        {'structure': 'bcc',  'label': 'β-Pm', 'a': 4.10,  'T_min': 1163, 'T_max': 1315},
    ], 'Tm': 1315, 'rt_structure': 'dhcp'},

    'Sm': {'symbol': 'Sm', 'z': 62, 'phases': [
        {'structure': 'rhombohedral', 'label': 'α-Sm', 'a': 3.629, 'c': 26.207, 'T_min': 0, 'T_max': 1007, 'note': '9R Sm-type'},
        {'structure': 'bcc',          'label': 'β-Sm', 'a': 4.07,  'T_min': 1007, 'T_max': 1345},
    ], 'Tm': 1345, 'rt_structure': 'rhombohedral'},

    'Eu': {'symbol': 'Eu', 'z': 63, 'phases': [
        {'structure': 'bcc', 'label': 'Eu', 'a': 4.581, 'T_min': 0, 'T_max': 1099},
    ], 'Tm': 1099, 'rt_structure': 'bcc'},

    'Gd': {'symbol': 'Gd', 'z': 64, 'phases': [
        {'structure': 'hcp', 'label': 'α-Gd', 'a': 3.634, 'c': 5.781, 'T_min': 0, 'T_max': 1508},
        {'structure': 'bcc', 'label': 'β-Gd', 'a': 4.06,  'T_min': 1508, 'T_max': 1585},
    ], 'Tm': 1585, 'rt_structure': 'hcp', 'magnetic': {'type': 'ferromagnetic', 'Tc': 293}},

    'Tb': {'symbol': 'Tb', 'z': 65, 'phases': [
        {'structure': 'hcp', 'label': 'α-Tb', 'a': 3.605, 'c': 5.696, 'T_min': 0, 'T_max': 1562},
        {'structure': 'bcc', 'label': 'β-Tb', 'a': 4.02,  'T_min': 1562, 'T_max': 1629},
    ], 'Tm': 1629, 'rt_structure': 'hcp'},

    'Dy': {'symbol': 'Dy', 'z': 66, 'phases': [
        {'structure': 'hcp', 'label': 'α-Dy', 'a': 3.590, 'c': 5.654, 'T_min': 0, 'T_max': 1654},
        {'structure': 'bcc', 'label': 'β-Dy', 'a': 4.03,  'T_min': 1654, 'T_max': 1680},
    ], 'Tm': 1680, 'rt_structure': 'hcp'},

    'Ho': {'symbol': 'Ho', 'z': 67, 'phases': [
        {'structure': 'hcp', 'label': 'α-Ho', 'a': 3.578, 'c': 5.618, 'T_min': 0, 'T_max': 1701},
        {'structure': 'bcc', 'label': 'β-Ho', 'a': 3.96,  'T_min': 1701, 'T_max': 1734},
    ], 'Tm': 1734, 'rt_structure': 'hcp'},

    'Er': {'symbol': 'Er', 'z': 68, 'phases': [
        {'structure': 'hcp', 'label': 'α-Er', 'a': 3.559, 'c': 5.587, 'T_min': 0, 'T_max': 1695},
        {'structure': 'bcc', 'label': 'β-Er', 'a': 3.94,  'T_min': 1695, 'T_max': 1802},
    ], 'Tm': 1802, 'rt_structure': 'hcp'},

    'Tm': {'symbol': 'Tm', 'z': 69, 'phases': [
        {'structure': 'hcp', 'label': 'α-Tm', 'a': 3.538, 'c': 5.554, 'T_min': 0, 'T_max': 1731},
        {'structure': 'bcc', 'label': 'β-Tm', 'a': 3.85,  'T_min': 1731, 'T_max': 1818},
    ], 'Tm': 1818, 'rt_structure': 'hcp'},

    'Yb': {'symbol': 'Yb', 'z': 70, 'phases': [
        {'structure': 'fcc', 'label': 'α-Yb', 'a': 5.486, 'T_min': 0, 'T_max': 1033},
        {'structure': 'bcc', 'label': 'β-Yb', 'a': 4.45,  'T_min': 1033, 'T_max': 1097},
    ], 'Tm': 1097, 'rt_structure': 'fcc'},

    'Lu': {'symbol': 'Lu', 'z': 71, 'phases': [
        {'structure': 'hcp', 'label': 'Lu', 'a': 3.503, 'c': 5.551, 'T_min': 0, 'T_max': 1925},
    ], 'Tm': 1925, 'rt_structure': 'hcp'},

    'Hf': {'symbol': 'Hf', 'z': 72, 'phases': [
        {'structure': 'hcp', 'label': 'α-Hf', 'a': 3.196, 'c': 5.051, 'T_min': 0, 'T_max': 2020},
        {'structure': 'bcc', 'label': 'β-Hf', 'a': 3.60,  'T_min': 2020, 'T_max': 2506},
    ], 'Tm': 2506, 'rt_structure': 'hcp'},

    'Ta': {'symbol': 'Ta', 'z': 73, 'phases': [
        {'structure': 'bcc', 'label': 'Ta', 'a': 3.301, 'T_min': 0, 'T_max': 3290},
    ], 'Tm': 3290, 'rt_structure': 'bcc'},

    'W': {'symbol': 'W', 'z': 74, 'phases': [
        {'structure': 'bcc', 'label': 'W', 'a': 3.165, 'T_min': 0, 'T_max': 3695},
    ], 'Tm': 3695, 'rt_structure': 'bcc'},

    'Re': {'symbol': 'Re', 'z': 75, 'phases': [
        {'structure': 'hcp', 'label': 'Re', 'a': 2.760, 'c': 4.458, 'T_min': 0, 'T_max': 3459},
    ], 'Tm': 3459, 'rt_structure': 'hcp'},

    'Os': {'symbol': 'Os', 'z': 76, 'phases': [
        {'structure': 'hcp', 'label': 'Os', 'a': 2.734, 'c': 4.320, 'T_min': 0, 'T_max': 3306},
    ], 'Tm': 3306, 'rt_structure': 'hcp'},

    'Ir': {'symbol': 'Ir', 'z': 77, 'phases': [
        {'structure': 'fcc', 'label': 'Ir', 'a': 3.839, 'T_min': 0, 'T_max': 2719},
    ], 'Tm': 2719, 'rt_structure': 'fcc'},

    'Pt': {'symbol': 'Pt', 'z': 78, 'phases': [
        {'structure': 'fcc', 'label': 'Pt', 'a': 3.924, 'T_min': 0, 'T_max': 2041},
    ], 'Tm': 2041, 'rt_structure': 'fcc'},

    'Au': {'symbol': 'Au', 'z': 79, 'phases': [
        {'structure': 'fcc', 'label': 'Au', 'a': 4.078, 'T_min': 0, 'T_max': 1337},
    ], 'Tm': 1337, 'rt_structure': 'fcc'},

    'Hg': {'symbol': 'Hg', 'z': 80, 'phases': [{'structure': 'rhombohedral', 'label': 'Hg',
            'a': 3.005, 'c': 6.704, 'T_min': 0, 'T_max': 234}], 'Tm': 234},

    'Tl': {'symbol': 'Tl', 'z': 81, 'phases': [
        {'structure': 'hcp', 'label': 'α-Tl', 'a': 3.457, 'c': 5.525, 'T_min': 0, 'T_max': 507},
        {'structure': 'bcc', 'label': 'β-Tl', 'a': 3.88,  'T_min': 507, 'T_max': 577},
    ], 'Tm': 577, 'rt_structure': 'hcp'},

    'Pb': {'symbol': 'Pb', 'z': 82, 'phases': [
        {'structure': 'fcc', 'label': 'Pb', 'a': 4.951, 'T_min': 0, 'T_max': 601},
    ], 'Tm': 601, 'rt_structure': 'fcc'},

    'Bi': {'symbol': 'Bi', 'z': 83, 'phases': [{'structure': 'rhombohedral', 'label': 'Bi',
            'a': 4.746, 'c': 11.862, 'T_min': 0, 'T_max': 545}], 'Tm': 545},

    'Po': {'symbol': 'Po', 'z': 84, 'phases': [
        {'structure': 'sc', 'label': 'α-Po', 'a': 3.359, 'T_min': 0, 'T_max': 309},
        {'structure': 'rhombohedral', 'label': 'β-Po', 'a': 3.366, 'T_min': 309, 'T_max': 527},
    ], 'Tm': 527, 'rt_structure': 'sc'},

    'At': {'symbol': 'At', 'z': 85, 'phases': [], 'Tm': 575},

    'Rn': {'symbol': 'Rn', 'z': 86, 'phases': [], 'Tm': 202},

    # ================================================================
    #  Period 7
    # ================================================================
    'Fr': {'symbol': 'Fr', 'z': 87, 'phases': [
        {'structure': 'bcc', 'label': 'Fr', 'a': 5.64, 'T_min': 0, 'T_max': 300},
    ], 'Tm': 300, 'rt_structure': 'bcc'},

    'Ra': {'symbol': 'Ra', 'z': 88, 'phases': [
        {'structure': 'bcc', 'label': 'Ra', 'a': 5.148, 'T_min': 0, 'T_max': 973},
    ], 'Tm': 973, 'rt_structure': 'bcc'},

    'Ac': {'symbol': 'Ac', 'z': 89, 'phases': [
        {'structure': 'fcc', 'label': 'Ac', 'a': 5.311, 'T_min': 0, 'T_max': 1323},
    ], 'Tm': 1323, 'rt_structure': 'fcc'},

    'Th': {'symbol': 'Th', 'z': 90, 'phases': [
        {'structure': 'fcc', 'label': 'α-Th', 'a': 5.084, 'T_min': 0, 'T_max': 1633},
        {'structure': 'bcc', 'label': 'β-Th', 'a': 4.11,  'T_min': 1633, 'T_max': 2023},
    ], 'Tm': 2023, 'rt_structure': 'fcc'},

    'Pa': {'symbol': 'Pa', 'z': 91, 'phases': [
        {'structure': 'tetragonal', 'label': 'α-Pa', 'a': 3.929, 'c': 3.238, 'T_min': 0, 'T_max': 1443},
        {'structure': 'bcc',       'label': 'β-Pa', 'a': 3.81,  'T_min': 1443, 'T_max': 1841},
    ], 'Tm': 1841, 'rt_structure': 'tetragonal'},

    'U': {'symbol': 'U', 'z': 92, 'phases': [
        {'structure': 'orthorhombic', 'label': 'α-U', 'a': 2.854, 'b': 5.870, 'c': 4.955, 'T_min': 0, 'T_max': 941},
        {'structure': 'tetragonal',   'label': 'β-U', 'a': 10.758, 'c': 5.656, 'T_min': 941, 'T_max': 1049, 'note': '30 atoms/cell'},
        {'structure': 'bcc',          'label': 'γ-U', 'a': 3.524, 'T_min': 1049, 'T_max': 1408},
    ], 'Tm': 1408, 'rt_structure': 'orthorhombic'},

    'Np': {'symbol': 'Np', 'z': 93, 'phases': [
        {'structure': 'orthorhombic', 'label': 'α-Np', 'a': 6.663, 'b': 4.723, 'c': 4.887, 'T_min': 0, 'T_max': 553},
        {'structure': 'tetragonal',   'label': 'β-Np', 'a': 4.887, 'c': 3.389, 'T_min': 553, 'T_max': 849},
        {'structure': 'bcc',          'label': 'γ-Np', 'a': 3.52,  'T_min': 849, 'T_max': 917},
    ], 'Tm': 917, 'rt_structure': 'orthorhombic'},

    'Pu': {'symbol': 'Pu', 'z': 94, 'phases': [
        {'structure': 'monoclinic',    'label': 'α-Pu',  'a': 6.183, 'b': 4.822, 'c': 10.963, 'T_min': 0, 'T_max': 395, 'note': '16 atoms/cell'},
        {'structure': 'monoclinic',    'label': 'β-Pu',  'a': 9.284, 'b': 10.463, 'c': 7.859, 'T_min': 395, 'T_max': 479, 'note': '34 atoms/cell'},
        {'structure': 'orthorhombic',  'label': 'γ-Pu',  'a': 3.159, 'b': 5.768, 'c': 10.162, 'T_min': 479, 'T_max': 487},
        {'structure': 'fcc',           'label': 'δ-Pu',  'a': 4.637, 'T_min': 487, 'T_max': 593},
        {'structure': 'bct',           'label': 'δ-Pu', 'a': 3.340, 'c': 4.447, 'T_min': 593, 'T_max': 736},
        {'structure': 'bcc',           'label': 'ε-Pu',  'a': 3.638, 'T_min': 736, 'T_max': 913},
    ], 'Tm': 913, 'rt_structure': 'monoclinic', 'note': '6 allotropes, most of any element'},

    'Am': {'symbol': 'Am', 'z': 95, 'phases': [
        {'structure': 'dhcp', 'label': 'α-Am', 'a': 3.468, 'c': 11.248, 'T_min': 0, 'T_max': 1042},
        {'structure': 'fcc',  'label': 'β-Am', 'a': 4.89,  'T_min': 1042, 'T_max': 1449},
    ], 'Tm': 1449, 'rt_structure': 'dhcp'},

    'Cm': {'symbol': 'Cm', 'z': 96, 'phases': [
        {'structure': 'dhcp', 'label': 'α-Cm', 'a': 3.496, 'c': 11.331, 'T_min': 0, 'T_max': 1550},
        {'structure': 'fcc',  'label': 'β-Cm', 'a': 4.382, 'T_min': 1550, 'T_max': 1613},
    ], 'Tm': 1613, 'rt_structure': 'dhcp'},
}


# ================================================================
#  Legacy / derived: ELEMENTS_BY_STRUCTURE (room-temperature phase only)
# ================================================================

def _build_elements_by_structure():
    """Build ELEMENTS_BY_STRUCTURE from ELEMENT_PHASE_DATA using rt_structure."""
    result = {}
    for sym, data in ELEMENT_PHASE_DATA.items():
        rt = data.get('rt_structure')
        if rt and data.get('phases') and len(data['phases']) > 0:
            if rt not in result:
                result[rt] = []
            result[rt].append(sym)
    return result

ELEMENTS_BY_STRUCTURE = _build_elements_by_structure()


# ================================================================
#  Convenience lookup
# ================================================================

def get_element_phases(symbol):
    """Get all phases for an element."""
    data = ELEMENT_PHASE_DATA.get(symbol, {})
    return data.get('phases', [])

def get_rt_phase(symbol):
    """Get the room-temperature phase info."""
    data = ELEMENT_PHASE_DATA.get(symbol, {})
    if data.get('phases'):
        return data['phases'][0]

def get_element_lattice_constants(symbol, structure_type=None):
    """Get (a, c) for element at its standard structure."""
    data = ELEMENT_PHASE_DATA.get(symbol, {})
    if not data.get('phases'):
        return None
    if structure_type:
        for phase in data['phases']:
            if phase['structure'] == structure_type:
                return phase.get('a'), phase.get('c')
    return data['phases'][0].get('a'), data['phases'][0].get('c')

def get_melting_point(symbol):
    """Get melting point in K."""
    data = ELEMENT_PHASE_DATA.get(symbol, {})
    return data.get('Tm')

def has_phase_transition(symbol):
    """Check if element has solid-solid phase transition."""
    data = ELEMENT_PHASE_DATA.get(symbol, {})
    return len(data.get('phases', [])) > 1

def get_phase_at_temperature(symbol, T):
    """Get the stable phase at a given temperature T (K)."""
    data = ELEMENT_PHASE_DATA.get(symbol, {})
    for phase in data.get('phases', []):
        if phase['T_min'] <= T < phase['T_max']:
            return phase
    return None


# ================================================================
#  Structure generation functions
# ================================================================

LATTICE_CONSTANTS = {
    ('Sn', 'diamond'): (6.489,),
    ('In', 'bct'): (3.252, 4.946),
}

ELEMENT_PRIMITIVE_ATOMS = {
    'fcc': 4,
    'bcc': 2,
    'hcp': 2,
    'dhcp': 4,
    'diamond': 8,
    'sc': 1,
    'bct': 2,
}


def generate_element_structure(element,
                               output_dir="structures",
                               target_atoms=100,
                               structure_type=None,
                               verbose=True):
    """
    Generate POSCAR file for element with automatic supercell expansion.

    Parameters
    ----------
    element : str
    output_dir : str
    target_atoms : int
    structure_type : str or None
        One of: 'fcc', 'bcc', 'hcp', 'diamond', 'sc', 'bct', 'dhcp'
        If None, auto-detect from ELEMENT_PHASE_DATA rt_structure
    verbose : bool

    Returns
    -------
    atoms : ase.Atoms
    poscar_path : str
    """
    if structure_type is None:
        # Use rt_structure from ELEMENT_PHASE_DATA
        rt = ELEMENT_PHASE_DATA.get(element, {}).get('rt_structure')
        if rt:
            structure_type = rt
        else:
            structure_type = 'fcc'

    lattice_key = (element, structure_type)
    has_custom_lattice = lattice_key in LATTICE_CONSTANTS

    if structure_type == 'dhcp':
        # dhcp: ABAC stacking, 4 atoms per unit cell
        # Use hcp primitive and double c-axis
        if has_custom_lattice:
            a, c4 = LATTICE_CONSTANTS[lattice_key]
            c = c4 / 2  # per layer
        else:
            # Get from ELEMENT_PHASE_DATA
            rt_phase = ELEMENT_PHASE_DATA.get(element, {}).get('phases', [{}])[0]
            a = rt_phase.get('a', 3.5)
            c4 = rt_phase.get('c', 11.8)
            c = c4 / 2

        atoms = bulk(element, 'hcp', a=a, c=c)

        # Build supercell with 4-layer dhcp stacking
        # hcp uses 2-layer AB, we need 4-layer ABAC
        # Simplification: just double the z-axis to get 4 atoms/cell
        supercell_matrix = np.eye(3, dtype=int)
        supercell_matrix[2, 2] = 2
        atoms = make_supercell(atoms, supercell_matrix)

        # Adjust positions for ABAC stacking
        pos = atoms.get_positions()
        cell = atoms.get_cell()
        frac = pos @ np.linalg.inv(cell)
        # For atoms originally at z~0.25 (layer B), shift to z~0.75 (layer C)
        for i in range(len(atoms)):
            if 0.4 < frac[i, 2] < 0.6:
                pass  # keep A and C layers
        atoms.set_positions(frac @ cell)

        # Target atom count
        primitive_atoms = 4  # ABAC = 4 layers
        n = int(np.floor((target_atoms / primitive_atoms) ** (1 / 3)))
        n = max(1, n)
        supercell_matrix = np.array([[n, 0, 0], [0, n, 0], [0, 0, max(1, n // 2)]], dtype=int)
        atoms = make_supercell(atoms, supercell_matrix)

    elif structure_type in ('fcc', 'bcc', 'diamond', 'sc'):
        if has_custom_lattice:
            a = LATTICE_CONSTANTS[lattice_key][0]
            atoms = bulk(element, structure_type, a=a, cubic=True)
        else:
            atoms = bulk(element, structure_type, cubic=True)
    elif structure_type == 'hcp':
        if has_custom_lattice:
            a, c = LATTICE_CONSTANTS[lattice_key]
            atoms = bulk(element, 'hcp', a=a, c=c)
        else:
            atoms = bulk(element, 'hcp')
    elif structure_type == 'bct':
        if has_custom_lattice:
            a, c = LATTICE_CONSTANTS[lattice_key]
        else:
            a, c = 3.252, 4.946
        positions = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])
        cell = np.array([[a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, c]])
        cart_positions = positions @ cell
        atoms = Atoms(element * 2, positions=cart_positions, cell=cell, pbc=True)
    else:
        raise ValueError(f"Unsupported structure type: {structure_type}")

    primitive_atoms = ELEMENT_PRIMITIVE_ATOMS.get(structure_type, 2)

    n = int(np.floor((target_atoms / primitive_atoms) ** (1 / 3)))
    n = max(1, n)
    supercell_matrix = np.eye(3, dtype=int) * n

    supercell = make_supercell(atoms, supercell_matrix)
    n_atoms = len(supercell)

    if n_atoms > target_atoms and structure_type in ('fcc', 'bcc', 'diamond', 'sc', 'bct'):
        n = int(np.floor((target_atoms / primitive_atoms) ** (1 / 3)))
        n = max(1, n)
        supercell_matrix = np.eye(3, dtype=int) * n
        supercell = make_supercell(atoms, supercell_matrix)
        n_atoms = len(supercell)

    os.makedirs(output_dir, exist_ok=True)
    poscar_path = os.path.join(output_dir, f"{element}-{structure_type.upper()}.POSCAR")
    write(poscar_path, supercell, format='vasp', direct=True)

    if verbose:
        print(f"Generated {element} ({structure_type.upper()}): {n_atoms} atoms -> {poscar_path}")

    return supercell, poscar_path


def batch_generate_structures(elements, output_dir="structures", target_atoms=100, structure_type=None):
    """Generate POSCAR files for multiple elements."""
    results = {}
    for element in elements:
        try:
            atoms, poscar_path = generate_element_structure(
                element, output_dir, target_atoms,
                structure_type=structure_type, verbose=True
            )
            results[element] = {'atoms': atoms, 'path': poscar_path, 'n_atoms': len(atoms)}
        except Exception as e:
            print(f"Error generating structure for {element}: {e}")
            results[element] = None
    return results


def generate_all_typical_elements(output_dir="structures", target_atoms=100):
    """Generate POSCAR files for all typical metallic elements."""
    results = {}
    total_count = 0
    success_count = 0

    # Only generate for elements with 'rt_structure' set
    elements_to_generate = {}
    for sym, data in ELEMENT_PHASE_DATA.items():
        rt = data.get('rt_structure')
        if rt and rt in ('fcc', 'bcc', 'hcp', 'diamond', 'sc', 'bct', 'dhcp'):
            if rt not in elements_to_generate:
                elements_to_generate[rt] = []
            elements_to_generate[rt].append(sym)

    print(f"Generating POSCAR files for typical metallic elements (target: {target_atoms} atoms)...")
    print("=" * 70)

    for structure_type, elements_list in elements_to_generate.items():
        print(f"\n{structure_type.upper()} structures:")
        for element in elements_list:
            total_count += 1
            try:
                atoms, poscar_path = generate_element_structure(
                    element, output_dir=output_dir, target_atoms=target_atoms,
                    structure_type=structure_type, verbose=True
                )
                results[element] = {
                    'structure': structure_type,
                    'atoms': atoms, 'path': poscar_path, 'n_atoms': len(atoms)
                }
                success_count += 1
            except Exception as e:
                print(f"  ERROR: Failed to generate {element} ({structure_type}): {e}")
                results[element] = None

    print("\n" + "=" * 70)
    print(f"Generation complete: {success_count}/{total_count} successful")
    print(f"Output directory: {output_dir}")
    return results

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



def _get_lattice_from_data(symbol, structure_type):
    """Look up lattice constants from ELEMENT_PHASE_DATA."""
    data = ELEMENT_PHASE_DATA.get(symbol, {})
    for phase in data.get("phases", []):
        if phase["structure"] == structure_type:
            a = phase.get("a")
            c = phase.get("c")
            if c:
                return (a, c)
            return (a,)
    raise ValueError(f"{symbol} {structure_type} not in ELEMENT_PHASE_DATA")
def _get_mass_map(type_map):

    mass_map_ref = {
        'H': 1.00794, 'He': 4.0026,

        'Li': 6.941,  'Be': 9.0122,  'B': 10.811,  'C': 12.0107, 
        'N':  14.0067, 'O': 15.9994, 'F': 18.9984, 'Ne': 20.1797, 

        'Na': 22.9897, 'Mg': 24.305, 'Al': 26.9815, 'Si': 28.0855, 
        'P':  30.9738, 'S': 32.065, 'Cl': 35.453, 'Ar': 39.948,

        'K': 39.0983, 'Ca': 40.078, 'Sc': 44.9559, 'Ti': 47.867, 
        'V': 50.9415, 'Cr': 51.9961, 'Mn': 54.938, 'Fe': 55.845,
        'Co': 58.9332, 'Ni': 58.6934, 'Cu': 63.546, 'Zn': 65.38, 
        'Ga': 69.723, 'Ge': 72.64, 'As': 74.9216, 'Se': 78.96,
        'Br': 79.904, 'Kr': 83.798, 'Rb': 85.4678, 'Sr': 87.62, 
        'Y': 88.9059, 'Zr': 91.224, 'Nb': 92.9064, 'Mo': 95.94,
        'Tc': 98.9063, 'Ru': 101.07, 'Rh': 102.9055, 'Pd': 106.42, 
        'Ag': 107.8682, 'Cd': 112.411, 'In': 114.818, 'Sn': 118.710, 
        'Sb': 121.760, 'Te': 127.60, 'I': 126.904, 'Xe': 131.293, 
        
        'Cs': 132.9055, 'Ba': 137.327, 'La': 138.9055, 'Ce': 140.116, 
        'Pr': 140.9077, 'Nd': 144.242, 'Pm': 145.00, 'Sm': 150.36, 
        'Eu': 151.964,  'Gd': 157.25, 'Tb': 158.9253, 'Dy': 162.50, 
        'Ho': 164.9303, 'Er': 167.259, 'Tm': 168.9342, 'Yb': 173.04,
        'Lu': 174.967,  'Hf': 178.49, 'Ta': 180.9479, 'W': 183.84, 
        'Re': 186.207,  'Os': 190.23, 'Ir': 192.217,  'Pt': 195.084, 
        'Au': 196.9665, 'Hg': 200.59, 'Tl': 204.3833, 'Pb': 207.2, 
        'Bi': 208.9804, 'Po': 209.0,  'At': 210.0, 'Rn': 222.0, 
        'Fr': 223.0,    'Ra': 226.0,  'Ac': 227.0, 'Th': 232.0381, 
        'Pa': 231.03588, 'U': 238.02891, 'Np': 237.0482, 'Pu': 244.0642, 
        'Am': 243.06138, 'Cm': 247.0703, 'Bk': 247.0703, 'Cf': 251.0796, 
        'Es': 252.0829, 'Fm': 257.0951, 'Md': 258.0984, 'No': 259.1010, 
        'Lr': 262.1096, 'Rf': 261.1086, 'Db': 262.1138, 'Sg': 266.1218, 
        'Bh': 264.1254, 'Hs': 265.1388, 'Mt': 266.1516, 'Ds': 269.1649, 
        'Rg': 272.174, 'Cn': 277.1926, 'Nh': 278.1884, 'Fl': 281.1972, 
        'Mc': 282.1935, 'Lv': 285.2078, 'Ts': 286.2102, 'Og': 289.2184
    }

    mass_map = []
    for key in type_map:
        mass_map.append(mass_map_ref[key])

    return mass_map

def _get_lattice_str(lattice_type, lattice_param):

    if lattice_type == 'fcc':
        strs = 'lattice             fcc %.6f \n'%(lattice_param)
        
    elif lattice_type == 'bcc':
        strs = 'lattice             bcc %.6f \n'%(lattice_param)
    elif lattice_type == 'sc':
        strs = 'lattice             sc %.6f \n'%(lattice_param)
    elif lattice_type == 'hcp':
        strs = 'lattice             hcp %.6f %.6f \n'%(lattice_param[0], lattice_param[1])
    elif lattice_type == 'diamond':
        strs = 'lattice             diamond %.6f \n'%(lattice_param)
    
    elif lattice_type == 'Mg2SiO4-fo1':

        strs =  'variable        a equal 4.75430  \n'
        strs += 'variable        b equal 10.20100   \n'
        strs += 'variable        c equal 5.98190 \n'

        strs += 'lattice               custom        1.0                                &  \n'
        strs += '                      a1            $a         0.0        0.0          &  \n'
        strs += '                      a2            0.0        $b         0.0          &  \n'
        strs += '                      a3            0.0        0.0        $c           &  \n'
        strs += '                      basis         0.00000    0.00000    0.00000      &  \n'
        strs += '                      basis         0.00000    0.00000    0.00000      &  \n'   
        strs += '                      basis         0.50000    0.50000    0.00000      &  \n'  
        strs += '                      basis         0.00000    0.00000    0.50000      &  \n'
        strs += '                      basis         0.50000    0.50000    0.50000      &  \n'
        strs += '                      basis         0.98450    0.27420    0.25000      &  \n'
        strs += '                      basis         0.01550    0.72580    0.75000      &  \n'
        strs += '                      basis         0.48450    0.22580    0.75000      &  \n'
        strs += '                      basis         0.51550    0.77420    0.25000      &  \n'
        strs += '                      basis         0.42500    0.09640    0.25000      &  \n'
        strs += '                      basis         0.57500    0.90360    0.75000      &  \n'
        strs += '                      basis         0.92500    0.40360    0.75000      &  \n'
        strs += '                      basis         0.07500    0.59640    0.25000      &  \n'
        strs += '                      basis         0.76750    0.08990    0.25000      &  \n'
        strs += '                      basis         0.23250    0.91010    0.75000      &  \n'
        strs += '                      basis         0.26750    0.41010    0.75000      &  \n'
        strs += '                      basis         0.73250    0.58990    0.25000      &  \n'
        strs += '                      basis         0.22710    0.44070    0.25000      &  \n'
        strs += '                      basis         0.77290    0.55930    0.75000      &  \n'
        strs += '                      basis         0.72710    0.05930    0.75000      &  \n'
        strs += '                      basis         0.27290    0.94070    0.25000      &  \n'
        strs += '                      basis         0.26380    0.16920    0.02330      &  \n'
        strs += '                      basis         0.73620    0.83080    0.97670      &  \n'
        strs += '                      basis         0.76380    0.33080    0.97670      &  \n'
        strs += '                      basis         0.23620    0.66920    0.02330      &  \n'
        strs += '                      basis         0.73620    0.83080    0.52330      &  \n'
        strs += '                      basis         0.26380    0.16920    0.47670      &  \n'
        strs += '                      basis         0.23620    0.66920    0.47670      &  \n'
        strs += '                      basis         0.76380    0.33080    0.52330         \n\n'

        strs += 'region                box block 0 ${Nx} 0 ${Ny} 0 ${Nz} units lattice \n\n'
        strs += 'create_box            3 box \n'
        strs += 'create_atoms          3 box  & \n\n'

        strs += '                     basis   1   1   basis   2  1   basis   3  1   basis   4  1   & \n'
        strs += '                     basis   5   1   basis   6  1   basis   7  1   basis   8  1   & \n\n'
        strs += '                     basis   9   2   basis  10  2   basis  11  2   basis  12  2   & \n'
        strs += '                     basis  13   3   basis  14  3   basis  15  3   basis  16  3   & \n'
        strs += '                     basis  17   3   basis  18  3   basis  19  3   basis  20  3   & \n'
        strs += '                     basis  21   3   basis  22  3   basis  23  3   basis  24  3   & \n'
        strs += '                     basis  25   3   basis  26  3   basis  27  3   basis  28  3   \n\n'


    elif lattice_type == 'Mg2SiO4-fo2':

        strs =  'lattice         custom    1.0     & \n'
        strs += '                a1    4.6830000877      0.0000000000         0.0000000000  & \n'
        strs += '                a2   -1.2499409133      9.1247875820         0.0000000000  & \n'
        strs += '                a3   -1.5492150373     -0.4930850674         5.0623401655  & \n'
        strs += '                basis   0.9714  0.0966  0.6415  & \n'
        strs += '                basis   0.0419  0.9354  0.0567  & \n'
        strs += '                basis   0.0026  0.6656  0.7362  & \n'
        strs += '                basis   0.5065  0.516   0.8488  & \n'
        strs += '                basis   0.5386  0.7905  0.1388  & \n'
        strs += '                basis   0.0105  0.3666  0.9614  & \n'
        strs += '                basis   0.4746  0.2414  0.5588  & \n'
        strs += '                basis   0.5065  0.5161  0.3489  & \n'
        strs += '                basis   0.0048  0.3875  0.4564  & \n'
        strs += '                basis   0.0083  0.6445  0.2412  & \n'
        strs += '                basis   0.6136  0.1239  0.0951  & \n'
        strs += '                basis   0.3995  0.9081  0.6027  & \n'
        strs += '                basis   0.2702  0.843   0.8193  & \n'
        strs += '                basis   0.7865  0.6995  0.968   & \n'
        strs += '                basis   0.7628  0.2287  0.3717  & \n'
        strs += '                basis   0.2756  0.0474  0.4928  & \n'
        strs += '                basis   0.7429  0.1891  0.8786  & \n'
        strs += '                basis   0.2266  0.3325  0.7297  & \n'
        strs += '                basis   0.2506  0.8032  0.3262  & \n'
        strs += '                basis   0.7833  0.6947  0.4029  & \n'
        strs += '                basis   0.7924  0.4622  0.1735  & \n'
        strs += '                basis   0.7372  0.9845  0.2048  & \n'
        strs += '                basis   0.2204  0.5698  0.5243  & \n'
        strs += '                basis   0.776   0.4475  0.645   & \n'
        strs += '                basis   0.2685  0.1109  0.0031  & \n'
        strs += '                basis   0.237   0.5847  0.0528  & \n'
        strs += '                basis   0.23    0.3374  0.2947  & \n'
        strs += '                basis   0.7446  0.9211  0.695 \n\n'

        strs += 'region          box block 0 ${Nx} 0 ${Ny} 0 ${Nz} units lattice \n\n'

        strs += 'create_box      3 box \n'
        strs += 'create_atoms    3 box  & \n'
        strs += '                basis   1   1   basis   2  1   basis   3  1   basis   4  1   & \n'
        strs += '                basis   5   1   basis   6  1   basis   7  1   basis   8  1   & \n'
        strs += '                basis   9   2   basis  10  2   basis  11  2   basis  12  2   & \n'
        strs += '                basis  13   3   basis  14  3   basis  15  3   basis  16  3   & \n'
        strs += '                basis  17   3   basis  18  3   basis  19  3   basis  20  3   & \n'
        strs += '                basis  21   3   basis  22  3   basis  23  3   basis  24  3   & \n'
        strs += '                basis  25   3   basis  26  3   basis  27  3   basis  28  3   \n\n'

    elif lattice_type == 'Mg2SiO4-fo3':

        strs =  'variable        a equal 2.6400 \n'
        strs += 'variable        b equal 8.5960   \n'
        strs += 'variable        c equal 9.0400 \n'

        strs += 'lattice         custom    1.0     & \n'
        strs += '                a1      $a       0.0     0.0     & \n'
        strs += '                a2      0.0      $b      0.0     & \n'
        strs += '                a3      0.0      0.0     $c      & \n'
        strs += '                basis         0.50000    0.37300    0.35000  &  \n'
        strs += '                basis         0.50000    0.62700    0.85000  & \n'
        strs += '                basis         0.00000    0.87300    0.35000  & \n'
        strs += '                basis         0.00000    0.12700    0.85000  & \n'
        strs += '                basis         0.50000    0.88900    0.66200  &  \n'
        strs += '                basis         0.50000    0.11100    0.16200  &  \n'
        strs += '                basis         0.00000    0.38900    0.66200  &  \n'
        strs += '                basis         0.00000    0.61100    0.16200  &  \n'
        strs += '                basis         0.00000    0.85700    0.00000  &   \n'
        strs += '                basis         0.00000    0.14300    0.50000  &  \n'
        strs += '                basis         0.50000    0.35700    0.00000  &  \n'
        strs += '                basis         0.50000    0.64300    0.50000  &  \n'
        strs += '                basis         0.50000    0.98300    0.46000  &   \n'
        strs += '                basis         0.50000    0.01700    0.96000  &  \n'
        strs += '                basis         0.00000    0.48300    0.46000  &  \n'
        strs += '                basis         0.00000    0.51700    0.96000  &  \n'
        strs += '                basis         0.50000    0.23100    0.57200  &  \n'
        strs += '                basis         0.50000    0.76900    0.07200  &  \n'
        strs += '                basis         0.00000    0.73100    0.57200  &  \n'
        strs += '                basis         0.00000    0.26900    0.07200  &  \n'
        strs += '                basis         0.50000    0.70600    0.27100  &  \n'
        strs += '                basis         0.50000    0.29400    0.77100  &  \n'
        strs += '                basis         0.00000    0.20600    0.27100  &  \n'
        strs += '                basis         0.00000    0.79400    0.77100  &  \n'
        strs += '                basis         0.50000    0.46600    0.18100  & \n'
        strs += '                basis         0.50000    0.53400    0.68100  & \n'
        strs += '                basis         0.00000    0.96600    0.18100  &  \n'
        strs += '                basis         0.00000    0.03400    0.68100    \n\n'


        strs += 'region          box block 0 ${Nx} 0 ${Ny} 0 ${Nz} units lattice \n\n'

        strs += 'create_box      3 box \n'
        strs += 'create_atoms    3 box  & \n'
        strs += '                basis   1   1   basis   2  1   basis   3  1   basis   4  1   & \n'
        strs += '                basis   5   1   basis   6  1   basis   7  1   basis   8  1   & \n'
        strs += '                basis   9   2   basis  10  2   basis  11  2   basis  12  2   & \n'
        strs += '                basis  13   3   basis  14  3   basis  15  3   basis  16  3   & \n'
        strs += '                basis  17   3   basis  18  3   basis  19  3   basis  20  3   & \n'
        strs += '                basis  21   3   basis  22  3   basis  23  3   basis  24  3   & \n'
        strs += '                basis  25   3   basis  26  3   basis  27  3   basis  28  3   \n\n'

    else:
        raise ValueError('Lattice type %s not supported'%(lattice_type))
    
    return strs


# ================================================================
#  Phase segment splitting (for automated DPGEN exploration)
# ================================================================

SKIP_MAGNETIC = {'Fe', 'Co', 'Ni', 'Cr', 'Mn',
                 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm'}
SKIP_RADIOACTIVE = {'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac',
                    'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm'}
SKIP_GAS = {'H', 'He', 'N', 'O', 'F', 'Ne', 'Cl', 'Ar',
            'Br', 'Kr', 'I', 'Xe', 'Rn', 'P', 'S', 'Se', 'Te'}
SKIP_LANTHANIDE = {str(el) for el in
                   ['La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu',
                    'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu']}
SKIP_ACTINIDE = {'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm'}
SUPPORTED_STRUCTURES = {'fcc', 'bcc', 'hcp', 'dhcp',
                        'diamond', 'sc', 'bct'}


def get_viable_elements():
    result = []
    for sym, data in ELEMENT_PHASE_DATA.items():
        rt = data.get('rt_structure')
        if rt is None or rt not in SUPPORTED_STRUCTURES:
            continue
        if data.get('Tm') is None:
            continue
        if sym in SKIP_MAGNETIC:
            continue
        if sym in SKIP_RADIOACTIVE:
            continue
        if sym in SKIP_GAS:
            continue
        if sym in SKIP_LANTHANIDE:
            continue
        if len(data.get('phases', [])) == 0:
            continue
        result.append(sym)
    return sorted(result, key=lambda s: ELEMENT_PHASE_DATA[s].get('z', 0))


def get_phase_segments(element, T_range=(300, None), overlap_rule='auto',
                       drop_below_T=200):
    data = ELEMENT_PHASE_DATA.get(element)
    if data is None:
        raise ValueError(f"Element '{element}' not in ELEMENT_PHASE_DATA")
    Tm = data.get('Tm')
    if Tm is None:
        raise ValueError(f"Element '{element}' has no Tm (sublimes)")
    T_min_range, T_max_range = T_range
    if T_max_range is None:
        T_max_range = 2 * Tm
    phases = data['phases']
    if not phases:
        raise ValueError(f"Element '{element}' has no phases data")
    segs = []
    for phase in phases:
        st = phase['structure']
        if st not in SUPPORTED_STRUCTURES:
            continue
        if phase['T_max'] < drop_below_T:
            continue
        core_low = max(phase['T_min'], T_min_range)
        core_high = min(phase['T_max'], T_max_range)
        if core_low >= core_high:
            continue
        segs.append({
            'label': element + '-' + st.upper(),
            'structure': st,
            'T_core': (core_low, core_high),
        })
    if not segs:
        raise ValueError(
            f"Element '{element}': no phases in T_range "
            f"{T_min_range}-{T_max_range}K after filtering")
    _deduplicate_labels(segs)
    for i in range(len(segs)):
        if i > 0:
            T_lower_cut = segs[i]['T_core'][0]
            delta_lower = _compute_overlap(T_lower_cut,
                                           is_melting=False,
                                           overlap_rule=overlap_rule)
        else:
            delta_lower = 0
        if i < len(segs) - 1:
            T_upper_cut = segs[i]['T_core'][1]
            delta_upper = _compute_overlap(T_upper_cut,
                                           is_melting=(T_upper_cut >= Tm * 0.99),
                                           overlap_rule=overlap_rule)
        else:
            delta_upper = _compute_overlap(Tm, is_melting=True,
                                           overlap_rule=overlap_rule)
        segs[i]['T_explore'] = (
            max(segs[i]['T_core'][0] - delta_lower, T_min_range),
            min(segs[i]['T_core'][1] + delta_upper, T_max_range),
        )
    last_structure = segs[-1]['structure']
    segs.append({
        'label': element + '-LIQ',
        'structure': last_structure,
        'T_core': (Tm, T_max_range),
        'T_explore': (max(Tm - _compute_overlap(Tm, True, overlap_rule),
                         T_min_range),
                       T_max_range),
    })
    return segs


def _deduplicate_labels(segs):
    seen = {}
    for seg in segs:
        lbl = seg['label']
        if lbl not in seen:
            seen[lbl] = 0
        else:
            seen[lbl] += 1
            seg['label'] = f"{lbl}-{seen[lbl] + 1}"


def _compute_overlap(T_ref, is_melting, overlap_rule):
    if overlap_rule == 'auto':
        return max(0.2 * T_ref, 100)
    elif callable(overlap_rule):
        return overlap_rule(T_ref, is_melting)
    else:
        return float(overlap_rule)
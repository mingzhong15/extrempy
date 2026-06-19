import os
import glob
import shutil
import copy

POTCAR_MAP = {
    'PBE54': {
        # Alkali metals
        'Li': {'variant': '', 'ZVAL': 0},
        'Na': {'variant': '', 'ZVAL': 0},
        'K':  {'variant': '', 'ZVAL': 0},
        'Rb': {'variant': '', 'ZVAL': 0},
        'Cs': {'variant': '', 'ZVAL': 0},
        # Alkaline earth
        'Be': {'variant': '', 'ZVAL': 0},
        'Mg': {'variant': '', 'ZVAL': 0},
        'Ca': {'variant': '', 'ZVAL': 0},
        'Sr': {'variant': '', 'ZVAL': 0},
        'Ba': {'variant': '_sv', 'ZVAL': 0},
        # Post-transition / simple metals
        'Al': {'variant': '', 'ZVAL': 0},
        # Ga/Sb/Bi excluded from get_viable_elements (no rt_structure for
        # orthorhombic/rhombohedral — complex crystal, not ASE-bulk-friendly)
        'Ga': {'variant': '', 'ZVAL': 0},
        'In': {'variant': '', 'ZVAL': 0},
        'Sn': {'variant': '', 'ZVAL': 0},
        'Tl': {'variant': '_sv', 'ZVAL': 0},
        'Pb': {'variant': '_sv', 'ZVAL': 0},
        'Bi': {'variant': '_sv', 'ZVAL': 0},
        # Metalloids / semiconductors
        'Si': {'variant': '', 'ZVAL': 0},
        'Ge': {'variant': '', 'ZVAL': 0},
        'Sb': {'variant': '', 'ZVAL': 0},
        # Transition metals (3d)
        'Sc': {'variant': '_pv', 'ZVAL': 0},
        'Ti': {'variant': '_pv', 'ZVAL': 0},
        'V':  {'variant': '_pv', 'ZVAL': 0},
        'Cu': {'variant': '_pv', 'ZVAL': 0},
        'Zn': {'variant': '_pv', 'ZVAL': 0},
        # Transition metals (4d)
        'Y':  {'variant': '_pv', 'ZVAL': 0},
        'Zr': {'variant': '_pv', 'ZVAL': 0},
        'Nb': {'variant': '_pv', 'ZVAL': 0},
        'Mo': {'variant': '_pv', 'ZVAL': 0},
        'Tc': {'variant': '_pv', 'ZVAL': 0},
        'Ru': {'variant': '_pv', 'ZVAL': 0},
        'Rh': {'variant': '_pv', 'ZVAL': 0},
        'Pd': {'variant': '_pv', 'ZVAL': 0},
        'Ag': {'variant': '_pv', 'ZVAL': 0},
        'Cd': {'variant': '_pv', 'ZVAL': 0},
        # Transition metals (5d)
        'Hf': {'variant': '_pv', 'ZVAL': 0},
        'Ta': {'variant': '_pv', 'ZVAL': 0},
        'W':  {'variant': '_pv', 'ZVAL': 0},
        'Re': {'variant': '_pv', 'ZVAL': 0},
        'Os': {'variant': '_pv', 'ZVAL': 0},
        'Ir': {'variant': '_pv', 'ZVAL': 0},
        'Pt': {'variant': '_pv', 'ZVAL': 0},
        'Au': {'variant': '_pv', 'ZVAL': 0},
    },
    'PBE52': {},
}

FALLBACK_VARIANTS = ['', '_pv', '_sv', '_d', '_h']


class PotcarMap:
    def __init__(self, set_name='PBE54', potcar_dir='.'):
        if set_name not in POTCAR_MAP:
            raise KeyError(f"POTCAR_MAP set '{set_name}' not found. "
                           f"Available: {list(POTCAR_MAP.keys())}")
        self.set_name = set_name
        self.potcar_dir = potcar_dir
        self.map = copy.deepcopy(POTCAR_MAP[set_name])

    def _resolve_path(self, element):
        entry = self.map.get(element)
        if entry is None:
            raise KeyError(f"Element '{element}' not in POTCAR_MAP['{self.set_name}']. "
                           f"Add it or use a different set.")
        preferred = entry.get('variant', '')
        candidates = [preferred] if preferred else ['']
        candidates += [v for v in FALLBACK_VARIANTS if v not in candidates]
        for variant in candidates:
            name = element + variant
            potcar_path = os.path.join(self.potcar_dir, name, 'POTCAR')
            if os.path.exists(potcar_path):
                return potcar_path, variant
        raise FileNotFoundError(
            f"No POTCAR found for {element} in {self.potcar_dir}. "
            f"Tried variants: {candidates}")

    @staticmethod
    def _parse_zval(potcar_path):
        with open(potcar_path, 'r') as f:
            for line in f:
                if 'ZVAL' in line:
                    val_str = line.split('=')[-1].strip().split()[0]
                    return float(val_str)
        raise ValueError(f"No ZVAL found in {potcar_path}")

    @staticmethod
    def _parse_titel(potcar_path):
        with open(potcar_path, 'r') as f:
            for line in f:
                if 'TITEL' in line:
                    return line.strip()
        raise ValueError(f"No TITEL found in {potcar_path}")

    def build(self, elements):
        potcar_paths = []
        zval_list = []
        for element in elements:
            path, variant_used = self._resolve_path(element)
            zval = self._parse_zval(path)
            titel = self._parse_titel(path)
            entry = self.map[element]
            if entry['ZVAL'] != 0:
                if abs(entry['ZVAL'] - zval) > 0.01:
                    raise ValueError(
                        f"ZVAL mismatch for {element}: "
                        f"POTCAR_MAP expects {entry['ZVAL']}, "
                        f"but {path} has {zval}")
            else:
                entry['ZVAL'] = zval
            if variant_used != entry.get('variant', '') and entry.get('variant') is not None:
                print(f"NOTE: {element} prefers variant '{entry.get('variant')}', "
                      f"but using '{variant_used}' (auto fallback)")
            zval_list.append(zval)
            potcar_paths.append(path)
            print(f"  {element}: ZVAL={zval:g}, TITEL={titel}")
        potcar_text = b''
        for path in potcar_paths:
            with open(path, 'rb') as f:
                potcar_text += f.read()
        return potcar_text, zval_list, potcar_paths

    def write_potcar(self, elements, output_path):
        text, zvals, paths = self.build(elements)
        with open(output_path, 'wb') as f:
            f.write(text)
        return zvals, paths

from typing import Dict, List, Tuple, Optional, Union

class OpenMXReader:
    """
    Handler for OpenMX input files (in.dat).

    This class can read, modify, and generate OpenMX input files.
    """

    def __init__(self, input_path: Optional[str] = None):
        """
        Initialize OpenMX input handler.

        Args:
            input_path: Path to existing OpenMX input file. If None, creates empty handler.
        """
        self.input_path = input_path
        self.lines: List[str] = []
        self.variables: Dict[str, str] = {}  # Key-value pairs for simple variables
        
        # Structure sections (indices in self.lines)
        self.atoms_number_idx: Optional[int] = None
        self.species_coord_unit_idx: Optional[int] = None
        self.structure_start_idx: Optional[int] = None
        self.structure_end_idx: Optional[int] = None
        self.unit_vectors_unit_idx: Optional[int] = None
        self.unit_vectors_start_idx: Optional[int] = None
        self.unit_vectors_end_idx: Optional[int] = None
        
        if input_path:
            self.read(input_path)

    def read(self, input_path: str):
        """
        Read and parse OpenMX input file.

        Args:
            input_path: Path to OpenMX input file (in.dat)
        """
        self.input_path = input_path
        with open(input_path, 'r', encoding='utf-8') as f:
            self.lines = [line.rstrip('\n\r') for line in f]
        
        self._parse()

    def _parse(self):
        """Parse the loaded file to extract variables and structure sections."""
        self.variables = {}
        
        # Find structure section boundaries
        # We need to find complete structure section including:
        # - Atoms.Number line
        # - Atoms.SpeciesAndCoordinates.Unit line
        # - <Atoms.SpeciesAndCoordinates ... Atoms.SpeciesAndCoordinates>
        # - Atoms.UnitVectors.Unit line
        # - <Atoms.UnitVectors ... Atoms.UnitVectors>
        
        self.atoms_number_idx = None
        self.species_coord_unit_idx = None
        self.structure_start_idx = None
        self.structure_end_idx = None
        self.unit_vectors_unit_idx = None
        self.unit_vectors_start_idx = None
        self.unit_vectors_end_idx = None
        
        for i, line in enumerate(self.lines):
            # Parse simple key-value pairs (with optional comments)
            stripped = line.strip()
            if not stripped or stripped.startswith('#'):
                continue
            
            # Check for Atoms.Number
            if stripped.startswith('Atoms.Number'):
                self.atoms_number_idx = i
            # Check for Atoms.SpeciesAndCoordinates.Unit
            elif stripped.startswith('Atoms.SpeciesAndCoordinates.Unit'):
                self.species_coord_unit_idx = i
            # Check for structure sections
            elif '<Atoms.SpeciesAndCoordinates' in line:
                self.structure_start_idx = i
            elif 'Atoms.SpeciesAndCoordinates>' in line:
                self.structure_end_idx = i
            # Check for Atoms.UnitVectors.Unit
            elif stripped.startswith('Atoms.UnitVectors.Unit'):
                self.unit_vectors_unit_idx = i
            elif '<Atoms.UnitVectors' in line:
                self.unit_vectors_start_idx = i
            elif 'Atoms.UnitVectors>' in line:
                self.unit_vectors_end_idx = i
            
            # Parse simple variables (key value # comment)
            # Match patterns like "key value" or "key value # comment"
            parts = stripped.split('#', 1)
            main_part = parts[0].strip()
            
            if not main_part or main_part.startswith('<') or main_part.endswith('>'):
                continue
            
            # Try to parse as key-value pair
            words = main_part.split()
            if len(words) >= 2:
                key = words[0]
                value = ' '.join(words[1:])
                self.variables[key] = value

    def set_variable(self, key: str, value: Union[str, int, float]):
        """
        Set or modify a configuration variable.

        Args:
            key: Variable name (e.g., 'System.Name', 'scf.energycutoff')
            value: Value to set (will be converted to string)
        """
        value_str = str(value)
        self.variables[key] = value_str
        
        # Update in lines if it exists
        found = False
        for i, line in enumerate(self.lines):
            stripped = line.strip()
            if not stripped or stripped.startswith('#'):
                continue
            
            parts = stripped.split('#', 1)
            main_part = parts[0].strip()
            words = main_part.split()
            
            if len(words) >= 2 and words[0] == key:
                # Update the line, preserve comment if exists
                comment = ' # ' + parts[1] if len(parts) > 1 else ''
                # Try to preserve original spacing
                indent = len(line) - len(line.lstrip())
                new_line = ' ' * indent + key + '  ' * (30 - len(key)) + value_str + comment
                self.lines[i] = new_line.rstrip()
                found = True
                break
        
        # If not found, add to the end before structure section or at the end
        if not found:
            insert_idx = self.structure_start_idx if self.structure_start_idx is not None else len(self.lines)
            self.lines.insert(insert_idx, f"{key:<20} {value_str}")
            # Re-parse to update indices
            self._parse()

    def get_variable(self, key: str, default: Optional[str] = None) -> Optional[str]:
        """
        Get a configuration variable value.

        Args:
            key: Variable name
            default: Default value if key not found

        Returns:
            Variable value as string, or default if not found
        """
        return self.variables.get(key, default)

    def write(self, output_path: str):
        """
        Write OpenMX input file.

        Args:
            output_path: Path to output file
        """
        with open(output_path, 'w', encoding='utf-8') as f:
            for line in self.lines:
                f.write(line + '\n')

    def __str__(self) -> str:
        """Return string representation of the input file."""
        return '\n'.join(self.lines)


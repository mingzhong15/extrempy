import numpy as np
import os
import matplotlib.pyplot as plt
import glob
import json
import copy
import shutil
import subprocess
from typing import List, Tuple, Dict, Optional, Union, Set, Any

# Color mapping for different elements
element_colors = {
    'H': '#FFFFFF', 'He': '#D9FFFF', 'Li': '#CC80FF', 'Be': '#C2FF00',
    'B': '#FFB5B5', 'C': '#909090', 'N': '#3050F8', 'O': '#FF0D0D',
    'F': '#90E050', 'Ne': '#B3E3F5', 'Na': '#AB5CF2', 'Mg': '#8AFF00',
    'Al': '#BFA6A6', 'Si': '#F0C8A0', 'P': '#FF8000', 'S': '#FFFF30',
    'Cl': '#1FF01F', 'Ar': '#80D1E3', 'K': '#8F40D4', 'Ca': '#3DFF00',
    'Sc': '#E6E6E6', 'Ti': '#BFC2C7', 'V': '#A6A6AB', 'Cr': '#8A99C7',
    'Mn': '#9C7AC7', 'Fe': '#E06633', 'Co': '#F090A0', 'Ni': '#50D050',
    'Cu': '#C88033', 'Zn': '#7D80B0', 'Ga': '#C28F8F', 'Ge': '#668F8F',
    'As': '#BD80E3', 'Se': '#FFA100', 'Br': '#A62929', 'Kr': '#5CB8D1',
    'Rb': '#702EB0', 'Sr': '#00FF00', 'Y': '#94FFFF', 'Zr': '#94E0E0',
    'Nb': '#73C2C9', 'Mo': '#54B5B5', 'Tc': '#3B9E9E', 'Ru': '#248F8F',
    'Rh': '#0A7D8C', 'Pd': '#006985', 'Ag': '#C0C0C0', 'Cd': '#FFD98F',
    'In': '#A67573', 'Sn': '#668080', 'Sb': '#9E63B5', 'Te': '#D47A00',
    'I': '#940094', 'Xe': '#429EB0', 'Cs': '#57178F', 'Ba': '#00C900',
    'La': '#70D4FF', 'Ce': '#FFFFC7', 'Pr': '#D9FFC7', 'Nd': '#C7FFC7',
    'Pm': '#A3FFC7', 'Sm': '#8FFFC7', 'Eu': '#61FFC7', 'Gd': '#45FFC7',
    'Tb': '#30FFC7', 'Dy': '#1FFFC7', 'Ho': '#00FF9C', 'Er': '#00E675',
    'Tm': '#00D452', 'Yb': '#00BF38', 'Lu': '#00AB24', 'Hf': '#4DC2FF',
    'Ta': '#4DA6FF', 'W': '#2194D6', 'Re': '#267DAB', 'Os': '#266696',
    'Ir': '#175487', 'Pt': '#D0D0E0', 'Au': '#FFD123', 'Hg': '#B8B8D0',
    'Tl': '#A6544D', 'Pb': '#575961', 'Bi': '#9E4FB5', 'Po': '#AB5C00',
    'At': '#754F45', 'Rn': '#428296', 'Fr': '#420066', 'Ra': '#007D00',
    'Ac': '#70ABFA', 'Th': '#00BAFF', 'Pa': '#00A1FF', 'U': '#008FFF',
    'Np': '#0080FF', 'Pu': '#006BFF', 'Am': '#545CF2', 'Cm': '#785CE3',
    'Bk': '#8A4FE3', 'Cf': '#A136D4', 'Es': '#B31FD4', 'Fm': '#B31FBA',
    'Md': '#B30DA6', 'No': '#BD0D87', 'Lr': '#C70066'
}

"""
Comprehensive POSCAR processing classes for 2D materials
Author: Assistant
Date: 2024

Structure:
- POSCARReader: Base class for reading and parsing POSCAR files
- POSCARProcessor: Extended class for processing, manipulation, and visualization

Features:
- Read and parse POSCAR files
- Lattice operations (rotation, scaling, shearing, transformation)
- Create bilayer structures
- Apply sliding operations along different crystal axes
- Random perturbations of atomic positions and layer spacing
- Export modified structures
- 3D visualization with lattice vectors and enhanced plotting
"""


class POSCARReader:
    """
    Base class for reading and parsing POSCAR files
    Provides core functionality for loading and accessing POSCAR data
    """
    
    def __init__(self, poscar_file: Optional[str] = None, verbose: bool = True):
        """
        Initialize POSCARReader
        
        Args:
            poscar_file: Path to POSCAR file (optional)
            verbose: Whether to print loading messages (default: True)
        """
        self.poscar_file = poscar_file
        self.verbose = verbose
        self.scale = 1.0
        self.lattice = np.zeros((3, 3))
        self.elements = []
        self.n_atoms = []
        self.coord_type = "cartesian"
        self.coords = np.array([])
        self.atom_types = []
        self.total_atoms = 0
        self.comment = ""  # Initialize comment field
        
        if poscar_file:
            self.load_poscar(poscar_file)

    def load_poscar(self, filename: str) -> None:
        """
        Load POSCAR file and parse all information with automatic coordinate conversion
        
        Args:
            filename: Path to POSCAR file
        """
        if not os.path.exists(filename):
            raise FileNotFoundError(f"POSCAR file {filename} not found!")
        
        with open(filename, 'r') as f:
            lines = f.readlines()
        
        # Read comment (first line)
        if len(lines) > 0:
            self.comment = lines[0].strip()
        else:
            self.comment = ""
        
        # Read scaling factor
        self.scale = float(lines[1].strip())
        
        # Read lattice vectors
        self.lattice = np.zeros((3, 3))
        for i in range(3):
            self.lattice[i] = [float(x) for x in lines[2+i].split()]
        
        # Read element names and merge duplicates
        raw_elements = lines[5].split()
        raw_n_atoms = [int(x) for x in lines[6].split()]
        
        # Merge duplicate elements
        self.elements = []
        self.n_atoms = []
        element_counts = {}
        
        for element, count in zip(raw_elements, raw_n_atoms):
            if element in element_counts:
                element_counts[element] += count
            else:
                element_counts[element] = count
        
        # Convert to lists maintaining order
        for element in raw_elements:
            if element not in self.elements:
                self.elements.append(element)
                self.n_atoms.append(element_counts[element])
        
        # Read coordinate type
        self.coord_type = lines[7].strip().lower()
        
        # Read atomic coordinates
        self.coords = []
        self.atom_types = []
        line_idx = 8
        
        for i, n in enumerate(self.n_atoms):
            for j in range(n):
                coord = [float(x) for x in lines[line_idx].split()[:3]]
                self.coords.append(coord)
                self.atom_types.append(self.elements[i])
                line_idx += 1
        
        self.coords = np.array(self.coords)
        self.total_atoms = len(self.coords)
        
        # Convert fractional coordinates to cartesian if needed
        if self.coord_type == 'direct':
            # Convert fractional coordinates to cartesian coordinates
            # Cartesian = Fractional * Lattice_matrix
            self.coords = np.dot(self.coords, self.lattice)
            # Update coordinate type to reflect the conversion
            self.coord_type = 'cartesian'
            if self.verbose:
                print(f"Converted fractional coordinates to cartesian coordinates")
        
        if self.verbose:
            print(f"Successfully loaded POSCAR: {self.total_atoms} atoms, {len(self.elements)} elements")
    
    def get_structure_info(self) -> Dict:
        """
        Get comprehensive structure information
        
        Returns:
            Dictionary containing all structure information
        """
        return {
            'scale': self.scale,
            'lattice': self.lattice.copy(),
            'elements': self.elements.copy(),
            'n_atoms': self.n_atoms.copy(),
            'coord_type': self.coord_type,
            'coords': self.coords.copy(),
            'atom_types': self.atom_types.copy(),
            'total_atoms': self.total_atoms,
            'lattice_parameters': self.get_lattice_parameters(),
            'layer_thickness': self.get_layer_thickness(),
            'bond_length': self.get_bond_length()
        }
    
    def get_lattice_parameters(self) -> Dict:
        """
        Calculate lattice parameters (a, b, c, alpha, beta, gamma)
        
        Returns:
            Dictionary with lattice parameters
        """
        a = np.linalg.norm(self.lattice[0])
        b = np.linalg.norm(self.lattice[1])
        c = np.linalg.norm(self.lattice[2])
        
        alpha = np.arccos(np.dot(self.lattice[1], self.lattice[2]) / (b * c)) * 180 / np.pi
        beta = np.arccos(np.dot(self.lattice[0], self.lattice[2]) / (a * c)) * 180 / np.pi
        gamma = np.arccos(np.dot(self.lattice[0], self.lattice[1]) / (a * b)) * 180 / np.pi
        
        return {
            'a': a, 'b': b, 'c': c,
            'alpha': alpha, 'beta': beta, 'gamma': gamma
        }
    
    def get_layer_thickness(self) -> float:
        """
        Calculate layer thickness along c-axis
        
        Returns:
            Layer thickness in Angstrom
        """
        if len(self.coords) == 0:
            return 0.0
        z_coords = self.coords[:, 2]
        return np.max(z_coords) - np.min(z_coords)

    def _get_coords_with_pbc(self) -> np.ndarray:
        """
        Get coordinates with periodic boundary conditions
        """
        coords = self.coords.copy()
        coords[:, 0] = np.mod(coords[:, 0], self.lattice[0, 0])
        coords[:, 1] = np.mod(coords[:, 1], self.lattice[1, 1])
        coords[:, 2] = np.mod(coords[:, 2], self.lattice[2, 2])
        return coords

    def get_bond_length(self) -> float:
        """
        Calculate bond length based on element types
        
        Returns:
            Bond length in Angstrom
        """
        if len(self.coords) == 0:
            return 0.0
        
        x_coords = self.coords[:, 0]
        y_coords = self.coords[:, 1]
        z_coords = self.coords[:, 2]

        extend_coords = self._get_coords_with_pbc()

        xij = x_coords.reshape(-1,1) - extend_coords[:,0].reshape(1,-1)
        yij = y_coords.reshape(-1,1) - extend_coords[:,1].reshape(1,-1)
        zij = z_coords.reshape(-1,1) - extend_coords[:,2].reshape(1,-1)

        rij = np.sqrt(xij**2 + yij**2 + zij**2)

        self.bond_length = np.min(rij[rij>0])

        return self.bond_length

    def print_summary(self) -> None:
        """
        Print structure summary
        """
        print("\n" + "="*50)
        print("POSCAR Structure Summary")
        print("="*50)
        print(f"Comment: {self.comment}")
        print(f"Total atoms: {self.total_atoms}")
        print(f"Elements: {', '.join(self.elements)}")
        print(f"Atom counts: {dict(zip(self.elements, self.n_atoms))}")
        
        lattice_params = self.get_lattice_parameters()
        print(f"\nLattice parameters:")
        print(f"  a = {lattice_params['a']:.3f} Å")
        print(f"  b = {lattice_params['b']:.3f} Å") 
        print(f"  c = {lattice_params['c']:.3f} Å")
        print(f"  α = {lattice_params['alpha']:.1f}°")
        print(f"  β = {lattice_params['beta']:.1f}°")
        print(f"  γ = {lattice_params['gamma']:.1f}°")
        
        print(f"Layer thickness: {self.get_layer_thickness():.3f} Å")
        print(f"Coordinate type: {self.coord_type}")
        print(f"Bond length: {self.get_bond_length():.3f} Å")
        print("="*50)


class POSCARProcessor(POSCARReader):
    """
    Extended class for processing, manipulating, and visualizing POSCAR structures
    Inherits from POSCARReader and adds processing capabilities including:
    - Structure manipulation (bilayer creation, centering)
    - Lattice operations (rotation, scaling, transformation)
    - Visualization and plotting
    - VASP file export
    """
    
    def __init__(self, poscar_file: Optional[str] = None, verbose: bool = True):
        """
        Initialize POSCARProcessor
        
        Args:
            poscar_file: Path to POSCAR file (optional)
            verbose: Whether to print loading messages (default: True)
        """
        # Call parent class initializer
        super().__init__(poscar_file=poscar_file, verbose=verbose)
    
    def create_bilayer(self, stacking: str = 'AA', vacuum: float = 15.0, 
                      interlayer_distance: Optional[float] = None) -> 'POSCARProcessor':
        """
        Create bilayer structure from monolayer
        
        Args:
            stacking: 'AA' or 'AB' stacking
            vacuum: Vacuum layer thickness in Angstrom
            interlayer_distance: Distance between layers in Angstrom
            
        Returns:
            New POSCARProcessor instance with bilayer structure
        """
        if interlayer_distance is None:
            interlayer_distance = np.max([3.5, 2 * self.get_bond_length()])  # Default interlayer distance
            print(f"  interlayer_distance: {interlayer_distance:.4f} Å")
        
        # Create new instance
        bilayer = copy.deepcopy(self)
        
        # Calculate new c-vector length
        layer_thickness = self.get_layer_thickness()
        new_c_length = 2 * layer_thickness + interlayer_distance + vacuum
        
        # Update lattice vector
        bilayer.lattice[2, 2] = new_c_length
        
        # Create second layer coordinates
        second_layer_coords = self.coords.copy()
        
        if stacking == 'AA':
            # AA stacking: second layer directly above first layer
            second_layer_coords[:, 2] += layer_thickness + interlayer_distance
        elif stacking == 'AB':
            # AB stacking: second layer shifted by (a/2, b/2)
            second_layer_coords[:, 0] += self.lattice[0, 0] / 2
            second_layer_coords[:, 1] += self.lattice[1, 1] / 2
            second_layer_coords[:, 2] += layer_thickness + interlayer_distance
        else:
            raise ValueError("Stacking must be 'AA' or 'AB'")
        
        # Combine coordinates and atom types - ensure proper element preservation
        bilayer.coords = np.vstack([self.coords, second_layer_coords])
        
        # Explicitly preserve element types for each atom
        bilayer.atom_types = list(self.atom_types) + list(self.atom_types)
        
        # Update number of atoms - ensure elements list is preserved
        bilayer.n_atoms = [2 * n for n in self.n_atoms]
        bilayer.total_atoms = len(bilayer.coords)
        
        # Ensure elements list is properly preserved
        #bilayer.elements = self.elements.copy()
        
        # Update comment
        bilayer.comment = f"{self.comment} (bilayer {stacking} stacking)"
        
        
        return bilayer

    def center_structure(self) -> None:
        """
        Center the structure in the unit cell (only z-direction)
        """
        if len(self.coords) == 0:
            return
        
        # Only center in z-direction
        z_center = np.mean(self.coords[:, 2])
        self.coords[:, 2] -= z_center
        self.coords[:, 2] += self.lattice[2, 2] / 2
        
        print("Centered structure in z-direction only")
    
    def export_poscar(self, filename: str) -> None:
        """
        Export structure to POSCAR file with proper element ordering and layer structure
        
        Args:
            filename: Output filename
        """
        with open(filename, 'w') as f:
            # Write comment
            f.write(f"{self.comment}\n")
            
            # Write scaling factor
            f.write(f"{self.scale:20.16f}\n")
            
            # Write lattice vectors
            for i in range(3):
                f.write(f"{self.lattice[i,0]:20.16f}{self.lattice[i,1]:20.16f}{self.lattice[i,2]:20.16f}\n")
            
            # Write element names
            f.write(" ".join(self.elements) + "\n")
            
            # Write number of atoms
            f.write(" ".join(map(str, self.n_atoms)) + "\n")
            
            # Write coordinate type
            f.write(f"{self.coord_type}\n")
            
            # For bilayer structures, write atoms by element but preserve layer structure
            if "bilayer" in self.comment.lower():
                # Identify layers based on z-coordinates
                z_coords = self.coords[:, 2]
                z_center = np.mean(z_coords)
                bottom_indices = np.where(z_coords <= z_center)[0]
                top_indices = np.where(z_coords > z_center)[0]
                
                # For each element, write bottom layer atoms first, then top layer atoms
                for element in self.elements:
                    # Find all atoms of this element in bottom layer
                    bottom_atoms = [idx for idx in bottom_indices if self.atom_types[idx] == element]
                    # Find all atoms of this element in top layer
                    top_atoms = [idx for idx in top_indices if self.atom_types[idx] == element]
                    
                    # Write bottom layer atoms of this element
                    for idx in bottom_atoms:
                        coord = self.coords[idx]
                        f.write(f"{coord[0]:20.16f}{coord[1]:20.16f}{coord[2]:20.16f} {element}\n")
                    
                    # Write top layer atoms of this element
                    for idx in top_atoms:
                        coord = self.coords[idx]
                        f.write(f"{coord[0]:20.16f}{coord[1]:20.16f}{coord[2]:20.16f} {element}\n")
            else:
                # For non-bilayer structures, use standard order
                for i, coord in enumerate(self.coords):
                    f.write(f"{coord[0]:20.16f}{coord[1]:20.16f}{coord[2]:20.16f} {self.atom_types[i]}\n")
        
        print(f"Structure exported to {filename}")


    def visualize_clean_projections(self, 
                                figsize: tuple = (15, 5),
                                atom_size: float = 120,
                                save_path: str = None,
                                original_atom_count: int = None) -> None:
        """
        Create clean three-view projections (XY, XZ, YZ) with bonds and unit cell only
        
        Args:
            figsize: Figure size
            atom_size: Size of atoms in the plot
            save_path: Path to save the figure
            original_atom_count: Number of original atoms (for distinguishing original vs new atoms)
        """
        if len(self.coords) == 0:
            raise ValueError("No coordinates loaded!")
        
        fig, axes = plt.subplots(1, 3, figsize=figsize)
        
        # Use the merged projection function for all three views
        self._plot_projection_clean(axes[0], "Top View (XY)", 'xy', original_atom_count)
        self._plot_projection_clean(axes[1], "Side View (XZ)", 'xz', original_atom_count)
        self._plot_projection_clean(axes[2], "Front View (YZ)", 'yz', original_atom_count)
        
        # Overall title
        fig.suptitle(f"{self.comment}", fontsize=16, fontweight='bold')
        
        # Save figure if path provided
        if save_path:
            plt.savefig(save_path, dpi=200, bbox_inches='tight', facecolor='white')
            print(f"Clean projections figure saved to {save_path}")
        
        plt.tight_layout()
        plt.show()

    def _plot_projection_clean(self, ax, title, projection_type, original_atom_count=None):
        """
        Plot structure projection with atom distinction
        
        Args:
            ax: matplotlib axis object
            title: plot title
            projection_type: 'xy', 'xz', or 'yz' - determines which coordinates to plot
            original_atom_count: Number of original atoms (for distinguishing original vs new atoms)
        """
        # Define coordinate indices and labels based on projection type
        coord_mapping = {
            'xy': ([0, 1], ['X (Å)', 'Y (Å)'], self.lattice[:2, :2]),
            'xz': ([0, 2], ['X (Å)', 'Z (Å)'], np.array([[self.lattice[0, 0], 0], [0, self.lattice[2, 2]]])),
            'yz': ([1, 2], ['Y (Å)', 'Z (Å)'], np.array([[self.lattice[1, 1], 0], [0, self.lattice[2, 2]]]))
        }
        
        if projection_type not in coord_mapping:
            raise ValueError(f"Invalid projection_type: {projection_type}. Must be 'xy', 'xz', or 'yz'")
        
        coord_indices, axis_labels, lattice_matrix = coord_mapping[projection_type]
        
        # Plot atoms with distinction
        unique_elements = list(set(self.atom_types))
        for element in unique_elements:
            element_indices = [i for i, atom_type in enumerate(self.atom_types) if atom_type == element]
            element_coords = self.coords[element_indices]
            
            color = element_colors.get(element, '#808080')
            
            # Determine if atoms are original or new based on index
            if original_atom_count is not None:
                original_indices = [i for i in element_indices if i < original_atom_count]
                new_indices = [i for i in element_indices if i >= original_atom_count]
                
                # Plot original atoms with solid edges
                if original_indices:
                    original_coords = self.coords[original_indices]
                    ax.scatter(original_coords[:, coord_indices[0]], original_coords[:, coord_indices[1]],
                              c=color, s=120, label=f'{element} (original)', alpha=0.9, 
                              edgecolors='black', linewidth=2.0, linestyle='-')
                
                # Plot new atoms with dashed edges
                if new_indices:
                    new_coords = self.coords[new_indices]
                    ax.scatter(new_coords[:, coord_indices[0]], new_coords[:, coord_indices[1]],
                              c=color, s=120, label=f'{element} (new)', alpha=0.9, 
                              edgecolors='black', linewidth=2.0, linestyle='--')
            else:
                # If no distinction needed, plot all atoms with solid edges
                ax.scatter(element_coords[:, coord_indices[0]], element_coords[:, coord_indices[1]],
                          c=color, s=120, label=element, alpha=0.9, 
                          edgecolors='black', linewidth=2.0, linestyle='-')
        
        # Plot unit cell projection (simplified - just the box)
        vertices = np.array([
            [0, 0], [1, 0], [1, 1], [0, 1]
        ])
        cart_vertices = np.dot(vertices, lattice_matrix)
        
        # Plot unit cell edges
        edges = [[0, 1], [1, 2], [2, 3], [3, 0]]
        for edge in edges:
            points = cart_vertices[edge]
            ax.plot(points[:, 0], points[:, 1], 'k-', linewidth=2, alpha=0.6)
        
        ax.set_xlabel(axis_labels[0], fontsize=12)
        ax.set_ylabel(axis_labels[1], fontsize=12)
        ax.set_title(title, fontsize=12, fontweight='bold')
        
        # Only show legend for XY projection to avoid clutter
        if projection_type == 'xy':
            ax.legend(loc=1, bbox_to_anchor=(1.0, -0.4), frameon=False)
        
        ax.set_aspect('equal')
        
        # Remove borders and grid
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.set_xticks([])
        ax.set_yticks([])


    def export_vasp(self, output_dir: str, poscar_filename: str = "POSCAR", 
                            potcar_dir: str = None, incar_file: str = None,
                            submit_file: str = None, platform: str = None) -> None:
        """
        Export complete VASP input files including POSCAR, POTCAR, INCAR, and job submission
        
        Args:
            output_dir: Output directory for VASP files
            poscar_filename: Name for POSCAR file (default: "POSCAR")
            potcar_dir: Directory containing VASP pseudopotential files
            incar_file: Path to INCAR file to copy (optional)
            submit_file: Path to job submission file (optional)
            platform: Platform name for job submission (e.g., 'bh', 'slurm', 'pbs')
        """
        
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        
        # Get the folder name (file_name) from output_dir
        file_name = os.path.basename(output_dir.rstrip('/'))
        
        # Export POSCAR file
        poscar_path = os.path.join(output_dir, poscar_filename)
        self.export_poscar(poscar_path)
        print(f"POSCAR exported to: {poscar_path}")
        
        # Handle POTCAR files if potcar_dir is provided
        if potcar_dir:
            if not os.path.exists(potcar_dir):
                raise FileNotFoundError(f"POTCAR directory {potcar_dir} not found!")
            
            potcar_path = os.path.join(output_dir, "POTCAR")
            self._create_potcar(potcar_dir, potcar_path)
            print(f"POTCAR created at: {potcar_path}")
        
        # Handle INCAR file if provided
        if incar_file:
            if not os.path.exists(incar_file):
                raise FileNotFoundError(f"INCAR file {incar_file} not found!")
            
            incar_dest = os.path.join(output_dir, "INCAR")
            shutil.copy2(incar_file, incar_dest)
            print(f"INCAR copied to: {incar_dest}")
        
        # Handle job submission file if provided
        if submit_file:
            if not os.path.exists(submit_file):
                raise FileNotFoundError(f"Submit file {submit_file} not found!")
            
            submit_dest = os.path.join(output_dir, os.path.basename(submit_file))
            
            # If it's a JSON file (for BH platform), modify the job_name
            if submit_file.endswith('.json'):
                self._copy_and_modify_job_json(submit_file, submit_dest, file_name)
            else:
                # For other platforms, just copy the file
                shutil.copy2(submit_file, submit_dest)
            
            print(f"Submit file copied to: {submit_dest}")
            
            # Handle platform-specific job submission
            if platform:
                self._submit_job(output_dir, submit_dest, platform)
        
        print(f"VASP input files successfully exported to: {output_dir}")

    def _create_potcar(self, potcar_dir: str, output_potcar: str) -> None:
        """
        Create POTCAR file by concatenating individual element POTCAR files
        with enhanced variant search for elements like Ba_sv, Ba_h, etc.
        
        Args:
            potcar_dir: Directory containing VASP pseudopotential files (folders named by elements)
            output_potcar: Output path for the combined POTCAR file
        """
        
        found_potcars = []
        missing_elements = []
        temp_potcar_files = []
        
        # Common VASP pseudopotential variants
        common_variants = ['_sv', '_h', '_pv', '_d', '_f', '_p', '_s']
        
        # Search for POTCAR files for each element
        for element in self.elements:
            element_found = False
            potcar_file = None
            
            # First try the exact element name
            element_dir = os.path.join(potcar_dir, element)
            exact_potcar = os.path.join(element_dir, "POTCAR")
            
            if os.path.exists(exact_potcar):
                potcar_file = exact_potcar
                element_found = True
                print(f"Found POTCAR for {element}: {exact_potcar}")
            else:
                # Try variants if exact match not found
                print(f"Exact match not found for {element}, searching variants...")
                
                for variant in common_variants:
                    variant_dir = os.path.join(potcar_dir, f"{element}{variant}")
                    variant_potcar = os.path.join(variant_dir, "POTCAR")
                    
                    if os.path.exists(variant_potcar):
                        potcar_file = variant_potcar
                        element_found = True
                        print(f"Found POTCAR for {element} (variant {variant}): {variant_potcar}")
                        break
                
                # If still not found, try to find any directory starting with the element name
                if not element_found:
                    pattern = os.path.join(potcar_dir, f"{element}*")
                    matching_dirs = glob.glob(pattern)
                    
                    for dir_path in matching_dirs:
                        if os.path.isdir(dir_path):
                            potential_potcar = os.path.join(dir_path, "POTCAR")
                            if os.path.exists(potential_potcar):
                                potcar_file = potential_potcar
                                element_found = True
                                variant_name = os.path.basename(dir_path)
                                print(f"Found POTCAR for {element} (variant {variant_name}): {potential_potcar}")
                                break
            
            if element_found and potcar_file:
                # Copy POTCAR to temporary file with element suffix
                temp_potcar = os.path.join(os.path.dirname(output_potcar), f"POTCAR.{element}")
                shutil.copy2(potcar_file, temp_potcar)
                found_potcars.append((element, temp_potcar))
                temp_potcar_files.append(temp_potcar)
            else:
                missing_elements.append(element)
                print(f"Warning: POTCAR not found for element {element} in any variant")
        
        # Check if all elements were found
        if missing_elements:
            print(f"Warning: POTCAR files not found for elements: {missing_elements}")
            print("Available element directories in POTCAR directory:")
            for item in os.listdir(potcar_dir):
                item_path = os.path.join(potcar_dir, item)
                if os.path.isdir(item_path):
                    potcar_path = os.path.join(item_path, "POTCAR")
                    status = "✓" if os.path.exists(potcar_path) else "✗"
                    print(f"  {status} {item}/POTCAR")
        
        # Create combined POTCAR file by concatenating in element order
        with open(output_potcar, 'wb') as outfile:
            for element, temp_potcar in found_potcars:
                try:
                    with open(temp_potcar, 'rb') as infile:
                        shutil.copyfileobj(infile, outfile)
                    print(f"Added POTCAR for {element} to combined file")
                except Exception as e:
                    print(f"Error reading POTCAR for {element} ({temp_potcar}): {e}")
        
        # Clean up temporary POTCAR files
        for temp_file in temp_potcar_files:
            try:
                os.remove(temp_file)
            except Exception as e:
                print(f"Warning: Could not remove temporary file {temp_file}: {e}")
        
        if found_potcars:
            print(f"Combined POTCAR created with {len(found_potcars)} elements in order: {[elem for elem, _ in found_potcars]}")
        else:
            print("Warning: No POTCAR files were found and combined")

    def _copy_and_modify_job_json(self, source_json: str, dest_json: str, job_name: str) -> None:
        """
        Copy JSON job file and modify the job_name field
        
        Args:
            source_json: Source JSON file path
            dest_json: Destination JSON file path
            job_name: New job name to set
        """
        
        try:
            # Read the original JSON file
            with open(source_json, 'r') as f:
                job_config = json.load(f)
            
            # Modify the job_name
            job_config['job_name'] = job_name
            
            # Write the modified JSON file
            with open(dest_json, 'w') as f:
                json.dump(job_config, f, indent=2)
            
            print(f"Job JSON modified: job_name set to '{job_name}'")
            
        except json.JSONDecodeError as e:
            print(f"Error parsing JSON file {source_json}: {e}")
            # Fallback: just copy the file without modification
            shutil.copy2(source_json, dest_json)
        except Exception as e:
            print(f"Error modifying JSON file: {e}")
            # Fallback: just copy the file without modification
            shutil.copy2(source_json, dest_json)

    def _submit_job(self, output_dir: str, submit_file: str, platform: str) -> None:
        """
        Submit job based on platform
        
        Args:
            output_dir: Output directory containing job files
            submit_file: Path to job submission file
            platform: Platform name for job submission
        """
        
        submit_filename = os.path.basename(submit_file)
        submit_path = os.path.join(output_dir, submit_filename)
        
        if platform.lower() == 'bh':
            # BH platform: use lbg job submit command
            if submit_filename.endswith('.json'):
                try:
                    # Change to output directory and submit job
                    original_cwd = os.getcwd()
                    os.chdir(output_dir)
                    
                    # Submit job using lbg command
                    cmd = ['lbg', 'job', 'submit', '-i', submit_filename, '-p', './','-r', './']
                    print(f"Submitting job with command: {' '.join(cmd)}")
                    
                    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
                    print(f"Job submitted successfully!")
                    print(f"Output: {result.stdout}")
                    
                    if result.stderr:
                        print(f"Warning: {result.stderr}")
                        
                except subprocess.CalledProcessError as e:
                    print(f"Error submitting job: {e}")
                    print(f"Error output: {e.stderr}")
                except FileNotFoundError:
                    print("Error: 'lbg' command not found. Please ensure lbg is installed and in PATH.")
                finally:
                    os.chdir(original_cwd)
            else:
                print(f"Warning: BH platform requires .json submit file, got {submit_filename}")
        
        elif platform.lower() == 'slurm':
            # SLURM platform: use sbatch command
            try:
                original_cwd = os.getcwd()
                os.chdir(output_dir)
                
                cmd = ['sbatch', submit_filename]
                print(f"Submitting SLURM job with command: {' '.join(cmd)}")
                
                result = subprocess.run(cmd, capture_output=True, text=True, check=True)
                print(f"SLURM job submitted successfully!")
                print(f"Output: {result.stdout}")
                
            except subprocess.CalledProcessError as e:
                print(f"Error submitting SLURM job: {e}")
                print(f"Error output: {e.stderr}")
            except FileNotFoundError:
                print("Error: 'sbatch' command not found. Please ensure SLURM is installed.")
            finally:
                os.chdir(original_cwd)
        
        elif platform.lower() == 'pbs':
            # PBS platform: use qsub command
            try:
                original_cwd = os.getcwd()
                os.chdir(output_dir)
                
                cmd = ['qsub', submit_filename]
                print(f"Submitting PBS job with command: {' '.join(cmd)}")
                
                result = subprocess.run(cmd, capture_output=True, text=True, check=True)
                print(f"PBS job submitted successfully!")
                print(f"Output: {result.stdout}")
                
            except subprocess.CalledProcessError as e:
                print(f"Error submitting PBS job: {e}")
                print(f"Error output: {e.stderr}")
            except FileNotFoundError:
                print("Error: 'qsub' command not found. Please ensure PBS is installed.")
            finally:
                os.chdir(original_cwd)
        
        else:
            print(f"Warning: Unsupported platform '{platform}'. Supported platforms: bh, slurm, pbs")


import os
from tqdm import tqdm

import numpy as np
import matplotlib.pyplot as plt

from extrempy.constant import *


class SEDCalc:
    def __init__(self, position: np.ndarray, velocity: np.ndarray, type_array: np.ndarray,
                 cells: np.ndarray, dt: float = 1.0, SKIP: int = 0, INTERVAL: int = 1):
        """Initialize SEDCalc for spectral energy density calculation.

        Args:
            position (np.ndarray): Position array with shape (nframes, natoms, 3)
            velocity (np.ndarray): Velocity array with shape (nframes, natoms, 3) 
            type_array (np.ndarray): Atom type array with shape (natoms,)
            cells (np.ndarray): Cell vectors array with shape (3, 3)
            dt (float, optional): Time step in fs. Defaults to 1.0.
            SKIP (int, optional): Number of frames to skip. Defaults to 0.
            INTERVAL (int, optional): Interval between frames. Defaults to 1.
        """
        # Store input parameters
        self.position = position
        self.velocity = velocity
        self.type = type_array
        self.cells = cells
        self.dt = dt
        self.SKIP = SKIP
        self.INTERVAL = INTERVAL

        # Calculate basic parameters
        self.numb_frames = position.shape[0]
        self.numb_atoms = position.shape[1]
        self.type_list = np.unique(type_array)
        self.numb_types = len(self.type_list)

        # Calculate frames to use
        self.dump_idx = np.arange(self.SKIP, self.numb_frames, self.INTERVAL)
        self.numb_frames_used = self.dump_idx.shape[0]

        # Calculate delta k
        self.delta_k = 2 * np.pi / np.diag(self.cells)

    def _init_jkt(self):
        """Initialize arrays for current-current correlation function."""
        if self.numb_types == 1:
            self.jkt_re = np.zeros((self.numb_frames_used, 3))
            self.jkt_im = np.zeros((self.numb_frames_used, 3))
        else:
            self.jkt_re = np.zeros(
                (self.numb_frames_used, 3 * (self.numb_types + 1)))
            self.jkt_im = np.zeros(
                (self.numb_frames_used, 3 * (self.numb_types + 1)))

    def _calc_theta(self, k_vec: np.ndarray):
        """Calculate phase factor theta = k·r."""
        self.theta = np.dot(self.position[self.current_frame], k_vec)

    def _calc_jkt(self, tdx: int):
        """Calculate current-current correlation function j(k,t)."""
        current_vel = self.velocity[self.current_frame]

        # Calculate for all atoms
        for idx in range(3):
            self.jkt_re[tdx, idx] = np.sum(
                current_vel[:, idx] * np.cos(self.theta))
            self.jkt_im[tdx,
                        idx] = np.sum(-current_vel[:, idx] * np.sin(self.theta))

        # Calculate for each element type if multiple types exist
        if self.numb_types > 1:
            for edx in range(self.numb_types):
                cri = self.type == edx + 1
                for idx in range(3):
                    self.jkt_re[tdx, idx + 3 * (edx + 1)] = np.sum(
                        current_vel[:, idx][cri] * np.cos(self.theta[cri]))
                    self.jkt_im[tdx, idx + 3 * (edx + 1)] = np.sum(
                        -current_vel[:, idx][cri] * np.sin(self.theta[cri]))

    def _calc_ckw(self, dt: float, jdx: int):
        """Calculate spectral energy density C(k,ω)."""
        freq = np.fft.fftfreq(self.numb_frames_used, d=dt)

        tmp = np.transpose(
            np.vstack((self.jkt_re[:, jdx], self.jkt_im[:, jdx])))
        j_kt = np.zeros(self.numb_frames_used, dtype=np.complex128)
        j_kt.real = tmp[:, 0]
        j_kt.imag = tmp[:, 1]

        j_kw = np.fft.fft(j_kt)
        j_kw_conj = np.fft.fft(j_kt.conjugate())

        ckw = np.abs((j_kw + j_kw_conj) / 2)**2 / \
            (self.numb_frames_used * self.numb_atoms)

        out = np.transpose(np.vstack((freq, ckw)))
        return out[np.lexsort(out[:, ::-1].T)]

    def calc_sed(self, k_vec_tmp: np.ndarray, nk: int = 1, save_dir: str = None, suffix: str = None, tqdm: bool = False):
        """Calculate spectral energy density for given k-vector.

        Args:
            k_vec_tmp (np.ndarray): k-vector direction
            nk (int, optional): k-vector multiplier. Defaults to 1.
            save_dir (str, optional): Directory to save results. Defaults to None.
            suffix (str, optional): Suffix for output files. Defaults to None.
            tqdm (bool, optional): Whether to use tqdm. Defaults to False.

        Returns:
            np.ndarray: Combined SED for all directions
        """
        self._init_jkt()

        unit_k_vec = k_vec_tmp / np.linalg.norm(k_vec_tmp)
        if suffix is None:
            suffix = '%.3f%.3f%.3f' % (
                unit_k_vec[0], unit_k_vec[1], unit_k_vec[2])

        k_vec = self.delta_k * unit_k_vec * nk

        if tqdm:
            tqdm_bar = tqdm(range(self.numb_frames_used),
                            desc=f'Calculating SED for k={nk}')
        else:
            tqdm_bar = range(self.numb_frames_used)

        for tdx in tqdm_bar:
            self.current_frame = self.dump_idx[tdx]
            self._calc_theta(k_vec)
            self._calc_jkt(tdx)

        labels = ['vx', 'vy', 'vz']
        total_spectrum = np.zeros((self.numb_frames_used, 2))

        for jdx in range(3):
            ckw = self._calc_ckw(self.dt * self.INTERVAL, jdx)

            if jdx == 0:
                total_spectrum = np.array(ckw)
            else:
                total_spectrum[:, 1] += ckw[:, 1]

            if save_dir:
                os.makedirs(save_dir, exist_ok=True)
                np.savetxt(
                    os.path.join(save_dir, f'ckw_{labels[jdx]}_{nk}k_{suffix}.txt'),
                    ckw,
                    header=f'component-{labels[jdx]} nk={nk} frequency (1/fs) spectrum (arb. unit.)',
                    comments='# '
                )

        if save_dir:
            np.savetxt(
                os.path.join(save_dir, f'ckw_tot_{nk}k_{suffix}.txt'),
                total_spectrum,
                header=f'nk={nk} frequency (1/fs) spectrum (arb. unit.)',
                comments='# '
            )

        return total_spectrum

    def plot_sed_1d(self, data: np.ndarray, k_value: float, ax: plt.Axes = None) -> tuple[plt.Figure, plt.Axes]:
        """Plot 1D spectral energy density.

        Args:
            data (np.ndarray): SED data array with shape (nfreq, 2)
            k_value (float): k-vector value
            ax (plt.Axes, optional): Matplotlib axes to plot on. Defaults to None.

        Returns:
            tuple[plt.Figure, plt.Axes]: Figure and axes objects
        """
        if ax is None:
            fig, ax = plt.subplots(figsize=(6, 4), dpi=200)
        else:
            fig = ax.figure

        ax.plot(data[:, 0], data[:, 1], '-o', ms=2,
                label=f'k={k_value}' + r'$\rm{\mathring A^{-1}}$')

        ax.set_xlim(0, )
        ax.set_xlabel(r'$\omega$ ($\rm{fs^{-1}}$)')
        ax.set_ylabel(r'$C(k,\omega)$')
        ax.legend(fontsize=6)
        plt.tight_layout()
        plt.show()

        return fig, ax

    def generate_sed_2d(self, save_dir: str, k_vec_tmp: np.ndarray, nk_range: tuple[int, int], suffix: str = None) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Generate 2D SED data by calculating multiple k-vectors.

        Args:
            save_dir (str): Directory to save results
            k_vec_tmp (np.ndarray): k-vector direction
            nk_range (tuple[int, int]): Range of nk values (start, end)
            suffix (str, optional): Suffix for output files. Defaults to None.

        Returns:
            tuple[np.ndarray, np.ndarray, np.ndarray]: K-mesh, Omega-mesh and SED data
        """
        print("\n=== Starting 2D SED Generation ===")

        k_list = []
        data_list = []
        freq = None

        unit_k_vec = k_vec_tmp / np.linalg.norm(k_vec_tmp)
        progress_bar = tqdm(range(*nk_range))

        for nk in progress_bar:
            try:
                progress_bar.set_description(f"Calculating 2D SED for k={nk}")
                k_value = np.linalg.norm(self.delta_k * unit_k_vec * nk)

                # Calculate SED for this k-vector
                spectrum = self.calc_sed(k_vec_tmp=k_vec_tmp, nk=nk,
                                         save_dir=save_dir, suffix=suffix, tqdm=False)

                if freq is None:
                    freq = spectrum[:, 0]

                k_list.append(k_value)
                data_list.append(spectrum[:, 1])

            except Exception as e:
                print(f"Error processing nk={nk}: {str(e)}")
                continue

        if not data_list:
            raise ValueError(
                "No data was collected. Check if the SED calculations were successful.")

        # Convert to numpy arrays
        k_list = np.array(k_list)
        data_array = np.column_stack(data_list)

        # Create meshgrid
        K, Omega = np.meshgrid(k_list, freq)

        # Save data
        if save_dir:
            np.savetxt(os.path.join(
                save_dir, f'sed_2d_{suffix}.dat'), data_array)
            np.savetxt(os.path.join(save_dir, f'knum_2d_{suffix}.dat'), K)
            np.savetxt(os.path.join(save_dir, f'freq_2d_{suffix}.dat'), Omega)
            print(f"\nSaved 2D SED data to {save_dir}")

        print("=== 2D SED Generation Complete ===\n")
        return K, Omega, data_array

    def plot_sed_2d(self, K: np.ndarray, Omega: np.ndarray, Z: np.ndarray,
                    k_max: float = None, ax: plt.Axes = None) -> tuple[plt.Figure, plt.Axes]:
        """Plot 2D spectral energy density.

        Args:
            K (np.ndarray): K-mesh grid
            Omega (np.ndarray): Omega-mesh grid
            Z (np.ndarray): SED data array
            k_max (float, optional): Maximum k value to plot. Defaults to None.
            ax (plt.Axes, optional): Matplotlib axes to plot on. Defaults to None.

        Returns:
            tuple[plt.Figure, plt.Axes]: Figure and axes objects
        """
        if ax is None:
            fig, ax = plt.subplots(figsize=(6, 4), dpi=200)
        else:
            fig = ax.figure

        contour = ax.contourf(K, Omega, Z, levels=100, cmap='Blues')

        ax.set_ylim(0, )
        if k_max is not None:
            ax.set_xlim(0, k_max)
        else:
            ax.set_xlim(0, )

        ax.set_xlabel(r'$k$ ($\rm{\mathring A^{-1}}$)')
        ax.set_ylabel(r'$\omega$ ($\rm{fs^{-1}}$)')

        plt.colorbar(contour, ax=ax, label=r'$C(k,\omega)$')
        plt.tight_layout()
        plt.show()

        return fig, ax

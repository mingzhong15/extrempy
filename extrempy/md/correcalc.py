from tqdm import tqdm
from typing import Callable

import numpy as np
import matplotlib.pyplot as plt

try
    from extrempy.constant import *
except  
    from ..constant import *


class TimeCorrelationCalc:
    def __init__(self, position: np.ndarray, velocity: np.ndarray, dt: float = 1.0, time_corre_max: int = 0):
        """Initialize TimeCorrelationCalc for time correlation function calculations.

        Args:
            position (np.ndarray): Position array with shape (nframes, natoms, 3)
            velocity (np.ndarray): Velocity array with shape (nframes, natoms, 3)
            dt (float, optional): Time step in fs. Defaults to 1.0.
            time_corre_max (int, optional): Maximum correlation time. If 0, uses nframes/2. Defaults to 0.
        """
        # Store input parameters
        self.position = position
        self.velocity = velocity
        self.dt = dt

        # Calculate basic parameters
        self.numb_frames = position.shape[0]
        self.numb_atoms = position.shape[1]

        # Set up correlation time parameters
        if time_corre_max == 0:
            self.time_corre_max = int(self.numb_frames / 2)
        else:
            self.time_corre_max = int(time_corre_max / self.dt)

        self.time = np.arange(self.time_corre_max) * self.dt

        strgs = '='*70 + '\n'
        strgs += ' Time Correlation Calculator is initialized \n'
        strgs += '='*70 + '\n'
        print(strgs)

    def _calc_autocorrelation(self, quantity_func: Callable, skip: int = 0, interval: int = 0, normalize: bool = True):
        """Calculate the time correlation function for a given quantity.

        Args:
            quantity_func (callable): Function that returns quantity array for given frame
            skip (int, optional): Number of initial frames to skip. Defaults to 0.
            interval (int, optional): Interval between reference frames. Defaults to 0.
            normalize (bool, optional): Whether to normalize correlation. Defaults to True.

        Returns:
            np.ndarray: Time correlation function
        """
        if interval == 0:
            skip = 0
            interval = self.time_corre_max

        correlation = np.zeros(self.time_corre_max)
        count = 0

        print(f"Starting autocorrelation calculation...")
        print(f"Total frames to process: {(self.time_corre_max - skip) // interval}")

        for init_idx in tqdm(range(skip, self.time_corre_max, interval),
                             desc="Calculating correlation"):
            q0 = quantity_func(init_idx)  # Get initial quantity

            for tau_idx in range(self.time_corre_max):
                if init_idx + tau_idx >= self.numb_frames:
                    break

                # Get quantity at time t
                qt = quantity_func(init_idx + tau_idx)
                correlation[tau_idx] += np.sum(qt * q0)

            count += 1

        correlation /= count
        print("Correlation calculation completed.")

        if normalize:
            correlation /= correlation[0]  # Normalize by initial value
            print("Correlation normalized.")

        return correlation

    def calc_vacf(self, skip: int = 0, interval: int = 0):
        """Calculate velocity autocorrelation function (VACF).

        Args:
            skip (int, optional): Frames to skip. Defaults to 0.
            interval (int, optional): Frame interval. Defaults to 0.

        Returns:
            np.ndarray: VACF array
        """
        def velocity_func(idx):
            return self.velocity[idx]

        vacf = self._calc_autocorrelation(
            velocity_func, skip, interval, normalize=True)
        return vacf

    def calc_pacf(self, skip: int = 0, interval: int = 0):
        """Calculate position autocorrelation function (PACF).

        Args:
            skip (int, optional): Frames to skip. Defaults to 0.
            interval (int, optional): Frame interval. Defaults to 0.

        Returns:
            np.ndarray: PACF array
        """
        def position_func(idx):
            return self.position[idx]

        pacf = self._calc_autocorrelation(
            position_func, skip, interval, normalize=True)
        return pacf

    def calc_power_spectrum(self, func_correlation: np.ndarray):
        """Calculate power spectrum from correlation function.

        Args:
            func_correlation (np.ndarray): Correlation function array

        Returns:
            tuple: (frequencies, spectrum)
        """
        freq = np.fft.fftfreq(self.time.shape[0], d=self.dt)
        spectra = np.abs(np.fft.fft(func_correlation))

        # Only return positive frequencies
        positive_freq_mask = freq >= 0
        return freq[positive_freq_mask], spectra[positive_freq_mask]

    def plot_time_correlation(self, results: dict, save_dir: str = None):
        """Calculate and plot time correlation functions and VDOS.

        Args:
            results (dict): Dictionary containing:
                - 'time': Time array
                - 'vacf': Velocity autocorrelation function
                - 'pacf': Position autocorrelation function if calc_pacf is True
                - 'freq': Frequency array for VDOS
                - 'vdos': Vibrational density of states
            save_dir (str, optional): Directory to save plots. Defaults to None.
        """

        # Create figure with subplots
        if 'pacf' in results.keys():
            fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))
        else:
            fig, (ax1, ax3) = plt.subplots(1, 2, figsize=(10, 5))

        # Plot VACF
        ax1.plot(results['time'], results['vacf'], '-')
        ax1.set_xlabel('Time (fs)')
        ax1.set_ylabel('VACF')
        ax1.set_title('Velocity Autocorrelation')

        # Plot PACF
        if 'pacf' in results.keys():
            ax2.plot(results['time'], results['pacf'], '-')
            ax2.set_xlabel('Time (fs)')
            ax2.set_ylabel('PACF')
            ax2.set_title('Position Autocorrelation')

        # Plot VDOS
        ax3.scatter(results['freq'], results['vdos'], s=2)
        ax3.set_xlabel('Frequency (THz)')
        ax3.set_ylabel('VDOS')
        ax3.set_title('Vibrational Density of States')

        plt.tight_layout()

        # Save plots if directory is provided
        if save_dir:
            os.makedirs(save_dir, exist_ok=True)
            if 'pacf' in results.keys():
                np.savetxt(os.path.join(save_dir, 'time_correlation.txt'),
                           np.column_stack(
                               [results['time'], results['vacf'], results['pacf']]),
                           header='Time(fs) VACF PACF')
            else:
                np.savetxt(os.path.join(save_dir, 'time_correlation.txt'),
                           np.column_stack(
                               [results['time'], results['vacf']]),
                           header='Time(fs) VACF')
            np.savetxt(os.path.join(save_dir, 'vdos.txt'),
                       np.column_stack([results['freq'], results['vdos']]),
                       header='Frequency(THz) VDOS')
            fig.savefig(os.path.join(save_dir, 'time_correlation.png'),
                        dpi=300, bbox_inches='tight')
        plt.show()
        if 'pacf' in results.keys():
            return [fig, ax1, ax2, ax3]
        else:
            return [fig, ax1, ax3]

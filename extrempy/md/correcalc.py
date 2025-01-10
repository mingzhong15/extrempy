from extrempy.constant import *
from extrempy.md.base import MDSys

class TimeCorrelationCalc(MDSys):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.time_corre_max = int(self.numb_frames / 2)
        self.time = np.arange(self.time_corre_max) * self.dump_freq * self.dt

        strgs = '# ================================= # \n'
        strgs += ' Time Correlation Calculator is initialized \n'
        strgs += '# ================================= # \n'
        print(strgs)

    def _calc_correlation(self, quantity_func, skip=0, interval=0, normalize=True):
        """
        Calculate the time correlation function for a given quantity.

        Parameters:
        -----------
        quantity_func : function
            A function that returns the quantity of interest (e.g., velocity, position) for a given frame index.
        skip : int, optional
            The number of initial frames to skip. Default is 0.
        interval : int, optional
            The interval between reference frames. Default is 0, which means only using initial frame.
        normalize : bool, optional
            Whether to normalize the correlation function. Default is True.

        Returns:
        --------
        correlation : np.ndarray
            The computed time correlation function.
        """
        if interval == 0:
            skip = 0
            interval = self.time_corre_max
        else:
            skip = skip
            interval = interval

        correlation = np.zeros(self.time_corre_max)
        count = 0

        for init_idx in range(skip, self.time_corre_max, interval):
            self._read_dump(idx=init_idx * self.dump_freq)
            q0 = quantity_func()  # Get the initial quantity

            for tau_idx in range(self.time_corre_max):
                self._read_dump(idx=(init_idx + tau_idx) * self.dump_freq)
                qt = quantity_func()  # Get the quantity at time t

                correlation[tau_idx] += np.sum(qt * q0)  # Compute the correlation

            count += 1

        correlation /= count

        if normalize:
            correlation /= np.sum(q0 * q0)  # Normalize by the initial value

        return correlation

    def calc_vacf(self, skip=0, interval=0):
        """
        Calculate the velocity autocorrelation function (VACF).

        Parameters:
        -----------
        skip : int, optional
            The number of initial frames to skip. Default is 0.
        interval : int, optional
            The interval between reference frames. Default is 0, which means only using initial frame.

        Returns:
        --------
        vacf : np.ndarray
            The computed velocity autocorrelation function.
        """
        def velocity_func():
            return np.array(self.velocity)

        vacf = self._calc_correlation(velocity_func, skip, interval, normalize=True)
        return vacf

    def calc_pacf(self, skip=0, interval=0):
        """
        Calculate the position autocorrelation function (PACF).

        Parameters:
        -----------
        skip : int, optional
            The number of initial frames to skip. Default is 0.
        interval : int, optional
            The interval between reference frames. Default is 0, which means using all frames.

        Returns:
        --------
        pacf : np.ndarray
            The computed position autocorrelation function.
        """
        def position_func():
            return np.array(self.position)

        pacf = self._calc_correlation(position_func, skip, interval, normalize=True)
        return pacf

    # Add more correlation functions as needed...
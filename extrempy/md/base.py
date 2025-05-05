import re
from tqdm import tqdm

import polars as pl
import numpy as np

try:
    from extrempy.constant import *
    from extrempy.md.sedcalc import SEDCalc
    from extrempy.md.correcalc import TimeCorrelationCalc
except:
    from ..constant import *
    from .sedcalc import SEDCalc
    from .correcalc import TimeCorrelationCalc


class System():
    """The molecular dynamic postprocessing system (for single trajectory file) of Extrempy.
    """

    def __init__(self, file_name: str, dt: float = 1, type_name: list[str] = None, fmt: str = None) -> None:
        """Initialize the System class.

        Args:
            file_name (str): The name of the file to read.
            dt (float, optional): The time step. Defaults to 1.
            type_name (list[str], optional): The type name. Defaults to None.
            fmt (str, optional): The format of the file. Defaults to None.
        """
        self.__file_name = file_name
        self.__timestep = dt
        self.__type_name = type_name
        self.__fmt = fmt

        # read the file
        self.__data, self.__box, self.__boundary = self._read_file_by_format()

        # update the position and velocity
        self.update_pos()
        if "vx" in self.__data.columns:
            self.update_vel()

        # If the type name is provided, add the type name to self.data
        if self.__type_name is not None:
            assert self.__data["type"].max() <= len(self.__type_name)

            type_dict = {str(i): self.__type_name[i - 1]
                         for i in self.__data["type"].unique()}

            self.__data = self.__data.with_columns(
                pl.col("type")
                .cast(pl.Utf8)
                .replace_strict(type_dict)
                .alias("type_name")
            )
        # check the 'id' column
        assert "id" in self.__data.columns
        # sort the data by 'id'
        self.__data = self.__data.sort("id")

    @property
    def fmt(self) -> str:
        """The format of the file.

        Returns:
            str: The format of the file.
        """
        if self.__fmt is None:
            return self._get_format_()
        else:
            return self.__fmt

    @property
    def file_name(self) -> str:
        """The name of the file.

        Returns:
            str: The name of the file.
        """
        return self.__file_name

    @property
    def timestep(self) -> float:
        """The time step.

        Returns:
            float: The time step.
        """
        return self.__timestep

    @property
    def type_name(self) -> list[str]:
        """The type name.

        Returns:
            list[str]: The type name.
        """
        return self.__type_name

    @property
    def data(self) -> pl.DataFrame:
        """All information contained in the trajectory file.

        Returns:
            pl.DataFrame: The datasets.
        """
        return self.__data

    @property
    def box(self) -> np.ndarray:
        """The simulation box.

        Returns:
            np.ndarray: The simulation box.
        """
        return self.__box

    @property
    def cells(self) -> np.ndarray:
        """The simulation cells.

        Returns:
            np.ndarray: The simulation box.
        """
        return self.__box[:-1]

    @property
    def boundary(self) -> np.ndarray:
        """The simulation boundary.

        Returns:
            np.ndarray: The boundary.
        """
        return self.__boundary

    def update_pos(self):
        """Call it only if you modify the positions information by modify the data."""
        assert "xu" in self.__data.columns, "Must contains the position information."
        self.__pos = np.c_[self.__data["xu"],
                           self.__data["yu"], self.__data["zu"]]
        self.__pos.flags.writeable = False

        # Apply periodic boundary conditions if 'x' column is not present
        if "x" not in self.__data.columns:
            self.periodic_boundary_condition()

    def periodic_boundary_condition(self):
        """Apply periodic boundary conditions to the positions, and the 'x', 'y', 'z' will be added to self.data"""
        box_lengths = np.array([
            np.linalg.norm(self.__box[0]),  # x length
            np.linalg.norm(self.__box[1]),  # y length
            np.linalg.norm(self.__box[2])   # z length
        ])

        # Apply PBC and store as new columns
        pbc_pos = self.__pos - self.__box[3]  # Subtract origin
        pbc_pos = pbc_pos - np.floor(pbc_pos / box_lengths) * box_lengths

        # Add the new columns to self.data
        self.__data = self.__data.with_columns([
            pl.Series("x", pbc_pos[:, 0]),
            pl.Series("y", pbc_pos[:, 1]),
            pl.Series("z", pbc_pos[:, 2])
        ])

    def update_vel(self):
        """Call it only if you modify the velocities information by modify the data."""
        assert "vx" in self.__data.columns, "Must contains the velocity information."
        self.__vel = np.c_[self.__data["vx"],
                           self.__data["vy"], self.__data["vz"]]
        self.__vel.flags.writeable = False

    @property
    def pos(self):
        """particle position information. Do not change it directly.
        If you want to modify the positions, modify the data and call self.update_pos()

        Returns:
            np.ndarray: position information.
        """

        return self.__pos

    @property
    def vel(self):
        """particle velocity information. Do not change it directly.
        If you want to modify the velocities, modify the data and call self.update_vel()

        Returns:
            np.ndarray: velocity information.
        """
        if "vx" in self.__data.columns:
            return self.__vel
        else:
            raise "No Velocity found."

    @property
    def N(self) -> int:
        """particle number.

        Returns:
            int: particle number.
        """
        return self.__data.shape[0]

    @property
    def vol(self) -> float:
        """system volume.

        Returns:
            float: system volume.
        """
        return np.inner(self.__box[0], np.cross(self.__box[1], self.__box[2]))

    @property
    def rho(self) -> float:
        """system number density.

        Returns:
            float: system number density.
        """
        return self.N / self.vol

    def _get_format_(self) -> str:
        """Attempt to infer the file format from the filename, supporting both prefix and suffix formats.

        Returns:
            str: The format of the file.
        """
        # Assuming the format name is in the suffix
        try:
            postfix = os.path.split(self.__file_name)[-1].split(".")[-1]
            fmt = postfix.lower()
            assert fmt in [
                "dump"], "Failed to determine the format at the end of the filename"
            return fmt
        # Assuming the format name is in the prefix
        except:
            prefix = os.path.split(self.__file_name)[-1].split(".")[0]
            fmt = prefix.lower()
            assert fmt in [
                "dump"], "Filename must start or end with dump. (Allow uppercase letters, e.g. 1.dump, dump.1, 1.DUMP, DUMP.1)"
            return fmt

    def _read_file_by_format(self) -> tuple[pl.DataFrame, np.ndarray, np.ndarray]:
        """Read the file based on the format. Currently only support dump format.

        Returns:
            tuple[pl.DataFrame, np.ndarray, np.ndarray]: The datasets, the simulation box, and the boundary.
        """
        if self.fmt == 'dump':
            return self._read_dump_()
        ####################### TODO: add other formats #######################
        else:
            raise ValueError(f"Unsupported format: {self.fmt}")

    def _read_dump_(self) -> tuple[pl.DataFrame, np.ndarray, np.ndarray]:
        """Read the dump file.

        Returns:
            tuple[pl.DataFrame, np.ndarray, np.ndarray]: The datasets, the simulation box, and the boundary.
        """
        assert self.fmt in ["dump"], "Only support dump format."
        dump_head = []
        if self.fmt == "dump":
            with open(self.file_name) as op:
                for _ in range(9):
                    dump_head.append(op.readline())

        line = dump_head[4].split()
        boundary = [1 if i == "pp" else 0 for i in line[-3:]]
        if "xy" in line:
            xlo_bound, xhi_bound, xy = np.array(dump_head[5].split(), float)
            ylo_bound, yhi_bound, xz = np.array(dump_head[6].split(), float)
            zlo_bound, zhi_bound, yz = np.array(dump_head[7].split(), float)
            xlo = xlo_bound - min(0.0, xy, xz, xy + xz)
            xhi = xhi_bound - max(0.0, xy, xz, xy + xz)
            ylo = ylo_bound - min(0.0, yz)
            yhi = yhi_bound - max(0.0, yz)
            zlo = zlo_bound
            zhi = zhi_bound
            box = np.array(
                [
                    [xhi - xlo, 0, 0],
                    [xy, yhi - ylo, 0],
                    [xz, yz, zhi - zlo],
                    [xlo, ylo, zlo],
                ]
            )
        else:
            box = np.array([i.split()[:2]
                           for i in dump_head[5:8]]).astype(float)
            xlo, xhi = np.array(dump_head[5].split(), float)
            ylo, yhi = np.array(dump_head[6].split(), float)
            zlo, zhi = np.array(dump_head[7].split(), float)
            box = np.array(
                [
                    [xhi - xlo, 0, 0],
                    [0, yhi - ylo, 0],
                    [0, 0, zhi - zlo],
                    [xlo, ylo, zlo],
                ]
            )
        col_names = dump_head[8].split()[2:]

        try:

            data = pl.read_csv(
                self.file_name,
                separator=" ",
                skip_rows=9,
                new_columns=col_names,
                columns=range(len(col_names)),
                has_header=False,
                truncate_ragged_lines=True,
            )
        except Exception:
            data = pl.read_csv(
                self.file_name,
                separator=" ",
                skip_rows=9,
                new_columns=col_names,
                columns=range(len(col_names)),
                has_header=False,
                truncate_ragged_lines=True,
                infer_schema_length=None,
            )

        if "xs" in data.columns:
            pos = data.select("xs", "ys", "zs").to_numpy() @ box[:-1]

            data = data.with_columns(x=pos[:, 0], y=pos[:, 1], z=pos[:, 2]).select(
                pl.all().exclude("xs", "ys", "zs")
            )

        return data, box, boundary

    def __repr__(self) -> str:
        """The representation of the System class.

        Returns:
            str: The representation of the System class.
        """
        return f"Filename: {self.__file_name}\nAtom Number: {self.N}\nSimulation Box:\n{self.__box}\nTimeStep: {self.__timestep} fs\nBoundary: {self.__boundary}\nParticle Information:\n{self.__data}"


class MDSys(list):
    """The molecular dynamic postprocessing system (for multiple trajectory files) of Extrempy.
    """

    def __init__(self, root_dir: str, dt: float = 1.0, traj_dir: str = 'traj', format: str = 'dump', type_name: list[str] = None, is_printf: bool = True):
        """Initialize the MDSys class.

        Args:
            root_dir (str): The root directory of the trajectory files.
            dt (float, optional): The time step. Defaults to 1.0.
            traj_dir (str, optional): The trajectory directory. Defaults to 'traj'.
            format (str, optional): The format of the trajectory files. Defaults to 'dump'.
            type_name (list[str], optional): The type name. Defaults to None.
            is_printf (bool, optional): Whether to print the initialization information. Defaults to True.
        """

        strgs = '='*70 + '\n'
        strgs += 'Molecular Dynamic PostProcessing System of Extrempy is initialized \n'
        strgs += '='*70 + '\n'
        print(strgs)
        assert 'dump' in format, "Currently only support dump format."
        self.__root_dir = root_dir
        self.__dt = dt
        self.__traj_dir = traj_dir
        self.__format = format
        self.__type_name = type_name
        self.__is_printf = is_printf

        # read and sort the trajectory file name
        dump_list = self.read_and_sort_files()

        # calculate the dump frequency
        self.dump_freq = self.calculate_dump_freq(dump_list)

        # if the format is dump, dump.*, *.dump etc., then the format is 'dump'
        if 'dump' in self.__format:
            __fmt = 'dump'
        else:
            raise ValueError(f"Unsupported format: {self.__format}")

        # read the trajectory files, ever single trajectory file is a System instance
        progress_bar = tqdm(dump_list)
        for filename in progress_bar:
            progress_bar.set_description(f"Reading {filename}")
            system = System(file_name=filename,
                            dt=self.__dt, fmt=__fmt, type_name=self.__type_name)
            # The MDSys instance is a list of System instances, and each System instance will be added to the MDSys instance
            self.append(system)

        # Extract data from each System instance
        pos_list, box_list, inverse_box_list, vel_list, type_list = [], [], [], [], []
        for system in self:
            pos_list.append(system.pos)
            box_list.append(system.box)
            inverse_box_list.append(np.linalg.inv(system.box[:-1]))
            vel_list.append(system.vel)
            type_list.append(system.data["type"])

        self.position = np.array(pos_list)
        self.box_list = np.array(box_list)
        self.inverse_box_list = np.array(inverse_box_list)
        self.velocity = np.array(vel_list)
        self.type = np.array(type_list)
        self.type_list = np.unique(self.type)

        self.numb_frames = self.position.shape[0]
        self.numb_atoms = self.position.shape[1]
        self.numb_types = len(self.type_list)
        self.cells = self.box_list[0][:-1]
        self.delta_k = 2 * np.pi / np.diag(self.cells)

    @property
    def root_dir(self) -> str:
        """The root directory of the trajectory files.

        Returns:
            str: The root directory of the trajectory files.
        """
        return self.__root_dir

    @property
    def dt(self) -> float:
        """The time step (unit: fs) of the dump file.

        Returns:
            float: The time step (unit: fs).
        """
        return self.__dt * self.dump_freq

    @property
    def format(self) -> str:
        """The format of the trajectory files.

        Returns:
            str: The format of the trajectory files.
        """
        return self.__format

    @property
    def traj_dir(self) -> str:
        """The trajectory directory.

        Returns:
            str: The trajectory directory.
        """
        return os.path.join(self.__root_dir, self.__traj_dir)

    @property
    def type_name(self) -> list[str]:
        """The type name.

        Returns:
            list[str]: The type name.
        """
        return self.__type_name

    @property
    def is_printf(self) -> bool:
        """Whether to print the initialization information.

        Returns:
            bool: Whether to print the initialization information.
        """
        return self.__is_printf

    def read_and_sort_files(self) -> list[str]:
        """Read all files in the directory and sort them.

        Returns:
            list[str]: The list of trajectory files.
        """
        # read all files in the directory
        if '*' not in self.__format:
            try:
                file_list = glob.glob(os.path.join(
                    self.traj_dir, '*.', self.__format))
                if len(file_list) == 0:
                    raise FileNotFoundError(
                        f"No files found in the directory {self.traj_dir} with the format {os.path.join('*.', self.__format)}.")
            except FileNotFoundError:
                file_list = glob.glob(os.path.join(
                    self.traj_dir, self.__format, '.*'))
                if len(file_list) == 0:
                    raise FileNotFoundError(
                        f"No files found in the directory {self.traj_dir} with the format {os.path.join(self.__format, '.*')}.")
        else:
            file_list = glob.glob(os.path.join(
                self.traj_dir, self.__format))
            if len(file_list) == 0:
                raise FileNotFoundError(
                    f"No files found in the directory {self.traj_dir} with the format {self.__format}.")

        # extract the number in the file name and sort
        file_list.sort(key=lambda x: int(
            re.findall(r'\d+', os.path.basename(x))[0]))

        return file_list

    @staticmethod
    def calculate_dump_freq(file_list: list[str]) -> int:
        """Calculate the dump frequency.

        Args:
            file_list (list[str]): The list of trajectory files.

        Returns:
            int: The dump frequency.
        """
        numbers = [int(re.findall(r'\d+', os.path.basename(f))[0])
                   for f in file_list]
        # calculate the difference between adjacent numbers
        intervals = [j - i for i, j in zip(numbers[:-1], numbers[1:])]

        # return the minimum interval
        return min(intervals)

    def _calc_sed_from_traj(self, save_dir: str, k_vec_tmp: np.ndarray, nk: int = 1,
                            SKIP: int = 0, INTERVAL: int = 1, suffix: str = None, plot: bool = True, loglocator: bool = True) -> np.ndarray:
        """Calculate spectral energy density (SED) from trajectory data.

        Args:
            save_dir (str): Directory to save the results
            k_vec_tmp (np.ndarray): k-vector direction
            nk (int, optional): k-vector multiplier. Defaults to 1.
            SKIP (int, optional): Number of frames to skip. Defaults to 0.
            INTERVAL (int, optional): Interval between frames. Defaults to 1.
            suffix (str, optional): Suffix for output files. Defaults to None.
            plot (bool, optional): Whether to plot the results. Defaults to True.
            loglocator (bool, optional): Whether to use log locator. Defaults to True.

        Returns:
            np.ndarray: Combined SED for all directions
        """
        # Create SEDCalc instance with trajectory data
        sed_calc = SEDCalc(
            position=self.position,    # (nframes, natoms, 3)
            velocity=self.velocity,    # (nframes, natoms, 3)
            type_array=self.type[0],   # Use first frame's type info
            cells=self.cells,          # (3, 3)
            dt=self.dt,
            SKIP=SKIP,
            INTERVAL=INTERVAL
        )

        # Calculate SED
        spectrum = sed_calc.calc_sed(
            k_vec_tmp=k_vec_tmp,
            nk=nk,
            save_dir=save_dir,
            suffix=suffix
        )
        if plot:
            fig, ax = sed_calc.plot_sed_1d(spectrum, nk, ax=None, loglocator=loglocator)
            return fig, ax
        return None, None

    def _calc_sed_2d_from_traj(self, save_dir: str, k_vec_tmp: np.ndarray, nk_range: tuple[float, float, float] = (1, 21, 0.5), SKIP: int = 0, INTERVAL: int = 1, suffix: str = None, plot: bool = True, k_max: float = None, loglocator: bool = True) -> tuple[plt.Figure, plt.Axes]:
        """Calculate and plot 2D spectral energy density.

        Args:
            save_dir (str): Directory to save results
            k_vec_tmp (np.ndarray): k-vector direction
            nk_range (tuple[float, float, float], optional): Range of nk values (start, end, step). Defaults to (1, 21, 0.5).
            SKIP (int, optional): Number of frames to skip. Defaults to 0.
            INTERVAL (int, optional): Interval between frames. Defaults to 1.
            suffix (str, optional): Suffix for output files. Defaults to None.
            plot (bool, optional): Whether to plot the results. Defaults to True.
            k_max (float, optional): Maximum k value. Defaults to None.
            loglocator (bool, optional): Whether to use log locator. Defaults to True.

        Returns:
            tuple[plt.Figure, plt.Axes]: Figure and axes objects if plot is True
        """
        # Create SEDCalc instance
        sed_calc = SEDCalc(
            position=self.position,
            velocity=self.velocity,
            type_array=self.type[0],
            cells=self.cells,
            dt=self.dt,
            SKIP=SKIP,
            INTERVAL=INTERVAL
        )

        # Generate 2D SED data
        K, Omega, Z = sed_calc.generate_sed_2d(
            save_dir=save_dir,
            k_vec_tmp=k_vec_tmp,
            nk_range=nk_range,
            suffix=suffix
        )

        if plot:
            # Plot 2D SED
            fig, ax = sed_calc.plot_sed_2d(K, Omega, Z, k_max=k_max, loglocator=loglocator)

            # Save plot if directory is provided
            if save_dir:
                fig.savefig(os.path.join(save_dir, f'sed_2d_{suffix}.png'),
                            dpi=300, bbox_inches='tight')

            return fig, ax
        return None, None

    def _calc_time_correlation(self, time_corre_max: int = 0, skip: int = 0, interval: int = 0, calc_pacf: bool = False, plot: bool = True, save_dir: str = None) -> dict:
        """Calculate time correlation functions (VACF, PACF) and VDOS.

        Args:
            time_corre_max (int, optional): Maximum correlation time. If 0, uses nframes/2. Defaults to 0.
            skip (int, optional): Number of frames to skip. Defaults to 0.
            interval (int, optional): Interval between frames. Defaults to 0.
            calc_pacf (bool, optional): Whether to calculate Position Autocorrelation Function (PACF). Defaults to False.
            plot (bool, optional): Whether to plot the results. Defaults to True.
            save_dir (str, optional): Directory to save the results. Defaults to None.

        Returns:
            dict: Dictionary containing:
                - 'time': Time array
                - 'vacf': Velocity autocorrelation function
                - 'pacf': Position autocorrelation function if calc_pacf is True
                - 'freq': Frequency array for VDOS
                - 'vdos': Vibrational density of states
        """

        # Create TimeCorrelationCalc instance
        tc_calc = TimeCorrelationCalc(
            position=self.position,    # (nframes, natoms, 3)
            velocity=self.velocity,    # (nframes, natoms, 3)
            dt=self.dt,             # timestep in fs
            time_corre_max=time_corre_max
        )

        # Calculate velocity autocorrelation function (VACF)
        vacf = tc_calc.calc_vacf(skip=skip, interval=interval)
        # Calculate position autocorrelation function (PACF) if calc_pacf is True
        if calc_pacf:
            pacf = tc_calc.calc_pacf(skip=skip, interval=interval)

        # Calculate vibrational density of states (VDOS)
        freq, vdos = tc_calc.calc_power_spectrum(vacf)

        # Store the results
        if calc_pacf:
            results = {
                'time': tc_calc.time,
                'vacf': vacf,
                'pacf': pacf,
                'freq': freq,
                'vdos': vdos
            }
        else:
            results = {
                'time': tc_calc.time,
                'vacf': vacf,
                'freq': freq,
                'vdos': vdos
            }

        # Plot the results if plot is True
        if plot:
            fig_ax = tc_calc.plot_time_correlation(results, save_dir)
            return results, fig_ax
        return results, None


if __name__ == "__main__":

    type_name = [
        "C",
    ]
    # sys = System(r"C:\Users\87627\Desktop\traj-1555\dump.0",
    #              type_name=type_name)
    # print(sys)
    TRAJ_DIR = r"C:\Users\87627\Desktop\traj-1555"
    OUT_DIR = TRAJ_DIR+r'\output'

    # timestep of MD (unit: fs)
    dt = 0.5

    sys = MDSys(TRAJ_DIR, format="dump.*", traj_dir='', type_name=type_name, dt=dt)

    # sys._calc_time_correlation(
    #     time_corre_max=1000,
    #     skip=0,
    #     interval=1,
    #     save_dir=OUT_DIR,
    #     calc_pacf=False,
    #     plot=True
    # )

    k_vec_tmp = np.array([1, 0, 0])  # GX方向

    fig, ax = sys._calc_sed_from_traj(
        save_dir=OUT_DIR,
        k_vec_tmp=k_vec_tmp,
        nk=10,
        suffix="GX",
        plot=True
    )

    fig, ax = sys._calc_sed_2d_from_traj(
        save_dir=OUT_DIR,
        k_vec_tmp=k_vec_tmp,
        nk_range=(1, 21),  # 计算nk从1到20
        suffix="GX",
        plot=True
    )
    

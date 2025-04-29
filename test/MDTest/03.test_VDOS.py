from extrempy.md import MDSys


# Name of the atom type, corresponding to the type in the dump file, e.g. ["C"], ["Mg", "O"], ["O", "H"], etc.
type_name = ["C",]

# Path to the directory containing the dump files
input_dir = r"C:\Users\87627\Desktop\traj-1555"

# Path to the output directory
output_dir = input_dir + r"\output"

# timestep of MD (unit: fs)
dt = 0.5

# Create a MDSys object
sys = MDSys(input_dir, format="dump.*", traj_dir='', type_name=type_name, dt=dt)

# Calculate the time correlation function (VACF and PACF(optional)), and the vibrational density of states (VDOS)
sys._calc_time_correlation(
    time_corre_max=1000,
    skip=0,
    interval=1,
    save_dir=output_dir,
    calc_pacf=False,
    plot=True
)
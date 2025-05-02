from extrempy.md import MDSys
import numpy as np


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

k0 = np.array([0, 0, 0])  # K0 direction
k1 = np.array([1, 0, 0])  # K1 direction
k_vec_tmp = k1 - k0  # GX direction

# Calculate the 1D SED for nk = 10 in the GX direction
fig, ax = sys._calc_sed_from_traj(
    save_dir=output_dir,
    k_vec_tmp=k_vec_tmp,
    nk=10, # calculate nk = 10 for example
    suffix="GX",
    plot=True
)

# Calculate the 2D SED for nk from 1 to 20 in the GX direction
fig, ax = sys._calc_sed_2d_from_traj(
    save_dir=output_dir,
    k_vec_tmp=k_vec_tmp,
    nk_range=(1, 21), 
    suffix="GX",
    plot=True
)
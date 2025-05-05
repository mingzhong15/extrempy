from extrempy.md import MDSys
import numpy as np
import matplotlib.pyplot as plt


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
supercell_size_of_k_vec_tmp = 15 # supercell size alone the k_vec_tmp direction
k_max = 2 * np.pi / (sys[0].cells[0, 0] / supercell_size_of_k_vec_tmp)

# Calculate the 1D SED for nk = 10 in the GX direction
fig, ax = sys._calc_sed_from_traj(
    save_dir=output_dir,
    k_vec_tmp=k_vec_tmp,
    nk=10, # calculate nk = 10 for example
    suffix="GX",
    plot=True,
    loglocator=True
)
# Any further settings of fig and ax can be added below
ax.set_xlim(0, 0.2)
plt.tight_layout()
plt.show()

# Calculate the 2D SED for nk from 0 to 20 with step 0.5 in the GX direction
fig, ax = sys._calc_sed_2d_from_traj(
    save_dir=output_dir,
    k_vec_tmp=k_vec_tmp,
    nk_range=(0, 21, 0.5), 
    suffix="GX",
    k_max=k_max,
    plot=True,
    loglocator=True
)
# Any further settings of fig and ax can be added below

plt.show()
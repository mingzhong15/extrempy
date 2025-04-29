from extrempy.md import MDSys


# Name of the atom type, corresponding to the type in the dump file, e.g. ["C"], ["Mg", "O"], ["O", "H"], etc.
type_name = ["C",]

# Path to the directory containing the dump files
input_dir = r"C:\Users\87627\Desktop\traj-1555"

# Path to the output directory
output_dir = input_dir + r"\output"

dt = 0.5
sys = MDSys(input_dir, format="dump.*", traj_dir='', type_name=type_name, dt=dt)
print(sys[0])
print(sys[-1])

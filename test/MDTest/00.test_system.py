from extrempy.md import System

# Name of the atom type, corresponding to the type in the dump file, e.g. ["C"], ["Mg", "O"], ["O", "H"], etc.
type_name = ["C",]

# Path to the single dump file
dump_file = r"C:\Users\87627\Desktop\traj-1555\dump.0"

# Create a System object
sys = System(dump_file, type_name=type_name)

# Print the System object
print(sys)

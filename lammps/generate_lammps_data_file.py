import os
from ase.io import read, write
from ase.data import atomic_numbers
# === Define element-to-type mapping ===
element_order = ['Cu', 'C', 'Cl', 'O', 'N', 'H', 'S']  
symbol_to_type = {sym: i + 1 for i, sym in enumerate(element_order)}
def convert_poscar_to_lammps(poscar_file="POSCAR", output_file="lammps2.data"):
    """Reads POSCAR from current folder and writes lammps.data with correct mapping."""
    if not os.path.isfile(poscar_file):
        print(f"POSCAR not found in current folder: {os.getcwd()}")
        return
    try:
        # Read POSCAR
        structure = read(poscar_file, format="vasp")

        # Reassign atomic numbers explicitly using ASE
        symbols = structure.get_chemical_symbols()
        structure.set_atomic_numbers([atomic_numbers[s] for s in symbols])

        # Write LAMMPS data file in the same folder
        write(output_file, structure, format="lammps-data")

        print(f"Successfully wrote {output_file} from {poscar_file}")
        print("Element type mapping:")
        for sym in element_order:
            print(f"  {sym} → {symbol_to_type[sym]}")

    except Exception as e:
        print(f" Error converting {poscar_file}: {e}")

# Run conversion in the current working directory
convert_poscar_to_lammps()


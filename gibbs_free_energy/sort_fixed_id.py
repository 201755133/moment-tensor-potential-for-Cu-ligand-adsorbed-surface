from ase.io import read
from ase.constraints import FixAtoms

# Load the structure from POSCAR
atoms = read('POSCAR', format='vasp')

# Initialize a list to store indices of atoms to be fixed
fixed_indices = []

# Open POSCAR file to check Selective Dynamics lines
with open('POSCAR', 'r') as f:
    lines = f.readlines()
    selective_dynamics = False
    # Find the start of the atomic positions section
    for i, line in enumerate(lines):
        if 'Selective dynamics' in line:
            selective_dynamics = True
        elif selective_dynamics and 'Direct' in line:
            # Start reading atom positions with F/F/F constraints after "Direct"
            position_start_index = i + 1
            break

    # Now, parse each line in the positions section for constraints
    for i, line in enumerate(lines[position_start_index:], start=position_start_index):
        parts = line.split()
        if len(parts) > 3:
            # Check if the constraint is `F F F` in the POSCAR format
            if parts[-3:] == ['F', 'F', 'F']:
                # Add atom index to the fixed list
                fixed_indices.append(i - position_start_index)

# Print all fixed atom indices
print("Fixed atom indices:", fixed_indices)

# Apply the FixAtoms constraint in ASE
constraint = FixAtoms(indices=fixed_indices)
atoms.set_constraint(constraint)

# Save the structure with constraints (optional)
atoms.write('CONSTRAINED_POSCAR', format='vasp')


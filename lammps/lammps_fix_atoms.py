import os
def extract_fixed_indices(poscar_path):
"""Extract 1-based atom indices marked F F F in POSCAR with Selective dynamics."""
    with open(poscar_path, "r") as f:
        lines = f.readlines()
    fixed_indices = []
    selective_dynamics = False
    position_start_index = None
    # Detect 'Selective dynamics' and locate start of atomic positions
    for i, line in enumerate(lines):
        if "selective dynamics" in line.lower():
            selective_dynamics = True
        elif selective_dynamics and ("direct" in line.lower() or "cartesian" in line.lower()):
            position_start_index = i + 1
            break
    if position_start_index is None:
        raise ValueError(f"No Selective dynamics block found in {poscar_path}")
    # Parse atomic positions and flags
    for i, line in enumerate(lines[position_start_index:], start=position_start_index):
        parts = line.split()
        if len(parts) >= 6 and parts[3:6] == ["F", "F", "F"]:
            fixed_indices.append(i - position_start_index + 1)  # 1-based index
    return fixed_indices
def update_lammps_file_in_current_dir():
    """Read POSCAR and insert fixed atom group line in in.lammps file in current folder."""
    poscar_path = "POSCAR"
    lammps_path = "in.lammps"
    # Check files
    if not os.path.isfile(poscar_path) or not os.path.isfile(lammps_path):
        print("Error: Missing POSCAR or in.lammps in current directory.")
        return
    try:
        fixed_ids = extract_fixed_indices(poscar_path)
    except ValueError as e:
        print(f"Error: {e}")
        return
    if not fixed_ids:
        print("Warning: No atoms marked as F F F in POSCAR.")
        return
    group_line = "group fixc id " + " ".join(map(str, fixed_ids)) + "\n"
    # Read LAMMPS input
    with open(lammps_path, "r") as f:
        lines = f.readlines()
    # Avoid duplicate insertion
    if any("group fixc id" in line for line in lines):
        print("Warning: group line already exists in in.lammps, skipping insertion.")
        return
    # Insert group line after 'read_data lammps.data'
    new_lines = []
    inserted = False
    for line in lines:
        new_lines.append(line)
        if not inserted and "read_data" in line and "lammps.data" in line:
            new_lines.append(group_line)
            inserted = True
    if inserted:
        with open(lammps_path, "w") as f:
            f.writelines(new_lines)
        print(f"Inserted group line into in.lammps: {group_line.strip()}")
    else:
        print("Warning: 'read_data lammps.data' not found, group line not inserted.")
if __name__ == "__main__":
    update_lammps_file_in_current_dir()


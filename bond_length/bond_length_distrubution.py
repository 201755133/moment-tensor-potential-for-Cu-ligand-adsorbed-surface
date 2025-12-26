import numpy as np
import matplotlib.pyplot as plt
from ase.io import read
from jarvis.core.atoms import Atoms as JarvisAtoms

# ===============================
# Step 1: Read POSCAR and Identify Atom Indices
# ===============================
# Read POSCAR using Jarvis Tools
jarvis_atoms = JarvisAtoms.from_poscar("POSCAR")
elements = jarvis_atoms.elements  # List of element symbols

# Global indices for Cu, O, and H atoms (based on POSCAR order)
cu_indices_global = [i for i, el in enumerate(elements) if el == "Cu"]
o_indices_global  = [i for i, el in enumerate(elements) if el == "O"]
h_indices_global  = [i for i, el in enumerate(elements) if el == "H"]

print("Total Cu atoms:", len(cu_indices_global))
print("Total O atoms:", len(o_indices_global))
print("Total H atoms:", len(h_indices_global))

# ===============================
# Step 2: Define Cutoff Distances
# ===============================
water_cutoff = 3.5  # Only consider Cu atoms that have an O or H within 3.5 Å
cutoff_CuO = 2.5    # For defining a Cu-O bond
cutoff_OH  = 3.0    # For defining a Cu-H bond
cutoff_CuCu = 3.0   # For defining a Cu-Cu bond

# ===============================
# Step 3: Read the XDATCAR Trajectory Using ASE
# ===============================
# Read all configurations from XDATCAR (ordering must match POSCAR)
traj = read("XDATCAR", index=":")
print("Total frames in trajectory:", len(traj))

# ===============================
# Step 4: Loop Over Frames and Compute Bond Lengths
# ===============================
# Lists to store bond lengths (over all frames)
CuO_lengths = []   # For Cu-O bonds
CuH_lengths = []   # For Cu-H bonds
CuCu_lengths = []  # For Cu-Cu bonds

# Combine water atoms indices (O and H) for filtering
water_indices_global = o_indices_global + h_indices_global

for frame in traj:
    positions = frame.get_positions()  # Array shape (N_atoms, 3)
    
    # Filter eligible Cu atoms: keep only those that have any water (O or H) within water_cutoff
    eligible_cu = []
    for i in cu_indices_global:
        pos_cu = positions[i]
        water_positions = positions[water_indices_global]
        if np.any(np.linalg.norm(water_positions - pos_cu, axis=1) < water_cutoff):
            eligible_cu.append(i)
    
    # --- Cu-O bond lengths ---
    # For each eligible Cu, compute distances to every O atom (global) that are within cutoff_CuO.
    for i in eligible_cu:
        pos_cu = positions[i]
        for j in o_indices_global:
            pos_o = positions[j]
            d = np.linalg.norm(pos_cu - pos_o)
            if d < cutoff_CuO:
                CuO_lengths.append(d)
                
    # --- Cu-H bond lengths ---
    # For each eligible Cu, compute distances to every H atom (global) that are within cutoff_OH.
    for i in eligible_cu:
        pos_cu = positions[i]
        for j in h_indices_global:
            pos_h = positions[j]
            d = np.linalg.norm(pos_cu - pos_h)
            if d < cutoff_OH:
                CuH_lengths.append(d)
                
    # --- Cu-Cu bond lengths ---
    # For each unique pair of eligible Cu atoms, compute distances that are within cutoff_CuCu.
    for idx, i in enumerate(eligible_cu):
        for j in eligible_cu[idx+1:]:
            d = np.linalg.norm(positions[i] - positions[j])
            if d < cutoff_CuCu:
                CuCu_lengths.append(d)

print("Total Cu-O bonds found:", len(CuO_lengths))
print("Total Cu-H bonds found:", len(CuH_lengths))
print("Total Cu-Cu bonds found:", len(CuCu_lengths))

# ===============================
# Step 5: Compute Histograms for Bond Length Distributions
# ===============================
# Define bin ranges for each bond type
bins_CuO = np.linspace(0, 4.0, 101)   # 0 to 4 Å, 100 bins for Cu-O
bins_CuH = np.linspace(0, 3.0, 101)    # 0 to 3 Å, 100 bins for Cu-H
bins_CuCu = np.linspace(0, 4.0, 101)   # 0 to 4 Å, 100 bins for Cu-Cu

# Compute histograms
hist_CuO, edges_CuO = np.histogram(CuO_lengths, bins=bins_CuO)
hist_CuH, edges_CuH = np.histogram(CuH_lengths, bins=bins_CuH)
hist_CuCu, edges_CuCu = np.histogram(CuCu_lengths, bins=bins_CuCu)

# Compute bin centers
centers_CuO = 0.5 * (edges_CuO[:-1] + edges_CuO[1:])
centers_CuH = 0.5 * (edges_CuH[:-1] + edges_CuH[1:])
centers_CuCu = 0.5 * (edges_CuCu[:-1] + edges_CuCu[1:])

# ===============================
# Step 6: Save Output Data Files (10 Decimal Places)
# ===============================
np.savetxt("bond_length_distribution_CuO.dat", 
           np.column_stack((centers_CuO, hist_CuO)), fmt="%.10f",
           header="r (Å)    Frequency for Cu-O", comments="")
np.savetxt("bond_length_distribution_CuH.dat", 
           np.column_stack((centers_CuH, hist_CuH)), fmt="%.10f",
           header="r (Å)    Frequency for Cu-H", comments="")
np.savetxt("bond_length_distribution_CuCu.dat", 
           np.column_stack((centers_CuCu, hist_CuCu)), fmt="%.10f",
           header="r (Å)    Frequency for Cu-Cu", comments="")

print("Bond length distribution data saved to files.")

# ===============================
# Step 7: Plot the Bond Length Distributions
# ===============================
plt.figure(figsize=(14, 6), dpi=300)

plt.subplot(1, 3, 1)
plt.plot(centers_CuO, hist_CuO, '-b', lw=3)
plt.xlabel("Cu-O bond length (Å)", fontsize=14)
plt.ylabel("Frequency", fontsize=14)
#plt.title("Bond Length Distribution: Cu-O", fontsize=16)

plt.subplot(1, 3, 2)
plt.plot(centers_CuH, hist_CuH, '-r', lw=3)
plt.xlabel("Cu-H bond length (Å)", fontsize=14)
plt.ylabel("Frequency", fontsize=14)
#plt.title("Bond Length Distribution: Cu-H", fontsize=16)

plt.subplot(1, 3, 3)
plt.plot(centers_CuCu, hist_CuCu, '-g', lw=3)
plt.xlabel("Cu-Cu bond length (Å)", fontsize=14)
plt.ylabel("Frequency", fontsize=14)
#plt.title("Bond Length Distribution: Cu-Cu", fontsize=16)

plt.tight_layout()
plt.savefig("bond_lengths.png", dpi=300, bbox_inches='tight') 
plt.show()


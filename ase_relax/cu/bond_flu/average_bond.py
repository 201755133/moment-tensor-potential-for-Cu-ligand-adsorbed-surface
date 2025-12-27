import numpy as np
import matplotlib.pyplot as plt
from ase.io import read
from jarvis.core.atoms import Atoms as JarvisAtoms

# ===============================
# Step 1: Read POSCAR and Identify Atom Indices
# ===============================
jarvis_atoms = JarvisAtoms.from_poscar("POSCAR")
elements = jarvis_atoms.elements

# Identify atomic indices
o_indices_global  = [i for i, el in enumerate(elements) if el == "O"]
h_indices_global  = [i for i, el in enumerate(elements) if el == "H"]
cu_indices_global = [i for i, el in enumerate(elements) if el == "Cu"]

# Define Cutoff Distances
cutoff_CuO = 2.5    # Cu-O bond
cutoff_CuH = 3.0    # Cu-H bond
cutoff_CuCu = 3.0   # Cu-Cu bond
cutoff_CuO_neigh = 3.5  # Cu-Cu atoms must be near O

# ===============================
# Step 2: Read Trajectory Data
# ===============================
traj = read("XDATCAR", index=":")
total_frames = len(traj)

# Storage for average bond lengths per frame
avg_CuO = []
avg_CuH = []
avg_CuCu = []

# ===============================
# Step 3: Loop Over Frames to Compute Average Bond Lengths
# ===============================
for frame in traj:
    positions = frame.get_positions()
    CuO_lengths, CuH_lengths, CuCu_lengths = [], [], []
    
    # Filter Cu atoms near O atoms
    eligible_cu = []
    for i in cu_indices_global:
        pos_cu = positions[i]
        if np.any(np.linalg.norm(positions[o_indices_global] - pos_cu, axis=1) < cutoff_CuO_neigh):
            eligible_cu.append(i)
    
    # Compute Cu-O bond lengths
    for i in eligible_cu:
        pos_cu = positions[i]
        for j in o_indices_global:
            d = np.linalg.norm(pos_cu - positions[j])
            if d < cutoff_CuO:
                CuO_lengths.append(d)
    
    # Compute Cu-H bond lengths
    for i in eligible_cu:
        pos_cu = positions[i]
        for j in h_indices_global:
            d = np.linalg.norm(pos_cu - positions[j])
            if d < cutoff_CuH:
                CuH_lengths.append(d)
    
    # Compute Cu-Cu bond lengths (among Cu atoms near O)
    for idx, i in enumerate(eligible_cu):
        for j in eligible_cu[idx+1:]:
            d = np.linalg.norm(positions[i] - positions[j])
            if d < cutoff_CuCu:
                CuCu_lengths.append(d)
    
    # Store average bond lengths for this frame
    avg_CuO.append(np.mean(CuO_lengths) if CuO_lengths else np.nan)
    avg_CuH.append(np.mean(CuH_lengths) if CuH_lengths else np.nan)
    avg_CuCu.append(np.mean(CuCu_lengths) if CuCu_lengths else np.nan)

# ===============================
# Step 4: Plot Average Bond Length vs. Number of Configurations
# ===============================
frames = np.arange(1, total_frames + 1)

plt.figure(figsize=(8, 6), dpi=300)
plt.plot(frames, avg_CuO, '-b', label='Cu-O')
plt.plot(frames, avg_CuH, '-r', label='Cu-H')
plt.plot(frames, avg_CuCu, '-g', label='Cu-Cu (near O)')
plt.xlabel("Number of Configurations", fontsize=24)
plt.ylabel("Average Bond Length (Å)", fontsize=24)
#plt.title("Average Bond Length per Configuration", fontsize=16)
plt.legend(fontsize=24)
plt.grid()
plt.savefig("average_bond_lengths.png", dpi=300, bbox_inches='tight')
plt.show()

# ===============================
# Step 5: Save Data
# ===============================
data = np.column_stack((frames, avg_CuO, avg_CuH, avg_CuCu))
np.savetxt("average_bond_lengths.dat", data, fmt="%.10f", header="Frame Cu-O Cu-H Cu-Cu", comments="")
print("Saved average bond lengths to file.")

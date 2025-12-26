import numpy as np
import matplotlib.pyplot as plt
from jarvis.core.atoms import Atoms as JarvisAtoms
from ase.io import read
from scipy.stats import gaussian_kde

# -------------------------------
# Step 1: Read POSCAR to Get Atom Indices
# -------------------------------
# Read the POSCAR using Jarvis Tools
jarvis_atoms = JarvisAtoms.from_poscar("POSCAR")
elements = jarvis_atoms.elements  # List of element symbols

# Extract indices for Cu, O, and H atoms
cu_indices = [i for i, el in enumerate(elements) if el == "Cu"]
o_indices  = [i for i, el in enumerate(elements) if el == "O"]
h_indices  = [i for i, el in enumerate(elements) if el == "H"]

print("Total Cu atoms:", len(cu_indices))
print("Total O atoms:", len(o_indices))
print("Total H atoms:", len(h_indices))

# -------------------------------
# Step 2: Define Bonding Cutoffs
# -------------------------------
cutoff_CuO = 2.5  # Cutoff distance for Cu-O bond (Å); adjust as needed
cutoff_OH  = 1.2  # Cutoff distance for O-H bond (Å); adjust as needed

# -------------------------------
# Step 3: Read XDATCAR Trajectory using ASE
# -------------------------------
# Read all configurations from XDATCAR (assumes the order is consistent with POSCAR)
traj = read("XDATCAR", index=":")  # This loads a list of ASE Atoms objects
print("Total configurations in XDATCAR:", len(traj))

# -------------------------------
# Step 4: Loop Over Configurations and Compute Cu-O-H Angles
# -------------------------------
all_angles = []     # List to store all computed angles (in degrees)
frame_angles = []   # List to store tuples: (frame_number, angle)

for frame_num, frame in enumerate(traj):
    positions = frame.get_positions()  # (n_atoms, 3) positions for the current frame

    # For each oxygen atom, find nearby Cu and H neighbors within the specified cutoffs
    for i in o_indices:
        pos_O = positions[i]

        # Find Cu neighbors for the current O
        cu_neighbors = []
        for j in cu_indices:
            pos_Cu = positions[j]
            if np.linalg.norm(pos_Cu - pos_O) < cutoff_CuO:
                cu_neighbors.append(j)

        # Find H neighbors for the current O
        h_neighbors = []
        for k in h_indices:
            pos_H = positions[k]
            if np.linalg.norm(pos_H - pos_O) < cutoff_OH:
                h_neighbors.append(k)

        # Compute the Cu-O-H angle (with O as the vertex) for every combination
        for cu in cu_neighbors:
            for h in h_neighbors:
                pos_Cu = positions[cu]
                pos_H  = positions[h]
                vec_Cu = pos_Cu - pos_O
                vec_H  = pos_H  - pos_O
                cos_angle = np.dot(vec_Cu, vec_H) / (np.linalg.norm(vec_Cu) * np.linalg.norm(vec_H))
                angle = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
                all_angles.append(angle)
                frame_angles.append((frame_num, angle))

print("Total Cu-O-H angles computed:", len(all_angles))

# -------------------------------
# Step 5: Write Angle Data to a File with 10 Decimal Places
# -------------------------------
output_filename = "CuOH_angles.dat"
with open(output_filename, "w") as f:
    # Write a header (optional)
    f.write("# Frame_number Cu-O-H_angle(degrees) (10 decimal places)\n")
    for frame_num, angle in frame_angles:
        f.write(f"{frame_num:5d} {angle:20.10f}\n")
print(f"Angle data written to {output_filename}")

# -------------------------------
# Step 6: Compute and Plot the Smooth Angle Distribution (KDE)
# -------------------------------
if len(all_angles) > 0:
    all_angles = np.array(all_angles)
    kde = gaussian_kde(all_angles)
    x_vals = np.linspace(0, 180, 400)  # Adjust range if necessary
    density = kde(x_vals)
    
    # Find the angle at maximum density
    max_idx = np.argmax(density)
    max_angle = x_vals[max_idx]
    # Save x values (angles) and probability densities to a file
    kde_output_filename = "CuOH_angle_distribution.dat"
    with open(kde_output_filename, "w") as f:
        f.write("# Cu-O-H_Angle (degrees) Probability_Density\n")
        for x, y in zip(x_vals, density):
            f.write(f"{x:10.4f} {y:15.10f}\n")

    print(f"KDE data written to {kde_output_filename}")

    # Plot the KDE (smooth distribution)
    plt.figure(figsize=(10, 6), dpi= 300)
    plt.plot(x_vals, density, color='blue', lw=3, label='Cu-O-H Angle Distribution')
    plt.axvline(x=max_angle, color='red', linestyle='--', lw=2, 
                label=f'Maximum at {max_angle:.1f}°')
    plt.xlabel("Cu-O-H Angle (degrees)", fontsize=16)
    plt.ylabel("Probability Density", fontsize=16)
    #plt.title("Smooth Distribution of Cu-O-H Angles", fontsize=18)
    plt.legend(fontsize=14)
    plt.tight_layout()
    plt.savefig("angle_distribution.png", dpi=300, bbox_inches='tight')  
    plt.show()
else:
    print("No angles were computed; check cutoff distances and structure.")


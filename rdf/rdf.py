import numpy as np
import matplotlib.pyplot as plt
from ase.io import read
from jarvis.core.atoms import Atoms as JarvisAtoms
from scipy.ndimage import gaussian_filter1d  # For smoothing

# -------------------------------
# Step 1: Read POSCAR to Get Atom Indices
# -------------------------------
# Use Jarvis Tools to read POSCAR and obtain element labels.
jarvis_atoms = JarvisAtoms.from_poscar("POSCAR") 
elements = jarvis_atoms.elements  # List of element symbols

# Global indices for each element type
cu_indices_global = [i for i, el in enumerate(elements) if el == "Cu"]
o_indices_global  = [i for i, el in enumerate(elements) if el == "O"]
h_indices_global  = [i for i, el in enumerate(elements) if el == "H"]

print("Total Cu atoms:", len(cu_indices_global))
print("Total O atoms:", len(o_indices_global))
print("Total H atoms:", len(h_indices_global))

# -------------------------------
# Step 2: Define Parameters for Filtering and RDF
# -------------------------------
# Only consider Cu atoms that are "near" water atoms.
water_cutoff = 3.5  # (Å) If a Cu has any O or H closer than this, it is included.

# RDF parameters
r_min = 0.0
r_max = 10.0
n_bins = 200
bins = np.linspace(r_min, r_max, n_bins+1)
bin_centers = 0.5 * (bins[:-1] + bins[1:])
dr = bins[1] - bins[0]

# Initialize histograms (will accumulate counts over frames)
hist_CuCu = np.zeros(n_bins)
hist_CuO  = np.zeros(n_bins)
hist_CuH  = np.zeros(n_bins)

# Normalization accumulators:
total_pairs_CuCu = 0.0  # Sum over frames of (number of filtered Cu pairs)
total_ref_CuO  = 0.0    # Sum over frames of (number of filtered Cu * number of O)
total_ref_CuH  = 0.0    # Sum over frames of (number of filtered Cu * number of H)

# -------------------------------
# Step 3: Read the XDATCAR Trajectory using ASE
# -------------------------------
# Read all configurations from XDATCAR (assumes the ordering of atoms is identical to POSCAR)
traj = read("XDATCAR", index=":")
N_frames = len(traj)
print("Total frames in trajectory:", N_frames)

# Get the cell from the first frame and calculate surface area (assuming z is non-periodic/vacuum)
cell = traj[0].get_cell()
# Calculate surface area (assuming xy plane is the surface)
A = np.linalg.norm(np.cross(cell[0], cell[1]))  # Area of the surface
print("Surface area:", A, "Å²")

# -------------------------------
# Step 4: Loop Over Frames and Accumulate RDF Data with PBC
# -------------------------------
def pbc_distance(pos1, pos2, cell):
    """Calculate distance between two points with PBC consideration"""
    delta = pos1 - pos2
    # Apply PBC to each coordinate
    for i in range(3):
        delta[i] -= round(delta[i]/cell[i,i]) * cell[i,i]
    return np.linalg.norm(delta)

for frame in traj:
    positions = frame.get_positions()  # shape (N_atoms, 3)
    cell = frame.get_cell()

    # Filter Cu: only include those Cu that have any water atom (O or H) within water_cutoff.
    water_indices = o_indices_global + h_indices_global  # Combined water atoms indices
    filtered_cu = []
    for i in cu_indices_global:
        pos_cu = positions[i]
        for j in water_indices:
            d = pbc_distance(pos_cu, positions[j], cell)
            if d < water_cutoff:
                filtered_cu.append(i)
                break  # Only need one water atom within cutoff
    
    N_filtered = len(filtered_cu)
    if N_filtered < 1:
        continue  # Skip frame if no Cu meets the criteria

    # -------------------------------
    # Cu-Cu: Compute distances for each unique pair in filtered Cu.
    for idx, i in enumerate(filtered_cu):
        for j in filtered_cu[idx+1:]:
            d = pbc_distance(positions[i], positions[j], cell)
            if d < r_max:
                bin_index = int((d - r_min) / dr)
                if 0 <= bin_index < n_bins:
                    hist_CuCu[bin_index] += 1
    total_pairs_CuCu += N_filtered * (N_filtered - 1) / 2.0

    # -------------------------------
    # Cu-O: For each filtered Cu, compute distances to every O atom (global).
    for i in filtered_cu:
        for j in o_indices_global:
            d = pbc_distance(positions[i], positions[j], cell)
            if d < r_max:
                bin_index = int((d - r_min) / dr)
                if 0 <= bin_index < n_bins:
                    hist_CuO[bin_index] += 1
    total_ref_CuO += N_filtered * len(o_indices_global)

    # -------------------------------
    # Cu-H: For each filtered Cu, compute distances to every H atom (global).
    for i in filtered_cu:
        for j in h_indices_global:
            d = pbc_distance(positions[i], positions[j], cell)
            if d < r_max:
                bin_index = int((d - r_min) / dr)
                if 0 <= bin_index < n_bins:
                    hist_CuH[bin_index] += 1
    total_ref_CuH += N_filtered * len(h_indices_global)

# -------------------------------
# Step 5: Normalize Histograms to Obtain g(r) - Using AREA for surface systems
# -------------------------------
g_CuCu = np.zeros(n_bins)
g_CuO  = np.zeros(n_bins)
g_CuH  = np.zeros(n_bins)

# Smoothing parameter (adjust as needed)
sigma = 1.0  # Higher sigma = more smoothing

for i in range(n_bins):
    # For surface systems, we use 2D normalization (area instead of volume)
    shell_area = 2 * np.pi * bin_centers[i] * dr  # Circumference of ring
    
    # For Cu-Cu, normalization with total pairs count
    if total_pairs_CuCu > 0:
        g_CuCu[i] = (hist_CuCu[i] / total_pairs_CuCu) / (shell_area / A)
    
    # For Cu-O, normalization with total reference pairs (filtered Cu * total O)
    if total_ref_CuO > 0:
        g_CuO[i] = (hist_CuO[i] / total_ref_CuO) / (shell_area / A)
    
    # For Cu-H, normalization with total reference pairs (filtered Cu * total H)
    if total_ref_CuH > 0:
        g_CuH[i] = (hist_CuH[i] / total_ref_CuH) / (shell_area / A)

# Apply Gaussian smoothing to the RDFs
g_CuCu_smooth = gaussian_filter1d(g_CuCu, sigma=sigma)
g_CuO_smooth = gaussian_filter1d(g_CuO, sigma=sigma)
g_CuH_smooth = gaussian_filter1d(g_CuH, sigma=sigma)

# -------------------------------
# Step 6: Write Output Data Files with 10-Decimal Precision
# -------------------------------
np.savetxt("rdf_CuCu.dat", np.column_stack((bin_centers, g_CuCu, g_CuCu_smooth)), fmt="%.10f",
           header="r (Å)    g(r) for Cu-Cu    g(r)_smooth", comments="")
np.savetxt("rdf_CuO.dat",  np.column_stack((bin_centers, g_CuO, g_CuO_smooth)),  fmt="%.10f",
           header="r (Å)    g(r) for Cu-O     g(r)_smooth", comments="")
np.savetxt("rdf_CuH.dat",  np.column_stack((bin_centers, g_CuH, g_CuH_smooth)),  fmt="%.10f",
           header="r (Å)    g(r) for Cu-H     g(r)_smooth", comments="")

print("RDF data written to 'rdf_CuCu.dat', 'rdf_CuO.dat', and 'rdf_CuH.dat'.")

# -------------------------------
# Step 7: Plot the RDF Curves
# -------------------------------
plt.figure(figsize=(12, 8))
plt.plot(bin_centers, g_CuCu_smooth, label="Cu-Cu RDF", lw=3)
plt.plot(bin_centers, g_CuO_smooth,  label="Cu-O RDF",  lw=3)
plt.plot(bin_centers, g_CuH_smooth,  label="Cu-H RDF",  lw=3)
plt.xlabel("r (Å)", fontsize=16)
plt.ylabel("g(r)", fontsize=16)
plt.legend(fontsize=14)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.grid(alpha=0.3)
plt.tight_layout()

# Save the figure as a high-resolution PNG file
plt.savefig("Radial_Distribution_Functions_Surface.png", dpi=300, bbox_inches='tight')
plt.show()

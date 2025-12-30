import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple

def read_chgdiff(filename: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Read CHGDIFF.vasp file and return:
    - charge density as 3D numpy array
    - lattice vectors as 3x3 numpy array
    - grid size (nx, ny, nz)
    """
    with open(filename, 'r') as f:
        f.readline()  # Title
        scale = float(f.readline())  # Universal scale factor

        # Read lattice vectors and scale them
        lattice = np.array([list(map(float, f.readline().split())) for _ in range(3)])
        lattice *= scale

        elements = f.readline().split()
        natoms = list(map(int, f.readline().split()))
        total_atoms = sum(natoms)

        f.readline()  # Coordinate type: Direct/Cartesian
        for _ in range(total_atoms):
            f.readline()

        # Read grid dimensions
        while True:
            line = f.readline()
            if not line:
                raise ValueError("Grid dimensions not found")
            parts = line.split()
            if len(parts) == 3:
                try:
                    nx, ny, nz = map(int, parts)
                    break
                except ValueError:
                    continue

        # Read charge density data (1D list)
        data = []
        while True:
            line = f.readline()
            if not line or "augmentation" in line.lower():
                break
            data.extend(map(float, line.split()))

        expected = nx * ny * nz
        if len(data) != expected:
            print(f"Warning: Expected {expected} values, got {len(data)}")
            if len(data) < expected:
                data.extend([0.0] * (expected - len(data)))
            else:
                data = data[:expected]

        # Convert to 3D array in Fortran order (z changes fastest)
        charge = np.array(data).reshape((nx, ny, nz), order='F')

        return charge, lattice, np.array([nx, ny, nz])

def calculate_plane_average(charge: np.ndarray, lattice: np.ndarray, axis: int = 2) -> np.ndarray:
    """
    Planar averaging:
    f(z) = (1/A) ∫∫ ρ(x, y, z) dx dy
    Returns: averaged 1D profile.
    """
    # Calculate area of unit cell perpendicular to chosen axis
    if axis == 0:
        area = np.linalg.norm(np.cross(lattice[1], lattice[2]))
        avg = np.mean(charge, axis=(1, 2))  # yz-plane
    elif axis == 1:
        area = np.linalg.norm(np.cross(lattice[0], lattice[2]))
        avg = np.mean(charge, axis=(0, 2))  # xz-plane
    else:
        area = np.linalg.norm(np.cross(lattice[0], lattice[1]))
        avg = np.mean(charge, axis=(0, 1))  # xy-plane

    return avg

def calculate_charge_transfer(charge_diff: np.ndarray, lattice: np.ndarray,
                              grid: np.ndarray, axis: int = 2) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Calculate:
    - Plane-averaged charge density difference Δρ(z)
    - Cumulative charge transfer Q(z) = ∫ Δρ(z) dz
    """
    avg_diff = calculate_plane_average(charge_diff, lattice, axis)
    dz = np.linalg.norm(lattice[axis]) / grid[axis]  # Grid spacing along z
    z_positions = np.arange(grid[axis]) * dz
    cumulative = np.cumsum(avg_diff) * dz  # e

    return z_positions, avg_diff, cumulative

def main():
    try:
        charge_diff, lattice, grid = read_chgdiff('CHGDIFF.vasp')
        print(f"Read CHGDIFF.vasp successfully: grid = {grid.tolist()}")

        # Compute along z-axis
        z, avg_rho, cumulative_Q = calculate_charge_transfer(charge_diff, lattice, grid, axis=2)

        # Get the z-lattice parameter (z-direction lattice vector)
        z_lattice = np.linalg.norm(lattice[2])

        # Divide by the z-lattice parameter
        avg_rho /= z_lattice  # Plane-averaged charge density in e/Å²
        cumulative_Q /= z_lattice  # Cumulative charge transfer in e/Å

        # Save the revised results
        np.savetxt('revised_charge_transfer.dat', np.column_stack((z, avg_rho, cumulative_Q)),
                   header='Position(Ang)  Plane_Averaged_Rho(e/A^3)  Cumulative_Charge(e)',
                   fmt='%12.6f  %15.8e  %15.8e')
        print("Saved revised results to revised_charge_transfer.dat")

        print(f"Total charge transfer (after dividing by z-lattice): {cumulative_Q[-1]:.4f} e")
        print(f"Charge density range (after dividing by z-lattice): {np.min(avg_rho):.2f} to {np.max(avg_rho):.2f} e/Å²")

        # Plot with dual y-axis
        fig, ax1 = plt.subplots(figsize=(10, 5))
        ax2 = ax1.twinx()

        ax1.plot(z, avg_rho, 'b-', label='Plane-Averaged Charge Density')
        ax2.plot(z, cumulative_Q, 'r--', label='Cumulative Charge Transfer')

        ax1.set_xlabel('z (Å)')
        ax1.set_ylabel('Δρ (e/Å²)', color='b')
        ax2.set_ylabel('Q (e/Å)', color='r')

        ax1.tick_params(axis='y', labelcolor='b')
        ax2.tick_params(axis='y', labelcolor='r')

        fig.tight_layout()
        plt.title("Plane-Averaged Charge Density & Cumulative Charge Transfer (Divided by z-lattice)")
        plt.savefig("charge_transfer_dual_axis_z_lattice.png")
        plt.show()

    except Exception as e:
        print(f"Error: {str(e)}")

if __name__ == "__main__":
    main()


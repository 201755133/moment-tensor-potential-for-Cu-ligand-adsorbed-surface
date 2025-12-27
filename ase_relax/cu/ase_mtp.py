import os
import subprocess
import numpy as np
from ase.calculators.calculator import Calculator, all_changes


# ============================================================
# Write ASE Atoms → MLIP-3 CFG
# ============================================================
def write_cfg(atoms, filename, symbol_to_type):
    with open(filename, "w") as f:
        f.write("BEGIN_CFG\n")

        f.write("Size\n")
        f.write(f"{len(atoms)}\n")

        f.write("Supercell\n")
        for v in atoms.cell:
            f.write(f"{v[0]:.10f} {v[1]:.10f} {v[2]:.10f}\n")

        f.write("AtomData:  id type       cartes_x      cartes_y      cartes_z\n")
        for i, atom in enumerate(atoms):
            t = symbol_to_type[atom.symbol]
            x, y, z = atom.position
            f.write(
                f"{i+1:6d} {t:4d} "
                f"{x:14.6f} {y:14.6f} {z:14.6f}\n"
            )

        f.write("END_CFG\n")


# ============================================================
# Read MLIP-3 output CFG (robust to wrapped lines)
# ============================================================
def read_cfg_results(filename):
    energy = None
    stress = None
    natoms = None

    in_atomdata = False
    atom_numbers = []

    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line == "Size":
                natoms = int(next(f).strip())

            elif line == "Energy":
                energy = float(next(f).strip())

            elif line.startswith("PlusStress"):
                stress = list(map(float, next(f).split()))

            elif line.startswith("AtomData"):
                in_atomdata = True
                atom_numbers.clear()
                continue

            elif line == "END_CFG":
                in_atomdata = False

            elif in_atomdata:
                for token in line.split():
                    try:
                        atom_numbers.append(float(token))
                    except ValueError:
                        pass

    if natoms is None or energy is None:
        raise RuntimeError("Failed to parse MLIP output CFG")

    expected = natoms * 8
    if len(atom_numbers) < expected:
        raise RuntimeError(
            f"AtomData incomplete: got {len(atom_numbers)}, expected {expected}"
        )

    forces = np.empty((natoms, 3))
    for i in range(natoms):
        base = i * 8
        forces[i, 0] = atom_numbers[base + 5]
        forces[i, 1] = atom_numbers[base + 6]
        forces[i, 2] = atom_numbers[base + 7]

    return energy, forces, stress


# ============================================================
# ASE Calculator Wrapper for MLIP-3 (MPI-enabled, efficient)
# ============================================================
class ASE_MTP(Calculator):
    implemented_properties = ["energy", "forces", "stress"]

    def __init__(
        self,
        pot_mtp,
        symbols,
        mlp_cmd="mlp",
        mpi_cmd=None,
        workdir="mlip_work",
        **kwargs,
    ):
        """
        pot_mtp : path to pot.mtp
        symbols : list of elements IN TRAINING ORDER
        mlp_cmd : mlp executable
        mpi_cmd : e.g. "mpirun -n 192" (None = serial)
        workdir : persistent directory for MLIP I/O
        """
        super().__init__(**kwargs)

        self.pot_mtp = pot_mtp
        self.mlp_cmd = mlp_cmd
        self.mpi_cmd = mpi_cmd
        self.symbol_to_type = {s: i for i, s in enumerate(symbols)}

        self.workdir = workdir
        os.makedirs(self.workdir, exist_ok=True)

        self.in_cfg = os.path.join(self.workdir, "input.cfg")
        self.out_cfg = os.path.join(self.workdir, "output.cfg")

    def calculate(self, atoms=None, properties=["energy"],
                  system_changes=all_changes):

        super().calculate(atoms, properties, system_changes)

        # Write CFG
        write_cfg(atoms, self.in_cfg, self.symbol_to_type)

        # Build command
        if self.mpi_cmd:
            cmd = f"{self.mpi_cmd} {self.mlp_cmd} calc-efs {self.pot_mtp} {self.in_cfg} {self.out_cfg}"
        else:
            cmd = f"{self.mlp_cmd} calc-efs {self.pot_mtp} {self.in_cfg} {self.out_cfg}"

        ret = subprocess.call(cmd, shell=True)
        if ret != 0:
            raise RuntimeError("mlp calc-efs failed")

        # Read results
        energy, forces, stress = read_cfg_results(self.out_cfg)

        self.results["energy"] = energy
        self.results["forces"] = forces

        if stress is not None:
            self.results["stress"] = np.array(stress)


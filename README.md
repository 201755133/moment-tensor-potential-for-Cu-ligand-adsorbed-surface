# Active Learning with Moment Tensor Potentials (MTP)
Solid–Liquid Interface Simulations using MLIP-3, VASP, LAMMPS, ASE, and JARVIS-Tools
# Overview

This repository provides a complete active learning workflow for developing and applying Moment Tensor Potentials (MTP) to solid–liquid interfaces, with a focus on electrochemical and catalytic systems.
We employ MLIP-3 interfaced with VASP (for DFT reference calculations) and LAMMPS (for large-scale molecular dynamics). The trained potentials are further integrated with ASE for structure optimization and atomistic simulations.
The workflow is designed to be computationally efficient, scalable, and reproducible, making it suitable for long-time MD simulations and interface sampling.
# Key Features
Active learning–based MTP training using MLIP-3
Initial low-cost training with MTP level = 8
Progressive refinement to MTP level = 16 using accumulated configurations
Solid–liquid interface sampling, DFT reference calculations using VASP,  Large-scale MD using LAMMPS,  structure relaxation and optimization are done via ASE. 
Post-processing and analysis using custom Python codes built on JARVIS-tools
# Software and Tools
MLIP-3 (Moment Tensor Potentials) https://mlip.skoltech.ru/, 
LAMMPS (MD simulations) https://www.lammps.org/, 
VASP (DFT reference data) https://www.vasp.at/, 
and ASE (Atomic Simulation Environment) https://wiki.fysik.dtu.dk/ase/, 
# Analysis and Post-processing
JARVIS-tools https://github.com/usnistgov/jarvis, 
and Python (NumPy, SciPy, Matplotlib, ASE, JARVIS)
# Methodology
# 1. Active Learning Strategy
We adopt an active learning framework to efficiently generate high-quality training data for solid–liquid interfaces.
Start with a low-complexity MTP (level = 8), which reduces initial computational cost, enabling fast exploration of configuration space
, Perform MD simulations using LAMMPS + MTP, Identify extrapolative or high-uncertainty configurations, Compute reference energies, forces, and stresses using VASP.
Retrain the MTP with the expanded dataset, and increase model expressiveness to MTP level = 16, and the final potential is used for production simulations
# 2. Moment Tensor Potential Training
Training performed using MLIP-3, Descriptor level progression:Level 8 → initial exploration, Level 16 → final high-accuracy model. 
Our potential training data includes: Solid–liquid interface structures, thermalized MD snapshots, water bulk structures, Cu bulk structures obtained from materials projects, Cu(100) surface binding with ligand, and isolated ligand structures in a cubic box, and the final potential file is provided in the Supporting Information of the article:
“Fluorine-free CO₂ Electrolyzer” 
# Correspondence
 All inquiries related to the machine learning methodology, potential training, or analysis workflows should be addressed to: Ghulam Abbas (Corresponding Author) Email: hgabbas71@gmail.com 

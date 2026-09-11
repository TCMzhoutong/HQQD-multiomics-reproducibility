# Molecular docking and molecular dynamics

The docking workflow is driven by `01_inputs/docking_pairs.csv` and the scripts in `13_scripts/`. Set `OBABEL_EXECUTABLE`, `VINA_EXECUTABLE`, and `PYMOL_EXECUTABLE` when these programs are not available on `PATH`.

The active design contains 11 ligands, four targets, and 44 ligand-target combinations, with three Vina repetitions per combination. Molecular-dynamics inputs and compact analysis outputs are retained for four systems. Large trajectories and checkpoints are not bundled and are listed in `../metadata/EXTERNAL_DATA_MANIFEST.tsv`.

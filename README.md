# PhD Scripts

A collection of scripts written for use with different computational biology/chemistry packages throughout my PhD.

---

## General Use

### [pdb_download.sh](pdb_download.sh) (2023/10/11)

Download a list of PDB files using a text file containing one 4-character PDB
code on each line.

### [align_pdb.py](align_pdb.py) (2024/02/02)

Align a list of PDB files to a reference structure using [PyMOL](https://pymol.org/2/) (Schrödinger, Inc.).

## Molecular Dynamics (AMBER)

### [split.py](split.py) (2023/01/23)

Determine the number of monovalent salt ions are required for an explicit solvent molecular dynamics (MD) simulation using the [SPLIT method](http://archive.ambermd.org/202002/0194.html). Method created by Dr. Matias Machado (2020).

Designed for use with LEaP ([AmberTools](https://ambermd.org/tutorials/basic/tutorial9/index.php)).

### [amber2csm.py](amber2csm.py) (2024/01/15)

Re-format pdb files from MD simulations using AMBER and GLYCAM force fields via cpptraj (AmberTools). 3-letter sugar residue codes are renamed to match those found on the [Protein Data Bank](https://www.rcsb.org) (PDB) for compatibility with [CSM-carbohydrate](https://biosig.lab.uq.edu.au/csm_carbohydrate/), a webserver predicting protein-carbohydrate binding affinity ([Nguyen _et al._ 2022](https://doi.org/10.1093/bib/bbab512)).

### [fixatomnum.py](fixatomnum.py) (2024/01/19)

Fix the atom numbering (columns 7-11) for PDB files that have had new atoms added to them using pdb4amber ([AmberTools](https://ambermd.org/tutorials/basic/tutorial9/index.php); [Case _et al._ 2023](https://doi.org/10.1021/acs.jcim.3c01153)), or [MolProbity](http://molprobity.biochem.duke.edu) ([Williams _et al._ 2017](https://doi.org/10.1002/pro.3330)).

---

## Molecular Docking (AutoDock Vina/Vina-Carb)

### [adv_config.py](adv_config.py) (2023/09/29)

Generates a template config file for use with AutoDock Vina/Vina-Carb using the grid dimensions from an grid AutoDockTools grid dimensions file.

---

## Flow Chemistry

### [read_pressure.py](read_pressure.py) (2023/12/15)

Measure and record the pressure within the current Manchester [Future BRH](https://futurebrh.com) flow chemistry setup using _Fluigent_ Microfluidic In-Line Pressure Sensor (IPS). More information about the IPS can be found here: <https://www.fluigent.com/research/instruments/sensors/pressure-unit/>

Requires the [Fluigent SDK](https://github.com/Fluigent/fgt-SDK).

---

# PhD Scripts

A collection of scripts written for use with different computational biology/chemistry packages throughout my PhD.

---

## General

### [pdb_download.sh](pdb_download.sh) (2023/10/11)

Download a list of PDB files using a text file containing one 4-character PDB
code on each line.

#### To-Do

- [ ] Parallel downloads

---

## AMBER

### [split.py](split.py) (2023/01/23)

Determine the number of monovalent salt ions are required for an explicit solvent molecular dynamics (MD) simulation using the SPLIT method. Method created by Dr. Matias Machado (2020). A more detailed explanation can be found here: <http://archive.ambermd.org/202002/0194.html>

Designed for use with LEaP (AmberTools).

---

## AutoDock Vina/Vina-Carb

### [adv_config.py](adv_config.py) (2023/09/29)

Generates a template config file for use with AutoDock Vina/Vina-Carb using the grid dimensions from an grid AutoDockTools grid dimensions file.

---

## Flow Chemistry

### [read_pressure.py](read_pressure.py) (2023/12/15)

Measure and record the pressure within the current Manchester [Future BRH](https://futurebrh.com) flow chemistry setup using _Fluigent_ Microfluidic In-Line Pressure Sensor (IPS). More information about the IPS can be found here: <https://www.fluigent.com/research/instruments/sensors/pressure-unit/>

Requires the Fluigent SDK: <https://github.com/Fluigent/fgt-SDK>

---

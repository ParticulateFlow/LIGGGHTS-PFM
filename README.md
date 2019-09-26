# LIGGGHTS

LIGGGHTS® - LAMMPS Improved for General Granular and Granular Heat Transfer Simulations - is a discrete element method (DEM) particle simulation software.
LIGGGHTS® is part of the [CFDEM®project](https://www.cfdem.com) and is based on the molecular dynamics simulation code [LAMMPS](https://lammps.sandia.gov/).

[![CircleCI](https://circleci.com/gh/ParticulateFlow/LIGGGHTS.svg?style=shield&circle-token=8905cdbf813717ce628dd05a454d0f7581110907)](https://circleci.com/gh/ParticulateFlow/LIGGGHTS)
[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

## Disclaimer

> This is an academic adaptation of the LIGGGHTS® software package, released by the
[Department of Particulate Flow Modelling at Johannes Kepler University in Linz, Austria.](https://www.jku.at/pfm)
> LIGGGHTS® and CFDEM® are registered trademarks, and this offering is not approved or
endorsed by DCS Computing GmbH, the official producer of the LIGGGHTS® and CFDEM®coupling software.

## Installation

This is a short summary of how to install LIGGGHTS on Linux. A more comprehensive guide can be found in the documentation.

### Install prerequisites

```bash
sudo apt-get install build-essential cmake openmpi-bin libopenmpi-dev python-dev
```

We recommend installing LIGGGHTS to a directory named `CFDEM`, especially when used with CFDEMcoupling.

```bash
cd
mkdir -p CFDEM
cd CFDEM
```

Clone or download the LIGGGHTS source from the repository.

### Build LIGGGHTS with CMake

```bash
cd LIGGGHTS
mkdir -p src-build
cd src-build
cmake ../src/
make
```

### Build LIGGGHTS with make

```bash
cd LIGGGHTS
mkdir -p src-build
cd src
make fedora
cp lmp_fedora ../src-build/liggghts
```

### Add an alias

You may want to create a permanent alias for the executable.

```bash
gedit ~/.bashrc &
alias liggghts='~/CFDEM/LIGGGHTS/src-build/liggghts'
source ~/.bashrc
```

## Getting Started

Navigate to the tutorials folder to run the chute_wear example case

```bash
cd ~/CFDEM/LIGGGHTS/examples/LIGGGHTS/Tutorials_public/chute_wear
```

Start the simulation by typing

```bash
liggghts -in in.chute_wear
```

## License

[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

- This software is distributed under the [GNU General Public License](https://opensource.org/licenses/GPL-2.0).
- Copyright © 2009-     JKU Linz
- Copyright © 2012-2015 DCS Computing GmbH, Linz
- Copyright © 2003      Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain rights in this software.


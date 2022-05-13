# Chemical Kinetcs Mechansms Reduction Framework

This project presents a framework under development to produce reduced chemical kinetics mechanisms from canonical models, more specifically for thermal autoignition (IDT) and laminar flame speed (LFS).

The project is developed within Labcet - UFSC.

## Current Features

This project is a Work In Progress (WIP). Currently, the following functionalities are implemented, albeit test cases are required for the entire project.

- IDT and LFS interface reactor models;
- DRG and DRGEP reduction methods;

## Guides

### IDT and LFS reactor models

The `reactors.py` under src provide an easy to use interface of both 0D, constant mass, constant pressure, adiabatic reactor and a free flame model (adapted from cantera). The implementation here focus on providing a similar way of running the cannonical models with minimun user interference, producing data ready to be used for the reduction.

A simple example of how to use those reactors can be found on `run_reactors.py` script.

### DRG and DRGEP reduction methods

The reduction implementation can be found on `reduction.py` under src. An example on how to use can be found on the script `main_app.py`.

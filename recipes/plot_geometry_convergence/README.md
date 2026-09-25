#cclib_plot_geometry_convergence.py
Given a computational chemistry logfile for a geometry optimization, plot how the convergence parameters change over the course of the optimization.

```
usage: plot_geometry_convergence.py [-h]
                                    [--scaling-energy-change SCALING_ENERGY_CHANGE]
                                    compchemfilename [compchemfilename ...]

positional arguments:
  compchemfilename

options:
  -h, --help            show this help message and exit
  --scaling-energy-change SCALING_ENERGY_CHANGE
                        Factor to scale the energy change by
```

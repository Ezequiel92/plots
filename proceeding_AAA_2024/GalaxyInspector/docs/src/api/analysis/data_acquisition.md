# Data acquisition functions

These functions are used internally and none are exported. 

These functions depend on the particulars of the simulation code (e.g. units).

These function read the different output files from the simulation (snapshots, FoF catalogs, sfr.txt, etc.), and load the data into memory as dictionaries.

---

```@autodocs
Modules = [GalaxyInspector]
Pages   = ["analysis/data_acquisition.jl"]
```

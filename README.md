# Random-Walk-Simulation
Simulation designed to probe collective decision making in ant colonies by modelling ants in a nest as particles in a box. Results are compared to experimental data to see what real life behavior can and cannot be explained by this model.

- Code is broken into two different geometries: Particles in a Box and Particles in a Disk.
- The bulk of the simulation is run via the c++ files (main.cpp, functions.cpp, & headers.cpp), while the animation is run through matplotlib in the jupyter notebook.The text files are OUTPUTS from the c++ code that are read in to the Jupyter notebook. To change parameters, edit "main.cpp", not "params.txt".
- Typically, elastic collisions have no effect on the simulation being run, but they are accessible via a flag in the main file.

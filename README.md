# Random Walk Simulation
![alt text](https://github.com/jakehanson/Random-Walk-Simulation/blob/new/GIFS/timing_visit_example.gif)

This is a general purpose brownian motion / random walk simulator designed to explore potential similarities between ants in their nests and particles in a disk. Collisions are typically handled using conservation of momentum and the 'hard body sphere' approximation. The simulation is event based which means timestep is dictated by the time until the next event (either collisions between particles or a particle and the wall). This avoids potential pitfalls due to simultaneous events and speeds up computation time, but one should be aware that visualizations of the data will either have to be smoothed in time or appear jumpy due to the non-uniform timestep. More information about this type of simulation can be found here: http://introcs.cs.princeton.edu/java/assignments/collisions.html.

## General Overview

The code has three primary files:
* **main.cpp** - sets up and runs simulation

* **functions.cpp** - contains the functions required by main.cpp

* **header.h** - standard header for function prototypes and libraries

There are also two subdirectories:
* **GIFS** - example animations from simulation output

* **MODELS** - various models that are derived from, but differ, from the source code


## Download and Run
To run the code:

* **Download**: the primary files main.cpp, functions.cpp, header.h (or the entire directory)
* **Compile:** the code must be compiled to make an executable
  * g++ -std=c++11 main.cpp functions.cpp -o run_sim.exe
* **Run Option 1**: the code can be run from the command line without any arguments
  * ./run_sim.exe
* **Run Option 2:** the code can also be running by accepting the number of ants as an integer argument
  *  ./run_sim.exe n_ants

As a first check that everything is working, ./run_sim.exe should produce to output files 'params.txt' and 'output.txt'.

## Initializing Parameters

Once you have the source code in your local directory you can modify the simulation parameters by opening main.cpp in a text editor and changing the variables and flags at the top.

The code is designed to model ants in a nest as particles in a disk. This being said, there are certain flags that violate the laws of physics in order to test different assumptions about ant behavior. The following descriptions provide information about the differences 

#### Conservation of Momentum vs Fixed Velocity
For most molecular dynamics, you will want collisions to conserve momentum. However, for simulating ants a more realistic assumption is that each ant maintains a fixed velocity. The split between these two options is controlled with the 'fixed_velo' flag. If false, the simulation will obey conservation of momentum, if true, particles reflect specularly but keep their same velocity at all times.

#### Starting Ant in Aperture Center
For ant simulations, we are usually interested in how long an ant spends inside a nest of a given geometry. Therefore there is a flag called 'start_in_center' to start an ant in the center of an aperture and a flag called 'exit_flag' to end the simulation when this ant leaves, thereby timing a single ants visit. For molecular dynamics, these should both be set to false to prevent biasing results.

#### Other parameters

* **num_ants** = number of ants if not passed as command line arg after ./run_sim.exe
* **R** = nest radius [units of length]
* **a** = aperture size [units of length]
* **velo** = starting velocity [units of speed]
* **r_enc** = radius for ants
* **max_init** = maximum number of tries to randomly initialize ants
* **max_steps** = maximum number of events before terminating simulation
* **single_particle** = set this to true if only simulating one particle

## Output

The primary output of running the simulation is a file called output.txt that contains the position and velocity of each ant as a function of time. In addition, the files contains the time the ant exited the nest as well a running number of collisions it has undergone.

The secondary output is params.txt which basically rewrites the parameters used to run the simulation. This file is especially relevant if simulation parameters need to be read into a plotting device.

## Animation
There is a jupyter notebook called Event Based Animation.ipynb that was used to generate animation gifs. While github can render the jupyter notebook in the browser, it must be run locally to create the animation. If you are unfamiliar with jupyter notebooks it is probably easier to read output.txt into your favorite plotting software.

## FEEL FREE TO EMAIL ANY QUESTIONS! 

jake.hanson@asu.edu

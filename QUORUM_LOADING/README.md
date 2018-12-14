# QUORUM LOADING

This directory examines quorum loading, which is essentially the process by which nests build up to quorum. We assume a fixed entry rate, and allow free outflow. We test the hypothesis that certain nest geometries preferentially hold recruiters, thereby allowing them to reach quorum.

Note, this branch is under active development. Several updates have been implemented, including careful handling of overlap upon nest entry and other resolutions to small glitches.

## Syntax

To compile, download this directory and enter the following:

g++ -std=c++11 -03 main.cpp functions.cpp -o run_sim.exe

To run, the syntax is either:

./run_sim.exe (this runs with the default params)

or 

./run_sim entry_rate max_ants sim_time

where entry_rate (ants/second), max_ants, and sim_time (seconds) are floating point values separated by a space.

jake.hanson@asu.edu

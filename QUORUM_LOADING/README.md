# QUORUM LOADING

This directory looks at how the encounter rate changes as a function of time for nests of varying quality. The hypothesis is that quorum is sensed via encounter rate, thus, one would expect that high quality nests are able to reach and maintain a high encounter rate while low quality nests cannot. We investigate this hypothesis within the framework of the random walk model by supplying two nests of different quality (different aperture size) with a constant influx of ants. The idea is that the good nest should naturally "load" and reach quorum while the mediocre nest should not.

## Syntax

To compile, download this directory and enter the following:

g++ -std=c++11 -03 main.cpp functions.cpp -o run_sim.exe

To run, the syntax is either:

./run_sim.exe (this runs with the default params)

or 

./run_sim entry_rate max_ants sim_time

where entry_rate (ants/second), max_ants, and sim_time (seconds) are all floating point values separated by a space.


jake.hanson@asu.edu

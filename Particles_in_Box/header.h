#pragma once
#include <vector>
#include <iostream>
#include <random>

/* Structure to hold ant colony */
struct Ants
{
	std::vector<int> ant_name;  // Vector to hold x-positions
	std::vector<double> x_positions;  // Vector to hold x-positions
	std::vector<double> y_positions;  // Vector to hold y-positions
	std::vector<double> x_velocities;  // Vector to hold x-velocites
	std::vector<double> y_velocities;  // Vector to hold y-velocites
	std::vector<double> event_time;  // Vector to hold timestep
	std::vector<double> exit_time;  // Vector to hold time of exit
	std::vector<double> collisions;  // Vector to hold the number of collisions
	std::vector<bool> nest_flag;  // Vector to hold nest flag (in nest == True)

	//Ants(size_t num_ants,double temp);  // signature for constructor. construct instance of class
	Ants(int num_ants);  // signature for constructor. construct instance of class
	void populate(double L, double H, double T,double r_enc,int max_init);  // populate function without passing a seed
	void populate(double L, double H, double T,double r_enc,int max_init, int seed);  // populate function with a seed
	void populate(double L, double H, double T,double r_enc,int max_init, std::random_device &rd);  // populate function with reference to rd

	/* Function to update particle locations */
	void update(double L, double H,double a,double t_min,int index1); // version 1 is for two particles colliding
	void update(double t_min,double r_enc,int index2,int index3); // version 2 is for a particle colliding with wall

};

/* Operator overload to write ant colony to output stream */
std::ostream &operator<<(std::ostream &out, Ants const &ants);

/* Function to return the time interval until the next event */
double get_t_wall(Ants ants,double L,double H,double r_enc, int &index1);
double get_t_ant(Ants ants,double r_enc,double t_wall, int &index2, int &index3);




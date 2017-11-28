#pragma once
#include <vector>
#include <iostream>
#include <random>

/* Structure to hold ant colony */
struct Ants
{
	std::vector<int> ant_name;  // Vector to hold x-positions
	std::vector<long double> x_positions;  // Vector to hold x-positions
	std::vector<long double> y_positions;  // Vector to hold y-positions
	std::vector<long double> x_velocities;  // Vector to hold x-velocites
	std::vector<long double> y_velocities;  // Vector to hold y-velocites
	std::vector<long double> event_time;  // Vector to hold timestep
	std::vector<long double> exit_time;  // Vector to hold time of exit
	std::vector<long double> collisions;  // Vector to hold the number of collisions
	std::vector<bool> nest_flag;  // Vector to hold nest flag (in nest == True)

	//Ants(size_t num_ants,long double temp);  // signature for constructor. construct instance of class
	Ants(int num_ants);  // signature for constructor. construct instance of class
	void populate(long double R, long double T,long double r_enc,int max_init,bool start_in_center);  // populate function without passing a seed
	void populate(long double R, long double T,long double r_enc,int max_init,int seed,bool start_in_center);  // populate function with a seed
	void populate(long double R, long double T,long double r_enc,int max_init, std::random_device &rd,bool start_in_center);  // populate function with reference to rd
	void populate(long double R, long double a, long double T, long double r_enc); // populate function for a single ant at top-center of nest
	
	/* Function to update particle locations */
	void update(long double R,long double r_enc,long double a,long double t_min,int index1); // version 1 is ant-to-ant collision w/ cons of momentum
	void update(long double t_min,long double r_enc,int index1,int index2,long double velo); // version 1 is ant-to-ant collision w/ fixed velo
	void update(long double t_min,long double r_enc,int index2,int index3); // version 3 is for a particle colliding with wall

};

/* Operator overload to write ant colony to output stream */
std::ostream &operator<<(std::ostream &out, Ants const &ants);

/* Function to return the time interval until the next event */
long double get_t_wall(Ants ants,long double R,long double r_enc, int &index1);
long double get_t_ant(Ants ants,long double r_enc,long double t_wall, int &index2, int &index3);




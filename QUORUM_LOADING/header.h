#pragma once
#include <vector>
#include <iostream>
#include <random>
#include <iomanip>      // std::setprecision

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
	void populate(void);
	int enter(long double R, long double a, long double velo, long double r_enc,double entry_rate);

	/* Function to update particle locations */
	void update(long double t_entry,long double R, long double r_enc); // ant entry
	void update(long double R,long double r_enc,long double a,long double t_min,int index1,long double machine_tol); // ant-to-wall collision
	void update(long double t_min,long double r_enc,int index1,int index2,long double velo); // ant-to-ant collision

};

/* Operator overload to write ant colony to output stream */
std::ostream &operator<<(std::ostream &out, Ants const &ants);

/* Function to return the time interval until the next event */
long double get_t_wall(Ants ants,long double R,long double r_enc,int &index1,long double t_entry,long double machine_tol);
long double get_t_ant(Ants ants,long double r_enc,long double t_wall, int &index2, int &index3,long double machine_tol);




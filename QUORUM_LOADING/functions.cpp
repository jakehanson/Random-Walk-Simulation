#include "header.h"
#include <stdio.h>
#include <cmath>
#include <stdexcept>

/* Define function to overload the operator << such that if the input is std:ostream &out and &ants it prints as follows*/
std::ostream &operator<<(std::ostream &out, Ants const &ants)
{
	for (size_t i=0; i< ants.x_positions.size();++i){
		out << ants.ant_name[i] << "\t";
		out << ants.x_positions[i] << "\t";
		out << ants.y_positions[i] << "\t";
		out << ants.x_velocities[i] << "\t";
		out << ants.y_velocities[i] << "\t";
		out << ants.event_time[i] << "\t";
		out << ants.exit_time[i] << "\t";
		out << ants.collisions[i] << "\t";
		out << ants.nest_flag[i] << std::endl;
	}
	return out;
}

/* Define our constructor. Pass Args to this to create an instance of Ant Class */
Ants::Ants(int num_ants)
{
	ant_name = std::vector<int>(num_ants);
	x_positions = std::vector<long double>(num_ants);
	y_positions = std::vector<long double>(num_ants);
	x_velocities = std::vector<long double>(num_ants);
	y_velocities = std::vector<long double>(num_ants);
	event_time = std::vector<long double>(num_ants);
	exit_time = std::vector<long double>(num_ants);
	collisions = std::vector<long double>(num_ants);
	nest_flag = std::vector<bool>(num_ants);
}

/* Define Method to Populate our arrays within Ant Class */
void Ants::populate(void)
{
	for(int i=0;i<x_positions.size();i++){
		ant_name.at(i) = i;
		x_positions.at(i) = 0.;
		y_positions.at(i) = 0.;
		x_velocities.at(i) = 0.;
		y_velocities.at(i) = 0.;
		event_time.at(i) = 0.;
		exit_time.at(i) = -1;
		collisions.at(i) = 0.;
		nest_flag.at(i) = false;
	}
}


/* Define method for adding an ant in the center of aperture */
int Ants::enter(long double R, long double a, long double velo, long double r_enc,double entry_rate)
{
	// First figure out which ant should enter
	int entry_ant = -1;
	for(int i=0;i<x_positions.size();i++){
		if(nest_flag.at(i) == false and exit_time.at(i) == -1){ // ant must not be in nest nor have previously left nest
			entry_ant = i;
			break;
		}
	}

	// If all ants have already entered, return -1
	if(entry_ant == -1){
		return entry_ant;
	}
	// Otherwise, initialize new ant entering from top center
	else{
		long double entry_time;
		entry_time = (entry_ant+1)/entry_rate;  // Get time of entry using the entry number
		
		// Init Positions
		x_positions.at(entry_ant) = 0.0; // centered
		y_positions.at(entry_ant) = R-r_enc; // located at top

		// Random Initial Velocity
		std::random_device rd; // Init random device
		long double seed = rd(); // random seed
		//long double seed = 500; // seed with specific value (creates reproducible simulation)
		std::mt19937 gen(seed); // mersenne-twister generator with seed
		std::uniform_real_distribution<long double> uni(0.,1); // uniform dist on (0,1]
		
		long double clearance = a/(4*R); // angle needed to clear edge
		long double velo_angle = (M_PI+clearance)+uni(gen)*(M_PI-2*clearance); // angle between pi+clearance and 2pi-clearance (downward)
		x_velocities.at(entry_ant) = velo*cos(velo_angle);
		y_velocities.at(entry_ant) = velo*sin(velo_angle);

		// Init all other parameters
		ant_name.at(entry_ant) = entry_ant;
		event_time.at(entry_ant) = entry_time;
		collisions.at(entry_ant) = 0.0;
		exit_time.at(entry_ant) = -1.;
		nest_flag.at(entry_ant) = true;
		return entry_ant;		
	}
}


/* Calculate time before a wall-ant collision */
long double get_t_wall(Ants ants,long double R,long double r_enc,int &index1,long double t_entry,long double machine_tol){
	int num_ants = ants.x_positions.size();
	long double x,y,v_x,v_y;  // variables used to find collision time
	long double t_min;
	long double t_collide;
	bool init = false;

	long double x_test,y_test;

	t_min = t_entry; // set minimum time before event to be the next entry time

	//Find t_min for all particles
	for(int i=0; i<num_ants; ++i){
		if(ants.nest_flag[i]){
			//std::cout << "CALCULTING FOR ANT " << i << std::endl;
			long double x = ants.x_positions[i];
			long double y = ants.y_positions[i];
			long double v_x = ants.x_velocities[i];
			long double v_y = ants.y_velocities[i];
			
			t_collide = (-1*(x*v_x+y*v_y)+
					sqrt(pow(x*v_x+y*v_y,2)-(pow(v_x,2)+pow(v_y,2))*(pow(x,2)+pow(y,2)-pow(R-r_enc,2))))
					/(pow(v_x,2)+pow(v_y,2));

			if(t_collide < t_min){
				t_min = t_collide;
				index1 = i;
			}
		}		
	}

	// If we never beat t_entry, just return a time greater than t_entry
	if(t_min == t_entry){
		return 2*t_entry;
	}
	// Otherwise, make sure we don't overstep due to numerical imprecision
	else{
		int counter = 0;
		int max_tries = 10;
		bool success_flag = false;
		long double x_step;
		long double y_step;
		long double distance;		
		while(counter < max_tries){
			x_step = ants.x_positions.at(index1)+ants.x_velocities.at(index1)*t_min;
			y_step = ants.y_positions.at(index1)+ants.y_velocities.at(index1)*t_min;
			distance = std::sqrt(x_step*x_step+y_step*y_step);
			// Check overstep
			if(distance > (R-r_enc)){
				//std::cout << "STEPPED OUTSIDE DOMAIN! ATTEMPTING TO SHRINK TIMESTEP BY MACHINE TOL..." << std::endl;
				t_min -= machine_tol; // reduce timestep by machine precision
				counter += 1; // we will try again
			}else{
				success_flag = true;
				break;		
			}
		}
		if(success_flag){
			// if(counter != 0){
			// 	std::cout << "-->ERROR RESOLVED IN " << counter << " ATTEMPT(S)." << std::endl;
			// }
			return t_min;
		}else{
			std::cout << "ATTEMPTS  = " << counter << std::endl;
			std::cout << "FINAL DISTANCE = " << std::setprecision(20) << distance << std::endl;
			std::cout << "R - r_enc  = " << std::setprecision(20) << R-r_enc << std::endl;
			throw std::runtime_error("UNABLE TO CALCULATE WALL TIME!");
			return -1; // fail
		}
	}
}


/* Calculate the time before the next ant-ant collision */
long double get_t_ant(Ants ants,long double r_enc,long double t_wall, int &index2, int &index3,long double machine_tol){
	int num_ants = ants.x_positions.size();
	long double sigma = 2*r_enc; // cross sectional distance
	long double v_dot_r;
	long double v_dot_v;
	long double r_dot_r;
	long double d;
	long double t_collide;

	long double t_min = 2.*t_wall; //Initialize minimum time as something more than wall collision time
	// Note: Do not change this to t_min = t_wall! Numerical imprecision will MESS SHIZ UP!!

	// For each ant in nest get the collision time to every other ant
	for(int i=0;i<num_ants;++i){
		if(ants.nest_flag[i]){
			for(int j=i+1;j<num_ants;++j){
				if(ants.nest_flag[j]){

					v_dot_r = (ants.x_velocities[j]-ants.x_velocities[i])
								*(ants.x_positions[j]-ants.x_positions[i])
							 	+ (ants.y_velocities[j]-ants.y_velocities[i])
							 	*(ants.y_positions[j]-ants.y_positions[i]);
					v_dot_v = pow(ants.x_velocities[j]-ants.x_velocities[i],2)
								+ pow(ants.y_velocities[j]-ants.y_velocities[i],2);
					r_dot_r = pow(ants.x_positions[j]-ants.x_positions[i],2)
								+pow(ants.y_positions[j]-ants.y_positions[i],2);
					d = pow(v_dot_r,2)-v_dot_v*(r_dot_r-pow(sigma,2));

					// If they collide, get time at which it occurs
					if(v_dot_r < 0. and d > 0. and r_dot_r > pow(sigma,2)){
						
						t_collide = -1.*(v_dot_r+sqrt(d))/v_dot_v;

						// Check for numerical errors (Just discard these runs for now...)
						if(d < 0){
							throw std::runtime_error("DISTANCE VIOLATION IN COLLISION CALCULATION.");
						}
						if (t_collide < 0.0){
							throw std::runtime_error("NEGATIVE TIME CALCULATED IN ANT COLLISION");
						}

						// Check if collision is next event
						if(t_collide < t_min){
							t_min = t_collide;
							index2 = i;
							index3 = j;
						}
					}
				}
			}
		}
	}

	// If we never beat t_wall, just return t_wall as the minimum time
	if(t_min >= t_wall){
		return t_min;
	}
	// Else we need to make sure that we don't overstep and have particles lock together
	else{
		int counter = 0;
		int max_tries = 10; // tries to resolve issue
		bool success_flag = false;
		long double x1;
		long double y1;
		long double x2;
		long double y2;
		long double separation;

		while(counter < max_tries){
			x1 = ants.x_positions.at(index2)+ants.x_velocities.at(index2)*t_min;
			y1 = ants.y_positions.at(index2)+ants.y_velocities.at(index2)*t_min;
			x2 = ants.x_positions.at(index3)+ants.x_velocities.at(index3)*t_min;
			y2 = ants.y_positions.at(index3)+ants.y_velocities.at(index3)*t_min;
			separation = std::sqrt(pow(x1-x2,2)+pow(y1-y2,2));
			// Check overstep
			if(separation < 2.*r_enc){
				// std::cout << "OVERLAP ISSUE!! ATTEMPTING TO SHRINK TIMESTEP BY MACHINE TOL..." << std::endl;
				t_min -= machine_tol; // reduce timestep by machine precision
				counter += 1; // we will try again
			}else{
				success_flag = true;
				break;		
			}
		}

		if(success_flag){
			// if(counter != 0){
			// 	std::cout << "--> ERROR RESOLVED IN " << counter << " TIMESTEPS." << std::endl;
			// }
			return t_min;
		}else{
			throw std::runtime_error("UNABLE TO RESOLVE ANT-TO-ANT COLLISION!");
			return -1; // fail	
		}
	}
}


/* This update method is for a new ant entering */
void Ants::update(long double t_entry,long double R, long double r_enc){
	long double distance;
	/* Update Particle Locations and velocities*/
	for(int i=0;i<x_positions.size();++i){
		if(nest_flag.at(i)){
			x_positions.at(i) = x_positions.at(i)+x_velocities.at(i)*t_entry;
			y_positions.at(i) = y_positions.at(i)+y_velocities.at(i)*t_entry;
			event_time.at(i) = event_time.at(i)+t_entry;

			// Make sure we won't initalize the new particle on top of an existing particle
			if(std::abs(x_positions.at(i))<r_enc and std::abs(y_positions.at(i)-(R-r_enc))<r_enc){
				throw std::runtime_error("ERROR! OVERLAP IN ANT ENTRY.");
			}

		}
	}
}


/* This update method is for wall collisions*/
void Ants::update(long double R,long double r_enc, long double a,long double t_min,int index1){
	int num_ants = x_positions.size();
	long double angle1; // angle to the location where the particle hits the wall
	long double angle2; // angle to the initial position of the particle
	long double final_angle; // angle of where the particle would be t_min after colliding with wall


	/* Update Particle Locations and velocities*/
	for(int i=0;i<num_ants;++i){

		// The particle that is hitting the wall switches its velocity component
		if(i == index1){
			long double x = x_positions[i];
			long double y = y_positions[i];
			long double v_x = x_velocities[i];
			long double v_y = y_velocities[i];
			long double x_new = x_positions[i]+x_velocities[i]*t_min;
			long double y_new = y_positions[i]+y_velocities[i]*t_min;
			long double x_temp,y_temp;

			// Reflect particle off wall and store results
			if(((x >= 0 and y >= 0) || (x >= 0 and y < 0))){
				angle1 = atan(y_new/x_new);
				angle2 = atan(y/x);
			}else{
				angle1 = M_PI+atan(y_new/x_new);
				angle2 = M_PI+atan(y/x);
			}
			final_angle = 2*angle1-angle2;
			x_temp = sqrt(pow(x,2)+pow(y,2))*cos(final_angle);  //where the particle would be t_min after collision
			y_temp = sqrt(pow(x,2)+pow(y,2))*sin(final_angle);  //where the particle would be t_min after collision

			x_velocities[i] = (x_temp-x_new)/t_min;
			y_velocities[i] = (y_temp-y_new)/t_min;
			x_positions[i] = x_new;
			y_positions[i] = y_new;
			event_time[i] = event_time[i]+t_min;				

			//Did we make it out?
			if(x_positions[i] >= 0 and y_positions[i] >= 0){
				//std::cout << "Where are we? " << x_new << y_new << sqrt(pow(x,2)+pow(y,2)) << std::endl;
				angle1 = atan(x_positions[i]/y_positions[i]);
				if(angle1 <= a/(2*R) ){
					//std::cout << "Where are we? " << x_new << "\t" << y_new << "\t" << sqrt(pow(x,2)+pow(y,2)) << std::endl;
					nest_flag[i] = false;
					exit_time[i] = event_time[i];
				}
			}
			if(x_positions[i] < 0 and y_positions[i] >= 0){
				//std::cout << "Where are we? " << x_new << y_new << sqrt(pow(x,2)+pow(y,2)) << std::endl;
				angle1 = atan(-1*x_positions[i]/y_positions[i]);
				if(angle1 < a/(2*R) ){
					//std::cout << "Where are we? " << x_new << "\t" << y_new << "\t" << sqrt(pow(x,2)+pow(y,2)) << std::endl;
					nest_flag[i] = false;
					exit_time[i] = event_time[i];
				}
			}

		// update particles not colliding with wall
		}else{
			if(nest_flag.at(i)){
				x_positions[i] = x_positions[i]+x_velocities[i]*t_min;
				y_positions[i] = y_positions[i]+y_velocities[i]*t_min;
				event_time[i] = event_time[i]+t_min;				
			}
		}
	}
}


/* This update method is for ant-to-ant collisions with fixed velocity*/
void Ants::update(long double t_min,long double r_enc,int index1,int index2,long double velo){
	int num_ants = x_positions.size();
	long double sigma = 2*r_enc;

	/* Update all particle locations */
	for(int i=0;i<num_ants;++i){
		if(nest_flag.at(i)){
			x_positions[i] = x_positions[i]+t_min*x_velocities[i];
			y_positions[i] = y_positions[i]+t_min*y_velocities[i];
			event_time[i] = event_time[i]+t_min;
		}
	}

	// Find the two particles that are colliding and update their velocities
	// Elastic, frictionless collision
	for (int j=0;j<num_ants;++j){
		if (j == index1){
			for (int k=0;k<num_ants;++k){
				if (k == index2){
					long double del_x, del_y, del_vx, del_vy;
					long double J, J_x, J_y;
					long double v_dot_r;
					long double vx1,vy1,vx2,vy2;

					// Initialize Variables
					del_x = x_positions[k]-x_positions[j];
					del_y = y_positions[k]-y_positions[j];
					del_vx = x_velocities[k]-x_velocities[j];
					del_vy = y_velocities[k]-y_velocities[j];
					v_dot_r = del_vx*del_x+del_vy*del_y;

					// Define impulse for hard body spheres (assuming unit mass)
					J = v_dot_r/sigma;
					J_x = J*del_x/sigma;
					J_y = J*del_y/sigma;

					// Use impulse to get new heading
					vx1 = x_velocities[j]+J_x;
					vy1 = y_velocities[j]+J_y;
					vx2 = x_velocities[k]-J_x;
					vy2 = y_velocities[k]-J_y;

					// Fix velocity and use heading for direction
					x_velocities[j] = velo*vx1/sqrt(pow(vx1,2)+pow(vy1,2));
					y_velocities[j] = velo*vy1/sqrt(pow(vx1,2)+pow(vy1,2));
					x_velocities[k] = velo*vx2/sqrt(pow(vx2,2)+pow(vy2,2));
					y_velocities[k] = velo*vy2/sqrt(pow(vx2,2)+pow(vy2,2));

					// Tick collision counter
					collisions[j]++;
					collisions[k]++;
				}
			}
		}
	}
}



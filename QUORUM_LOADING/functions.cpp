#include "header.h"
#include <stdio.h>
#include <cmath>
#include <stdexcept>

// THIS IS THE NEW COLLISION RESOLUTION!!

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


/* Define Method to check for overlap */
int Ants::check(double r_enc)
{
	double d;
	double x1,y1,x2,y2;

	// For each ant, check if there is overlap
	for(int i=0;i<x_positions.size();i++){
		if(nest_flag.at(i) == true){
			x1 = x_positions.at(i);
			y1 = y_positions.at(i);
			for(int j=i+1;j<x_positions.size();j++){
				if(nest_flag.at(j) == true){
					x2 = x_positions.at(j);
					y2 = y_positions.at(j);
					d = std::sqrt(pow(x2-x1,2)+pow(y2-y1,2));
					// check for penetration, allowing for numerical error
					if(d < 2.*r_enc){
						//std::cout << "INVALID d = " << d << std::endl;
						return 0;
					}
				}
			}

		}
	}
	return 1; // if there are no problems, return 1
}


/* Define method for adding an ant somewhere in the aperture */
int Ants::enter(long double R, long double a, long double velo, long double r_enc,double entry_rate,long double machine_tol)
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
	// Otherwise, initialize new ant entering from somewhere in the aperture
	else{
		long double entry_time;
		int n_tries = 50; // 50 tries to enter the nest
		int entry_flag = 0; // flag for trying to enter
		entry_time = (entry_ant+1)/entry_rate;  // Get time of entry using the entry number

		// Random Initial Position and velocity
		std::random_device rd; // Init random device
		long double seed = rd(); // random seed
		//long double seed = 2; // seed with specific value (creates reproducible simulation)
		std::mt19937 gen(seed); // mersenne-twister generator with seed
		std::uniform_real_distribution<long double> uni(0.,1); // uniform dist on (0,1]

		// Try to enter
		long double x,y;
		long double entry_size = (a-2*r_enc)/R; // angular size of aperture given the size of the ant
		long double entry_angle = (M_PI/2-entry_size/2.)+uni(gen)*(entry_size);  // random angle in aperture
		long double x1,y1,x2,y2; // clearance positions on left and right
		long double v_x,v_y; // projected velocities
		long double alpha,beta; // minimum angles needed to reach clearance positions
		long double launch_angle;
		long double t_cross; // time to crossover into nest

		for(int i=0;i<n_tries;i++){
			if(entry_flag == 0){
				// Initialize Position
				x = R*cos(entry_angle); 
				y = R*sin(entry_angle);

				// Initialize starting velocity
				x1 = -(R-r_enc)*sin(a/(2.*R));
				y1 = (R-r_enc)*cos(a/(2.*R));
				x2 = -x1;
				y2 = y1;
				alpha = atan((y-y1)/std::abs(x1-x));
				beta = atan((y-y2)/(x2-x));
				launch_angle = M_PI+alpha+uni(gen)*(M_PI-(alpha+beta)); // randomly puts us in [pi+alpha,2pi-beta]
				v_x = velo*cos(launch_angle);
				v_y = velo*sin(launch_angle);

				// Enter Nest (Note we must take negative root of time equation in order to cross properly)
				t_cross = (-(x*v_x+y*v_y)-std::sqrt(pow(x*v_x+y*v_y,2)-pow(velo,2)*(2*R*r_enc-pow(r_enc,2))))/pow(velo,2); // From R to R-r_enc
				t_cross = t_cross + machine_tol; // make sure we enter the nest beyond numerical precision
				if(t_cross < machine_tol){
					throw std::runtime_error("INVALID CROSSOVER TIME");
				}else{
					x_positions.at(entry_ant) = x;
					y_positions.at(entry_ant) = y;
					x_velocities.at(entry_ant) = v_x;
					y_velocities.at(entry_ant) = v_y;
					nest_flag.at(entry_ant) = true;
					update_ant(entry_ant,t_cross); // move ant into nest
					entry_flag = check(r_enc); // check for overlap
					if(entry_flag == 1){
						//std::cout << "GOOD! ENTRY IN " << i+1 << " TRIES!" << std::endl;
						break;
					}
				}
			}
			
		}
		// If we were able to enter, update information
		if(entry_flag == 0){
			throw std::runtime_error("UNABLE TO ENTER NEST");
		}else{
			ant_name.at(entry_ant) = entry_ant;
			event_time.at(entry_ant) = entry_time;
			collisions.at(entry_ant) = 0.0;
			exit_time.at(entry_ant) = -1.;
		}
	
	return entry_ant; // return index of entering ant
	}
}


/* Calculate time before a wall-ant collision */
long double get_t_wall(Ants ants,long double R,long double r_enc,int &index1,long double t_entry,long double machine_tol){
	int num_ants = ants.x_positions.size();
	long double x,y,v_x,v_y;  // variables used to find collision time
	long double t_min;
	long double t_collide;
	long double x_new,y_new;
	long double t1, t2; // roots of t

	t_min = 10.*t_entry; // set minimum time before event to be longer than the next entry time
	
	//Find t_min for all particles
	for(int i=0; i<num_ants; ++i){
		if(ants.nest_flag[i]){
			x = ants.x_positions[i];
			y = ants.y_positions[i];
			v_x = ants.x_velocities[i];
			v_y = ants.y_velocities[i];

			t_collide = (-1*(x*v_x+y*v_y)+
					sqrt(pow(x*v_x+y*v_y,2)-(pow(v_x,2)+pow(v_y,2))*(pow(x,2)+pow(y,2)-pow(R-r_enc,2))))
					/(pow(v_x,2)+pow(v_y,2));

			if(t_collide > 2.*R/0.1){
				std::cout << "ERROR! TIME TO HIT WALL IS LARGERS THAN LARGEST CHORD!" << std::endl;
			}

			// Check if collision is the next event
			if(t_collide < t_min){
				t_min = t_collide;
				t_min = t_min - machine_tol; // make sure wall collision won't take us out of bounds
				index1 = i;
			}
		}		
	}
	return t_min;
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

	long double t_min = 10.*t_wall; //Initialize minimum time as something more than wall collision time
	// Note: Do not change this to t_min = t_wall! Numerical imprecision will mess things up

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
					d = pow(v_dot_r,2.)-v_dot_v*(r_dot_r-pow(sigma,2.));

					// If they collide, get time at which it occurs
					if(v_dot_r < 0. and d > 0.){
						t_collide = -1.*(v_dot_r+sqrt(d))/v_dot_v;
						// Check if collision is next event
						if(t_collide < t_min){
							t_min = t_collide;
							t_min = t_min - machine_tol; // make sure particles won't penetrate due to numerical precision
							index2 = i; // store colliding particles
							index3 = j;
						}
					}
				}
			}
		}
	}
	if(t_min < 0){

		throw std::runtime_error("WHAT?");
	}
	return t_min;
}


/* This update method is for a new ant entering */
// Note, if we are running this, no ant is exiting because there isn't a wall collision
// There also won't be any overlap because an ant-to-ant collision isn't being run
void Ants::update_all(long double t){
	// Update Particle Locations and velocities
	for(int i=0;i<x_positions.size();++i){
		if(nest_flag.at(i)){
			update_ant(i,t); // update locations
			event_time.at(i) = event_time.at(i)+t; // update times
		}
	}
}


void Ants::update_ant(int index,long double t){
	x_positions.at(index) = x_positions.at(index)+x_velocities.at(index)*t;
	y_positions.at(index) = y_positions.at(index)+y_velocities.at(index)*t;
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

	/* Update all particle locations */
	for(int i=0;i<x_positions.size();++i){
		if(nest_flag.at(i)){
			x_positions[i] = x_positions[i]+t_min*x_velocities[i];
			y_positions[i] = y_positions[i]+t_min*y_velocities[i];
			event_time[i] = event_time[i]+t_min;
		}
	}

	// Find the two particles that are colliding and update their velocities
	// Specular reflection
	for (int j=0;j<num_ants;++j){
		if (j == index1){
			for (int k=0;k<num_ants;++k){
				if (k == index2){
					long double del_x, del_y;
					long double theta; // coordinate rotation
					double vx1,vy1,vx2,vy2; // initial velocities
					double tx1,ty1,tx2,ty2; // temporary velocities
					double ux1,uy1,ux2,uy2; // new velocities

					// Initialize Variables
					vx1 = x_velocities[j];
					vy1 = y_velocities[j];
					vx2 = x_velocities[k];
					vy2 = y_velocities[k];
					del_x = x_positions[k]-x_positions[j];
					del_y = y_positions[k]-y_positions[j];
					theta = atan(del_y/del_x);

					// First rotate into new coordinate system
					tx1 = vx1*cos(theta)+vy1*sin(theta);
					ty1 = -vx1*sin(theta)+vy1*cos(theta);
					tx2 = vx2*cos(theta)+vy2*sin(theta);
					ty2 = -vx2*sin(theta)+vy2*cos(theta);

					// Flip x-velocity
					tx1 = -tx1;
					tx2 = -tx2;

					// De-rotate
					ux1 = tx1*cos(theta)-ty1*sin(theta);
					uy1 = tx1*sin(theta)+ty1*cos(theta);
					ux2 = tx2*cos(theta)-ty2*sin(theta);
					uy2 = tx2*sin(theta)+ty2*cos(theta);

					// Update
					x_velocities[j] = ux1;			
					y_velocities[j] = uy1;			
					x_velocities[k] = ux2;			
					y_velocities[k] = uy2;			

					// Tick collision counter
					collisions[j]++;
					collisions[k]++;
				}
			}
		}
	}
}


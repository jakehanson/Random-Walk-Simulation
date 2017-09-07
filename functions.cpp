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
/* Note: Its a constructor because there's no return type and method name 'Ants' matches class name. */
//Ants::Ants(size_t num_ants, long double temp)
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
/*Overload it with ways of initializing*/
//Method 1: explicitly provide seed to populate
void Ants::populate(long double R, long double velo,long double r_enc,int max_init, int seed, bool start_in_center)
{
	std::mt19937 gen(seed); // create generator with default seed
	std::normal_distribution<long double> nd(0.,1.0); // init normal dist with mean 0 and stddev 1
	std::uniform_real_distribution<long double> uni(0.,1); // init uniform dist on (0,1]

	// Define polar coords for position and velocity
	long double radius;
	long double pos_angle; 
	long double velo_angle;

	// If start_in_center initialize ant 0 in center of aperture with downward velo
	int index = 0;
	if(start_in_center == true)
	{
		// init position
		x_positions.at(0) = 0.0; // centered
		y_positions.at(0) = R-r_enc; // located at top

		// init velocity
		velo_angle = M_PI+uni(gen)*M_PI; // angle between pi and 2 pi (downward)
		x_velocities.at(0) = velo*cos(velo_angle);
		y_velocities.at(0) = velo*sin(velo_angle);

		// Init all other parameters
		ant_name.at(0) = 0;
		event_time.at(0) = 0.0;
		collisions.at(0) = 0.0;
		exit_time.at(0) = -1.;
		nest_flag.at(0) = true;

		index = index + 1;
	}

	// Initialize all positions and velocities
	for (int i=index; i<x_positions.size();++i)
	{
		radius = uni(gen)*(R-r_enc);
		pos_angle = uni(gen)*2*M_PI; // position angle
		velo_angle = uni(gen)*2*M_PI; // velocity angle

		ant_name.at(i) = i;
		x_positions.at(i) = radius*cos(pos_angle); // x=rcos(theta): r in (0,R-r_enc) and theta in (0,2pi)		
		y_positions.at(i) = radius*sin(pos_angle); // y=rsin(theta)
		x_velocities.at(i) = velo*cos(velo_angle); // initialize velocity in x direction
		y_velocities.at(i) = velo*sin(velo_angle); // initialize velicity in y direction
		event_time.at(i) = 0.0;
		collisions.at(i) = 0.0;
		exit_time.at(i) = -1.;
		nest_flag.at(i) = true;
	}
	
	// Make sure particles don't overlap!
	int counter = 0; // keep track of tries to initialize
	for (int j=1; j<x_positions.size();++j){
	
		int tally=0; // keeps track of all the other particles
		bool all_clear=false; // must ensure the particle clears all neighbors
		long double r;

		while(all_clear==false){
			// as long as we haven't cleared everyone keep going through and redrawing
			if(x_positions.size() == 1){
				break; // if we only have one particle this is pointless
			}
			// check all particle pairs
			for(int k=0;k<x_positions.size();++k){
				if(k!=j){
					r = sqrt(pow(x_positions[k]-x_positions[j],2)+pow(y_positions[k]-y_positions[j],2));
					if(r<2.*r_enc){
						// if overlap, redraw angle and radius
						radius = uni(gen)*(R-r_enc);
						pos_angle = uni(gen)*2*M_PI;
						x_positions.at(j) = radius*cos(pos_angle);		
						y_positions.at(j) = radius*sin(pos_angle);
						tally = 0; // reset talley
						counter++;
						if(counter>max_init){
							all_clear=true;
							throw std::runtime_error("FAILED TO INITIALIZE PARTICLES! NOT ENOUGH NEST SPACE.");
						}
						break; //break out of for loop
					}else{
						tally++; // we have successfully dodged another particle
						if(tally == x_positions.size()-1){
							all_clear=true; // we can only get here if we dodge all particles
						}
					}
				}
			}
		}
	}

	/* Make sure we obeyed the start in center flag */
	if(start_in_center == true){
		if(x_positions.at(0) != 0.0 or y_positions.at(0) != (R-r_enc)){
			throw std::runtime_error("START IN CENTER FLAG VIOLATED DURING INITIALIZATION.");
		}
	}
}
//Method 2: Let a random device (rd) provide seed for you
void Ants::populate(long double R, long double velo,long double r_enc,int max_init,bool start_in_center)
{
    std::random_device rd;
	populate(R,velo,r_enc,max_init,rd,start_in_center); // this will call method 3
}
//Method 3: Pass a reference to a previously created random device (rd) and use that to seed populate
void Ants::populate(long double R, long double velo,long double r_enc,int max_init, std::random_device &rd,bool start_in_center)
{
	populate(R,velo,r_enc,max_init,rd(),start_in_center); // this will call method 1
}

//Method 4: Single ant starts at the top of the nest with random downward velocity
void Ants::populate(long double R, long double a, long double velo, long double r_enc)
{
	// Initialize positions at top center
	if(x_positions.size() != 1){
		throw std::runtime_error("INITIALIZATION FAILED! SINGLE ANT FLAG USED WITH MULTIPLE ANTS.");
	}
	else{
		// Init random device and seed it
		std::random_device rd; // random device
		long double seed = rd(); // seed
		std::mt19937 gen(seed); // create generator with seed
		std::uniform_real_distribution<long double> uni(0.,1); // init uniform dist on (0,1]
		
		// Init Positions
		x_positions.at(0) = 0.0; // centered
		y_positions.at(0) = R-r_enc; // located at top

		// Init Velocity
		long double velo_angle = M_PI+uni(gen)*M_PI; // angle between pi and 2 pi (downward)
		x_velocities.at(0) = velo*cos(velo_angle);
		y_velocities.at(0) = velo*sin(velo_angle);

		// Init all other parameters
		ant_name.at(0) = 0;
		event_time.at(0) = 0.0;
		collisions.at(0) = 0.0;
		exit_time.at(0) = -1.;
		nest_flag.at(0) = true;
	}
}


/* Calculate time before a wall-ant collision */
long double get_t_wall(Ants ants,long double R,long double r_enc,int &index1){
	int num_ants = ants.x_positions.size();
	long double x,y,v_x,v_y;  // variables used to find collision time
	long double t_min;
	long double t_collide;
	bool init = false;

	//Initialize t_min with the first particle still in nest
	while(init==false){
		for(int i=0; i<num_ants; ++i){
			if(ants.nest_flag[i]){
				long double x = ants.x_positions[i];
				long double y = ants.y_positions[i];
				long double v_x = ants.x_velocities[i];
				long double v_y = ants.y_velocities[i];

				t_min = (-1*(x*v_x+y*v_y)+
					sqrt(pow(x*v_x+y*v_y,2)-(pow(v_x,2)+pow(v_y,2))*(pow(x,2)+pow(y,2)-pow(R-r_enc,2))))
					/(pow(v_x,2)+pow(v_y,2));
				index1 = i;
				init = true; // we have initialized t_min
				break;
			}
		}
	}

	//Find t_min for all particles
	for(int i=0; i<num_ants; ++i){
		if(ants.nest_flag[i]){
			long double x = ants.x_positions[i];
			long double y = ants.y_positions[i];
			long double v_x = ants.x_velocities[i];
			long double v_y = ants.y_velocities[i];
			
			//Get time to hit wall
			t_collide = (-1*(x*v_x+y*v_y)+
					sqrt(pow(x*v_x+y*v_y,2)-(pow(v_x,2)+pow(v_y,2))*(pow(x,2)+pow(y,2)-pow(R-r_enc,2))))
					/(pow(v_x,2)+pow(v_y,2));

			if(t_collide < t_min and pow(x,2.)+pow(y,2)<pow(R-r_enc,2.)){
				t_min = t_collide;
				index1 = i;
			}
		}		
	}
	//std::cout << "wall time:\t" << t_min << "\n";
	if (t_min < 0.){
		throw std::runtime_error("NEGATIVE TIME TO HIT WALL!");
	}
	return(t_min);
}


/* Calculate the time before the next ant-ant collision */
long double get_t_ant(Ants ants,long double r_enc,long double t_wall, int &index2, int &index3){
	int num_ants = ants.x_positions.size();
	long double sigma = 2*r_enc; // cross sectional distance
	long double v_dot_r;
	long double v_dot_v;
	long double r_dot_r;
	long double d;
	long double t_collide;

	long double t_min = 2.*t_wall; //Initialize minimum time as something more than wall collision time
	// Note: Do not change this to t_min = t_wall!

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

						// Check for errors
						if(d < 0){
							throw std::runtime_error("Distance Violation for Ant-to-Ant Collision!");
						}
						if (t_collide < 0.0){
							throw std::runtime_error("Negative Time for Ant-to-Ant Collision!");
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
	//std::cout << "Ant time:\t" << t_min << "\n";
	return(t_min);
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
			x_positions[i] = x_positions[i]+x_velocities[i]*t_min;
			y_positions[i] = y_positions[i]+y_velocities[i]*t_min;
			event_time[i] = event_time[i]+t_min;
		}
	}
}

/* This update method is for ant-to-ant collisions using conservation of momentum*/
void Ants::update(long double t_min,long double r_enc,int index1,int index2){
	int num_ants = x_positions.size();
	long double sigma = 2*r_enc;

	/* Update all particle locations */
	for(int i=0;i<num_ants;++i){
		x_positions[i] = x_positions[i]+t_min*x_velocities[i];
		y_positions[i] = y_positions[i]+t_min*y_velocities[i];
		event_time[i] = event_time[i]+t_min;
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
					x_velocities[index1] = x_velocities[j]+J_x;
					y_velocities[index1] = y_velocities[j]+J_y;
					x_velocities[k] = x_velocities[k]-J_x;
					y_velocities[k] = y_velocities[k]-J_y;

					// Tick collision counter
					collisions[j]++;
					collisions[k]++;
				}
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
		x_positions[i] = x_positions[i]+t_min*x_velocities[i];
		y_positions[i] = y_positions[i]+t_min*y_velocities[i];
		event_time[i] = event_time[i]+t_min;
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



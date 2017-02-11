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
//Ants::Ants(size_t num_ants, double temp)
Ants::Ants(int num_ants)
{
	ant_name = std::vector<int>(num_ants);
	x_positions = std::vector<double>(num_ants);
	y_positions = std::vector<double>(num_ants);
	x_velocities = std::vector<double>(num_ants);
	y_velocities = std::vector<double>(num_ants);
	event_time = std::vector<double>(num_ants);
	exit_time = std::vector<double>(num_ants);
	collisions = std::vector<double>(num_ants);
	nest_flag = std::vector<bool>(num_ants);
}


/* Define Method to Populate our arrays within Ant Class */
/*Overload it with ways of initializing*/
//Method 1: explicitly provide seed to populate
void Ants::populate(double R, double T,double r_enc,int max_init, int seed)
{
	std::mt19937 gen(seed); // create generator with default seed
	std::normal_distribution<double> nd(0.,1.0); // init normal dist with mean 0 and stddev 1
	std::uniform_real_distribution<double> uni(0.,1); // init uniform dist on (0,1]

	// Initialize all positions and velocities
	int counter = 0;
	for (int i=0; i<x_positions.size();++i)
	{
		double radius = uni(gen)*(R-r_enc);
		double angle = uni(gen)*2*M_PI;
		ant_name.at(i) = i                                         ;
		x_positions.at(i) = radius*cos(angle); // x=rcos(theta): r in (0,R-r_enc) and theta in (0,2pi)
		y_positions.at(i) = radius*sin(angle); // y=rsin(theta)
		x_velocities.at(i) = sqrt(T)*nd(gen);
		y_velocities.at(i) = sqrt(T)*nd(gen);
		event_time.at(i) = 0.0;
		collisions.at(i) = 0.0;
		exit_time.at(i) = -1.;
		nest_flag.at(i) = true;
	}

	for (int j=0; j<x_positions.size();++j){
	
		int tally=0; // keeps track of  all the other particles
		bool all_clear=false; // must ensure the particle clears all neighbors
		double r;

		while(all_clear==false){
			// as long as we haven't cleared everyone keep going through and redrawing
			if(x_positions.size() == 1){
				break; // if we only have one particle this is pointless
			}
			// check all particle pairs
			for(int k=0;k<x_positions.size();++k){
				if(k!=j){
					r = sqrt(pow(x_positions[k]-x_positions[j],2)+pow(y_positions[k]-y_positions[j],2));
					if(r<2*r_enc){
						x_positions.at(j) = uni(gen)*(R-r_enc)*cos(uni(gen)*2.*M_PI); // redraw
						y_positions.at(j) = uni(gen)*(R-r_enc)*sin(uni(gen)*2.*M_PI);
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
}
//Method 2: Let a random device (rd) provide seed for you
void Ants::populate(double R, double T,double r_enc,int max_init)
{
    std::random_device rd;
	populate(R,T,r_enc,max_init,rd); // this will call method 3
}
//Method 3: Pass a reference to a previously created random device (rd) and use that to seed populate
void Ants::populate(double R, double T,double r_enc,int max_init, std::random_device &rd)
{
	populate(R,T,r_enc,max_init,rd()); // this will call method 1
}

//Method 4: Single ant starts at the top of the nest with random downward velocity
void Ants::populate(double R, double a, double T, double r_enc)
{
	// Initialize positions at top center
	if(x_positions.size() != 1){
		throw std::runtime_error("INITIALIZATION FAILED! SINGLE ANT FLAG USED WITH MULTIPLE ANTS.");
	}
	else{
		bool x_velo = false; // flags to insist velocity is downward
		bool y_velo = false; // flags to insist velocity is downward
		
		// Init random device and seed it
		std::random_device rd; // random device
		double seed = rd(); // seed
		std::mt19937 gen(seed); // create generator with seed
		std::normal_distribution<double> nd(0.,1.0); // init normal dist with mean 0 and stddev 1
		
		// Init Positions
		x_positions.at(0) = 0.0; // centered
		y_positions.at(0) = R-r_enc; // located at top
		
		// Init velocities
		while(x_velo == false){
			x_velocities.at(0) = sqrt(T)*nd(gen);
			if(x_velocities.at(0) < 0){
				x_velo = true;
			}
		}
		while(y_velo == false){
			y_velocities.at(0) = sqrt(T)*nd(gen);
			if(y_velocities.at(0) < 0){
				y_velo = true;
			}
		}

		event_time.at(0) = 0.0;
		collisions.at(0) = 0.0;
		exit_time.at(0) = -1.;
		nest_flag.at(0) = true;
	}
}


/* Calculate time before a wall-ant collision */
double get_t_wall(Ants ants,double R,double r_enc,int &index1){
	int num_ants = ants.x_positions.size();
	double x,y,v_x,v_y;  // variables used to find collision time
	double t_min;
	double t_collide;
	bool init = false;

	//Initialize t_min with the first particle in nest
	while(init==false){
		for(int i=0; i<num_ants; ++i){
			if(ants.nest_flag[i]){
				double x = ants.x_positions[i];
				double y = ants.y_positions[i];
				double v_x = ants.x_velocities[i];
				double v_y = ants.y_velocities[i];

				t_min = (-1*(x*v_x+y*v_y)+
					sqrt(pow(x*v_x+y*v_y,2)-(pow(v_x,2)+pow(v_y,2))*(pow(x,2)+pow(y,2)-pow(R-r_enc,2))))
					/(pow(v_x,2)+pow(v_y,2));

				// Check if t_min is nan
				if (t_min != t_min){
					throw std::runtime_error("FAILED TO RESOLVE WALL COLLISION! TOO MANY ANTS");
				}

				index1 = i;
				init = true; // we have initialized t_min
				break;
			}
		}
	}

	//Find t_min for all particles
	for(int i=0; i<num_ants; ++i){
		if(ants.nest_flag[i]){
			double x = ants.x_positions[i];
			double y = ants.y_positions[i];
			double v_x = ants.x_velocities[i];
			double v_y = ants.y_velocities[i];
			
			//Get time to hit wall
			t_collide = (-1*(x*v_x+y*v_y)+
					sqrt(pow(x*v_x+y*v_y,2)-(pow(v_x,2)+pow(v_y,2))*(pow(x,2)+pow(y,2)-pow(R-r_enc,2))))
					/(pow(v_x,2)+pow(v_y,2));
			
			// Check if t_collide is nan
			if (t_collide != t_collide){
				throw std::runtime_error("FAILED TO RESOLVE WALL COLLISION! TOO MANY ANTS");
			}
	
			if(t_collide < t_min){
				t_min = t_collide;
				index1 = i;
			}
		}		
	}
	return(t_min);
}


/* Calculate the time before the next ant-ant collision */
double get_t_ant(Ants ants,double r_enc,double t_wall, int &index2, int &index3){
	int num_ants = ants.x_positions.size();
	double sigma = 2*r_enc;
	double v_dot_r;
	double v_dot_v;
	double r_dot_r;
	double d;
	double t_collide;

	double t_min = 2*t_wall; //Initialize minimum time as anytime more than wall collision time

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
					if(v_dot_r < 0. and d >= 0.){
						t_collide = -1.*(v_dot_r+sqrt(d))/v_dot_v;
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
	return(t_min);
}


/* This update method is for wall collisions*/
void Ants::update(double R,double r_enc, double a,double t_min,int index1){
	int num_ants = x_positions.size();
	double angle1; // angle to the location where the particle hits the wall
	double angle2; // angle to the initial position of the particle
	double final_angle; // angle of where the particle would be t_min after colliding with wall


	/* Update Particle Locations amd velocities*/
	for(int i=0;i<num_ants;++i){

		// The particle that is hitting the wall switches its velocity component
		if(i == index1){
			double x = x_positions[i];
			double y = y_positions[i];
			double v_x = x_velocities[i];
			double v_y = y_velocities[i];
			double x_new = x_positions[i]+x_velocities[i]*t_min;
			double y_new = y_positions[i]+y_velocities[i]*t_min;
			double x_temp,y_temp;

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

/* This update method is for particle collisions*/
void Ants::update(double t_min,double r_enc,int index1,int index2){
	int num_ants = x_positions.size();
	double sigma = 2*r_enc;

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
					double del_x, del_y, del_vx, del_vy;
					double J, J_x, J_y;
					double v_dot_r;

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

					// Update velocity
					x_velocities[j] = x_velocities[j]+J_x;
					y_velocities[j] = y_velocities[j]+J_y;
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



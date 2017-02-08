#include "header.h"
#include <stdio.h>
#include <cmath>
#include <stdexcept>

//GENERATE STATISTICS (DON'T REALLY EVEN NEED COLLISION COUNTER FOR THE ABOVE)
//GET TIME SPENT AS A FUNCTION OF ENCOUNTER RATE

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
/*Overload the shit out of it with ways of seeding random number*/
//Method 1: explicitly provide seed to populate
void Ants::populate(double L, double H, double T,double r_enc,int max_init, int seed)
{
	std::mt19937 gen(seed); // create generator with default seed
	std::normal_distribution<double> nd(0.,1.0); // init normal dist with mean 0 and stddev 1
	std::uniform_real_distribution<double> uni(0.,1);

	// Initialize all positions and velocites
	int counter = 0;
	for (int i=0; i<x_positions.size();++i)
	{
		ant_name.at(i) = i;
		x_positions.at(i) = r_enc+uni(gen)*(L-2*r_enc);
		y_positions.at(i) = r_enc+uni(gen)*(H-2*r_enc);
		x_velocities.at(i) = sqrt(T)*nd(gen);
		y_velocities.at(i) = sqrt(T)*nd(gen);
		event_time.at(i) = 0.0;
		collisions.at(i) = 0.0;
		exit_time.at(i) = -1.;
		nest_flag.at(i) = true;
	}

	for (int j=0; j<x_positions.size();++j){
	
		int tally=0;
		bool all_clear=false; // must ensure the particle clears all neighbors
		double r;

		while(all_clear==false){
			// as long as we haven't cleared everyone keep going through and redrawing
			if(x_positions.size() == 1){
				break; // if we only have one particle this is pointless
			}
			for(int k=0;k<x_positions.size();++k){
				if(k!=j){
					r = sqrt(pow(x_positions[k]-x_positions[j],2)+pow(y_positions[k]-y_positions[j],2));
					if(r<2*r_enc){
						x_positions.at(j) = r_enc+uni(gen)*(L-2*r_enc);
						y_positions.at(j) = r_enc+uni(gen)*(H-2*r_enc);
						tally = 0; // reset talley
						counter++;
						if(counter>max_init){
							//std::cout << "FATAL ERROR!! FAILED TO INITIALIZE PARTICLES." << std::endl;
							all_clear=true;
							throw std::runtime_error("FAILED TO INITIALIZE PARTICLES.");
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
void Ants::populate(double L, double H, double T,double r_enc,int max_init)
{
    std::random_device rd;
	populate(L,H,T,r_enc,max_init,rd); // this will call method 3
}
//Method 3: Pass a reference to a previously created random device (rd) and use that to seed populate
void Ants::populate(double L, double H, double T,double r_enc,int max_init, std::random_device &rd)
{
	populate(L,H,T,r_enc,max_init,rd()); // this will call method 1
}


/* Calculate time before a wall-ant collision */
double get_t_wall(Ants ants,double L,double H,double r_enc,int &index1){
	int num_ants = ants.x_positions.size();
	double t_min;
	double t_collide;
	bool init = false;

	//Initialize t_min with the first particle in nest
	while(init==false){
		for(int i=0; i<num_ants; ++i){
			if(ants.nest_flag[i]){
				if(ants.x_velocities[i] > 0){
					t_min = (L-r_enc-ants.x_positions[i])/ants.x_velocities[i];
				}else{
					t_min = (ants.x_positions[i]-r_enc)/std::abs(ants.x_velocities[i]);
				}
				index1 = i;
				init = true; // we have initialized t_min
				break;
			}
		}
	}
	if(ants.x_velocities[0] > 0){
		t_min = (L-r_enc-ants.x_positions[0])/ants.x_velocities[0];
	}
	else{
		t_min = (ants.x_positions[0]-r_enc)/std::abs(ants.x_velocities[0]);
	}
	index1 = 0;


	//Find t_min
	for(int i=0; i<num_ants; ++i){

		if(ants.nest_flag[i]){
			//Get time to hit vertical wall
			if(ants.x_velocities[i] > 0){
				t_collide = (L-r_enc-ants.x_positions[i])/ants.x_velocities[i];
			}
			else{
				t_collide = (ants.x_positions[i]-r_enc)/std::abs(ants.x_velocities[i]);
			}
			if(t_collide < t_min){
				t_min = t_collide;
				index1 = i;
			}
			//Get time to hit horizontal wall
			if(ants.y_velocities[i] > 0){
				t_collide = (H-r_enc-ants.y_positions[i])/ants.y_velocities[i];
			}
			else{
				t_collide = (ants.y_positions[i]-r_enc)/std::abs(ants.y_velocities[i]);
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
void Ants::update(double L, double H,double a,double t_min,int index1){
	int num_ants = x_positions.size();

	/* Update Particle Locations */
	for(int i=0;i<num_ants;++i){

		x_positions[i] = x_positions[i]+t_min*x_velocities[i];
		y_positions[i] = y_positions[i]+t_min*y_velocities[i];
		event_time[i] = event_time[i]+t_min;

		// The particle that is hitting the wall switches its velocity component
		if(i == index1){

			// Figure out which wall we hit
			double left_dist = x_positions[i];
			double right_dist = L - x_positions[i];
			double bottom_dist = y_positions[i];
			double top_dist = H - y_positions[i];
			double vert_dist,horiz_dist;

			if(left_dist < right_dist){
				vert_dist = left_dist;
			}else{
				vert_dist = right_dist;}
			if(bottom_dist<top_dist){
				horiz_dist = bottom_dist;
			}else{
				horiz_dist = top_dist;
				// Did we make it out?
				if(L/2.-a/2.<x_positions[i] && x_positions[i] < L/2.+a/2.){
					nest_flag[i] = false;
					exit_time[i] = event_time[i];
				}
			}

			// Switch velocity accordingly
			if(vert_dist < horiz_dist){
				x_velocities[i] = -1*x_velocities[i];
			}else{
				y_velocities[i] = -1*y_velocities[i];
			}
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



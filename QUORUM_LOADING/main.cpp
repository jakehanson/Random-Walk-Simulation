#include "header.h"
#include <fstream>

// CAREFULLY READ ALL CODE... THROW ERRORS WHERE FAIL
// SAVE TO GIT
// GET DATA

/* Function to define geometry, initialize ants, and run/write simulation */
int main(int argc, char** argv)
{
	//First get variable params from either first arg of command line or hardcode
	int max_ants;
	double entry_rate;
	double sim_time; // we run the sim for a fixed length of time rather than number of steps
	if(argc == 1){
		entry_rate = 1/10.; // 1 ant enters every 30 seconds
		max_ants = 20;
		sim_time = 10000; // max time to run sim
	}else if(argc == 4){
		entry_rate = atof(argv[1]);
		max_ants = atoi(argv[2]);
		sim_time = atof(argv[3]);
	}else{
		std::cout << "ERROR! INCORRECT NUMBER OF INPUT VARIABLES!" << std::endl;
	}

	/* Define Fixed Params */
	long double R = 2;  // Box Radius (in cm)
	long double a = 2.0;  // aperture size (in cm)
	long double velo = .1; // initial ant velocity (in cm/s)
	long double r_enc = 0.1;  // encounter radius (cm)
	long double machine_tol = 1e-12; // corrects for numerical inaccuracies
	long double t = 0; // keeps track of simulation time
	long double t_wall,t_ant,t_entry; // variables to store the time until next event
	int index1,index2,index3; // indices for particles involved in collision
	double ants_in_nest; // keeps track of ants in nest
	int overlap_flag; // catches overlap issues

	// Write parameters to file
	std::ofstream params("params.txt");
	params << "R\ta\tvelocity\tr_enc\tentry_rate\tmax_ants" << std::endl;
	params << R << "\t" << a << "\t" << velo << "\t" << r_enc << "\t" << entry_rate << "\t" << max_ants << std::endl;
	params.close();

	// Open file to hold simulation data
	std::ofstream data_file("output.txt");
	data_file << "Name\t" << "x\t" << "y\t" << "v_x\t" << "v_y\t" << "event_time\t"
				<< "exit_time\t" << "collisions\t" << "in_nest" << std::endl;

	/* Check Params and Initialize */
	Ants ants = Ants(max_ants); // instance of Ants called ants
	try{
		ants.populate(); // initialize all values
		if(a/R >=(2/3.*M_PI)){
			throw std::runtime_error("APERTURE SIZE TOO LARGE."); 
		}else if(R < 0.1 or R > 20){
			throw std::runtime_error("INVALID NEST SIZE");
		}else if(velo < 0.01 or velo > 10){
			throw std::runtime_error("INVALID VELOCITY.");
		}else if(r_enc > a/2.){
			throw std::runtime_error("ENCOUNTER RADIUS TOO LARGE.");
		}else{
			//std::cout << "INPUT PARAMS VALID\n";
			//if all params are valid machine precision of 1e-12 will be safe
		}
	}
	catch(std::runtime_error &e){
		std::cerr << "RUNTIME ERROR: " << e.what() << std::endl;
		return 1;
	}
	data_file << ants; //prints all of ants to name.txt using overloaded operator <<


	/* Run Sim */
	double entry_counter; // keeps track of which ant is entering the nest
	while(t <= sim_time){

		ants_in_nest = std::accumulate(ants.nest_flag.begin(),ants.nest_flag.end(),0); // keep track of number of ants in nest
		std::cout << t << "\t" << ants_in_nest << std::endl;

		// If there are no ants in the nest, allow the next ant to enter
		if(std::accumulate(ants.nest_flag.begin(),ants.nest_flag.end(),0) == 0){
			try{
				entry_counter = ants.enter(R,a,velo,r_enc,entry_rate,machine_tol); // run enter method
			}catch(std::runtime_error &e){
				std::cerr << "RUNTIME ERROR: " << e.what() << std::endl;
				return 1;		
			}
			// Check if we are out of ants
			if(entry_counter == -1){
				std::cout << "MAXIMUM ANTS REACHED!" << std::endl;
				t = sim_time + 1; // breaks us out of loop
				break;
			// Else, update time
			}else{
				t = ants.event_time.at(entry_counter);
			}
			// Write to file
			data_file << ants;
		//If there are ants in the nest, see what event happens first
		}else{
			try{
				// Get time until next entry
				if(entry_counter+1 < max_ants){
					t_entry = (entry_counter+2)/entry_rate - t;
				}else{
					t_entry = 2*sim_time; // no ants left, so set entry time > than simulation time
				}
				t_wall = get_t_wall(ants,R,r_enc,index1,t_entry,machine_tol); // how long until an ant collides with the wall?
				t_ant = get_t_ant(ants,r_enc,t_wall,index2,index3,machine_tol); // how long until two ants collide?
				// Run correct update
				if(t_entry < t_wall and t_entry < t_ant){
					//std::cout << "RUN ENTRY" << std::endl;
					ants.update_all(t_entry); // move existing ants to new location
					entry_counter = ants.enter(R,a,velo,r_enc,entry_rate,machine_tol); // try to add a new ant
					t = t+t_entry; // update time
				}
				if(t_wall < t_entry and t_wall < t_ant){
					//std::cout << "RUN WALL" << std::endl;
					ants.update(R,r_enc,a,t_wall,index1);	// run ant-to-wall collisions for single particle
					overlap_flag = ants.check(r_enc); // check there isn't penetration
					if(overlap_flag == 0){
						data_file << ants; // write the overlapping data to file
						throw std::runtime_error("OVERLAP ERROR");
					}
					t = ants.event_time.at(index1); // update time

				}
				if(t_ant < t_wall and t_ant < t_entry){
					//std::cout << "RUN ANT" << std::endl;
					//std::cout << "\t t_ant = " << t_ant << std::endl;
					ants.update(t_ant,r_enc,index2,index3,velo); // run an ant-to-ant collision with fixed velocity
					overlap_flag = ants.check(r_enc); // check there isn't penetration
					if(overlap_flag == 0){
						data_file << ants; // write the overlapping data to file
						throw std::runtime_error("OVERLAP ERROR");
					}
					t = ants.event_time.at(index2);
				}
				if(t <= sim_time){
					data_file << ants; // write to file
				}else{
					std::cout << "MAXIMUM TIME REACHED!" << std::endl;
					break; // just in case...
				}
			}catch(std::runtime_error &e){
				std::cerr << "RUNTIME ERROR: " << e.what() << std::endl;
				return 1;		
			}
		}
	}

	// If we didn't catch an error, the treatment was good	
	std::cout << "Good\t" << std::accumulate(ants.nest_flag.begin(),ants.nest_flag.end(),0) << std::endl;
	data_file.close();
	//std::cout << "SIMULATION COMPLETE." << std::endl;

}

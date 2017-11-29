#include "header.h"
#include <fstream>

// NEED TO UPDATE HEADER
// NEED TO SOLVE ENTRY OVERLAP PROBLEM
// NEED TO OUTPUT ENCOUNTER RATE AS A FUNCTION OF TIME


/* Function to define geometry, initialize ants, and run/write simulation */
int main(int argc, char** argv)
{

	//First get variable params from either first arg of command line or hardcode
	int max_ants;
	double entry_rate;
	double sim_time; // we run the sim for a fixed length of time rather than number of steps
	if(argc == 1){
		entry_rate = 0.1; // 1 ant enters every 10 seconds
		max_ants = 100;
		sim_time = 100;
	}else if(argc == 4){
		entry_rate = atof(argv[1]);
		max_ants = atoi(argv[2]);
		sim_time = atof(argv[3]);
	}else{
		std::cout << "ERROR! INCORRECT NUMBER OF INPUT VARIABLES!" << std::endl;
	}

	/* Define Variables */
	long double R = 2;  // Box Radius (in cm)
	long double a = 1;  // aperture size (in cm)
	long double velo = .1; // initial ant velocity (in cm/s)
	long double r_enc = 0.1;  // encounter radius (cm)
	double t = 0; // keeps track of simulation time

	long double t_wall,t_ant,t_entry; // variables to store the time until next event
	int index1,index2,index3; // indices for particles involved in collision

	// Write parameters to file
	std::ofstream params("params.txt");
	params << "R\ta\tvelocity\tr_enc\tentry_rate\tmax_ants" << std::endl;
	params << R << "\t" << a << "\t" << velo << "\t" << r_enc << "\t" << entry_rate << "\t" << max_ants << std::endl;
	params.close();

	// Open file to hold simulation data
	std::ofstream data_file("output.txt");
	data_file << "Name\t" << "x\t" << "y\t" << "v_x\t" << "v_y\t" << "event_time\t"
				<< "exit_time\t" << "collisions\t" << "in_nest" << std::endl;

	/* Initialize ants */	 	
	Ants ants = Ants(max_ants); // instance of Ants called ants
	try{
		ants.populate(R,a,velo,r_enc); // initialize all values
	}
	catch(std::runtime_error &e){
		std::cerr << "RUNTIME ERROR: " << e.what() << std::endl;
		return 1;
	}
	
	data_file << ants; //prints all of ants to name.txt using overloaded operator <<

	/* Run Sim */
	double entry_counter = -1; // keeps track of which ant is to enter	
	while(t <= sim_time){

		std::cout << "Time = " << t << std::endl;

		// If there are no ants in the nest, allow the next ant to enter
		if(std::accumulate(ants.nest_flag.begin(),ants.nest_flag.end(),0) == 0){
			t = (entry_counter+1)/entry_rate; // so we are at the next entry time
			entry_counter = ants.entry(R,a,velo,r_enc,t); // have an ant enter
			if(entry_counter = -1){
				std::cout << "MAXIMUM ANTS REACHED!" << std::endl;
				t = sim_time + 1;
				break;
			}
			data_file << ants; // write to file

		// If there are ants in the nest, see what event happens first
		}else{
			t_entry = (entry_counter+1)/entry_rate - t; // how long until next ant enters?
			t_wall = get_t_wall(ants,R,r_enc,index1); // how long until colliding with wall?
			t_ant = get_t_ant(ants,r_enc,t_wall,index2,index3); // how long before colliding with another ant?
			if(t_entry < t_wall and t_entry < t_ant){
				t = t+t_entry; // update time
				ants.update(t_entry); // moves ants to new spots in nest
				entry_counter = ants.entry(R,a,velo,r_enc,t); // have an ant enter
				if(entry_counter = -1){
					std::cout << "MAXIMUM ANTS REACHED!" << std::endl;
					t = sim_time + 1;
					break;
				}
				data_file << ants;
			}
			else if(t_wall < t_entry and t_wall < t_ant){
				ants.update(R,r_enc,a,t_wall,index1);	// run ant-to-wall collisions for single particle
				t = ants.event_time.at(index1); // update time
			}
			else if(t_ant < t_wall and t_ant < t_entry){
				ants.update(t_ant,r_enc,index2,index3,velo); // run an ant-to-ant collision with fixed velocity
				t = ants.event_time.at(index2);
			}else{
				std::cout << "ERROR! NO CASE SATISFIED!" << std::endl;
			}
		}

		data_file << ants; //write to file
		ants_in_nest = std::accumulate(ants.nest_flag.begin(),ants.nest_flag.end(),0); // keep track of number of ants in nest
	}

	//std::cout << "Good\t" << num_ants << "\t" << ants.event_time[0] << "\t" << ants.collisions[0] << std::endl;


	data_file.close();
	std::cout << "SIMULATION COMPLETE." << std::endl;

}

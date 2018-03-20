#include "header.h"
#include <fstream>

// CHANGES:
// Get rid of anything to do with brood
// During Initialization of recruiter entry allow position to vary
// Switch up order of what is taken as an arg
// Not currently writing to data file or params, just enc_rate
// throw runtime error if max time reached
// Get rid of longest chord error (functions line 244 or so)


/* Function to define geometry, initialize ants, and run/write simulation */
int main(int argc, char** argv)
{
	//First get variable params from either first arg of command line or hardcode
	int N_recruiters = 50;
	long double recruiter_rate = 1/300.;
	long double sim_time = 1000000; // we run the sim for a fixed length of time rather than number of steps
	/* Define Fixed Params */
	long double R;  // Box Radius (in cm)
	long double a;  // aperture size (in cm)	
	if(argc == 1){
		R = 5.0;
		a = 0.5;
	}else if(argc == 3){
		R = atof(argv[1]);
		a = atof(argv[2]); // rate after sensing quorum
		R = 1.0 + R/5.; // have to rescale since BASH doesn't allow floats this puts us in (1,5) with 4 radii samples per integer
		a = a/0.5; // rescale a to be 2.0
	}else{
		std::cout << "ERROR! INCORRECT NUMBER OF INPUT VARIABLES!" << std::endl;
	}

	//std::cout << "R = " << R << "\t" << " a = " << a << std::endl;

	/* Define Fixed Params */
	long double velo = .1; // initial ant velocity (in cm/s)
	long double r_enc = 0.1;  // encounter radius (cm)
	long double machine_tol = 1e-12; // corrects for numerical inaccuracies
	long double t = 0; // keeps track of simulation time
	long double t_wall,t_ant,t_recruiter; // variables to store the time until next event
	int index1,index2,index3; // indices for particles involved in collision
	double ants_in_nest,N_transporters; // keeps track of ants in nest and the number of tranporters
	int overlap_flag; // catches overlap issues

	// Write parameters to file
	// std::ofstream params("params.txt");
	// params << "R\ta\tvelocity\tr_enc\trecruiter_rate\tN_recruiters" << std::endl;
	// params << R << "\t" << a << "\t" << velo << "\t" << r_enc << "\t" << recruiter_rate << "\t" << N_recruiters << std::endl;
	// params.close();

	// Open file to hold simulation data
	// std::ofstream data_file("output.txt");
	// data_file << "Name\t" << "x\t" << "y\t" << "v_x\t" << "v_y\t" << "entry_time\t" << "event_time\t"
	// 			<< "exit_time\t" << "quorum_flag\t" << "collisions\t" << "in_nest" << std::endl;

	/* Check Params and Initialize */
	Ants ants = Ants(N_recruiters); // instance of Ants called ants. Stores recruiter and info
	try{
		ants.populate(); // initialize all values
		if(a/R >=(2/3.*M_PI)){
			throw std::runtime_error("APERTURE SIZE TOO LARGE."); 
		}else if(R < 0.1 or R > 11){
			throw std::runtime_error("INVALID NEST SIZE");
		}else if(velo < 0.01 or velo > 5){
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
	//data_file << ants; //prints all of ants to name.txt using overloaded operator <<


	/* Run Sim */
	long double entry_counter; // keeps track of which ant is entering the nest
	long double recruiter_time; // time of next event
	int entry_ant = -1; // initialize which ant should enter next

	while(t <= sim_time){
		ants_in_nest = std::accumulate(ants.nest_flag.begin(),ants.nest_flag.end(),0); // keep track of number of ants in nest
		N_transporters = std::accumulate(ants.quorum_flag.begin(),ants.quorum_flag.end(),0); // get number of transporters
		//std::cout << t << "\t" << ants_in_nest << "\t" << N_transporters << std::endl;

		try{
			// Get next recruiter entering nest
			entry_ant = -1; // reset
			for(int i=0;i<N_recruiters;i++){
				if(ants.nest_flag.at(i) == false and ants.exit_time.at(i) == -1){ // recruiters must not be in nest nor have previously left nest
					entry_ant = i;
					break;
				}
			}
			// If there are no recruiters left check if there are still ants inside
			if(entry_ant == -1){
				if(ants_in_nest == 0){
					break; // done
				}else{
					t_recruiter = sim_time+1; // not done yet, but no more recruitment
				}
			}else{
				t_recruiter = (entry_ant+1)/recruiter_rate-t; // time until next recruitment event
			}

			// Get t_ant ant t_wall
			t_wall = get_t_wall(ants,R,r_enc,index1,t_recruiter,machine_tol); // how long until an ant collides with the wall?
			t_ant = get_t_ant(ants,r_enc,t_wall,index2,index3,machine_tol); // how long until two ants collide?
			
			// Run correct update
			if(t_recruiter < t_wall and t_recruiter < t_ant){
				if(t_recruiter < 0){
					throw std::runtime_error("ERROR NEGATIVE RECRUITER TIME!");
				}else{
					//std::cout << "RUN RECRUITER" << std::endl;
					ants.move_all(t_recruiter); // move existing ants to new location
					ants.recruiter_enter(entry_ant,R,a,velo,r_enc,recruiter_rate,machine_tol); // try to add a new recruiter
					t = t+t_recruiter; // update time
				}
			}
			if(t_wall < t_recruiter and t_wall < t_ant){
				//std::cout << "RUN WALL" << std::endl;
				ants.update(R,r_enc,a,t_wall,index1);	// run ant-to-wall collisions for single particle
				t = ants.event_time.at(index1); // update time
			}
			if(t_ant < t_wall and t_ant < t_recruiter){
				//std::cout << "RUN ANT" << std::endl;
				ants.update(t_ant,r_enc,index2,index3,velo); // run an ant-to-ant collision with fixed velocity
				t = ants.event_time.at(index2);
			}
			if(t <= sim_time){
				overlap_flag = ants.check(r_enc); // check there isn't penetration
				if(overlap_flag == 0){
					//data_file << ants; // write the overlapping data to file
					throw std::runtime_error("OVERLAP ERROR");
				}
				else{
					//data_file << ants; // write to file
				}
			}else{
				throw std::runtime_error("MAXIMUM TIME REACHED!");
				break; // just in case...
			}
		}catch(std::runtime_error &e){
			std::cerr << "RUNTIME ERROR: " << e.what() << std::endl;
			return 1;
		}
	}

	// If we didn't catch an error, the treatment was good	
	// std::cout << "Good\t" << std::accumulate(ants.nest_flag.begin(),ants.nest_flag.end(),0) << std::endl;
	// data_file.close();
	// std::cout << "SIMULATION COMPLETE." << std::endl;

	// Write Encounter rate info to file
	//std::ofstream enc_rate("enc_rate.txt");
	for(int i=0;i<ants.x_positions.size();i++){
		std::cout << a << "\t" << R << "\t" << 1/recruiter_rate << "\t" << N_recruiters << "\t" << ants.ant_name.at(i) << "\t" << ants.exit_time.at(i) << "\t" << ants.entry_time.at(i) << "\t" << ants.collisions.at(i) << "\t" << std::endl;
	}
	//enc_rate.close();
}

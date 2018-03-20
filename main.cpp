#include "header.h"
#include <fstream>


/* Function to define geometry, initialize ants, and run/write simulation */
int main(int argc, char** argv)
{

	//First get number of ants in sim either from first arg of command line or hardcode
	int num_ants;
	if(argc == 1){
		num_ants = 20; // number of ants in simulation
	}else{
		num_ants = atoi(argv[1]);
	}
	//std::cout << "Number of Ants:\t" << num_ants << std::endl;


	/* Define Variables */
	long double R = 2;  // Box Radius (in cm)
	long double a = 1;  // aperture size (in cm)
	long double velo = .1; // initial ant velocity (in cm/s)
	long double r_enc = 0.1;  // encounter radius (cm)
	bool single_particle = false; // true means we start a single particle starting IN APERTURE CENTER
	bool fixed_velo = true; // if true, ants will never change velocity, only heading
	bool start_in_center = true; // true means one ant 0 starts in center of aperture facing down
	bool exit_flag = true; // true means the trial ends when ant 0 leaves 
	int max_init = 500;  // max number of tries to initialize a given setup so no ants overlap 
	int max_steps= 50000; // max number of events in simulation (event is collision w/ wall or ant)

	int ants_in_nest = num_ants; // used to calculate how many ants are in nest
	bool collision_flag = true; // true means collisions are on
	bool all_clear = false;  // used to break out of main simulation
	int counter = 0; // counts simulation steps

	long double t_wall,t_ant; // variables to store the time until a collision	
	int index1,index2,index3; // indices for particles involved in collision
	//int seed = 8;  // pass the same seed to 'populate' method to initialize trial the same way

	// Write parameters to file
	std::ofstream params("params.txt");
	params << "R\ta\tvelocity\tr_enc\tnum_ants\tcollision_flag\tcenter_flag\texit_flag\tfixed_velo" << std::endl;
	params << R << "\t" << a << "\t" << velo << "\t" << r_enc << "\t" << num_ants << "\t"
			<< collision_flag << "\t" << start_in_center << "\t" << exit_flag << "\t" << fixed_velo << std::endl;
	params.close();

	// Open file to hold simulation data
	std::ofstream data_file("output.txt");
	data_file << "Name\t" << "x\t" << "y\t" << "v_x\t" << "v_y\t" << "event_time\t"
				<< "exit_time\t" << "collisions\t" << "in_nest" << std::endl;

	/* Initialize ants */	 	
	Ants ants = Ants(num_ants); // instance of Ants called ants
	try{
		if(single_particle){
			ants.populate(R,a,velo,r_enc); // initialize a single particle
		}
		else{
			ants.populate(R,velo,r_enc,max_init,start_in_center); // initialize positions and velocities
		}
	}
	catch(std::runtime_error &e){
		std::cerr << "RUNTIME ERROR: " << e.what() << std::endl;
		return 1;
	}
	
	data_file << ants; //prints all of ants to name.txt using overloaded operator <<

	/* Run Sim */	
	while(all_clear == false){
		
		try{
			t_wall = get_t_wall(ants,R,r_enc,index1); // how long double until colliding with wall?
		}
		catch(std::runtime_error &e){
			std::cerr << "RUNTIME ERROR: " << e.what() << std::endl;
			return 1;
		}

		if(collision_flag == true && single_particle == false){
			try{
				t_ant = get_t_ant(ants,r_enc,t_wall,index2,index3); // how long double before colliding with another ant?
			}
			catch(std::runtime_error &e){
				std::cerr << "RUNTIME ERROR: " << e.what() << std::endl;
				return 1;
			}
			if(t_ant < t_wall && fixed_velo == false){
				ants.update(t_ant,r_enc,index2,index3); // run an ant-to-ant collision with conservation of momentum
			}else if(t_ant < t_wall && fixed_velo == true){
				ants.update(t_ant,r_enc,index2,index3,velo); // run an ant-to-ant collision with fixed velocity
			}else{
				ants.update(R,r_enc,a,t_wall,index1);	// run ant-to-wall collisions for single particle
			}
		}else{
			ants.update(R,r_enc,a,t_wall,index1);	// run ant-to-wall collisions for single particle
		}

		data_file << ants; //write to file
		ants_in_nest = std::accumulate(ants.nest_flag.begin(),ants.nest_flag.end(),0); // keep track of number of ants in nest

		// Test to see if we've reached max iterations
		counter++;
		//std::cout << "Counter: " << counter << "/" << max_steps << std::endl;
		if(counter >= max_steps){
			//std::cout << "Maximum iterations reached -- trial terminated." << std::endl;
			std::cout << "Bad\t" << num_ants << "\t" << ants.event_time[0] << "\t" << ants.collisions[0] << std::endl;
			all_clear = true;
		}
		// If the trial ends when ant 0 leaves, check if ant 0 left
		if(exit_flag == true){
			if(ants.nest_flag.at(0)==0){
				//std::cout << "Ant 0 has left nest -- trial complete.\n";
				//std::cout << "Ants Left:\t" << std::accumulate(ants.nest_flag.begin(),ants.nest_flag.end(),0) << std::endl;
				std::cout << "Good\t" << num_ants << "\t" << ants.event_time[0] << "\t" << ants.collisions[0] << std::endl;
				all_clear = true;
			}
		// Else, see if all ants have left
		}else{
			//ants_in_nest = std::accumulate(ants.nest_flag.begin(),ants.nest_flag.end(),0);
			//std::cout << "Ants in Nest = " << ants_in_nest << std::endl;
			if(ants_in_nest == 0){
				std::cout << "Good\t" << num_ants << "\t" << ants.event_time[0] << "\t" << ants.collisions[0] << std::endl;
				//std::cout << "All ants have left nest -- trial complete.\n" << std::endl;
				all_clear = true;
			}			
		}
	}


	data_file.close();
	std::cout << "ANTS LEFT: " << ants_in_nest << std::endl;
	std::cout << "SIMULATION COMPLETE." << std::endl;

}

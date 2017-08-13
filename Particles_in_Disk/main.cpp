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
	bool single_particle = false; // true means we start a single particle starting IN APERTURE
	bool fixed_velo = true; // if true, ants will never change velocity, only heading
	bool start_in_center = true; // true means one ant 0 starts in center of aperture facing down
	bool exit_flag = true; // true means the trial ends when ant 0 leaves 
	int max_init = 500;  // max number of tries to initialize a given setup so no ants overlap 
	int max_steps= 10000; // max number of events in simulation (event is collision w/ wall or ant)

	int ants_in_nest; // used to calculate how many ants are in nest
	bool collision_flag = true; // true means collisions are on
	bool all_clear = false;  // used to break out of main simulation
	int counter = 0; // counts simulation steps

	long double t_wall,t_ant; // variables to store the time until a collision	
	int index1,index2,index3; // indices for particles involved in collision
	//int seed = 8;  // pass the same seed to 'populate' method to initialize trial the same way

	// Write parameters to file
	std::ofstream params("params.txt");
	params << "R\ta\tvelocity\tr_enc\tnum_ants\tcollision_flag\tcenter_flag\texit_flag" << std::endl;
	params << R << "\t" << a << "\t" << velo << "\t" << r_enc << "\t" << num_ants << "\t"
			<< collision_flag << "\t" << start_in_center << "\t" << exit_flag << std::endl;
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

		// Test to see if the trial is done
		if(exit_flag == true){
			if(ants.nest_flag.at(0)==0){
				//std::cout << "Ant 0 has left nest -- trial complete.\n";

				std::cout << "Good\t" << ants.event_time[0] << "\t" << ants.collisions[0] << std::endl;
				all_clear = true;
			}
		}else{
			//Calculate how many ants are still in nest
			ants_in_nest = std::accumulate(ants.nest_flag.begin(),ants.nest_flag.end(),0);
			if(ants_in_nest == 0){
				std::cout << "OK\t" << ants.event_time[0] << "\t" << ants.collisions[0] << std::endl;
				//std::cout << "All ants have left nest -- trial complete.\n" << std::endl;
				all_clear = true;
			}else{
				counter++;
				//std::cout << "Counter: " << counter << "/" << max_steps << std::endl;
				if(counter >= max_steps){
					//std::cout << "Maximum iterations reached -- trial terminated.\n" << std::endl;
					std::cout << "Bad\t" << ants.event_time[0] << "\t" << ants.collisions[0] << std::endl;
					all_clear = true;
				}
			}			
		}
	}


	data_file.close();
	// std::cout << "ANTS LEFT: " << ants_in_nest << std::endl;
	// std::cout << "SIMULATION COMPLETE." << std::endl;

}






/****************** NOTES *******************/
/*

g++ -std=c++11 main.cpp functions.cpp -o run_sim.exe

//prints ants to name.txt using operator overload
std::ofstream file("name.txt");
file << ants <<std::endl;
file.close();

//prints ants to cout (terminal) using operator overload
std::cout << ants;  //prints ants to terminal

// convert to string by passing to stringstream
std::stringstream strstream  //stringstream is a class strstream is an instance
//basically the above line does what std::ofstream or std::cout does
strstream << ants << std::endl;
std::string my_ant_string = strstream.str();
std::cout << my_ant_string;
"std::string" is the type for string. Like how "long double" is the type for long double

//ant update(ant J);
//ant initialize(ant J, int n, long double array[n]);

//for (auto const &x : t.x_positions){}
//drop the const and you can modify in place
*/
//	int L = 10;  // Box Length
//	int H = 10;  // Box Height
//	int n = 10;  // Number of ants
//	int counter;
//	long double value;

//	srand48(time(0));
	//printf("Time %ld\n",time(0));
	//srand48(1.0)
/*	
	for (counter=0;counter<=100;counter++){
		value = drand48();
		printf("%f\n",value);
	}
*/

//	ant J = {1.0,1.0,5.0,5.0,1,"jake","hanson","jake1","hanson1"};
//	ant K = {5.0,5.0,1.0,1.0,2,"nida","raja","nida1","raja1"};

//	ant ant_array[2];

//	long double init_pos[2] = {1.0,1.0};

//	initialize(J,2,init_pos);


/*
	printf("Before: %f\n",J.x);
	J = update(J);
	printf("After: %f\n",J.x);
*/	
//	printf("%f\t %f\t %f\t %f\t %s\t %s\n",J.x,J.y,J.v_x,J.v_y,J.name1,J.name2);
/*	if(strncmp("jake","jake",4) == 0){
		printf("YES!!!\n");
	};
*/
//	printf"
//	printf("%d\n",argc);
/*	for ( i=0; i<argc; i++){
		printf("%s\n", argv[i]);
	}
*/

/*
ant initialize(ant J, int n, long double array[n]){
	int i;

	for(i=0;i<n;i++){
		printf("%f\n",array[i]);
	}

	return J;
}
ant update(ant J){
	J.x = 1000;
	return J;
}*/


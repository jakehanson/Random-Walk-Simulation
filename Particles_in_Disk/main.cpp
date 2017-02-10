#include "header.h"
#include <fstream>


/* Function to define geometry, initialize ants, and run/write simulation */
int main(int argc, char** argv)
{
	/* Define Variables */
	double R = 2;  // Box Radius
	double a = 0.5;  // aperture size
	double temp = 5;  // temperature of box
	double r_enc = 0.1;  // encounter radius
	int num_ants = 30; // number of ants in simulation

	int max_init = 500;  // max number of tries to initialize a given setup 
	int max_steps= 1000; // max number of events in simulation (event is collision w/ wall or ant)

	int ants_in_nest; // used to calculate how many ants are in nest
	bool collision_flag = true; // true means collisions are on
	bool all_clear = false;  // used to break out of main simulation
	int counter = 0; // counts simulation steps

	double t_wall,t_ant; // variables to store the time until a collision	
	int index1,index2,index3; // indices for particles of interest
	//int seed = 8;  // pass the same seed to 'populate' method to initialize trial the same way

	// Write parameters to file
	std::ofstream params("params.txt");
	params << "R\ta\ttemp\tr_enc\tnum_ants\tcollision_flag" << std::endl;
	params << R << "\t" << a << "\t" << temp << "\t" << r_enc << "\t" << num_ants << "\t"
			<< collision_flag <<std::endl;
	params.close();

	// Open file to hold simulation data
	std::ofstream data_file("output.txt");
	data_file << "Name\t" << "x\t" << "y\t" << "v_x\t" << "v_y\t" << "event_time\t"
				<< "exit_time\t" << "collisions\t" << "in_nest" << std::endl;

	/* Initialize ants */	 	
	Ants ants = Ants(num_ants); // instance of Ants called ants
	try{
	ants.populate(R,temp,r_enc,max_init); // initialize positions and velocities
	}
	catch(std::runtime_error &e){
		std::cerr << "RUNTIME ERROR: " << e.what() << std::endl;
		return 1;
	}

	data_file << ants; //prints all of ants to name.txt using overloaded operator <<

	/* Run Sim */	
	while(all_clear == false){
		
		try{
			t_wall = get_t_wall(ants,R,r_enc,index1); // how long until colliding with wall?
			//std::cout << "Time: " << t_wall << std::endl;
		}
		catch(std::runtime_error &e){
			std::cerr << "RUNTIME ERROR: " << e.what() << std::endl;
			return 1;
		}
		if(collision_flag){
			t_ant = get_t_ant(ants,r_enc,t_wall,index2,index3); // how long before colliding with another ant?
		}
		if(collision_flag && t_ant < t_wall){
			ants.update(t_ant,r_enc,index2,index3); // run an ant collision
		}
		else{
			ants.update(R,r_enc,a,t_wall,index1);
		}

		//std::cout << "x-position of ant 1: " << ants.x_positions[0] << std::endl;		
		data_file << ants; //write to file

		//Calculate how many ants are still in nest
		ants_in_nest = std::accumulate(ants.nest_flag.begin(),ants.nest_flag.end(),0);
		if(ants_in_nest == 0){
			std::cout << "All Ants Have Left Nest" << std::endl;
			all_clear = true;
		}else{
			counter++;
			//std::cout << "Counter: " << counter << "/" << max_steps << std::endl;
			if(counter >= max_steps){
				std::cout << "Maximum Iterations Reached" << std::endl;
				all_clear = true;
			}
		}
	}

	// for(int a=0;a<ants.nest_flag.size();a++){
	// 	std::cout << ants.collisions[a] << std::endl;
	// }

	data_file.close();
	std::cout << "ANTS LEFT: " << ants_in_nest << std::endl;
	//std::cout << "SIMULATION COMPLETE." << std::endl;

}






/****************** NOTES *******************/
/*
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
"std::string" is the type for string. Like how "double" is the type for double

//ant update(ant J);
//ant initialize(ant J, int n, double array[n]);

//for (auto const &x : t.x_positions){}
//drop the const and you can modify in place
*/
//	int L = 10;  // Box Length
//	int H = 10;  // Box Height
//	int n = 10;  // Number of ants
//	int counter;
//	double value;

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

//	double init_pos[2] = {1.0,1.0};

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
ant initialize(ant J, int n, double array[n]){
	int i;

	for(i=0;i<n;i++){
		printf("%f\n",array[i]);
	}

	return J;
}
ant update(ant J){
	J.x = 1000;
	return J;
}

//#include <string>
//#include <fstream>
//#include <strstream>
//#include <iomanip>  //allows us to use cout::setprecision()

CAIRO Graphics for C or C++
output as svg for animation?
*/

/*
while (std::cin.good()){
	std::cin >> number;
	std::cout << (number*2) << std::endl;
}
*/
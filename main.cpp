#include <iostream>
#include <mpi.h>
#include <fstream>
#include <cstring>
#include "Grid.h"

double get_time(void) {		// Gets clock time with milliseconds accuracy and returns it as a double 
    struct timespec a;
    clock_gettime(CLOCK_MONOTONIC, &a);
    double t = (double) a.tv_nsec * 1e-6 + (double) a.tv_sec*1e3;
    return t;
}

void combine_output_files(MPI_Comm cartcomm, Grid &omega){
	int rank, size;
	int tag =0;
	MPI_Comm_rank(cartcomm, &rank);
	MPI_Comm_size(cartcomm, &size );
	if (rank == root){
		std::cout<< "------------------------------------------------"<<std::endl;
		std::cout<< "Writing output to file"<< std::flush;
		std::string local_input = omega.write_output();
		std::ofstream total_output("solution.txt");
		total_output << local_input;
		for (int source=0; source<size; ++source){
			if(source!=root){
				MPI_Status status;
				int buffsize=0;
				MPI_Probe(source, tag, cartcomm, &status);
				MPI_Get_count(&status, MPI_CHAR, &buffsize);
				char tmp_input[buffsize+1]={'\0'};
				MPI_Recv(&tmp_input, buffsize, MPI_CHAR, source, tag, cartcomm, &status);
				total_output << tmp_input;
			}
		}
		total_output.close();
		std::cout<< " done."<< std::endl;
	}
	else{
		std::string local_input = omega.write_output();
		int char_buff = local_input.length();
		MPI_Send(local_input.c_str(), char_buff, MPI_CHAR, root, tag, cartcomm);
	}
}
void print_usage(){
	std::cout<< "----------------------------------------------------------------------------------------"<<std::endl;
	std::cout<< "Usage: ./cgsolve nx ny iter np_x np_y [write_option]"<<std::endl;
	std::cout<< "At least 5 arguments are needed. The 6th argument is optional for writing file."<<std::endl;
	std::cout<< "np_x and np_y are number of processes in x and y direction."<<std::endl;
	std::cout<< "To write solution to a file, provide any arbitrary 6th argument."<<std::endl;
	std::cout<< "----------------------------------------------------------------------------------------"<<std::endl;
}

int main(int argc, char *argv[]){
	root = 0;
    double runtime;
    double start=0.0, end=0.0; 	//time
    int size, rank;
	MPI_Comm cartcomm;
    bool writefile=false;

    MPI_Init( &argc, &argv );
	MPI_Comm_size( MPI_COMM_WORLD, &size );
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );

	if(argc<6){
		if (rank==root){
			print_usage();
		}
		MPI_Finalize();
		return 0;
	}

	int nx_total = std::stoi(argv[1]); 	// nx global
	int ny_total = std::stoi(argv[2]);	// ny global
	int iter = std::stoi(argv[3]);	// number of iterations
	np_x = std::stoi(argv[4]);	// number of processes for x direction
	np_y = std::stoi(argv[5]);	// number of processes for y direction
	if(argv[6]){writefile = true;}

	if (np_x*np_y!=size){
		if (rank==root){
			std::cout<< "Invalid arguments: np_x*np_y must be equal to np (in this case np = "<<size<<")"<< std::endl;
			print_usage();
		}
		MPI_Finalize();
		return 0;
	}
	if (nx_total%np_x){
		if (rank==root){
			std::cout<< "Invalid arguments: nx must be divisible by np_x."<< std::endl;
			print_usage();
		}
		MPI_Finalize();
		return 0;
	}
	if (ny_total%np_y){
		if (rank==root){
			std::cout<< "Invalid arguments: ny must be divisible by np_y."<< std::endl;
			print_usage();	
		}
		MPI_Finalize();
		return 0;	
	}

	int nx_local = nx_total/np_x;	// local nx for each rank
	int ny_local = ny_total/np_y;	// local ny for each rank

	if (rank==root){
		std::cout<<"Number of processes: "<<np_x <<" x "<< np_y << std::endl;
		std::cout<<"Grid size per process: "<< nx_local << " x "<< ny_local << " = " << nx_local*ny_local << std::endl;
		std::cout<<"Mesh width: hx = "<< 1.0/(nx_total-1) << "\thy = "<< 1.0/(ny_total-1) << std::endl;
		std::cout<< "------------------------------------------------"<<std::endl;
	}

	int dims[2] = {np_x ,np_y};
	int coords[2];
	int periods[2] = {0,0};
	int reorder = 0;

	MPI_Cart_create( MPI_COMM_WORLD, 2, dims, periods, reorder, &cartcomm );
	MPI_Comm_rank( cartcomm, &rank );
	MPI_Cart_coords( cartcomm, rank, 2, coords );

	if (rank==root){
    	start = get_time();	// start clock
	}

	Grid omega(nx_local, ny_local, cartcomm, coords);	// create local Grid

	omega.set_BC_for_u(0.0);	// Boundary condition is set with zero initial guess for the interior
	omega.initialise_b();	// The right hand side matrix is created

	omega.cg(iter);		// CG is run 'iter' times
	MPI_Barrier(cartcomm);

	if(rank==root){		// cg iterations are done, so get end time
		end = get_time();	// end clock
   	 	runtime = end - start;		// Find elapsed time	
	}

	omega.set_solution(omega.b);
	omega.compute_error(omega.r,omega.b);
	double err_norm_local = omega.L2_norm(omega.r);
	double err_norm_global = 0.0;
	MPI_Reduce(&err_norm_local,&err_norm_global,1,MPI_DOUBLE,MPI_SUM,root,cartcomm);
	if(rank==root){
		std::cout << "Global error (L2_norm) = "<< err_norm_global << std::endl;
		std::cout << "Time to solution: " << runtime << " ms" <<std::endl;	
	}


	if(writefile){
		combine_output_files(cartcomm, omega);
	}

	MPI_Finalize();
	return 0;
}
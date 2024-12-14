#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <mpi.h>
#include "Grid.h"
#define PI 3.14159265358979
int np_x;
int np_y;
int root;

void Grid::set_BC_for_u(double init_val){	// Sets boundary condition and iterior of u
	init_guess(u, init_val);

	if (is_boundary[0]){for (int i=x_min; i<x_max; ++i){u(i-x_min,ny-1)=-cos(PI*hx*i);}}	// top boundary
	if (is_boundary[1]){for (int i=x_min; i<x_max; ++i){u(i-x_min,0)=cos(PI*hx*i);}}		// bottom boundary
	if (is_boundary[2]){for (int j=y_min; j<y_max; ++j){u(0,j-y_min)=cos(PI*hy*j);}}		// left boundary
	if (is_boundary[3]){for (int j=y_min; j<y_max; ++j){u(nx-1,j-y_min)=-cos(PI*hy*j);}}	// right boundary
}

void Grid::initialise_b(){		// Initialise the b array
	for (int i=x_min; i<x_max; ++i){
		for (int j=y_min; j<y_max; ++j){
			b(i-x_min,j-y_min)= 2*PI*PI*cos(PI*hx*i)*cos(PI*hy*j);
		}
	}
}

double Grid::product(const Data &arr1, const Data &arr2){	// vector-vector dot product
	double sum = 0;
	for (int i=1; i<nx-1; ++i){		// interior
		for (int j=1; j<ny-1; ++j){
			sum += arr1(i,j)*arr2(i,j);
		}
	}	
	if(!is_boundary[0]){for (int i=1; i<nx-1; ++i){sum += arr1(i,ny-1)*arr2(i,ny-1);}}	// interior of top boundary
	if(!is_boundary[1]){for (int i=1; i<nx-1; ++i){sum += arr1(i,0)*arr2(i,0);}}		// interior of bottom boundary
	if(!is_boundary[2]){for (int j=1; j<ny-1; ++j){sum += arr1(0,j)*arr2(0,j);}}		// interior of left boundary 
	if(!is_boundary[3]){for (int j=1; j<ny-1; ++j){sum += arr1(nx-1,j)*arr2(nx-1,j);}}	// interior of right boundary
	if(!is_boundary[0] && !is_boundary[2]){sum += arr1(0,ny-1)*arr2(0,ny-1);}			// top left corner
	if(!is_boundary[0] && !is_boundary[3]){sum += arr1(nx-1,ny-1)*arr2(nx-1,ny-1);}		// top right corner
	if(!is_boundary[1] && !is_boundary[2]){sum += arr1(0,0)*arr2(0,0);}					// bottom left corner
	if(!is_boundary[1] && !is_boundary[3]){sum += arr1(nx-1,0)*arr2(nx-1,0);}			// bottom right corner

	return sum;
}

Data Grid::product(double val, const Data &arr){	// scalar-vector product
	Data tmp(nx,ny);
	for (int i=1; i<nx-1; ++i){		// interior
		for (int j=1; j<ny-1; ++j){
			tmp(i,j) = val*arr(i,j);
		}
	}
	if(!is_boundary[0]){for (int i=1; i<nx-1; ++i){tmp(i,ny-1) = val*arr(i,ny-1);}}		// interior of top boundary
	if(!is_boundary[1]){for (int i=1; i<nx-1; ++i){tmp(i,0) = val*arr(i,0);}}			// interior of bottom boundary
	if(!is_boundary[2]){for (int j=1; j<ny-1; ++j){tmp(0,j) = val*arr(0,j);}}			// interior of left boundary 
	if(!is_boundary[3]){for (int j=1; j<ny-1; ++j){tmp(nx-1,j) = val*arr(nx-1,j);}}		// interior of right boundary
	if(!is_boundary[0] && !is_boundary[2]){tmp(0,ny-1) = val*arr(0,ny-1);}				// top left corner
	if(!is_boundary[0] && !is_boundary[3]){tmp(nx-1,ny-1) = val*arr(nx-1,ny-1);}		// top right corner
	if(!is_boundary[1] && !is_boundary[2]){tmp(0,0) = val*arr(0,0);}					// bottom left corner
	if(!is_boundary[1] && !is_boundary[3]){tmp(nx-1,0) = val*arr(nx-1,0);}				// bottom right corner
	return tmp;
}

Data Grid::stencil_mvm(const Data &arr){	// Multiplies matrix A with arr
	Data tmp(nx,ny);
	MPI_Request requests[4];

//------------------------------ set boundaries ---------------------------------
	if(!is_boundary[0]){for (int i=0; i<nx; ++i){top[i] = arr(i,ny-1);}}	// top
	if(!is_boundary[1]){for (int i=0; i<nx; ++i){bottom[i] = arr(i,0);}}	// bottom
	if(!is_boundary[2]){for (int j=0; j<ny; ++j){left[j] = arr(0,j);}}		// left
	if(!is_boundary[3]){for (int j=0; j<ny; ++j){right[j] = arr(nx-1,j);}}	// right

//------------------- send boundaries -------------------------------------------------------------------------------
	if(!is_boundary[0]){MPI_Isend(top, 1, Ghost_horizontal, nbr_rank[0], 0, cartcomm, &requests[0]);}		// send top boundary to top neighbor
	if(!is_boundary[1]){MPI_Isend(bottom, 1, Ghost_horizontal, nbr_rank[1], 1, cartcomm, &requests[1]);}	// send bottom boundary to bottom neighbor
	if(!is_boundary[2]){MPI_Isend(left, 1, Ghost_vertical, nbr_rank[2], 2, cartcomm, &requests[2]);}		// send left boundary to left neighbor
	if(!is_boundary[3]){MPI_Isend(right, 1, Ghost_vertical, nbr_rank[3], 3, cartcomm, &requests[3]);}		// send right boundary to right neighbor

//---------------------- recv ghosts ---------------------------------------------------------------------------------
	if(!is_boundary[0]){MPI_Irecv(top_ghost, 1, Ghost_horizontal, nbr_rank[0], 1, cartcomm, &requests[0]);}		// recv top from top neighbor
	if(!is_boundary[1]){MPI_Irecv(bottom_ghost, 1, Ghost_horizontal, nbr_rank[1], 0, cartcomm, &requests[1]);}	// recv bottom from bottom neighbor
	if(!is_boundary[2]){MPI_Irecv(left_ghost, 1, Ghost_vertical, nbr_rank[2], 3, cartcomm, &requests[2]);}		// recv left from left neighbor
	if(!is_boundary[3]){MPI_Irecv(right_ghost, 1, Ghost_vertical, nbr_rank[3], 2, cartcomm, &requests[3]);}		// recv right from right neighbor

//------------------------------- find iterior ------------------------------------
	for (int i=1; i<nx-1; ++i){
		for (int j=1; j<ny-1; ++j){
			tmp(i,j) = (- arr(i-1,j) + 2*arr(i,j) - arr(i+1,j))/(hx*hx) 
					 + (- arr(i,j-1) + 2*arr(i,j) - arr(i,j+1))/(hy*hy);
		}
	}
//------------------------------ wait for ghosts to arrive ------------------------------
	if(!is_boundary[0]){MPI_Wait(&requests[0], MPI_STATUSES_IGNORE);}
	if(!is_boundary[1]){MPI_Wait(&requests[1], MPI_STATUSES_IGNORE);}
	if(!is_boundary[2]){MPI_Wait(&requests[2], MPI_STATUSES_IGNORE);}
	if(!is_boundary[3]){MPI_Wait(&requests[3], MPI_STATUSES_IGNORE);}

//============================== find boundaries =============================================
// ------------------------------------------------------------ interior of top boundary
	if(!is_boundary[0]){
		for (int i=1; i<nx-1; ++i){
			tmp(i,ny-1) = (- arr(i-1,ny-1) + 2*arr(i,ny-1) - arr(i+1,ny-1))/(hx*hx) 
			 + (- arr(i,ny-2) + 2*arr(i,ny-1) - top_ghost[i])/(hy*hy);
		}
	}
//------------------------------------------------------------- interior of bottom boundary
	if(!is_boundary[1]){
		for (int i=1; i<nx-1; ++i){
				tmp(i,0) = (- arr(i-1,0) + 2*arr(i,0) - arr(i+1,0))/(hx*hx) 
						 + (- bottom_ghost[i] + 2*arr(i,0) - arr(i,1))/(hy*hy);
		}
	}
//------------------------------------------------------------- interior of left boundary	
	if(!is_boundary[2]){
		for (int j=1; j<ny-1; ++j){
			tmp(0,j) = (- left_ghost[j] + 2*arr(0,j) - arr(1,j))/(hx*hx) 
					 + (- arr(0,j-1) + 2*arr(0,j) - arr(0,j+1))/(hy*hy);		
		}
	}
//------------------------------------------------------------- interior of right boundary
	if(!is_boundary[3]){
		for (int j=1; j<ny-1; ++j){
			tmp(nx-1,j) = (- arr(nx-2,j) + 2*arr(nx-1,j) - right_ghost[j])/(hx*hx) 
			 + (- arr(nx-1,j-1) + 2*arr(nx-1,j) - arr(nx-1,j+1))/(hy*hy);
		}
	}
//------------------------------------------------------------- corners
// top left corner
	if(!is_boundary[0] && !is_boundary[2]){
		tmp(0,ny-1) = (- left_ghost[ny-1] + 2*arr(0,ny-1) - arr(1,ny-1))/(hx*hx)
		 + (- arr(0,ny-2) + 2*arr(0,ny-1) - top_ghost[0])/(hy*hy);
	}
// top right corner
	if(!is_boundary[0] && !is_boundary[3]){
		tmp(nx-1,ny-1) = (- arr(nx-2,ny-1) + 2*arr(nx-1,ny-1) - right_ghost[ny-1])/(hx*hx) 
		 + (- arr(nx-1,ny-2) + 2*arr(nx-1,ny-1) - top_ghost[nx-1])/(hy*hy);
	}
// bottom left corner
	if(!is_boundary[1] && !is_boundary[2]){
		tmp(0,0) = (- left_ghost[0] + 2*arr(0,0) - arr(1,0))/(hx*hx) 
			 + (- bottom_ghost[0] + 2*arr(0,0) - arr(0,1))/(hy*hy);
	}
// bottom right corner
	if(!is_boundary[1] && !is_boundary[3]){	 
		tmp(nx-1,0) = (- arr(nx-2,0) + 2*arr(nx-1,0) - right_ghost[0])/(hx*hx) 
				 + (- bottom_ghost[nx-1] + 2*arr(nx-1,0) - arr(nx-1,1))/(hy*hy);
	}

	return tmp;
}

void Grid::init_guess(Data &arr, double val){
	for (int i=0; i<nx; ++i){
		for (int j=0; j<ny; ++j){
			arr(i,j) = val;
		}
	}
}

void Grid::set_solution(Data &arr){		// Sets arr to the analyticial solution (for finding error)
	for (int i=x_min; i<x_max; ++i){
		for (int j=y_min; j<y_max; ++j){
			arr(i-x_min,j-y_min)=cos(PI*hx*i)*cos(PI*hy*j);
		}
	}
}

void Grid::compute_error(Data &err, const Data &sol){ // Computes the error at each node w.r.t the analytical solution (error is saved in r)
	for (int i=0; i<nx; ++i){
		for (int j=0; j<ny; ++j){
			err(i,j) = u(i,j) - sol(i,j);
		}
	}
}

double Grid::L2_norm(const Data &arr){ // returns the discrete L2 norm of arr as a double
	double square_sum=0.0;
	for (int i=0; i<nx; ++i){
		for (int j=0; j<ny; ++j){
			square_sum += arr(i,j)*arr(i,j);
		}
	}
	return sqrt(square_sum/(nx*np_x*ny*np_y));
}

void Grid::cg(int iter){
	double alpha, beta;
	double res=10, rr_new=0.0, rr_old, rr_new_local, dAd, dAd_local;
	int it;
	Data Ad(nx, ny);
	Ad.set_is_boundary(is_boundary);
	r = b - stencil_mvm(u);
	rr_new_local = product(r,r);
	MPI_Reduce(&rr_new_local,&rr_new,1,MPI_DOUBLE,MPI_SUM,root,cartcomm);
	MPI_Bcast(&rr_new,1,MPI_DOUBLE,root,cartcomm);
	d = r;
	if (rank==root){
		std::cout<< std::endl<< "#Iter\tGlobal Residual"<<std::endl;
		std::cout<< "------------------------------------------------"<<std::endl;
	}
	for (it=1; it<=iter; ++it){
		if (rank==root){
			res = sqrt(rr_new/(nx*np_x*ny*np_y)); // discrete L2 norm
			std::cout<<it << "\t" << res << std::endl;
		}
		Ad = stencil_mvm(d); 
		dAd_local = product(d,Ad);
		dAd=0.0;
		MPI_Reduce(&dAd_local,&dAd,1,MPI_DOUBLE,MPI_SUM,root,cartcomm);
		MPI_Bcast(&dAd,1,MPI_DOUBLE,root,cartcomm);
		alpha = rr_new / dAd;
		u = u + product(alpha,d);
		r = r - product(alpha,Ad);
		rr_old = rr_new;
		rr_new_local = product(r,r);
		rr_new=0.0;
		MPI_Reduce(&rr_new_local,&rr_new,1,MPI_DOUBLE,MPI_SUM,root,cartcomm);
		MPI_Bcast(&rr_new,1,MPI_DOUBLE,root,cartcomm);
		beta = rr_new / rr_old ;
		d = r + product(beta,d);
	}
	if(rank==root){
		std::cout<< "------------------------------------------------"<<std::endl;
		std::cout << "Global residual (L2-norm) = " << res << std::endl;
	}
}

std::string Grid::write_output() const{	// Function to write soluton u to a .txt file for gnuplot
	std::ofstream local_output;
	std::string filename = "solution_" + std::to_string(coords[0]) + "_" + std::to_string(coords[1]) + ".txt";
	local_output.open(filename);
	local_output << "#X\tY\tu \t from process ("<<coords[0]<<","<<coords[1]<<")"<< std::endl;
	local_output << std::fixed << std::setprecision(5);
	for (int i=x_min; i<x_max; ++i){
		for (int j=y_min; j<y_max; ++j){
			local_output << i*hx << "\t" << j*hy << "\t" << u(i-x_min,j-y_min) << std::endl;
		}
	}
	local_output.close();
	std::ostringstream local_input_string;
	std::ifstream local_input(filename);
	local_input_string << local_input.rdbuf();
	return local_input_string.str();
}

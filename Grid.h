#ifndef GRID_H
#define GRID_H
#include "Data.h"
#include <mpi.h>
extern int np_x;
extern int np_y;
extern int root;

struct Grid{	// The data for each process is stored in this class, so the grid data is local to every process

	Grid()= default;
	inline Grid(int Nx, int Ny, const MPI_Comm &cartcommin, const int coordin[]);
	~Grid(){
		delete[] top; delete[] top_ghost; delete[] bottom; delete[] bottom_ghost; 
		delete[] left; delete[] left_ghost; delete[] right; delete[] right_ghost; 
	}
	int nx, ny;	// local number of grid points
	int coords[2];
	int x_min, x_max; // indices of this local grid on the global grid
	int y_min, y_max;
	int grid_size = nx*ny;	// Total number of local grid points
	MPI_Comm cartcomm;
	int rank;
	double global_res;

	double hx = 1.0/(np_x*nx-1);	// Horizontal spacing between nodes
	double hy = 1.0/(np_y*ny-1); 	// Vertical spacing between nodes
	Data u, b, r, d;		// Data on the grid: u(solution), b (right hand side), r(residual), d(search direction)

	bool is_boundary[4]={false, false, false, false}; // {top, bottom, left, right} to check if this domain is in the boundary of global domain 
	int nbr_coords[4][2];	// cart coords of neighbor processes
	int nbr_rank[4];

	double *top, *top_ghost, *bottom, *bottom_ghost, *left, *left_ghost, *right, *right_ghost; // boundary nodes
	MPI_Datatype Ghost_horizontal;
	MPI_Datatype Ghost_vertical;

	void set_BC_for_u(double init_val);
	void initialise_b();
	void init_guess(Data &arr, double val);
	double product(const Data &arr1, const Data &arr2);
	Data product(double val, const Data &arr);
	Data stencil_mvm(const Data &arr);
	void cg(int iter);
	std::string write_output() const;
	void set_solution(Data &arr);
	void compute_error(Data &err, const Data &sol);
	double L2_norm(const Data &arr);
};

Grid::Grid(int Nx, int Ny, const MPI_Comm &cartcommin, const int coordin[]): nx(Nx), ny(Ny), cartcomm(cartcommin),
	u(nx, ny), 
	b(nx, ny), 
	r(nx, ny),
	d(nx, ny),
	top(new double[nx]()), top_ghost(new double[nx]()),
	bottom(new double[nx]()), bottom_ghost(new double[nx]()),
	left(new double[ny]()), left_ghost(new double[ny]()),
	right(new double[ny]()), right_ghost(new double[ny]())
	{
	coords[0]=coordin[0];
 	coords[1]=coordin[1];

 	// checking if this process is in the boundary of the gloabl domain
 	if (coords[1]==np_y-1){ // top
		is_boundary[0]=true;
	}
	if (coords[1]==0){ // bottom
		is_boundary[1]=true;
	}
	if (coords[0]==0){	// left
		is_boundary[2]=true;
	}
	if (coords[0]==np_x-1){	//right
		is_boundary[3]=true;
	}

	x_min = coords[0]*nx;
	x_max = x_min + nx ; // this is max but the loops dont run over this node
	y_min = coords[1]*ny;
	y_max = y_min + ny ;

	// finding the neighbor coordinates and ranks
	MPI_Cart_rank(cartcomm, coords, &rank);
	if(!is_boundary[0]){
		nbr_coords[0][0]=coords[0];
		nbr_coords[0][1]=coords[1]+1;
		MPI_Cart_rank(cartcomm, nbr_coords[0], &nbr_rank[0]);
	}
	if(!is_boundary[1]){
		nbr_coords[1][0]=coords[0];
		nbr_coords[1][1]=coords[1]-1;
		MPI_Cart_rank(cartcomm, nbr_coords[1], &nbr_rank[1]);
	}
	if(!is_boundary[2]){
		nbr_coords[2][0]=coords[0]-1;
		nbr_coords[2][1]=coords[1];
		MPI_Cart_rank(cartcomm, nbr_coords[2], &nbr_rank[2]);
	}
	if(!is_boundary[3]){
		nbr_coords[3][0]=coords[0]+1;
		nbr_coords[3][1]=coords[1];
		MPI_Cart_rank(cartcomm, nbr_coords[3], &nbr_rank[3]);
	}

	u.set_is_boundary(is_boundary);
	b.set_is_boundary(is_boundary);
	r.set_is_boundary(is_boundary);
	d.set_is_boundary(is_boundary);

	MPI_Type_contiguous(nx,MPI_DOUBLE,&Ghost_horizontal);	
	MPI_Type_contiguous(ny,MPI_DOUBLE,&Ghost_vertical);
	MPI_Type_commit(&Ghost_horizontal);
	MPI_Type_commit(&Ghost_vertical);
}

#endif
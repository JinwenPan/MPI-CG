#ifndef DATA_H
#define DATA_H
#include <memory>
#include <cassert>

struct Data{	

	Data()= default;
	Data(const int Nx, const int Ny) : nx(Nx), ny(Ny), arr(new double[nx*ny]()){}

	inline double& operator ()(int i, int j) {
		// assert(i<nx && j<ny);
		return arr[i*ny +j];
	} 		
	inline double operator ()(int i, int j) const {
		// assert(i<nx && j<ny); 
		return arr[i*ny +j];
	}
	Data operator +(const Data &other){
		Data tmp(nx, ny);
		for (int i=1; i<nx-1; ++i){
			for (int j=1; j<ny-1; ++j){
				tmp(i,j) = (*this)(i,j) + other(i,j);
			}
		}
		if(!is_boundary[0]){for (int i=1; i<nx-1; ++i){tmp(i,ny-1) = (*this)(i,ny-1) + other(i,ny-1);}}
		if(!is_boundary[1]){for (int i=1; i<nx-1; ++i){tmp(i,0) = (*this)(i,0) + other(i,0);}}
		if(!is_boundary[2]){for (int j=1; j<ny-1; ++j){tmp(0,j) = (*this)(0,j) + other(0,j);}}
		if(!is_boundary[3]){for (int j=1; j<ny-1; ++j){tmp(nx-1,j) = (*this)(nx-1,j) + other(nx-1,j);}}
		if(!is_boundary[0] && !is_boundary[2]){tmp(0,ny-1) = (*this)(0,ny-1) + other(0,ny-1);}
		if(!is_boundary[0] && !is_boundary[3]){tmp(nx-1,ny-1) = (*this)(nx-1,ny-1) + other(nx-1,ny-1);}
		if(!is_boundary[1] && !is_boundary[2]){tmp(0,0) = (*this)(0,0) + other(0,0);}
		if(!is_boundary[1] && !is_boundary[3]){tmp(nx-1,0) = (*this)(nx-1,0) + other(nx-1,0);}
		return tmp;
	}
	Data operator -(const Data &other){
		Data tmp(nx, ny);
		for (int i=1; i<nx-1; ++i){
			for (int j=1; j<ny-1; ++j){
				tmp(i,j) = (*this)(i,j) - other(i,j);
			}
		}
		if(!is_boundary[0]){for (int i=1; i<nx-1; ++i){tmp(i,ny-1) = (*this)(i,ny-1) - other(i,ny-1);}}
		if(!is_boundary[1]){for (int i=1; i<nx-1; ++i){tmp(i,0) = (*this)(i,0) - other(i,0);}}
		if(!is_boundary[2]){for (int j=1; j<ny-1; ++j){tmp(0,j) = (*this)(0,j) - other(0,j);}}
		if(!is_boundary[3]){for (int j=1; j<ny-1; ++j){tmp(nx-1,j) = (*this)(nx-1,j) - other(nx-1,j);}}
		if(!is_boundary[0] && !is_boundary[2]){tmp(0,ny-1) = (*this)(0,ny-1) - other(0,ny-1);}
		if(!is_boundary[0] && !is_boundary[3]){tmp(nx-1,ny-1) = (*this)(nx-1,ny-1) - other(nx-1,ny-1);}
		if(!is_boundary[1] && !is_boundary[2]){tmp(0,0) = (*this)(0,0) - other(0,0);}
		if(!is_boundary[1] && !is_boundary[3]){tmp(nx-1,0) = (*this)(nx-1,0) - other(nx-1,0);}
		return tmp;
	}
	Data& operator =(const Data &other){
		for (int i=1; i<nx-1; ++i){
			for (int j=1; j<ny-1; ++j){
				(*this)(i,j) = other(i,j);
			}
		}
		if(!is_boundary[0]){for (int i=1; i<nx-1; ++i){(*this)(i,ny-1) = other(i,ny-1);}}
		if(!is_boundary[1]){for (int i=1; i<nx-1; ++i){(*this)(i,0) = other(i,0);}}
		if(!is_boundary[2]){for (int j=1; j<ny-1; ++j){(*this)(0,j) = other(0,j);}}
		if(!is_boundary[3]){for (int j=1; j<ny-1; ++j){(*this)(nx-1,j) = other(nx-1,j);}}
		if(!is_boundary[0] && !is_boundary[2]){(*this)(0,ny-1) = other(0,ny-1);}
		if(!is_boundary[0] && !is_boundary[3]){(*this)(nx-1,ny-1) = other(nx-1,ny-1);}
		if(!is_boundary[1] && !is_boundary[2]){(*this)(0,0) = other(0,0);}
		if(!is_boundary[1] && !is_boundary[3]){(*this)(nx-1,0) = other(nx-1,0);}

		return *this;
	}
	void set_is_boundary(const bool isBoundary[]){
		for (int i=0; i<4; ++i){
			is_boundary[i]=isBoundary[i];
		}
	}
	bool is_boundary[4];

private:
	int nx, ny;
	std::shared_ptr<double[]> arr;
};

#endif
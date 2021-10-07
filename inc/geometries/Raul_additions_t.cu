#include <iostream>
#include <string>
#include <vector>
#include <functional>
#include <map>
#include <math.h>

#include "dg/algorithm.h"
#include "dg/geometries/geometries.h"

	double function(double x, double y, double z) {
	return sin(y);
	}
	
	double function_prime(double x, double y, double z) {
	return cos(y);
	}
	
	double vx(double x, double y, double z) {
	return sin(x);
	}
	double vy(double x, double y, double z) {
	return sin(y);
	}
	double vz(double x, double y, double z) {
	return sin(z);
	}
	double divergence(double x, double y, double z) {
	return cos(x)+cos(y);
	}

int main () {
	
	int nR=25;	//Size of grid in Radial direction
	int nZ=25;	//Size of grid in Z direction
	int nP=10;	//Size of grid in toroidal direction
	int n=1; //Order of the Legendre polynomia
	
	
	//for (int i=1; i<nR; i++){
	/*
	for (int i=1; i<10;i++)
	{
	double N=5*i;
	double p=1;
	double h=M_PI/N;
	
	dg::RealCartesianGrid3d<double> geom3d(0,M_PI,0,M_PI,0,M_PI, p, N, N, N, dg::PER,dg::PER,dg::PER);	
	dg::HVec weights=dg::create::volume(geom3d);
	
	
	
	const dg::HVec F=dg::construct<dg::DVec>(dg::evaluate(function, geom3d));
	const dg::HVec F_prime=dg::construct<dg::DVec>(dg::evaluate(function_prime, geom3d));
	dg::geo::Nablas nabla(geom3d);
	dg::HVec f=F;
	nabla.GradPerp_Z(F, f);
    
    dg::HVec error, abs_error;
    error=F_prime;   
    abs_error=error;
    dg::blas1::axpby(1,f,-1,error);
    dg::blas1::transform( error, abs_error, dg::ABS<double>());   
    double norm= sqrt(dg::blas2::dot(abs_error, weights, abs_error));
	std::cout<<"Norm ="<<norm<<" and N="<<N<<"\n";	
}

*/

	dg::RealCartesianGrid3d<double> geom3d(0,M_PI,0,M_PI,0,M_PI, n, nR, nZ, nP, dg::PER,dg::PER,dg::PER);	
	dg::DVec weights=dg::create::volume(geom3d);
	dg::DVec VX=dg::construct<dg::DVec>(dg::evaluate(vx, geom3d));
	dg::DVec VY=dg::construct<dg::DVec>(dg::evaluate(vy, geom3d));
	//dg::HVec VZ=dg::construct<dg::DVec>(dg::evaluate(vz, geom3d));
	dg::geo::Nablas<dg::RealCartesianGrid3d<double>, dg::DMatrix, dg::DVec> nabla(geom3d); //I WOULD NEED A DAMN MAGNETIC FIELD TO DEFINE IT, AND THIS TEST DOES NOT HAVE IT
	dg::DVec F;
	F=VX;
	dg::DVec F_theory=dg::construct<dg::DVec>(dg::evaluate(divergence, geom3d));
	nabla.div (VX, VY, F);
	
	/*
	 std::cout<<"Size of vectors: "<<F.size()<<"\n";
	for (int i=0; i<N*p; i++)
	{for (int j=0; j<N*p; j++)
		{for (int k=0; k<N; k++)
	{std::cout<<F[i*N*N+j*N+k]<<"  ";}
	std::cout<<"\n";
	}
	std::cout<<"\n"<<"\n";
}
*/

	dg::DVec error=F_theory;
	dg::DVec abs_error=error;
	dg::blas1::axpby(1,F,-1,error);
	dg::blas1::transform( error, abs_error, dg::ABS<double>());   
	double norm= sqrt(dg::blas2::dot(abs_error, weights, abs_error));
	std::cout<<"Norm:"<<norm<<"\n";
	std::cout<<"Finished"<<"\n";


}

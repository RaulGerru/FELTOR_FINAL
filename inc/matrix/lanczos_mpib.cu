// #define DG_DEBUG
#include <iostream>
#include <iomanip>
#include <mpi.h>

#include "lanczos.h"

const double lx = 2.*M_PI;
const double ly = 2.*M_PI;
dg::bc bcx = dg::DIR;
dg::bc bcy = dg::PER;
const double alpha = -0.5;
const double m = 4.;
const double n = 4.;

double lhs( double x, double y) {return sin(m*x)*sin(n*y);}
double rhs( double x, double y){ return (1.-(m*m+n*n)*alpha)*sin(m*x)*sin(n*y);}

using Matrix = dg::MDMatrix;
using Container = dg::MDVec;
using HDiaMatrix = cusp::dia_matrix<int, double, cusp::host_memory>;
using HCooMatrix = cusp::coo_matrix<int, double, cusp::host_memory>;

int main(int argc, char * argv[])
{
    MPI_Init(&argc, &argv);
    dg::Timer t;
    unsigned n, Nx, Ny;
//     std::cout << "# Type n, Nx and Ny\n";
//     std::cin >> n >> Nx >> Ny;
    MPI_Comm comm;
    dg::mpi_init2d( bcx, bcy, n, Nx, Ny, comm);
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
    if(rank==0)
    {
        std::cout <<"# You typed\n"
              <<"n:  "<<n<<"\n"
              <<"Nx: "<<Nx<<"\n"
              <<"Ny: "<<Ny<<std::endl;
    }
    unsigned max_iter;
    double eps = 1e-6;
    if(rank==0) std::cout << "# Type in max_iter and eps\n";
    if(rank==0) std::cin >> max_iter>> eps;
    if(rank==0) std::cout <<"# You typed\n"
              <<"max_iter:  "<<max_iter<<"\n"
              <<"eps: "<<eps <<std::endl;
    MPI_Bcast(  &max_iter,1 , dg::getMPIDataType<unsigned>(), 0, comm);
    MPI_Bcast(  &eps,1 , dg::getMPIDataType<double>(), 0, comm);

    dg::RealCartesianMPIGrid2d<double> grid( 0., lx, 0, ly, n, Nx, Ny, bcx, bcy, comm);

    const Container w2d = dg::create::weights( grid);

    Container x = dg::evaluate( lhs, grid), b(x), error(x), xexac(x);
    Container bexac = dg::evaluate( rhs, grid);
    dg::Helmholtz<dg::aRealMPIGeometry2d<double>, Matrix, Container> A( grid, alpha, dg::centered);
    {
        t.tic();
        dg::Lanczos< Container > lanczos(x, max_iter);
        lanczos.set_verbose(true);
        t.toc();
        if(rank==0) std::cout << "# Lanczos creation took "<< t.diff()<<"s   \n";

        HDiaMatrix T;
        if(rank==0) std::cout << "Lanczos:\n";

        t.tic();
        T = lanczos( A, x, b, 1., eps, true); 
        t.toc();

        if(rank==0) std::cout << "    iter: "<< lanczos.get_iter() << "\n";
        if(rank==0) std::cout << "    time: "<< t.diff()<<"s \n";
        dg::blas1::axpby(-1.0, bexac, 1.0, b,error);
        if(rank==0) std::cout << "    # Relative error between b=||x||_2 V^T T e_1 and b: \n";   
        double temp = dg::blas2::dot( w2d, error)/dg::blas2::dot( w2d, bexac);
        if(rank==0) std::cout << "    error: " << sqrt(temp) << " \n";   

        if(rank==0) std::cout << "\nM-Lanczos:\n";
        x = dg::evaluate( lhs, grid);
        t.tic();
        T = lanczos(A, x, b, w2d, eps, true); 
        t.toc();
        if(rank==0) std::cout << "    iter: "<< lanczos.get_iter() << "\n";
        if(rank==0) std::cout << "    time: "<< t.diff()<<"s \n";
        dg::blas1::axpby(-1.0, bexac, 1.0, b,error);
        if(rank==0) std::cout << "    # Relative error between b=||x||_M V^T T e_1 and b: \n";  
        temp = dg::blas2::dot( w2d, error)/dg::blas2::dot( w2d, bexac);
        if(rank==0) std::cout << "    error: " << sqrt(temp) << " \n";   

    } 
    {
        if(rank==0) std::cout << "\nM-CG: \n";
        t.tic();
        dg::MCG<Container> mcg(x, max_iter);
        t.toc();
        if(rank==0) std::cout << "#    M-CG creation took "<< t.diff()<<"s   \n";
        dg::blas1::scal(x, 0.0); //initialize with zero
        t.tic();
        HDiaMatrix T = mcg(A, x, b, w2d, eps, 1., true); 
        t.toc();

        dg::blas1::axpby(-1.0, xexac, 1.0, x, error);
        if(rank==0) std::cout << "    iter: "<< mcg.get_iter() << "\n";
        if(rank==0) std::cout << "    time: "<< t.diff()<<"s \n";
        if(rank==0) std::cout << "    # Relative error between x= R T^{-1} e_1 and x: \n";
        double temp = dg::blas2::dot( w2d, error)/dg::blas2::dot( w2d, bexac);
        if(rank==0) std::cout << "    error: " << sqrt(temp) << " \n";   
    }

    MPI_Finalize();
    return 0;
}

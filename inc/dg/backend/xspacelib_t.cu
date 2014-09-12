#include <iostream>

#include <cusp/print.h>

#include "file/read_input.h"
#include "xspacelib.cuh"
#include "cg.h"


double sine( double x, double y){ return 2.*sin(x)*sin(y);}
double solution( double x, double y){ return sin(x)*sin(y);}

using namespace std;
using namespace dg;

int main()
{
    const dg::Grid2d<double> grid( 0, 2.*M_PI, 0, 2.*M_PI, 3, 80, 80, dg::DIR, dg::PER);
    const double eps = 1e-6;

    DVec b = dg::evaluate( sine, grid), x( b.size(), 0);
    const DVec sol = evaluate( solution, grid);
    DVec w2d = create::weights( grid);
    DVec v2d = create::inv_weights( grid);

    Polarisation2dX<DVec> polarisation ( grid);
    DMatrix laplace = create::laplacianM( grid, not_normed);

    CG<DVec > cg( x, x.size());
    cout << "Test of w2d: "<<blas2::dot( w2d, b)<<endl;
    blas2::symv( w2d, b, b);
    cout << "Test of v2d: "<<blas2::dot( v2d, b)<<endl;

    cg( laplace, x, b, v2d, eps);

    blas1::axpby( 1., x, -1, sol, x);
    cout << "Norm "<< blas2::dot(w2d,x)<<endl;
    cout << "Rel error "<<sqrt( blas2::dot(w2d, x))/sqrt(  blas2::dot( w2d, sol) )<<endl;
    cout << "Bracket "<< blas2::dot(w2d,sol)<<endl;
    cout << "Rel error "<<(blas2::dot(w2d,sol)-M_PI*M_PI)/M_PI/M_PI<<endl;
    return 0;
}

#include <iostream>
#include <iomanip>
#include <functional>

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

#include "backend/typedefs.h"
#include "topology/evaluation.h"
#include "topology/grid.h"
#include "arakawa.h"
#include "runge_kutta.h"
#include "adaptive.h"


//![function]
void rhs(double t, const std::array<double,2>& y, std::array<double,2>& yp,
        double damping, double omega_0, double omega_drive){
    //damped driven harmonic oscillator
    // x -> y[0] , v -> y[1]
    yp[0] = y[1];
    yp[1] = -2.*damping*omega_0*y[1] - omega_0*omega_0*y[0] + sin(omega_drive*t);
}
//![function]

std::array<double, 2> solution( double t, double damping, double omega_0,
        double omega_drive)
{
    double tmp1 = (2.*omega_0*damping);
    double tmp2 = (omega_0*omega_0 - omega_drive*omega_drive)/omega_drive;
    double amp = 1./sqrt( tmp1*tmp1 + tmp2*tmp2);
    double phi = atan( 2.*omega_drive*omega_0*damping/
            (omega_drive*omega_drive-omega_0*omega_0));

    double x = amp*sin(omega_drive*t+phi)/omega_drive;
    double v = amp*cos(omega_drive*t+phi);
    return {x,v};
}

int main()
{
    std::cout << "Program to test correct implementation of adaptive methods in adaptive.h at the example of the damped driven harmonic oscillator. Errors should be small! \n";
    std::cout << std::scientific;
    //![doxygen]
    //... in main
    //set start and end time
    double t_start = 0., t_end = 1.;
    //set physical parameters and initial condition
    const double damping = 0.2, omega_0 = 1.0, omega_drive = 0.9;
    std::array<double,2> u_start = solution(t_start, damping, omega_0,
            omega_drive), u_end(u_start);
    //construct a functor with the right interface
    using namespace std::placeholders; //for _1, _2, _3
    auto functor = std::bind( rhs, _1, _2, _3, damping, omega_0, omega_drive);
    double dt= 0;
    //integration
    int counter = dg::integrateERK( "Dormand-Prince-7-4-5", functor, t_start,
            u_start, t_end, u_end, dt, dg::pid_control, dg::l2norm, 1e-6);
    //now compute error
    dg::blas1::axpby( 1., solution(t_end, damping, omega_0, omega_drive), -1.,
            u_end);
    std::cout << "With "<<counter<<"\t Dormand Prince steps norm of error is "
              << dg::l2norm( u_end)<<"\n";
    //![doxygen]
    std::cout << "Explicit Methods \n";
    std::vector<std::string> names{
        "Heun-Euler-2-1-2",
        "Bogacki-Shampine-4-2-3",
        "ARK-4-2-3 (explicit)",
        "Zonneveld-5-3-4",
        "ARK-6-3-4 (explicit)",
        "Sayfy-Aburub-6-3-4",
        "Cash-Karp-6-4-5",
        "Fehlberg-6-4-5",
        "Dormand-Prince-7-4-5",
        "ARK-8-4-5 (explicit)",
        "Verner-8-5-6",
        "Fehlberg-13-7-8",
        "Feagin-17-8-10"
    };
    for( auto name : names)
    {
        dt = 0;
        u_start = solution(t_start, damping, omega_0, omega_drive);
        counter = dg::integrateERK( name, functor, t_start, u_start, t_end,
                u_end, dt, dg::pid_control, dg::l2norm, 1e-6, 1e-10);

        std::array<double, 2> sol = solution(t_end, damping, omega_0, omega_drive);
        dg::blas1::axpby( 1.,sol  , -1., u_end);
        std::cout << "With "<<std::setw(6)<<counter<<" steps norm of error in "
                  <<std::setw(24)<<name<<"\t"<<dg::l2norm( u_end)<<"\n";
    }
    ///-------------------------------Implicit Methods----------------------//
    std::cout << "Implicit Methods \n";
    std::vector<std::string> implicit_names{
        "SDIRK-2-1-2",
        "Billington-3-3-2",
        "TRBDF2-3-3-2",
        "Kvaerno-4-2-3",
        "ARK-4-2-3 (implicit)",
        "Cash-5-2-4",
        "Cash-5-3-4",
        "SDIRK-5-3-4",
        "ARK-6-3-4 (implicit)",
        "Kvaerno-7-4-5",
        "ARK-8-4-5 (implicit)",
    };
    for( auto name : implicit_names)
    {
        dt = 0;
        u_start = solution(t_start, damping, omega_0, omega_drive);
        dg::Adaptive< dg::DIRKStep<
            std::array<double,2>, dg::FixedPointSolver<std::array<double,2> > > >
                pd( name, u_start, 100, 1e-14);
        counter = integrateAdaptive( pd, functor, t_start, u_start, t_end,
            u_end, dt, dg::pid_control, dg::l2norm, 1e-6, 1e-10);

        std::array<double, 2> sol = solution(t_end, damping, omega_0, omega_drive);
        dg::blas1::axpby( 1.,sol  , -1., u_end);
        std::cout << "With "<<std::setw(6)<<counter<<" steps norm of error in "
                  <<std::setw(24)<<name<<"\t"<<dg::l2norm( u_end)<<"\n";
    }
    ///---------------------------Test domain restriction-------------------//
    std::cout << "Test domain restriction \n";
    for( auto name : names)
    {
        double dt = 0;
        double t_start = 0;
        double t_end = 10;
        double u_start = 1.0, u_end;
        auto rhs = [](double t, double y, double& yp){
                yp = y;
        };
        unsigned counter = dg::integrateERK( name, rhs , t_start, u_start, t_end,
                u_end, dt, dg::pid_control, dg::l2norm, 1e-6, 1e-10,
                 dg::Grid1d( 0., 100., 1,1)  );
        double analytic = log( 100.);
        std::cout << "With "<<std::setw(6)<<counter<<" steps norm of error in "
                  <<std::setw(24)<<name<<"\t"<<fabs( t_end - analytic)<<"\n";
    }
    return 0;
}

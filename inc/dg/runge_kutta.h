#ifndef _DG_RK_
#define _DG_RK_

#include <cassert>
#include <array>

#include "backend/exceptions.h"
#include "tableau.h"
#include "blas1.h"
#include "implicit.h"

/*! @file
 * @brief Runge-Kutta explicit time-integrators
 */

namespace dg{

 /** @class hide_rhs
  * @tparam RHS The right hand side
        is a functor type with no return value (subroutine)
        of signature <tt> void operator()(value_type, const ContainerType&, ContainerType&)</tt>
        The first argument is the time, the second is the input vector, which the functor may \b not override, and the third is the output,
        i.e. y' = f(t, y) translates to f(t, y, y').
        The two ContainerType arguments never alias each other in calls to the functor.
  */
 /** @class hide_limiter
  * @tparam Limiter The filter or limiter class to use in connection with the time-stepper
        has a member function \c apply
        of signature <tt> void apply( const ContainerType&, ContainerType&)</tt>
        The first argument is the input vector, which the functor may \b not override, and the second is the output,
        i.e. y' = L( y) translates to L.apply( y, y').
        The two ContainerType arguments never alias each other in calls to the functor.
  */


/**
* @brief Embedded Runge Kutta explicit time-step with error estimate
* \f[
 \begin{align}
    k_i = f\left( t^n + c_i \Delta t, u^n + \Delta t \sum_{j=1}^{s-1} a_{ij} k_j\right) \\
    u^{n+1} = u^{n} + \Delta t\sum_{j=1}^s b_j k_j \\
    \tilde u^{n+1} = u^{n} + \Delta t\sum_{j=1}^s \tilde b_j k_j
 \end{align}
\f]

The method is defined by its (extended explicit) ButcherTableau, given by
the coefficients \c a, \c b and \c c,  and \c s is the number
of stages. The embedding is given by the coefficients \c bt (tilde b).

You can provide your own coefficients or use one of the embedded methods
in the following table:
@copydoc hide_explicit_butcher_tableaus

* @note The name of this class is in reminiscence of the ARKode library http://runge.math.smu.edu/arkode_dev/doc/guide/build/html/index.html
* @copydoc hide_ContainerType
* @ingroup time
*/
template< class ContainerType>
struct ERKStep
{
    using value_type = get_value_type<ContainerType>;//!< the value type of the time variable (float or double)
    using container_type = ContainerType; //!< the type of the vector class in use
    ///@copydoc RungeKutta::RungeKutta()
    ERKStep(){ m_k.resize(1); //this makes the copyable function work
    }
    ///@copydoc RungeKutta::construct()
    ERKStep( ConvertsToButcherTableau<value_type> tableau, const ContainerType& copyable): m_rk(tableau), m_k(m_rk.num_stages(), copyable)
        { }
    ///@copydoc RungeKutta::construct()
    void construct( ConvertsToButcherTableau<value_type> tableau, const ContainerType& copyable ){
        *this = ERKStep( tableau, copyable);
    }
    ///@copydoc RungeKutta::copyable()
    const ContainerType& copyable()const{ return m_k[0];}

    ///All subsequent calls to \c step method will ignore the first same as last property (useful if you want to implement an operator splitting)
    void ignore_fsal(){ m_ignore_fsal = true;}
    ///All subsequent calls to \c step method will enable the check for the first same as last property
    void enable_fsal(){ m_ignore_fsal = false;}

    ///@copydoc RungeKutta::step()
    ///@param delta Contains error estimate (u1 - tilde u1) on return (must have equal size as \c u0)
    template<class RHS>
    void step( RHS& rhs, value_type t0, const ContainerType& u0, value_type& t1, ContainerType& u1, value_type dt, ContainerType& delta);
    ///global order of the method given by the current Butcher Tableau
    unsigned order() const {
        return m_rk.order();
    }
    ///global order of the embedding given by the current Butcher Tableau
    unsigned embedded_order() const {
        return m_rk.embedded_order();
    }
    ///number of stages of the method given by the current Butcher Tableau
    unsigned num_stages() const{
        return m_rk.num_stages();
    }
  private:
    ButcherTableau<value_type> m_rk;
    std::vector<ContainerType> m_k;
    value_type m_t1 = 1e300;//remember the last timestep at which ERK is called
    bool m_ignore_fsal = false;
};

///@cond
template< class ContainerType>
template< class RHS>
void ERKStep<ContainerType>::step( RHS& f, value_type t0, const ContainerType& u0, value_type& t1, ContainerType& u1, value_type dt, ContainerType& delta)
{
    unsigned s = m_rk.num_stages();
    //0 stage: probe
    value_type tu = t0;
    if( t0 != m_t1 || m_ignore_fsal)
        f(t0, u0, m_k[0]); //freshly compute k_0
    //else take from last call
    //1 stage
    if( s>1) {
        tu = DG_FMA( m_rk.c(1),dt, t0);
        blas1::axpby( 1., u0, dt*m_rk.a(1,0), m_k[0], delta);
        f( tu, delta, m_k[1]);
    }
    //2 stage
    if( s>2) {
        tu = DG_FMA( m_rk.c(2),dt, t0);
        blas1::evaluate( delta, dg::equals(), PairSum(), 1., u0,
                            dt*m_rk.a(2,0),m_k[0],
                            dt*m_rk.a(2,1),m_k[1]);
        f( tu, delta, m_k[2]);
    }
    //3 stage
    if( s> 3){
        tu = DG_FMA( m_rk.c(3),dt, t0);
        blas1::evaluate( delta, dg::equals(), PairSum(), 1., u0,
                             dt*m_rk.a(3,0),m_k[0],
                             dt*m_rk.a(3,1),m_k[1],
                             dt*m_rk.a(3,2),m_k[2]);
        f( tu, delta, m_k[3]);
    }
    //4 stage
    if( s>4){
        tu = DG_FMA( m_rk.c(4),dt, t0);
        blas1::evaluate( delta, dg::equals(), PairSum(), 1.        , u0,
                             dt*m_rk.a(4,0),m_k[0],  dt*m_rk.a(4,1),m_k[1],
                             dt*m_rk.a(4,2),m_k[2],  dt*m_rk.a(4,3),m_k[3]);
        f( tu, delta, m_k[4]);
    }
    //5 stage
    if( s>5) {
        tu = DG_FMA( m_rk.c(5),dt, t0);
        blas1::evaluate( delta, dg::equals(), PairSum(), 1., u0,
                 dt*m_rk.a(5,0),m_k[0], dt*m_rk.a(5,1),m_k[1],
                 dt*m_rk.a(5,2),m_k[2], dt*m_rk.a(5,3),m_k[3],
                 dt*m_rk.a(5,4),m_k[4]);
        f( tu, delta, m_k[5]);
    }
    //6 stage
    if( s>6)
    {
        tu = DG_FMA( m_rk.c(6),dt, t0);
        blas1::evaluate( delta, dg::equals(), PairSum(), 1., u0,
                           dt*m_rk.a(6,0),m_k[0], dt*m_rk.a(6,1),m_k[1],
                           dt*m_rk.a(6,2),m_k[2], dt*m_rk.a(6,3),m_k[3],
                           dt*m_rk.a(6,4),m_k[4], dt*m_rk.a(6,5),m_k[5]);
        f( tu, delta, m_k[6]);
        for ( unsigned i=7; i<s; i++)
        {
            blas1::axpby( 1.,u0, dt*m_rk.a(i,0),m_k[0], delta); //l=0
            tu = DG_FMA( dt,m_rk.c(i),t0); //l=0
            for( unsigned l=1; l<i; l++)
                blas1::axpby( dt*m_rk.a(i,l), m_k[l],1., delta);
            f( tu, delta, m_k[i]);
        }
    }
    //Now add everything up to get solution and error estimate
    switch( s)
    {
        //the first is for Euler
        case 1:
                blas1::subroutine( dg::EmbeddedPairSum(),
                            u1, delta,
                            1., 0., u0,
                            dt*m_rk.b(0), dt*m_rk.d(0), m_k[0]); break;
        case 2: blas1::subroutine( dg::EmbeddedPairSum(),
                            u1, delta,
                            1., 0., u0,
                            dt*m_rk.b(0), dt*m_rk.d(0), m_k[0],
                            dt*m_rk.b(1), dt*m_rk.d(1), m_k[1]); break;
        case 3: blas1::subroutine( dg::EmbeddedPairSum(),
                            u1, delta,
                            1., 0., u0,
                            dt*m_rk.b(0), dt*m_rk.d(0), m_k[0],
                            dt*m_rk.b(1), dt*m_rk.d(1), m_k[1],
                            dt*m_rk.b(2), dt*m_rk.d(2), m_k[2]); break;
        case 4: blas1::subroutine( dg::EmbeddedPairSum(),
                            u1, delta,
                            1., 0., u0,
                            dt*m_rk.b(0), dt*m_rk.d(0), m_k[0],
                            dt*m_rk.b(1), dt*m_rk.d(1), m_k[1],
                            dt*m_rk.b(2), dt*m_rk.d(2), m_k[2],
                            dt*m_rk.b(3), dt*m_rk.d(3), m_k[3]); break;
        case 5: blas1::subroutine( dg::EmbeddedPairSum(),
                            u1, delta,
                            1., 0., u0,
                            dt*m_rk.b(0), dt*m_rk.d(0), m_k[0],
                            dt*m_rk.b(1), dt*m_rk.d(1), m_k[1],
                            dt*m_rk.b(2), dt*m_rk.d(2), m_k[2],
                            dt*m_rk.b(3), dt*m_rk.d(3), m_k[3],
                            dt*m_rk.b(4), dt*m_rk.d(4), m_k[4]); break;
        case 6: blas1::subroutine( dg::EmbeddedPairSum(),
                            u1, delta,
                            1., 0., u0,
                            dt*m_rk.b(0), dt*m_rk.d(0), m_k[0],
                            dt*m_rk.b(1), dt*m_rk.d(1), m_k[1],
                            dt*m_rk.b(2), dt*m_rk.d(2), m_k[2],
                            dt*m_rk.b(3), dt*m_rk.d(3), m_k[3],
                            dt*m_rk.b(4), dt*m_rk.d(4), m_k[4],
                            dt*m_rk.b(5), dt*m_rk.d(5), m_k[5]); break;
        default: blas1::subroutine( dg::EmbeddedPairSum(),
                            u1, delta,
                            1., 0., u0,
                            dt*m_rk.b(0), dt*m_rk.d(0), m_k[0],
                            dt*m_rk.b(1), dt*m_rk.d(1), m_k[1],
                            dt*m_rk.b(2), dt*m_rk.d(2), m_k[2],
                            dt*m_rk.b(3), dt*m_rk.d(3), m_k[3],
                            dt*m_rk.b(4), dt*m_rk.d(4), m_k[4],
                            dt*m_rk.b(5), dt*m_rk.d(5), m_k[5],
                            dt*m_rk.b(6), dt*m_rk.d(6), m_k[6]);
            //sum the rest
            for( unsigned i=7; i<s; i++)
            {
                dg::blas1::axpby( dt*m_rk.b(i), m_k[i], 1., u1);
                dg::blas1::axpby( dt*m_rk.d(i), m_k[i], 1., delta);
            }
    }
    //make sure (t1,u1) is the last call to f
    m_t1 = t1 = t0 + dt;
    if(!m_rk.isFsal() )
        f(t1,u1,m_k[0]);
    else
    {
        using std::swap;
        swap( m_k[0], m_k[s-1]); //enable free swap functions
    }
}
///@endcond


//MW: if ever we want to change the SolverType at runtime (with an input parameter e.g.) make it a new parameter in the solve method (either abstract type or template like RHS)


/*!
 * @brief Additive Runge Kutta (semi-implicit) time-step with error estimate
 * following
 * <a href="http://runge.math.smu.edu/arkode_dev/doc/guide/build/html/Mathematics.html#arkstep-additive-runge-kutta-methods">The ARKode library</a>
 *
 * Currently, the possible Butcher Tableaus for a fully implicit-explicit scheme
 * are the "ARK-4-2-3", "ARK-6-3-4" and "ARK-8-4-5" combinations.
 * So far we did not implement the use of a mass matrix \c M.
 * @attention When you use the ARKStep in combination with the Adaptive time
 * step algorithm pay attention to solve the implicit part with sufficient
 * accuracy. Else, the error propagates into the time controller, which will
 * then choose the timestep as if the implicit part was explicit i.e. far too
 * small (don't really know fully why though).
 *
 * @copydoc hide_SolverType
 * @copydoc hide_ContainerType
 * @ingroup time
 */
template<class ContainerType, class SolverType = dg::DefaultSolver<ContainerType>>
struct ARKStep
{
    using value_type = get_value_type<ContainerType>;//!< the value type of the time variable (float or double)
    using container_type = ContainerType; //!< the type of the vector class in use
    ///@copydoc RungeKutta::RungeKutta()
    ARKStep(){ }
    /*!@brief Construct with given name
     * @param name Currently, one of "ARK-4-2-3", "ARK-6-3-4" or "ARK-8-4-5"
     * @param ps Parameters that
     * are forwarded to the constructor of \c SolverType
     * @tparam SolverParams Type of parameters (deduced by the compiler)
     */
    template<class ...SolverParams>
    ARKStep( std::string name, SolverParams&& ...ps) :
         m_solver( std::forward<SolverParams>(ps)...),
         m_rhs( m_solver.copyable())
    {
        if( name == "ARK-4-2-3" )
        {
            m_rkE = ConvertsToButcherTableau<value_type>( "ARK-4-2-3 (explicit)");
            m_rkI = ConvertsToButcherTableau<value_type>( "ARK-4-2-3 (implicit)");
        }
        else if( name == "ARK-6-3-4" )
        {
            m_rkE = ConvertsToButcherTableau<value_type>( "ARK-6-3-4 (explicit)");
            m_rkI = ConvertsToButcherTableau<value_type>( "ARK-6-3-4 (implicit)");
        }
        else if( name == "ARK-8-4-5" )
        {
            m_rkE = ConvertsToButcherTableau<value_type>( "ARK-8-4-5 (explicit)");
            m_rkI = ConvertsToButcherTableau<value_type>( "ARK-8-4-5 (implicit)");
        }
        else
            throw dg::Error( dg::Message()<<"Unknown name");
        assert( m_rkE.num_stages() == m_rkI.num_stages());
        m_kE.assign(m_rkE.num_stages(), m_rhs);
        m_kI.assign(m_rkI.num_stages(), m_rhs);
    }
    ///@copydoc construct()
    template<class ...SolverParams>
    ARKStep( ConvertsToButcherTableau<value_type> ex_tableau,
             ConvertsToButcherTableau<value_type> im_tableau,
             SolverParams&& ...ps
             ):
         m_solver( std::forward<SolverParams>(ps)...),
         m_rhs( m_solver.copyable()),
         m_rkE(ex_tableau),
         m_rkI(im_tableau),
         m_kE(m_rkE.num_stages(), m_rhs),
         m_kI(m_rkI.num_stages(), m_rhs)
    {
        assert( m_rkE.num_stages() == m_rkI.num_stages());
    }
    /*!@brief Construct with two Butcher Tableaus
     *
     * The two Butcher Tableaus represent the parameters for the explicit
     * and implicit parts respectively. If both the explicit and implicit part
     * of your equations are nontrivial, they must be one of the "ARK-X-X-X" methods
     * listed in \c ConvertsToButcherTableau. Or you have your own tableaus of
     * course but both tableaus must have the same number of steps.
     *
     * @param ex_tableau Tableau for the explicit part
     * @param im_tableau Tableau for the implicit part (must have the same number of stages as \c ex_tableau )
     * @param ps Parameters that
     * are forwarded to the constructor of \c SolverType
     * @tparam SolverParams Type of parameters (deduced by the compiler)
     */
    template<class ...SolverParams>
    void construct(
             ConvertsToButcherTableau<value_type> ex_tableau,
             ConvertsToButcherTableau<value_type> im_tableau,
             SolverParams&& ...ps
             )
    {
        *this = ARKStep( ex_tableau, im_tableau, std::forward<SolverParams>(ps)...);
    }
    ///@brief Return an object of same size as the object used for construction
    ///@return A copyable object; what it contains is undefined, its size is important
    const ContainerType& copyable()const{ return m_rhs;}

    ///Write access to the internal solver for the implicit part
    SolverType& solver() { return m_solver;}
    ///Read access to the internal solver for the implicit part
    const SolverType& solver() const { return m_solver;}

    /**
    * @brief Advance one step
    *
    * @copydoc hide_explicit_implicit
    * @param t0 start time
    * @param u0 value at \c t0
    * @param t1 (write only) end time ( equals \c t0+dt on return, may alias \c t0)
    * @param u1 (write only) contains result on return (may alias u0)
    * @param dt timestep
    * @param delta Contains error estimate (u1 - tilde u1) on return (must have equal size as \c u0)
    * @note the implementation is such that on return the last call is the
    * explicit part \c ex at the new \c (t1,u1).
    * This is useful if \c ex holds
    * state, which is then updated to the new timestep and/or if \c im changes
    * the state of \c ex through the friend construct.
    */
    template< class Explicit, class Implicit>
    void step( Explicit& ex, Implicit& im, value_type t0, const ContainerType& u0, value_type& t1, ContainerType& u1, value_type dt, ContainerType& delta);
    ///@copydoc ERKStep::order()
    unsigned order() const {
        return m_rkE.order();
    }
    ///@copydoc ERKStep::embedded_order()
    unsigned embedded_order() const {
        return m_rkE.order();
    }
    ///@copydoc ERKStep::num_stages()
    unsigned num_stages() const{
        return m_rkE.num_stages();
    }
    private:
    SolverType m_solver;
    ContainerType m_rhs;
    ButcherTableau<value_type> m_rkE, m_rkI;
    std::vector<ContainerType> m_kE, m_kI;
    value_type m_t1 = 1e300;
};

///@cond
template<class ContainerType, class SolverType>
template< class Explicit, class Implicit>
void ARKStep<ContainerType, SolverType>::step( Explicit& ex, Implicit& im, value_type t0, const ContainerType& u0, value_type& t1, ContainerType& u1, value_type dt, ContainerType& delta)
{
    unsigned s = m_rkE.num_stages();
    value_type tu = t0;
    //0 stage
    //a^E_00 = a^I_00 = 0
    if( t0 != m_t1)
        ex(t0, u0, m_kE[0]); //freshly compute k_0
    im(t0, u0, m_kI[0]);

    //1 stage
    blas1::evaluate( m_rhs, dg::equals(), PairSum(), 1., u0,
            dt*m_rkE.a(1,0), m_kE[0],
            dt*m_rkI.a(1,0), m_kI[0]);
    tu = DG_FMA( m_rkI.c(1),dt, t0);
    //store solution in delta, init with last solution
    //blas1::copy( u0, delta); //if Implicit is zero this leads to unnecessary iterations
    blas1::copy( m_rhs, delta); //better init with rhs
    m_solver.solve( -dt*m_rkI.a(1,1), im, tu, delta, m_rhs);
    ex(tu, delta, m_kE[1]);
    im(tu, delta, m_kI[1]);

    //2 stage
    blas1::evaluate( m_rhs, dg::equals(), PairSum(), 1., u0,
             dt*m_rkE.a(2,0), m_kE[0],
             dt*m_rkE.a(2,1), m_kE[1],
             dt*m_rkI.a(2,0), m_kI[0],
             dt*m_rkI.a(2,1), m_kI[1]);
    tu = DG_FMA( m_rkI.c(2),dt, t0);
    //just take last solution as init
    blas1::copy( m_rhs, delta); //better init with rhs
    m_solver.solve( -dt*m_rkI.a(2,2), im, tu, delta, m_rhs);
    ex(tu, delta, m_kE[2]);
    im(tu, delta, m_kI[2]);
    //3 stage
    blas1::evaluate( m_rhs, dg::equals(), PairSum(), 1., u0,
             dt*m_rkE.a(3,0), m_kE[0],
             dt*m_rkE.a(3,1), m_kE[1],
             dt*m_rkE.a(3,2), m_kE[2],
             dt*m_rkI.a(3,0), m_kI[0],
             dt*m_rkI.a(3,1), m_kI[1],
             dt*m_rkI.a(3,2), m_kI[2]);
    tu = DG_FMA( m_rkI.c(3),dt, t0);
    blas1::copy( m_rhs, delta); //better init with rhs
    m_solver.solve( -dt*m_rkI.a(3,3), im, tu, delta, m_rhs);
    ex(tu, delta, m_kE[3]);
    im(tu, delta, m_kI[3]);
    //higher stages
    for( unsigned i=4; i<s; i++)
    {
        dg::blas1::copy( u0, m_rhs);
        for( unsigned j=0; j<s; j++)
            dg::blas1::axpbypgz( dt*m_rkE.a(i,j), m_kE[j],
                                 dt*m_rkI.a(i,j), m_kI[j], 1., m_rhs);
        tu = DG_FMA( m_rkI.c(i),dt, t0);
        blas1::copy( m_rhs, delta); //better init with rhs
        m_solver.solve( -dt*m_rkI.a(i,i), im, tu, delta, m_rhs);
        ex(tu, delta, m_kE[i]);
        im(tu, delta, m_kI[i]);
    }
    m_t1 = t1 = tu;
    // do up to 8 stages for ARK-8-4-5
    //Now compute result and error estimate
    blas1::subroutine( dg::EmbeddedPairSum(),
            u1, delta,
             1., 0., u0,
            dt*m_rkE.b(0), dt*m_rkE.d(0),m_kE[0],
            dt*m_rkE.b(1), dt*m_rkE.d(1),m_kE[1],
            dt*m_rkE.b(2), dt*m_rkE.d(2),m_kE[2],
            dt*m_rkE.b(3), dt*m_rkE.d(3),m_kE[3],
            dt*m_rkI.b(0), dt*m_rkI.d(0),m_kI[0],
            dt*m_rkI.b(1), dt*m_rkI.d(1),m_kI[1],
            dt*m_rkI.b(2), dt*m_rkI.d(2),m_kI[2],
            dt*m_rkI.b(3), dt*m_rkI.d(3),m_kI[3]);
    //sum the rest
    for( unsigned i=4; i<s; i++)
    {
        dg::blas1::axpbypgz( dt*m_rkE.b(i), m_kE[i],
                             dt*m_rkI.b(i), m_kI[i], 1., u1);
        dg::blas1::axpbypgz( dt*m_rkE.d(i), m_kE[i],
                             dt*m_rkI.d(i), m_kI[i], 1., delta);
    }
    //make sure (t1,u1) is the last call to ex
    ex(t1,u1,m_kE[0]);
}
///@endcond

/**
* @brief Runge-Kutta fixed-step explicit time-integration
* \f[
 \begin{align}
    k_i = f\left( t^n + c_i \Delta t, u^n + \Delta t \sum_{j=1}^{s-1} a_{ij} k_j\right) \\
    u^{n+1} = u^{n} + \Delta t\sum_{j=1}^s b_j k_j
 \end{align}
\f]

The method is defined by its (explicit) ButcherTableau, given by
the coefficients \c a, \c b and \c c,  and \c s is the number
of stages.

You can provide your own coefficients or use one of our predefined methods:
@copydoc hide_explicit_butcher_tableaus
The following code snippet demonstrates how to use the class for the integration of
the harmonic oscillator:

@snippet runge_kutta_t.cu function
@snippet runge_kutta_t.cu doxygen
* @ingroup time
*
* @note Uses only \c dg::blas1 routines to integrate one step.
* @copydoc hide_ContainerType
*/
template<class ContainerType>
struct RungeKutta
{
    using value_type = get_value_type<ContainerType>;//!< the value type of the time variable (float or double)
    using container_type = ContainerType; //!< the type of the vector class in use
    ///@brief No memory allocation, Call \c construct before using the object
    RungeKutta(){}
    ///@copydoc construct()
    RungeKutta( ConvertsToButcherTableau<value_type> tableau, const ContainerType& copyable): m_erk( tableau, copyable), m_delta( copyable)
        { }
    /**
    * @brief Reserve internal workspace for the integration
    *
    * @param tableau Tableau, name or identifier that \c ConvertsToButcherTableau
    * @param copyable vector of the size that is later used in \c step (
     it does not matter what values \c copyable contains, but its size is important;
     the \c step method can only be called with vectors of the same size)
    */
    void construct(ConvertsToButcherTableau<value_type> tableau, const ContainerType& copyable){
        *this = RungeKutta( tableau, copyable);
    }
    ///@brief Return an object of same size as the object used for construction
    ///@return A copyable object; what it contains is undefined, its size is important
    const ContainerType& copyable()const{ return m_delta;}
    /**
    * @brief Advance one step
    *
    * @copydoc hide_rhs
    * @param rhs right hand side subroutine
    * @param t0 start time
    * @param u0 value at \c t0
    * @param t1 (write only) end time ( equals \c t0+dt on return, may alias \c t0)
    * @param u1 (write only) contains result on return (may alias u0)
    * @param dt timestep
    * @note on return \c rhs(t1, u1) will be the last call to \c rhs (this is useful if \c RHS holds state, which is then updated to the current timestep)
    * @note About the first same as last property (fsal): Some Butcher tableaus
    * (e.g. Dormand-Prince or Bogacki-Shampine) have the property that the last value k_s of a
    * timestep is the same as the first value k_0 of the next timestep. This
    * means that we can save one call to the right hand side. This property is
    * automatically activated if \c tableau.isFsal() returns \c true and \c t0
    * equals \c t1 of the last call to \c step. You can deactivate it by
    * calling the \c ignore_fsal() method, which is useful for splitting methods
    * but increases the number of rhs calls by 1.
    */
    template<class RHS>
    void step( RHS& rhs, value_type t0, const ContainerType& u0, value_type& t1, ContainerType& u1, value_type dt){
        m_erk.step( rhs, t0, u0, t1, u1, dt, m_delta);
    }
    ///All subsequent calls to \c step method will ignore the first same as last property (useful if you want to implement an operator splitting)
    void ignore_fsal(){ m_erk.ignore_fsal();}
    ///All subsequent calls to \c step method will enable the check for the first same as last property
    void enable_fsal(){ m_erk.enable_fsal();}
    ///@copydoc ERKStep::order
    unsigned order() const {
        return m_erk.order();
    }
    ///@copydoc ERKStep::num_stages()
    unsigned num_stages() const{
        return m_erk.num_stages();
    }
  private:
    ERKStep<ContainerType> m_erk;
    ContainerType m_delta;
};

/**
 * @brief A filter that does nothing
 */
struct IdentityFilter
{
    /**
     * @brief copy in to out
     *
     * @copydoc hide_ContainerType
     * @param in (input)
     * @param out (copied version of in)
     */
    template<class ContainerType0, class ContainerType1>
    void apply( const ContainerType0& in, ContainerType1& out) const{
        dg::blas1::copy( in, out);
    }

};
/**
* @brief Shu-Osher fixed-step explicit time-integration with Slope Limiter / Filter
* \f[
 \begin{align}
    u_0 &= u_n \\
    u_i &= \Lambda\Pi \left(\sum_{j=0}^{i-1}\left[ \alpha_{ij} u_j + \Delta t \beta_{ij} f( t_j, u_j)\right]\right)\\
    u^{n+1} &= u_s
 \end{align}
\f]

where \f$ \Lambda\Pi\f$ is the limiter, \c i=1,...,s and \c s is the number of stages (i.e. the number of times the right hand side is evaluated.

The method is defined by its (explicit) ShuOsherTableau, given by
the coefficients \c alpha and \c beta,  and \c s is the number
of stages.
@note the original reference for the scheme is
 * <a href="https://doi.org/10.1016/0021-9991(88)90177-5">
 Chi-Wang Shu, Stanley Osher,
Efficient implementation of essentially non-oscillatory shock-capturing schemes,
Journal of Computational Physics,
Volume 77, Issue 2,
1988,
Pages 439-471</a>
@note This struct can be used to implement the RKDG methods with slope-limiter described in
<a href ="https://doi.org/10.1023/A:1012873910884">Cockburn, B., Shu, CW. Runge???Kutta Discontinuous Galerkin Methods for Convection-Dominated Problems. Journal of Scientific Computing 16, 173???261 (2001) </a>

You can use one of our predefined methods (only the ones that are marked with "Shu-Osher-Form"):
@copydoc hide_explicit_butcher_tableaus
* @ingroup time
*
* @note Uses only \c dg::blas1 routines to integrate one step.
* @copydoc hide_ContainerType
*/
template<class ContainerType>
struct ShuOsher
{
    using value_type = get_value_type<ContainerType>;//!< the value type of the time variable (float or double)
    using container_type = ContainerType; //!< the type of the vector class in use
    ///@copydoc RungeKutta::RungeKutta()
    ShuOsher(){}
    ///@copydoc RungeKutta::construct()
    ShuOsher( dg::ConvertsToShuOsherTableau<value_type> tableau, const ContainerType& copyable): m_t( tableau), m_u(  m_t.num_stages(), copyable), m_k(m_u), m_temp(copyable)
        { }
    ///@copydoc RungeKutta::construct()
    void construct(dg::ConvertsToShuOsherTableau<value_type> tableau, const ContainerType& copyable){
        *this = ShuOsher( tableau, copyable);
    }
    ///@copydoc RungeKutta::copyable()
    const ContainerType& copyable()const{ return m_temp;}

    /**
    * @brief Advance one step
    *
    * @copydoc hide_rhs
    * @copydoc hide_limiter
    * @param limiter the filter or limiter to use
    * @param rhs right hand side subroutine
    * @param t0 start time
    * @param u0 value at \c t0
    * @param t1 (write only) end time ( equals \c t0+dt on return, may alias \c t0)
    * @param u1 (write only) contains result on return (may alias u0)
    * @param dt timestep
    * @note on return \c rhs(t1, u1) will be the last call to \c rhs (this is useful if \c RHS holds state, which is then updated to the current timestep)
    */
    template<class RHS, class Limiter>
    void step( RHS& rhs, Limiter& limiter, value_type t0, const ContainerType& u0, value_type& t1, ContainerType& u1, value_type dt){
        unsigned s = m_t.num_stages();
        std::vector<value_type> ts( m_t.num_stages()+1);
        ts[0] = t0;
        if( t0 != m_t1 ) //this is the first time we call step
        {
            limiter.apply( u0, m_u[0]);
            rhs(ts[0], m_u[0], m_k[0]); //freshly compute k_0
        }
        else
            dg::blas1::copy(u0, m_u[0]);
        for( unsigned i=1; i<=s; i++)
        {

            dg::blas1::axpbypgz( m_t.alpha(i-1,0), m_u[0], dt*m_t.beta(i-1,0), m_k[0], 0., m_temp);
            ts[i] = m_t.alpha(i-1,0)*ts[0] + dt*m_t.beta(i-1,0);
            for( unsigned j=1; j<i; j++)
            {
                //about the i-1: it is unclear to me how the ShuOsher tableau makes implicit schemes
                dg::blas1::axpbypgz( m_t.alpha(i-1,j), m_u[j], dt*m_t.beta(i-1,j), m_k[j], 1., m_temp);
                ts[i] += m_t.alpha(i-1,j)*ts[j] + dt*m_t.beta(i-1,j);

            }
            if(i!=s)
            {
                limiter.apply( m_temp, m_u[i]);
                rhs(ts[i], m_u[i], m_k[i]);
            }
            else{
                limiter.apply( m_temp, u1);
                //make sure (t1,u1) is the last call to f
                rhs(ts[i], u1, m_k[0]);
            }
        }
        m_t1 = t1 = ts[s];
    }
    ///@copydoc ERKStep::order
    unsigned order() const {
        return m_t.order();
    }
    ///@copydoc ERKStep::num_stages()
    unsigned num_stages() const{
        return m_t.num_stages();
    }
  private:
    ShuOsherTableau<value_type> m_t;
    std::vector<ContainerType> m_u, m_k;
    ContainerType m_temp;
    value_type m_t1 = 1e300;
};

/*!
 * @brief diagonally implicit Runge Kutta time-step with error estimate
* \f[
 \begin{align}
    k_i = f\left( t^n + c_i \Delta t, u^n + \Delta t \sum_{j=1}^{s} a_{ij} k_j\right) \\
    u^{n+1} = u^{n} + \Delta t\sum_{j=1}^s b_j k_j \\
    \tilde u^{n+1} = u^{n} + \Delta t\sum_{j=1}^s \tilde b_j k_j
 \end{align}
\f]
 *
 * So far we did not implement the use of a mass matrix \c M.
 * You can provide your own coefficients or use one of the methods
 * in the following table:
 * @copydoc hide_implicit_butcher_tableaus
 *
 * @copydoc hide_SolverType
 * @copydoc hide_ContainerType
 * @ingroup time
 */
template<class ContainerType, class SolverType = dg::DefaultSolver<ContainerType>>
struct DIRKStep
{
    using value_type = get_value_type<ContainerType>;//!< the value type of the time variable (float or double)
    using container_type = ContainerType; //!< the type of the vector class in use
    ///@copydoc RungeKutta::RungeKutta()
    DIRKStep(){ }

    ///@copydoc construct()
    template<class ...SolverParams>
    DIRKStep( ConvertsToButcherTableau<value_type> im_tableau,
             SolverParams&& ...ps
             ):
         m_solver( std::forward<SolverParams>(ps)...),
         m_rhs( m_solver.copyable()),
         m_rkI(im_tableau),
         m_kI(m_rkI.num_stages(), m_rhs)
    {
    }

    /*!@brief Construct with a diagonally implicit Butcher Tableau
     *
     * The tableau may be one of the implict methods listed in
     * \c ConvertsToButcherTableau, or you provide your own tableau.
     *
     * @param im_tableau diagonally implicit tableau, name or identifier that \c ConvertsToButcherTableau
     * @param ps Parameters that
     * are forwarded to the constructor of \c SolverType
     * @tparam SolverParams Type of parameters (deduced by the compiler)
     */
    template<class ...SolverParams>
    void construct(
             ConvertsToButcherTableau<value_type> im_tableau,
             SolverParams&& ...ps
             )
    {
        *this = DIRKStep( im_tableau, std::forward<SolverParams>(ps)...);
    }
    ///@brief Return an object of same size as the object used for construction
    ///@return A copyable object; what it contains is undefined, its size is important
    const ContainerType& copyable()const{ return m_rhs;}

    ///Write access to the internal solver for the implicit part
    SolverType& solver() { return m_solver;}
    ///Read access to the internal solver for the implicit part
    const SolverType& solver() const { return m_solver;}

    /**
    * @brief Advance one step
    *
    * @copydoc hide_rhs
    * @param rhs right hand side subroutine
    * @param t0 start time
    * @param u0 value at \c t0
    * @param t1 (write only) end time ( equals \c t0+dt on return
    *   may alias \c t0)
    * @param u1 (write only) contains result on return (may alias u0)
    * @param dt timestep
    * @param delta Contains error estimate (u1 - tilde u1) on return (must have equal size as \c u0)
    */
    template< class RHS>
    void step( RHS& rhs, value_type t0, const ContainerType& u0, value_type& t1, ContainerType& u1, value_type dt, ContainerType& delta);
    ///@copydoc ERKStep::order()
    unsigned order() const {
        return m_rkI.order();
    }
    ///@copydoc ERKStep::embedded_order()
    unsigned embedded_order() const {
        return m_rkI.order();
    }
    ///@copydoc ERKStep::num_stages()
    unsigned num_stages() const{
        return m_rkI.num_stages();
    }

    private:
    SolverType m_solver;
    ContainerType m_rhs;
    ButcherTableau<value_type> m_rkI;
    std::vector<ContainerType> m_kI;
};

///@cond
template<class ContainerType, class SolverType>
template< class RHS>
void DIRKStep<ContainerType, SolverType>::step( RHS& rhs, value_type t0, const ContainerType& u0, value_type& t1, ContainerType& u1, value_type dt, ContainerType& delta)
{
    unsigned s = m_rkI.num_stages();
    value_type tu = t0;
    //0 stage
    //rhs = u0
    tu = DG_FMA( m_rkI.c(0),dt, t0);
    blas1::copy( u0, delta); //better init with rhs
    if( !(m_rkI.a(0,0)==0) )
        m_solver.solve( -dt*m_rkI.a(0,0), rhs, tu, delta, u0);
    rhs(tu, delta, m_kI[0]);

    //1 stage
    if( s>1){
        blas1::evaluate( m_rhs, dg::equals(), PairSum(), 1., u0,
                dt*m_rkI.a(1,0), m_kI[0]);
        tu = DG_FMA( m_rkI.c(1),dt, t0);
        //store solution in delta, init with last solution
        blas1::copy( m_rhs, delta); //better init with rhs
        m_solver.solve( -dt*m_rkI.a(1,1), rhs, tu, delta, m_rhs);
        rhs(tu, delta, m_kI[1]);
    }
    //2 stage
    if( s>2){
        blas1::evaluate( m_rhs, dg::equals(), PairSum(), 1., u0,
                 dt*m_rkI.a(2,0), m_kI[0],
                 dt*m_rkI.a(2,1), m_kI[1]);
        tu = DG_FMA( m_rkI.c(2),dt, t0);
        //just take last solution as init
        blas1::copy( m_rhs, delta); //better init with rhs
        m_solver.solve( -dt*m_rkI.a(2,2), rhs, tu, delta, m_rhs);
        rhs(tu, delta, m_kI[2]);
    }
    //3 stage and higher
    if( s>3){
        blas1::evaluate( m_rhs, dg::equals(), PairSum(), 1., u0,
                 dt*m_rkI.a(3,0), m_kI[0],
                 dt*m_rkI.a(3,1), m_kI[1],
                 dt*m_rkI.a(3,2), m_kI[2]);
        tu = DG_FMA( m_rkI.c(3),dt, t0);
        blas1::copy( m_rhs, delta); //better init with rhs
        m_solver.solve( -dt*m_rkI.a(3,3), rhs, tu, delta, m_rhs);
        rhs(tu, delta, m_kI[3]);
        for( unsigned i=4; i<s; i++)
        {
            dg::blas1::copy( u0, m_rhs);
            for( unsigned j=0; j<i; j++)
                dg::blas1::axpby( dt*m_rkI.a(i,j), m_kI[j], 1., m_rhs);
            tu = DG_FMA( m_rkI.c(i),dt, t0);
            blas1::copy( m_rhs, delta); //better init with rhs
            m_solver.solve( -dt*m_rkI.a(i,i), rhs, tu, delta, m_rhs);
            rhs(tu, delta, m_kI[i]);
        }
    }
    t1 = t0 + dt;
    //Now compute result and error estimate
    switch( s)
    {
        case 1: blas1::subroutine( dg::EmbeddedPairSum(),
                            u1, delta,
                            1., 0., u0,
                            dt*m_rkI.b(0), dt*m_rkI.d(0), m_kI[0]); break;
        case 2: blas1::subroutine( dg::EmbeddedPairSum(),
                            u1, delta,
                            1., 0., u0,
                            dt*m_rkI.b(0), dt*m_rkI.d(0), m_kI[0],
                            dt*m_rkI.b(1), dt*m_rkI.d(1), m_kI[1]); break;
        case 3: blas1::subroutine( dg::EmbeddedPairSum(),
                            u1, delta,
                            1., 0., u0,
                            dt*m_rkI.b(0), dt*m_rkI.d(0), m_kI[0],
                            dt*m_rkI.b(1), dt*m_rkI.d(1), m_kI[1],
                            dt*m_rkI.b(2), dt*m_rkI.d(2), m_kI[2]); break;
        default: blas1::subroutine( dg::EmbeddedPairSum(),
                            u1, delta,
                            1., 0., u0,
                            dt*m_rkI.b(0), dt*m_rkI.d(0), m_kI[0],
                            dt*m_rkI.b(1), dt*m_rkI.d(1), m_kI[1],
                            dt*m_rkI.b(2), dt*m_rkI.d(2), m_kI[2],
                            dt*m_rkI.b(3), dt*m_rkI.d(3), m_kI[3]);
            //sum the rest
            for( unsigned i=4; i<s; i++)
            {
                dg::blas1::axpby( dt*m_rkI.b(i), m_kI[i], 1., u1);
                dg::blas1::axpby( dt*m_rkI.d(i), m_kI[i], 1., delta);
            }
    }
}
///@endcond
/**
* @brief Runge-Kutta fixed-step implicit time-integration
* \f[
 \begin{align}
    k_i = f\left( t^n + c_i \Delta t, u^n + \Delta t \sum_{j=1}^{s} a_{ij} k_j\right) \\
    u^{n+1} = u^{n} + \Delta t\sum_{j=1}^s b_j k_j
 \end{align}
\f]

The method is defined by its (implicit) ButcherTableau, given by
the coefficients \c a, \c b and \c c,  and \c s is the number
of stages.

You can provide your own coefficients or use one of our predefined methods:
@copydoc hide_implicit_butcher_tableaus
* @ingroup time
*
* @note Uses only \c dg::blas1 routines to integrate one step.
* @copydoc hide_ContainerType
* @copydoc hide_SolverType
*/
template<class ContainerType, class SolverType = dg::DefaultSolver<ContainerType>>
struct ImplicitRungeKutta
{
    using value_type = get_value_type<ContainerType>;//!< the value type of the time variable (float or double)
    using container_type = ContainerType; //!< the type of the vector class in use
    ///@brief No memory allocation, Call \c construct before using the object
    ImplicitRungeKutta(){}

    ///@copydoc DIRKStep::construct()
    template<class ...SolverParams>
    ImplicitRungeKutta( ConvertsToButcherTableau<value_type> im_tableau,
             SolverParams&& ...ps
             ): m_dirk( im_tableau, std::forward<SolverParams>(ps)...), m_delta(m_dirk.copyable())
             {}
    ///@copydoc DIRKStep::construct()
    template<class ...SolverParams>
    void construct( ConvertsToButcherTableau<value_type> im_tableau,
             SolverParams&& ...ps
             )
    {
        *this = ImplicitRungeKutta( im_tableau, std::forward<SolverParams>(ps)...);
    }
    ///@brief Return an object of same size as the object used for construction
    ///@return A copyable object; what it contains is undefined, its size is important
    const ContainerType& copyable()const{ return m_delta;}
    ///Write access to the internal solver for the implicit part
    SolverType& solver() { return m_dirk.solver();}
    ///Read access to the internal solver for the implicit part
    const SolverType& solver() const { return m_dirk.solver();}
    /**
    * @brief Advance one step
    *
    * @copydoc hide_rhs
    * @param rhs right hand side subroutine
    * @param t0 start time
    * @param u0 value at \c t0
    * @param t1 (write only) end time ( equals \c t0+dt on return
    *   may alias \c t0)
    * @param u1 (write only) contains result on return (may alias u0)
    * @param dt timestep
    */
    template<class RHS>
    void step( RHS& rhs, value_type t0, const ContainerType& u0, value_type& t1, ContainerType& u1, value_type dt){
        m_dirk.step( rhs, t0, u0, t1, u1, dt, m_delta);
    }
    ///@copydoc ERKStep::order
    unsigned order() const {
        return m_dirk.order();
    }
    ///@copydoc ERKStep::num_stages()
    unsigned num_stages() const{
        return m_dirk.num_stages();
    }
  private:
    DIRKStep<ContainerType, SolverType> m_dirk;
    ContainerType m_delta;
};

///@addtogroup time
///@{

/**
 * @brief Integrate differential equation with an explicit Runge-Kutta scheme and a fixed number of steps
 *
 * @copydoc hide_rhs
 * @copydoc hide_ContainerType
 * @param tableau Tableau, name or identifier that \c ConvertsToButcherTableau
 * @param rhs The right-hand-side
 * @param t_begin initial time
 * @param begin initial condition
 * @param t_end final time
 * @param end (write-only) contains solution at \c t_end on return (may alias begin)
 * @param N number of steps
 */
template< class RHS, class ContainerType>
void stepperRK(ConvertsToButcherTableau<get_value_type<ContainerType>> tableau,
        RHS& rhs, get_value_type<ContainerType>  t_begin, const ContainerType&
        begin, get_value_type<ContainerType> t_end, ContainerType& end,
        unsigned N )
{
    using value_type = get_value_type<ContainerType>;
    RungeKutta<ContainerType > rk( tableau, begin);
    if( t_end == t_begin){ end = begin; return;}
    const value_type dt = (t_end-t_begin)/(value_type)N;
    dg::blas1::copy( begin, end);
    value_type t0 = t_begin;
    for( unsigned i=0; i<N; i++)
        rk.step( rhs, t0, end, t0, end, dt);
}


///@}

} //namespace dg

#endif //_DG_RK_

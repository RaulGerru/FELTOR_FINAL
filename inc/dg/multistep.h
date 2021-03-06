#pragma once

#include <map>
#include "implicit.h"
#include "runge_kutta.h"
#include "multistep_tableau.h"

//MW: if ever we want to change the SolverType at runtime (with an input parameter e.g.) make it a new parameter in the solve method (either abstract type or template like RHS)

/*! @file
  @brief contains multistep explicit& implicit time-integrators
  */
namespace dg{


/*! @class hide_explicit_implicit
* @tparam Explicit The explicit part of the right hand side
    is a functor type with no return value (subroutine)
    of signature <tt> void operator()(value_type, const ContainerType&, ContainerType&)</tt>
    The first argument is the time, the second is the input vector, which the functor may \b not override, and the third is the output,
    i.e. y' = f(t, y) translates to f(t, y, y').
        The two ContainerType arguments never alias each other in calls to the functor.
 * @tparam Implicit The implicit part of the right hand side
        is a functor type with no return value (subroutine)
        of signature <tt> void operator()(value_type, const ContainerType&, ContainerType&)</tt>
        The first argument is the time, the second is the input vector, which the functor may \b not override, and the third is the output,
        i.e. y' = f(t, y) translates to f(t, y, y').
        The two ContainerType arguments never alias each other in calls to the functor.
    Furthermore, if the \c DefaultSolver is used, the routines %weights(), %inv_weights() and %precond() must be callable
    and return diagonal weights, inverse weights and the preconditioner for the conjugate gradient method.
    The return type of these member functions must be useable in blas2 functions together with the ContainerType type.
   @note If Explicit is a class then the suggested way of implementing Implicit is as a **friend** to the
   Explicit class. This is because the only reason to write the implicit part separate from the
   explicit part is the interface to the timestepper. The friend construct helps to reduce
   duplicate code and memory consumption by making Implicit essentially an extension of Explicit.
 * @param ex explic part
 * @param im implicit part ( must be linear in its second argument and symmetric up to weights)
 */
/*!@class hide_note_multistep
* @note Uses only \c blas1::axpby routines to integrate one step.
* @note The difference between a multistep and a single step method like RungeKutta
* is that the multistep only takes one right-hand-side evaluation per step.
* This is advantageous if the right hand side is expensive to evaluate.
* Even though Runge Kutta methods can have a larger absolute timestep, if
* the effective timestep per rhs evaluation is compared, multistep methods
* generally win.
* @note a disadvantage of multistep is that timestep adaption is not easily done.
*/

/**
 * @brief Semi-implicit multistep time-integration
 * \f[
 * \begin{align}
     v^{n+1} = \sum_{q=0}^{s-1} a_q v^{n-q} + \Delta t\left[\left(\sum_{q=0}^{s-1}b_q  \hat E(t^{n}-q\Delta t, v^{n-q}) + \sum_{q=1}^{s} c_q \hat I( t^n - q\Delta t, v^{n-q})\right) + c_0\hat I(t^{n}+\Delta t, v^{n+1})\right]
     \end{align}
     \f]
     which discretizes
     \f[
     \frac{\partial v}{\partial t} = \hat E(t,v) + \hat I(t,v)
     \f]
    where \f$ \hat E \f$ contains the explicit and \f$ \hat I \f$ the implicit part of the equations.
    As an example, the coefficients for the 3-step, 3rd order "Karniadakis" scheme are
    \f[
    a_0 = \frac{18}{11}\ a_1 = -\frac{9}{11}\  a_2 = \frac{2}{11} \\
    b_0 = \frac{18}{11}\ b_1 = -\frac{18}{11}\ b_2 = \frac{6}{11} \\
    c_0 = \frac{6}{11}\quad c_1 = c_2 = c_3 = 0 \\
\f]
    You can use your own coefficients defined as a \c dg::MultistepTableau
    or use one of the predefined coefficients in
    @copydoc hide_imex_multistep_tableaus
 *
 * The necessary Inversion in the implicit part is provided by the \c SolverType class.
 * Per Default, a conjugate gradient method is used (therefore \f$ \hat I(t,v)\f$ must be linear in \f$ v\f$).
 *
 * The following code example demonstrates how to implement the method of manufactured solutions on a 2d partial differential equation with the dg library:
 * @snippet multistep_t.cu function
 * In the main function:
 * @snippet multistep_t.cu karniadakis
 * @note In our experience the implicit treatment of diffusive or hyperdiffusive
terms may significantly reduce the required number of time steps. This
outweighs the increased computational cost of the additional matrix inversions.
However, each PDE is different and general statements like this one should be
treated with care.
* @copydoc hide_note_multistep
* @copydoc hide_SolverType
* @copydoc hide_ContainerType
* @ingroup time
*/
template<class ContainerType, class SolverType = dg::DefaultSolver<ContainerType>>
struct ImExMultistep
{
    using value_type = get_value_type<ContainerType>;//!< the value type of the time variable (float or double)
    using container_type = ContainerType; //!< the type of the vector class in use
    ///@copydoc RungeKutta::RungeKutta()
    ImExMultistep(){}

    /*! @brief Reserve memory for integration and construct Solver
     *
     * @param tableau Tableau, name or identifier that \c ConvertsToMultistepTableau
     * @param ps Parameters that are forwarded to the constructor of \c SolverType
     * @tparam SolverParams Type of parameters (deduced by the compiler)
     */
    template<class ...SolverParams>
    ImExMultistep( ConvertsToMultistepTableau<value_type> tableau,
            SolverParams&& ...ps):
        m_t(tableau),
        m_solver( std::forward<SolverParams>(ps)...)
    {
        //only store implicit part if needed
        unsigned size_f = 0;
        for( unsigned i=0; i<m_t.steps(); i++ )
        {
            if( m_t.im( i+1) != 0 )
                size_f = i+1;
        }
        m_im.assign( size_f, m_solver.copyable());

        m_u.assign( m_t.steps(), m_solver.copyable());
        m_ex.assign( m_t.steps(), m_solver.copyable());
        m_tmp = m_solver.copyable();
        m_counter = 0;
    }
    /**
    * @brief Perfect forward parameters to one of the constructors
    *
    * @tparam Params deduced by the compiler
    * @param ps parameters forwarded to constructors
    */
    template<class ...Params>
    void construct( Params&& ...ps)
    {
        //construct and swap
        *this = ImExMultistep( std::forward<Params>( ps)...);
    }
    ///@brief Return an object of same size as the object used for construction
    ///@return A copyable object; what it contains is undefined, its size is important
    const ContainerType& copyable()const{ return m_tmp;}

    ///Write access to the internal solver for the implicit part
    SolverType& solver() { return m_solver;}
    ///Read access to the internal solver for the implicit part
    const SolverType& solver() const { return m_solver;}

    /**
     * @brief Initialize timestepper. Call before using the step function.
     *
     * This routine has to be called before the first timestep is made.
     * @copydoc hide_explicit_implicit
     * @param t0 The intital time corresponding to u0
     * @param u0 The initial value of the integration
     * @param dt The timestep saved for later use
     * @note the implementation is such that on return the last call is the
     * explicit part \c ex at \c (t0,u0).  This is useful if \c ex holds
     * state, which is then updated to that timestep and/or if \c im changes
     * the state of \c ex through the friend construct.
     * This might be interesting if the call to \c ex changes its state.
     */
    template< class Explicit, class Implicit>
    void init( Explicit& ex, Implicit& im, value_type t0, const ContainerType& u0, value_type dt);

    /**
    * @brief Advance one timestep
    *
    * @copydoc hide_explicit_implicit
    * @param t (write-only), contains timestep corresponding to \c u on return
    * @param u (write-only), contains next step of time-integration on return
    * @note the implementation is such that on return the last call is the
    * explicit part \c ex at the new \c (t,u).  This is useful if \c ex holds
    * state, which is then updated to the new timestep and/or if \c im changes
    * the state of \c ex through the friend construct.
    * @attention The first few steps after the call to the init function are performed with a semi-implicit Runge-Kutta method to initialize the multistepper
    */
    template< class Explicit, class Implicit>
    void step( Explicit& ex, Implicit& im, value_type& t, ContainerType& u);

  private:
    dg::MultistepTableau<value_type> m_t;
    SolverType m_solver;
    std::vector<ContainerType> m_u, m_ex, m_im;
    ContainerType m_tmp;
    value_type m_tu, m_dt;
    unsigned m_counter; //counts how often step has been called after init
};

///@cond
template< class ContainerType, class SolverType>
template< class RHS, class Diffusion>
void ImExMultistep<ContainerType, SolverType>::init( RHS& f, Diffusion& diff, value_type t0, const ContainerType& u0, value_type dt)
{
    m_tu = t0, m_dt = dt;
    unsigned s = m_t.steps();
    blas1::copy(  u0, m_u[s-1]);
    m_counter = 0;
    if( s-1-m_counter < m_im.size())
        diff( m_tu, m_u[s-1-m_counter], m_im[s-1-m_counter]);
    f( t0, u0, m_ex[s-1]); //f may not destroy u0
}

template<class ContainerType, class SolverType>
template< class RHS, class Diffusion>
void ImExMultistep<ContainerType, SolverType>::step( RHS& f, Diffusion& diff, value_type& t, ContainerType& u)
{
    unsigned s = m_t.steps();
    if( m_counter < s - 1)
    {
        std::map<unsigned, std::string> order2method{
            {1, "ARK-4-2-3"},
            {2, "ARK-4-2-3"},
            {3, "ARK-4-2-3"},
            {4, "ARK-6-3-4"},
            {5, "ARK-8-4-5"},
            {6, "ARK-8-4-5"},
            {7, "ARK-8-4-5"}
        };
        ARKStep<ContainerType, SolverType> ark( order2method.at( m_t.order()), m_solver);
        ContainerType tmp ( u);
        ark.step( f, diff, t, u, t, u, m_dt, tmp);
        m_counter++;
        m_tu = t;
        dg::blas1::copy( u, m_u[s-1-m_counter]);
        //only assign to f if we actually need to store it
        if( s-1-m_counter < m_im.size())
            diff( m_tu, m_u[s-1-m_counter], m_im[s-1-m_counter]);
        f( m_tu, m_u[s-1-m_counter], m_ex[s-1-m_counter]);
        m_solver = ark.solver(); // store the state of the solver
        return;
    }
    //compute right hand side of inversion equation
    dg::blas1::axpbypgz( m_t.a(0), m_u[0], m_dt*m_t.ex(0), m_ex[0], 0., m_tmp);
    for (unsigned i = 1; i < s; i++)
        dg::blas1::axpbypgz( m_t.a(i), m_u[i], m_dt*m_t.ex(i), m_ex[i], 1., m_tmp);
    for (unsigned i = 0; i < m_im.size(); i++)
        dg::blas1::axpby( m_dt*m_t.im(i+1), m_im[i], 1., m_tmp);
    t = m_tu = m_tu + m_dt;

    value_type alpha[2] = {2., -1.};
    //value_type alpha[2] = {1., 0.};
    if( s > 1 ) //everything higher than Euler
        dg::blas1::axpby( alpha[0], m_u[0], alpha[1],  m_u[1], u);
    else
        dg::blas1::copy( m_u[0], u);

    //Rotate 1 to the right (note the reverse iterator here!)
    std::rotate( m_u.rbegin(), m_u.rbegin() + 1, m_u.rend());
    std::rotate( m_ex.rbegin(), m_ex.rbegin() + 1, m_ex.rend());
    if( !m_im.empty())
        std::rotate( m_im.rbegin(), m_im.rbegin() + 1, m_im.rend());
    //compute implicit part
    m_solver.solve( -m_dt*m_t.im(0), diff, t, u, m_tmp);

    blas1::copy( u, m_u[0]); //store result
    if( 0 < m_im.size())
        diff( m_tu, m_u[0], m_im[0]); //call diff on new point
    f(m_tu, m_u[0], m_ex[0]); //call f on new point (AFTER diff!)

}
///@endcond


/**
* @brief EXPERIMENTAL: Implicit multistep time-integration with Limiter/Filter
* \f[
* \begin{align}
    \tilde v &= \sum_{i=0}^{s-1} a_i v^{n-i} + \Delta t \sum_{i=1}^{s} c_i\hat I(t^{n+1-i}, v^{n+1-i}) + \Delta t c_{0} \hat I (t + \Delta t, \tilde v) \\
    v^{n+1} &= \Lambda\Pi\left(\tilde v\right)
    \end{align}
    \f]

    which discretizes
    \f[
    \frac{\partial v}{\partial t} = \hat I(t,v)
    \f]
    where \f$ \hat I \f$ represents the right hand side of the equations.
    You can use your own coefficients defined as a \c dg::MultistepTableau
    or use one of the predefined coefficients in
    @copydoc hide_implicit_multistep_tableaus
    and (any imex tableau can be used in an implicit scheme, disregarding the explicit coefficients)
    @copydoc hide_imex_multistep_tableaus
*
* The necessary Inversion in the implicit part is provided by the \c SolverType class.
* Per Default, a conjugate gradient method is used (therefore \f$ \hat I(t,v)\f$ must be linear in \f$ v\f$). For nonlinear right hand side we recommend the AndersonSolver
*
* @note In our experience the implicit treatment of diffusive or hyperdiffusive
terms can significantly reduce the required number of time steps. This
outweighs the increased computational cost of the additional inversions.
However, each PDE is different and general statements like this one should be
treated with care.
* @copydoc hide_note_multistep
* @copydoc hide_SolverType
* @copydoc hide_ContainerType
* @ingroup time
* @attention The filter function inside the Implicit Multistep method is a
* somewhat experimental feature, so use this class over
* \c dg::ImplicitMultistep at your own risk
*/
template<class ContainerType, class SolverType = dg::DefaultSolver<ContainerType>>
struct FilteredImplicitMultistep
{

    using value_type = get_value_type<ContainerType>;//!< the value type of the time variable (float or double)
    using container_type = ContainerType; //!< the type of the vector class in use
    ///@copydoc RungeKutta::RungeKutta()
    FilteredImplicitMultistep(){}

    /*! @brief Reserve memory for integration and construct Solver
     *
     * @param tableau Tableau, name or identifier that \c ConvertsToMultistepTableau
     * @param ps Parameters that are forwarded to the constructor of \c SolverType
     * @tparam SolverParams Type of parameters (deduced by the compiler)
     */
    template<class ...SolverParams>
    FilteredImplicitMultistep( ConvertsToMultistepTableau<value_type> tableau,
        SolverParams&& ...ps):
        m_t( tableau),
        m_solver( std::forward<SolverParams>(ps)...)
    {
        unsigned size_f = 0;
        for( unsigned i=0; i<m_t.steps(); i++ )
        {
            if( m_t.im( i+1) != 0 )
                size_f = i+1;
        }
        m_f.assign( size_f, m_solver.copyable());
        m_u.assign( m_t.steps(), m_solver.copyable());
        m_tmp = m_solver.copyable();
        m_counter = 0;
    }

    /**
    * @brief Perfect forward parameters to one of the constructors
    *
    * @tparam Params deduced by the compiler
    * @param ps parameters forwarded to constructors
    */
    template<class ...Params>
    void construct(Params&& ...ps)
    {
        //construct and swap
        *this = FilteredImplicitMultistep(  std::forward<Params>(ps)...);
    }
    ///@brief Return an object of same size as the object used for construction
    ///@return A copyable object; what it contains is undefined, its size is important
    const ContainerType& copyable()const{ return m_tmp;}
    ///Write access to the internal solver for the implicit part
    SolverType& solver() { return m_solver;}
    ///Read access to the internal solver for the implicit part
    const SolverType& solver() const { return m_solver;}

    ///@copydoc FilteredExplicitMultistep::init()
    template<class RHS, class Limiter>
    void init(RHS& rhs, Limiter& limiter, value_type t0, const ContainerType& u0, value_type dt);

    ///@copydoc FilteredExplicitMultistep::step()
    template<class RHS, class Limiter>
    void step(RHS& rhs, Limiter& limiter, value_type& t, container_type& u);
    private:
    dg::MultistepTableau<value_type> m_t;
    SolverType m_solver;
    value_type m_tu, m_dt;
    std::vector<ContainerType> m_u;
    std::vector<ContainerType> m_f;
    ContainerType m_tmp;
    unsigned m_counter = 0; //counts how often step has been called after init
};

///@cond
template< class ContainerType, class SolverType>
template<class RHS, class Limiter>
void FilteredImplicitMultistep<ContainerType, SolverType>::init(RHS& rhs, Limiter& l, value_type t0,
    const ContainerType& u0, value_type dt)
{
    m_tu = t0, m_dt = dt;
    l.apply( u0, m_u[m_u.size()-1]);
    m_counter = 0;
    //only assign to f if we actually need to store it
    unsigned s = m_t.steps();
    if( s-1-m_counter < m_f.size())
        rhs( m_tu, m_u[s-1-m_counter], m_f[s-1-m_counter]);
}

template< class ContainerType, class SolverType>
template<class RHS, class Limiter>
void FilteredImplicitMultistep<ContainerType, SolverType>::step(RHS& rhs, Limiter& l, value_type& t, container_type& u)
{
    unsigned s = m_t.steps();
        //max( m_u.size(), m_f.size())
    if( m_counter < s - 1)
    {
        std::map<unsigned, enum tableau_identifier> order2method{
            {1, IMPLICIT_EULER_1_1},
            {2, TRAPEZOIDAL_2_2},
            {3, KVAERNO_4_2_3},
            {4, SDIRK_5_3_4},
            {5, KVAERNO_7_4_5},
            {6, KVAERNO_7_4_5},
            {7, KVAERNO_7_4_5}
        };
        ImplicitRungeKutta<ContainerType, SolverType> dirk(
                order2method.at(m_t.order()), m_solver);
        dirk.step( rhs, t, u, t, u, m_dt);
        m_counter++;
        m_tu = t;
        l.apply( u, m_u[s-1-m_counter]);
        dg::blas1::copy(  m_u[s-1-m_counter], u);
        //only assign to f if we actually need to store it
        if( s-1-m_counter < m_f.size())
            rhs( m_tu, m_u[s-1-m_counter], m_f[s-1-m_counter]);
        m_solver = dirk.solver(); // store the state of the solver
        return;
    }
    //compute right hand side of inversion equation
    dg::blas1::axpby( m_t.a(0), m_u[0], 0., m_tmp);
    for (unsigned i = 1; i < s; i++)
        dg::blas1::axpby( m_t.a(i), m_u[i], 1., m_tmp);
    for (unsigned i = 0; i < m_f.size(); i++)
        dg::blas1::axpby( m_dt*m_t.im(i+1), m_f[i], 1., m_tmp);
    t = m_tu = m_tu + m_dt;

    value_type alpha[2] = {2., -1.};
    //value_type alpha[2] = {1., 0.};
    if( s > 1 ) //everything higher than Euler
        dg::blas1::axpby( alpha[0], m_u[0], alpha[1],  m_u[1], u);
    else
        dg::blas1::copy( m_u[0], u);

    //Rotate 1 to the right (note the reverse iterator here!)
    std::rotate(m_u.rbegin(), m_u.rbegin() + 1, m_u.rend());
    if( !m_f.empty())
        std::rotate(m_f.rbegin(), m_f.rbegin() + 1, m_f.rend());
    m_solver.solve( -m_dt*m_t.im(0), rhs, t, u, m_tmp);

    l.apply( u, m_u[0]);
    if( 0 < m_f.size())
        rhs( m_tu, m_u[0], m_f[0]);
    dg::blas1::copy(  m_u[0], u);
}
///@endcond

/**
* @brief Implicit multistep time-integration
* \f[
* \begin{align}
    v^{n+1} &= \sum_{i=0}^{s-1} a_i v^{n-i} + \Delta t \sum_{i=1}^{s} c_i\hat I(t^{n+1-i}, v^{n+1-i}) + \Delta t c_{0} \hat I (t + \Delta t, v^{n+1}) \\
    \end{align}
    \f]

    which discretizes
    \f[
    \frac{\partial v}{\partial t} = \hat I(t,v)
    \f]
    where \f$ \hat I \f$ represents the right hand side of the equations.
    You can use your own coefficients defined as a \c dg::MultistepTableau
    or use one of the predefined coefficients in
    @copydoc hide_implicit_multistep_tableaus
    and (any imex tableau can be used in an implicit scheme, disregarding the explicit coefficients)
    @copydoc hide_imex_multistep_tableaus
*
* The necessary Inversion in the implicit part is provided by the \c SolverType class.
* Per Default, a conjugate gradient method is used (therefore \f$ \hat I(t,v)\f$ must be linear in \f$ v\f$). For nonlinear right hand side we recommend the AndersonSolver
*
* @note In our experience the implicit treatment of diffusive or hyperdiffusive
terms can significantly reduce the required number of time steps. This
outweighs the increased computational cost of the additional inversions.
However, each PDE is different and general statements like this one should be
treated with care.
* @copydoc hide_note_multistep
* @copydoc hide_SolverType
* @copydoc hide_ContainerType
* @ingroup time
*/
template<class ContainerType, class SolverType = dg::DefaultSolver<ContainerType>>
struct ImplicitMultistep
{

    using value_type = get_value_type<ContainerType>;//!< the value type of the time variable (float or double)
    using container_type = ContainerType; //!< the type of the vector class in use
    ///@copydoc RungeKutta::RungeKutta()
    ImplicitMultistep(){}

    /*! @brief Reserve memory for integration and construct Solver
     *
     * @param tableau Tableau, name or identifier that \c ConvertsToMultistepTableau
     * @param ps Parameters that are forwarded to the constructor of \c SolverType
     * @tparam SolverParams Type of parameters (deduced by the compiler)
     */
    template<class ...SolverParams>
    ImplicitMultistep( ConvertsToMultistepTableau<value_type> tableau, SolverParams&& ...ps): m_bdf( tableau,
            std::forward<SolverParams>(ps)...) {}

    /**
    * @brief Perfect forward parameters to one of the constructors
    *
    * @tparam Params deduced by the compiler
    * @param ps parameters forwarded to constructors
    */
    template<class ...Params>
    void construct(Params&& ...ps)
    {
        //construct and swap
        *this = ImplicitMultistep(  std::forward<Params>(ps)...);
    }
    ///@copydoc RungeKutta::copyable()
    const ContainerType& copyable()const{ return m_bdf.copyable();}
    ///Write access to the internal solver for the implicit part
    SolverType& solver() { return m_bdf.solver();}
    ///Read access to the internal solver for the implicit part
    const SolverType& solver() const { return m_bdf.solver();}

    ///@copydoc ExplicitMultistep::init()
    template<class RHS>
    void init(RHS& rhs, value_type t0, const ContainerType& u0, value_type dt){
        dg::IdentityFilter id;
        m_bdf.init( rhs, id, t0, u0, dt);
    }

    ///@copydoc ExplicitMultistep::step()
    template<class RHS>
    void step(RHS& rhs, value_type& t, container_type& u){
        dg::IdentityFilter id;
        m_bdf.step( rhs, id, t, u);
    }
    private:
    FilteredImplicitMultistep<ContainerType, SolverType> m_bdf;
};



/**
* @brief EXPERIMENTAL: General explicit linear multistep time-integration with Limiter / Filter
* \f[
* \begin{align}
    \tilde v &= \sum_{j=0}^{s-1} a_j v^{n-j} + \Delta t\left(\sum_{j=0}^{s-1}b_j  \hat f\left(t^{n}-j\Delta t, v^{n-j}\right)\right) \\
    v^{n+1} &= \Lambda\Pi \left( \tilde v\right)
    \end{align}
    \f]

    where \f$ \Lambda\Pi\f$ is the limiter, which discretizes
    \f[
    \frac{\partial v}{\partial t} = \hat f(t,v)
    \f]
    where \f$ f \f$ contains the equations.
    The coefficients for an order 3 "eBDF" scheme are given as an example:
    \f[
    a_0 = \frac{18}{11}\ a_1 = -\frac{9}{11}\ a_2 = \frac{2}{11} \\
    b_0 = \frac{18}{11}\ b_1 = -\frac{18}{11}\ b_2 = \frac{6}{11}
\f]
    You can use your own coefficients defined as a \c dg::MultistepTableau
    or use one of the predefined coefficients in
    @copydoc hide_explicit_multistep_tableaus

@note This scheme is the same as ExplicitMultistep with the additional option to use a filter
*
* @copydoc hide_note_multistep
* @copydoc hide_ContainerType
* @ingroup time
* @attention The filter function inside the Explicit Multistep method is a
* somewhat experimental feature, so use this class over
* \c dg::ExplicitMultistep at your own risk
*/
template<class ContainerType>
struct FilteredExplicitMultistep
{
    using value_type = get_value_type<ContainerType>;//!< the value type of the time variable (float or double)
    using container_type = ContainerType; //!< the type of the vector class in use
    ///@copydoc RungeKutta::RungeKutta()
    FilteredExplicitMultistep(){ m_u.resize(1); //this makes the copyable function work
    }

    /**
     * @brief Reserve memory for the integration
     *
     * Set the coefficients \f$ a_i,\ b_i\f$
     * @param tableau Tableau, name or identifier that \c ConvertsToMultistepTableau
     * @param copyable ContainerType of the size that is used in \c step
     * @note it does not matter what values \c copyable contains, but its size is important
     */
    FilteredExplicitMultistep( ConvertsToMultistepTableau<value_type> tableau,
            const ContainerType& copyable): m_t(tableau)
    {
        m_f.assign( m_t.steps(), copyable);
        m_u.assign( m_t.steps(), copyable);
        m_counter = 0;
    }
    /**
    * @brief Perfect forward parameters to one of the constructors
    *
    * @tparam Params deduced by the compiler
    * @param ps parameters forwarded to constructors
    */
    template<class ...Params>
    void construct( Params&& ...ps)
    {
        //construct and swap
        *this = FilteredExplicitMultistep( std::forward<Params>( ps)...);
    }
    ///@brief Return an object of same size as the object used for construction
    ///@return A copyable object; what it contains is undefined, its size is important
    const ContainerType& copyable()const{ return m_u[0];}

    /**
     * @brief Initialize timestepper. Call before using the step function.
     *
     * This routine has to be called before the first timestep is made.
     * @copydoc hide_rhs
     * @copydoc hide_limiter
     * @param rhs The rhs functor
     * @param limiter The limiter or filter to use
     * @param t0 The intital time corresponding to u0
     * @param u0 The initial value of the integration
     * @param dt The timestep saved for later use
     * @note the implementation is such that on return the last call to the explicit part \c ex is at \c (t0,u0).
     * This might be interesting if the call to \c ex changes its state.
     */
    template< class RHS, class Limiter>
    void init( RHS& rhs, Limiter& limiter, value_type t0, const ContainerType& u0, value_type dt);

    /**
    * @brief Advance one timestep
    *
    * @copydoc hide_rhs
    * @copydoc hide_limiter
    * @param rhs The rhs functor
    * @param limiter The limiter or filter to use
    * @param t (write-only), contains timestep corresponding to \c u on return
    * @param u (write-only), contains next step of time-integration on return
    * @note the implementation is such that on return the last call to the explicit part \c ex is at the new \c (t,u).
    * This might be interesting if the call to \c ex changes its state.
    * @attention The first few steps after the call to the init function are performed with a Runge-Kutta method
    */
    template< class RHS, class Limiter>
    void step( RHS& rhs, Limiter& limiter, value_type& t, ContainerType& u);

  private:
    dg::MultistepTableau<value_type> m_t;
    std::vector<ContainerType> m_u, m_f;
    value_type m_tu, m_dt;
    unsigned m_counter; //counts how often step has been called after init
};
///@cond
template< class ContainerType>
template< class RHS, class Limiter>
void FilteredExplicitMultistep<ContainerType>::init( RHS& f, Limiter& l, value_type t0, const ContainerType& u0, value_type dt)
{
    m_tu = t0, m_dt = dt;
    unsigned s = m_t.steps();
    l.apply( u0, m_u[s-1]);
    f(m_tu, m_u[s-1], m_f[s-1]); //call f on new point
    m_counter = 0;
}

template<class ContainerType>
template<class RHS, class Limiter>
void FilteredExplicitMultistep<ContainerType>::step(RHS& f, Limiter& l, value_type& t, ContainerType& u)
{
    unsigned s = m_t.steps();
    if( m_counter < s-1)
    {
        std::map<unsigned, enum tableau_identifier> order2method{
            {1, SSPRK_2_2},
            {2, SSPRK_2_2},
            {3, SSPRK_3_3},
            {4, SSPRK_5_4},
            {5, SSPRK_5_4},
            {6, SSPRK_5_4},
            {7, SSPRK_5_4}
        };
        ShuOsher<ContainerType> rk( order2method.at(m_t.order()), u);
        rk.step( f, l, t, u, t, u, m_dt);
        m_counter++;
        m_tu = t;
        blas1::copy(  u, m_u[s-1-m_counter]);
        f( m_tu, m_u[s-1-m_counter], m_f[s-1-m_counter]);
        return;
    }
    //compute new t,u
    t = m_tu = m_tu + m_dt;
    dg::blas1::axpby( m_t.a(0), m_u[0], m_dt*m_t.ex(0), m_f[0], u);
    for (unsigned i = 1; i < s; i++){
        dg::blas1::axpbypgz( m_t.a(i), m_u[i], m_dt*m_t.ex(i), m_f[i], 1., u);
    }
    //permute m_f[s-1], m_u[s-1]  to be the new m_f[0], m_u[0]
    std::rotate( m_f.rbegin(), m_f.rbegin()+1, m_f.rend());
    std::rotate( m_u.rbegin(), m_u.rbegin()+1, m_u.rend());
    //apply limiter
    l.apply( u, m_u[0]);
    blas1::copy( m_u[0], u); //store result
    f(m_tu, m_u[0], m_f[0]); //call f on new point
}
///@endcond

/**
* @brief General explicit linear multistep time-integration
* \f[
* \begin{align}
    v^{n+1} = \sum_{j=0}^{s-1} a_j v^{n-j} + \Delta t\left(\sum_{j=0}^{s-1}b_j  \hat f\left(t^{n}-j\Delta t, v^{n-j}\right)\right)
    \end{align}
    \f]

    which discretizes
    \f[
    \frac{\partial v}{\partial t} = \hat f(t,v)
    \f]
    where \f$ f \f$ contains the equations.
    The coefficients for an order 3 "eBDF" scheme are given as an example:
    \f[
    a_0 = \frac{18}{11}\ a_1 = -\frac{9}{11}\ a_2 = \frac{2}{11} \\
    b_0 = \frac{18}{11}\ b_1 = -\frac{18}{11}\ b_2 = \frac{6}{11}
\f]
    You can use your own coefficients defined as a \c dg::MultistepTableau
    or use one of the predefined coefficients in
    @copydoc hide_explicit_multistep_tableaus
*
* @copydoc hide_note_multistep
* @copydoc hide_ContainerType
 @ingroup time
*/
template<class ContainerType>
struct ExplicitMultistep
{
    using value_type = get_value_type<ContainerType>;//!< the value type of the time variable (float or double)
    using container_type = ContainerType; //!< the type of the vector class in use
    ///@copydoc RungeKutta::RungeKutta()
    ExplicitMultistep(){}
    ///@copydoc FilteredExplicitMultistep::FilteredExplicitMultistep(ConvertsToMultistepTableau<value_type>,const ContainerType&)
    ExplicitMultistep( ConvertsToMultistepTableau<value_type> tableau, const ContainerType& copyable): m_fem( tableau, copyable){ }
    /**
    * @brief Perfect forward parameters to one of the constructors
    *
    * @tparam Params deduced by the compiler
    * @param ps parameters forwarded to constructors
    */
    template<class ...Params>
    void construct(Params&& ...ps)
    {
        //construct and swap
        *this = ExplicitMultistep(  std::forward<Params>(ps)...);
    }
    ///@brief Return an object of same size as the object used for construction
    ///@return A copyable object; what it contains is undefined, its size is important
    const ContainerType& copyable()const{ return m_fem.copyable();}

    /**
     * @brief Initialize timestepper. Call before using the step function.
     *
     * This routine has to be called before the first timestep is made.
     * @copydoc hide_rhs
     * @param rhs The rhs functor
     * @param t0 The intital time corresponding to u0
     * @param u0 The initial value of the integration
     * @param dt The timestep saved for later use
     * @note the implementation is such that on return the last call to the explicit part \c ex is at \c (t0,u0).
     * This might be interesting if the call to \c ex changes its state.
     */
    template< class RHS>
    void init( RHS& rhs, value_type t0, const ContainerType& u0, value_type dt){
        dg::IdentityFilter id;
        m_fem.init( rhs, id, t0, u0, dt);
    }

    /**
    * @brief Advance one timestep
    *
    * @copydoc hide_rhs
    * @param rhs The rhs functor
    * @param t (write-only), contains timestep corresponding to \c u on return
    * @param u (write-only), contains next step of time-integration on return
    * @note the implementation is such that on return the last call to the explicit part \c ex is at the new \c (t,u).
    * This might be interesting if the call to \c ex changes its state.
    * @attention The first few steps after the call to the init function are performed with a Runge-Kutta method (of the same order) to initialize the multistepper
    */
    template< class RHS>
    void step( RHS& rhs, value_type& t, ContainerType& u){
        dg::IdentityFilter id;
        m_fem.step( rhs, id, t, u);
    }

  private:
    FilteredExplicitMultistep<ContainerType> m_fem;
};

/** @brief DEPRECATED  (use ImExMultistep and select "Karniadakis" from the multistep tableaus)
* @ingroup time
* @sa dg::ImExMultistep
*/
template<class ContainerType, class SolverType = dg::DefaultSolver<ContainerType>>
struct Karniadakis
{
    using value_type = get_value_type<ContainerType>;
    using container_type = ContainerType;
    Karniadakis(){}
    template<class ...SolverParams>
    Karniadakis( SolverParams&& ...ps): m_imex( "Karniadakis", std::forward<SolverParams> (ps)...) { }
    template<class ...Params>
    void construct( Params&& ...ps)
    {
        *this = Karniadakis( std::forward<Params>( ps)...);
    }
    const ContainerType& copyable()const{ return m_imex.copyable();}
    SolverType& solver() { return m_imex.solver();}
    const SolverType& solver() const { return m_imex.solver();}
    template< class Explicit, class Implicit>
    void init( Explicit& ex, Implicit& im, value_type t0, const ContainerType& u0, value_type dt){
        m_imex.init( ex, im, t0, u0, dt);
    }
    template< class Explicit, class Implicit>
    void step( Explicit& ex, Implicit& im, value_type& t, ContainerType& u){
        m_imex.step( ex, im, t, u);
    }
  private:
    ImExMultistep<ContainerType, SolverType> m_imex;
};


} //namespace dg

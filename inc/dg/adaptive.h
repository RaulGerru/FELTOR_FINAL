#pragma once

#include "implicit.h"
#include "runge_kutta.h"

namespace dg
{

///@addtogroup time
///@{
/*! @brief Compute \f[ \sqrt{\sum_i x_i^2}\f] using \c dg::blas1::dot
 *
 * The intention of this function is to used in the \c Adaptive timestepping class.
 * @param x Vector to take the norm of
 * @return \c sqrt(dg::blas1::dot(x,x))
* @copydoc hide_ContainerType
 */
template <class ContainerType>
get_value_type<ContainerType> l2norm( const ContainerType& x)
{
    return sqrt( dg::blas1::dot( x,x));
}
///\f[ h'= h \epsilon_n^{-0.58/p}\epsilon_{n-1}^{0.21/p}\epsilon_{n-2}^{-0.1/p}\f]
template<class real_type>
real_type pid_control( real_type dt_old, real_type eps_0, real_type eps_1, real_type eps_2, unsigned embedded_order, unsigned order)
{
    real_type m_k1 = -0.58, m_k2 = 0.21, m_k3 = -0.1;
    real_type factor = pow( eps_0, m_k1/(real_type)order)
                     * pow( eps_1, m_k2/(real_type)order)
                     * pow( eps_2, m_k3/(real_type)order);
    //std::cout <<" control "<< eps_0<<" "<<eps_1<<" "<<eps_2<<" "<<factor<<"\n";
    return dt_old*factor;
}
///\f[ h'= h \epsilon_n^{-0.8/p}\epsilon_{n-1}^{0.31/p}\f]
template<class real_type>
real_type pi_control( real_type dt_old, real_type eps_0, real_type eps_1, real_type eps_2, unsigned embedded_order, unsigned order)
{
    real_type m_k1 = -0.8, m_k2 = 0.31;
    real_type factor = pow( eps_0, m_k1/(real_type)order)
                     * pow( eps_1, m_k2/(real_type)order);
    return dt_old*factor;
}
///\f[ h'= h \epsilon_n^{-1/p}\f]
template<class real_type>
real_type i_control( real_type dt_old, real_type eps_0, real_type eps_1, real_type eps_2, unsigned embedded_order, unsigned order)
{
    real_type m_k1 = -1.;
    real_type factor = pow( eps_0, m_k1/(real_type)order);
    return dt_old*factor;
}
///@}

///@cond
template<class real_type>
struct PIDController
{
    PIDController( ){}
    real_type operator()( real_type dt_old, real_type eps_n, real_type eps_n1, real_type eps_n2, unsigned embedded_order, unsigned order)const
    {
        real_type factor = pow( eps_n,  m_k1/(real_type)order)
                         * pow( eps_n1, m_k2/(real_type)order)
                         * pow( eps_n2, m_k3/(real_type)order);
        real_type dt_new = dt_old*std::max( m_lower_limit, std::min( m_upper_limit, factor) );
        return dt_new;
    }
    void set_lower_limit( real_type lower_limit) {
        m_lower_limit = lower_limit;
    }
    void set_upper_limit( real_type upper_limit) {
        m_upper_limit = upper_limit;
    }
    private:
    real_type m_k1 = -0.58, m_k2 = 0.21, m_k3 = -0.1;
    real_type m_lower_limit = 0, m_upper_limit = 1e300;
};
namespace detail{
template<class real_type>
struct Tolerance
{
    Tolerance( real_type rtol, real_type atol, real_type size):m_rtol(rtol*sqrt(size)), m_atol( atol*sqrt(size)){}
    DG_DEVICE
    void operator()( real_type previous, real_type& delta) const{
        delta = delta/ ( m_rtol*fabs(previous) + m_atol);
    }
    private:
    real_type m_rtol, m_atol;
};
} //namespace detail
///@endcond
/*!@class hide_stepper
 *
 * @tparam Stepper A timestepper class that computes the actual timestep
 * and an error estimate, for example an embedded Runge Kutta method
 * \c dg::ERKStep or the additive method \c dg::ARKStep. But really,
 * you can also choose to use your own timestepper class. The requirement
 * is that there is a \c step member function that is called as
 * \b stepper.step( rhs, t0, u0, t1, u1, dt, delta)
 * or \b stepper.step( ex, im, t0, u0, t1, u1, dt, delta)
 * depending on whether a purely explicit/implicit or a semi-implicit stepper
 * is used.
 * Here, t0, t1 and dt are of type \b Stepper::real_type, u0,u1 and delta
 * are vector types of type \b Stepper::container& and rhs, ex and im are
 * functors implementing the equations that are forwarded from the caller.
 * The parameters t1, u1 and delta are output parameters and must be updated by
 * the stepper.
 * The \c Stepper must have a default constructor and a constructor that takes
 * \c Stepper::container& as the first parameter.
 * Also, it must have the \c order() and \c embedded_order() member functions that
 * return the (global) order of the method and its error estimate.
 */

//%%%%%%%%%%%%%%%%%%%Adaptive%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
/*!@brief Driver class for adaptive timestep integration
 *
 * In order to build an adaptive Time integrator you basically need three
 * ingredients: a \c Stepper, a \c ControlFunction and an \c ErrorNorm.
 * The \c Stepper does the actual computation and advances the solution one step further
 * with a given timestep \c dt. Furthermore, it has to come up with an estimate
 * of the error of the solution and indicate the order of that error.
 * With the \c ErrorNorm the error estimate can be converted to a scalar that
 * can be compared to given relative and absolute error tolerances \c rtol and \c atol.
 * Based on the comparison the step is either accepted or rejected. In both cases
 * the \c ControlFunction then comes up with an adapted
 * suggestion for the timestep in the next step, however, if the step was
 * rejected, we make the stepsize decrease by at least 10\%.
 * For more information on these concepts we recommend
 * <a href="http://runge.math.smu.edu/arkode_dev/doc/guide/build/html/Mathematics.html#">the mathematical primer</a> of the ARKode library.
 *
 * For an example on how to use this class in a practical example consider the following code snippet:
 * @snippet multistep_t.cu adaptive
 * @copydoc hide_stepper
 * @note On step rejection, choosing timesteps and introducing restrictions on the controller: here is a quote from professor G. Söderlind (the master of control functions) from a private e-mail conversation:
 *
@note "The issue is that most controllers are best left to do their work without interference. Each time one interferes with the control loop, one creates a transient that has to be dealt with by the controller.

@note It is indeed necessary to reject steps every once in a while, but I usually try to avoid it to the greatest possible extent. In my opinion, the fear of having a too large error on some steps is vastly exaggerated. You see, in some steps the error is too large, and in others it is too small, and all controllers I constructed are designed to be “expectation value correct” in the sense that if the errors are random, the too large and too small errors basically cancel in the long run.

@note Even so, there are times when the error is way out of proportion. But I usually accept an error that is up to, say 2*TOL, which typically won’t cause any problems. Of course, if one hits a sharp change in the solution, the error may be much larger, and the step recomputed. But then one must “reset" the controller, i.e., the controller keeps back information, and it is better to restart the controller to avoid back information to affect the next step."
@attention Should you use this class instead of a fixed stepsize Multistep say?
The thing to consider especially when solving
partial differential equations, is that the right hand side might be very costly
to evaluate. An adaptive stepper (especially one-step methods) usually calls this right hand side more often than a multistep (only one call per step). The additional
computation of the error norm in the Adaptive step might also be important since
the norm requires global communication in a parallel program (bad for scaling
to many nodes).
So does the Adaptive timestepper make up for its increased cost throught the
adaption of the stepsize? In some cases the answer might be a sound Yes.
Especially when
there are velocity peaks in the solution the multistep timestep might be restricted
too much. In other cases when the timestep does not need to be adapted much, a multistep method can be faster.
In any case the big advantage of Adaptive is that it usually just works (even though
it is not fool-proof) and you do not have to spend time finding a suitable timestep
like in the multistep method.
@ingroup time
 */
template<class Stepper>
struct Adaptive
{
    using container = typename Stepper::container; //!< the type of the vector class in use by \c Stepper
    using real_type = typename Stepper::real_type; //!< the value type of the time variable defined by \c Stepper (float or double)
    ///@brief No memory allocation,
    ///requires \c Stepper to have default constructor
    Adaptive(){}
    /*!@brief Allocate workspace and construct stepper
     * @param copyable vector of the size that is later used in \c step (
      it does not matter what values \c copyable contains, but its size is important;
      the \c step method can only be called with vectors of the same size)
     * @param ps Parameters that, together with \c copyable as the first parameter,
     * are forwarded to the constructor of \c Stepper
     * @tparam StepperParams Type of parameters (deduced by the compiler)
     */
    template<class ...StepperParams>
    Adaptive( const container& copyable, StepperParams&& ...ps): m_stepper(copyable, std::forward<StepperParams>(ps)...) , m_next(copyable), m_delta(copyable)
    {
        dg::blas1::copy( 1., m_next);
        m_size = dg::blas1::dot( m_next, 1.);
    }
    /*!@brief Guess an initial stepsize
     *
     * If you have wondered what stepsize you should choose in the u0ning,
     * don't freak out about it. Really, the initial stepsize is not that
     * important, the stepper does not even have to succeed. Usually the
     * control function will very(!) quickly adapt the stepsize in just one or
     * two steps (even if it's several orders of magnitude off in the u0ning).
     *
     * Currently, this function won't do much better than if you just choose a
     * smallish number yourself, but it's there for future improvements.
     */
    template<class Explicit, class ErrorNorm = real_type(const container&)>
    real_type guess_stepsize( Explicit& ex, real_type t0, const container& u0, enum direction dir, ErrorNorm& norm, real_type rtol, real_type atol);

    /*!@brief Explicit or Implicit adaptive step
     *
     * @param rhs The right hand side of the equation to integrate
     * @copydoc hide_adaptive_params
     * @copydoc hide_rhs
     * @copydoc hide_control_error
     */
    template< class RHS,
              class ControlFunction = real_type (real_type, real_type, real_type, real_type, unsigned, unsigned),
              class ErrorNorm = real_type( const container&)>
    void step( RHS& rhs,
              real_type t0,
              const container& u0,
              real_type& t1,
              container& u1,
              real_type& dt,
              ControlFunction& control,
              ErrorNorm& norm,
              real_type rtol,
              real_type atol
              )
    {
        m_stepper.step( rhs, t0, u0, m_t_next, m_next, dt, m_delta);
        return update( t0, u0, t1, u1, dt, control, norm , rtol, atol);
    }
    /*!@brief Semi-implicit adaptive step
     *
     * @copydoc hide_adaptive_params
     * @copydoc hide_explicit_implicit
     * @copydoc hide_control_error
     */
    template< class Explicit,
              class Implicit,
              class ControlFunction = real_type (real_type, real_type, real_type, real_type, unsigned, unsigned),
              class ErrorNorm = real_type( const container&)>
    void step( Explicit& ex,
              Implicit& im,
              real_type t0,
              const container& u0,
              real_type& t1,
              container& u1,
              real_type& dt,
              ControlFunction& control,
              ErrorNorm& norm,
              real_type rtol,
              real_type atol)
    {
        m_stepper.step( ex, im, t0, u0, m_t_next, m_next, dt, m_delta);
        return update( t0, u0, t1, u1, dt, control, norm , rtol, atol);
    }
    ///Return true if the last stepsize in step was rejected
    bool hasFailed() const {
        return m_failed;
    }
    private:
    template<   class ControlFunction = real_type (real_type, real_type, real_type, real_type, unsigned, unsigned),
                class ErrorNorm = real_type( const container&)>
    void update( real_type t0,
                const container& u0,
                real_type& t1,
                container& u1,
                real_type& dt,
                ControlFunction& control,
                ErrorNorm& norm,
                real_type rtol,
                real_type atol
              )
    {
        //std::cout << "Try stepsize "<<dt;
        dg::blas1::subroutine( detail::Tolerance<real_type>( rtol, atol, m_size), u0, m_delta);
        real_type eps0 = norm(m_delta);
        //std::cout << " error "<<eps0;
        if( eps0 > m_reject_limit || std::isnan( eps0) )
        {
            real_type dt_old = dt;
            dt = control( dt, eps0, m_eps1, m_eps2, m_stepper.embedded_order(), m_stepper.order());
            if( fabs( dt) > 0.9*fabs(dt_old))
                dt = 0.9*dt_old;
            //0.9*dt_old is a safety limit
            //that prevents an increase of the timestep in case the stepper fails
            m_failed = true;
            dg::blas1::copy( u0, u1);
            t1 = t0;
            //std::cout << " Failed! New stepsize: "<<dt;
        }
        else
        {
            dt = control( dt, eps0, m_eps1, m_eps2, m_stepper.embedded_order(), m_stepper.order());
            m_eps2 = m_eps1;
            m_eps1 = eps0;
            dg::blas1::copy( m_next, u1);
            t1 = m_t_next;
            m_failed = false;
            //std::cout << " Success " << t1<<" "<<u1[0]<<" "<<u1[1];
            //std::cout << " New stepsize "<<dt;
        }
        //std::cout << std::endl;
    }
    bool m_failed = false;
    Stepper m_stepper;
    container m_next, m_delta;
    real_type m_reject_limit = 2;
    real_type m_size, m_eps1=1, m_eps2=1;
    real_type m_t_next = 0;
};
template<class Stepper>
template<class Explicit, class ErrorNorm>
typename Adaptive<Stepper>::real_type Adaptive<Stepper>::guess_stepsize( Explicit& ex, real_type t0, const container& u0, enum direction dir, ErrorNorm& tol, real_type rtol, real_type atol)
{
    real_type desired_accuracy = rtol*tol(u0) + atol;
    ex( t0, u0, m_next);
    real_type dt = pow(desired_accuracy, 1./(real_type)m_stepper.order())/tol(m_next);
    if( dir != forward)
        dt*=-1.;
    return dt;
}
/*!@class hide_adaptive_params
 * @param t0 initial time
 * @param u0 value at \c t0
 * @param t1 (write only) end time ( equals \c t0+dt on output if the step was accepted, otherwise equals \c t0, may alias \c t0)
 * @param u1 (write only) contains the updated result on output if the step was accepted, otherwise a copy of \c u0 (may alias \c u0)
 * @param dt on input: timestep to try out (see dg::Adaptive::guess_stepsize() for an initial stepsize).
 * On output: stepsize proposed by the controller that can be used to continue the integration in the next step.
 * @param control The control function. Usually \c dg::pid_control is a good choice
 * @param norm The error norm. Usually \c dg::l2norm is a good choice, but for
 * very small vector sizes the time for the binary reproducible dot product might become
 * a performance bottleneck. Then it's time for your own implementation.
 * @param rtol the desired relative accuracy. Usually 1e-5 is a good choice.
 * @param atol the desired absolute accuracy. Usually 1e-10 is a good choice.
 * @note Try not to mess with dt. The controller is best left alone and it does a very good job choosing timesteps. But how do I output my solution at certain (equidistant) timesteps? First, think about if you really, really need that. Why is it so bad to have
 * output at non-equidistant timesteps? If you are still firm, then consider
 * using an interpolation scheme (cf. \c dg::Extrapolation). Let choosing the timestep
 * yourself be the very last option if the others are not viable
 * @note For partial differential equations the exact value of \c rtol and \c atol might
 * not be important. Due to the CFL condition there might be a sharp barrier in the
 * range of possible stepsizes and the controller usually does a good job finding
 * it and keeping the timestep "just right". However, don't make \c rtol too small, \c 1e-1 say, since then the controller might
 * get too close to the CFL barrier. The timestepper is still able
 * to crash, mind, even though the chances of that happening are
 * somewhat lower than in a fixed stepsize method.
 */

/*!@class hide_control_error
 *
 * @tparam ControlFunction function or Functor called as dt' = control( dt, eps0, eps1, eps2, order, embedded_order), where all parameters are of type real_type except the last two, which are unsigned
 * @tparam ErrorNorm function or Functor of type real_type( const ContainerType&)
 */

///@addtogroup time
///@{

/**
 * @brief Integrates a differential equation using a one-step explicit Timestepper, with adaptive stepsize-control and monitoring the sanity of integration
 *
 * @param adaptive An instance of the Adaptive class
 * @param rhs The right-hand-side
 * @copydoc hide_adaptive_params
 * @return number of steps
 * @copydoc hide_rhs
 * @copydoc hide_ContainerType
 */
template< class Adaptive,
          class RHS,
          class ContainerType,
          class ErrorNorm = get_value_type<ContainerType>( const ContainerType&),
          class ControlFunction = get_value_type<ContainerType> (get_value_type<ContainerType>, get_value_type<ContainerType>, get_value_type<ContainerType>, get_value_type<ContainerType>, unsigned, unsigned)>
int integrateAdaptive(Adaptive& adaptive,
                      RHS& rhs,
                      get_value_type<ContainerType> t0,
                      const ContainerType& u0,
                      get_value_type<ContainerType> t1,
                      ContainerType& u1,
                      get_value_type<ContainerType> dt,
                      ControlFunction control,
                      ErrorNorm norm,
                      get_value_type<ContainerType> rtol,
                      get_value_type<ContainerType> atol=1e-10
                      )
{
    using  real_type = get_value_type<ContainerType>;
    real_type t_current = t0, dt_current = dt;
    blas1::copy( u0, u1 );
    ContainerType& current(u1);
    if( t1 == t0)
        return 0;
    bool forward = (t1 - t0 > 0);
    if( dt == 0)
        dt_current = adaptive.guess_stepsize( rhs, t0, u0, forward ? dg::forward:dg::backward, norm, rtol, atol);

    int counter =0;
    while( (forward && t_current < t1) || (!forward && t_current > t1))
    {
        dt = dt_current;
        if( (forward && t_current+dt_current > t1) || (!forward && t_current + dt_current < t1) )
            dt_current = t1-t_current;
        // Compute a step and error
        adaptive.step( rhs, t_current, current, t_current, current, dt_current, control, norm, rtol, atol);
        counter++;
    }
    return counter;
}

///@brief Shortcut for \c dg::integrateAdaptive with an embedded ERK class as timestepper
///@snippet adaptive_t.cu function
///@snippet adaptive_t.cu doxygen
///@param name name of an embedded method that \c ConvertsToButcherTableau
///@param rhs The right-hand-side
///@copydoc hide_adaptive_params
///@return number of steps
///@copydoc hide_rhs
///@copydoc hide_ContainerType
template<class RHS, class ContainerType, class ErrorNorm = get_value_type<ContainerType>( const ContainerType&),
             class ControlFunction = get_value_type<ContainerType> (get_value_type<ContainerType>, get_value_type<ContainerType>, get_value_type<ContainerType>, get_value_type<ContainerType>, unsigned, unsigned)>
int integrateERK( std::string name,
                  RHS& rhs,
                  get_value_type<ContainerType> t0,
                  const ContainerType& u0,
                  get_value_type<ContainerType> t1,
                  ContainerType& u1,
                  get_value_type<ContainerType> dt,
                  ControlFunction control,
                  ErrorNorm norm,
                  get_value_type<ContainerType> rtol,
                  get_value_type<ContainerType> atol=1e-10
              )
{
    dg::Adaptive<dg::ERKStep<ContainerType>> pd( u0, name);
    return integrateAdaptive( pd, rhs, t0, u0, t1, u1, dt, control, norm, rtol, atol);
}
///@}
}//namespace dg

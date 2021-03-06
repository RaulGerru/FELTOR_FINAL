#ifndef _DG_EVE_
#define _DG_EVE_

#include <cmath>
#include "blas.h"
#include "functors.h"

/*! @file
 * EVE adds an estimator for the largest Eigenvalue of the generalized Eigenvalue problem
 *  @author Eduard Reiter and Matthias Wiesenberger
 */


namespace dg
{

/*! @brief The Eigen-Value-Estimator (EVE) finds largest Eigenvalue of \f[ M^{-1}Ax = \lambda_\max x\f]
 *
 * Estimate largest Eigenvalue
 * of a symmetric positive definite matrix \f$ A\f$ with possible preconditioner \f$ M^{-1}\f$ using the conjugate gradient (CG) method.
 * The unpreconditioned version is the algorithm suggested in <a href="http://www.iam.fmph.uniba.sk/amuc/ojs/index.php/algoritmy/article/view/421">Tichy, On error estimation in the conjugate gradient method: Normwise backward error, Proceedings of the Conference Algoritmy, 323-332, 2016 </a>.
 * The reason this works also for the preconditioned CG method is because
 * preconditioned CG is equivalent to applying the
 * unpreconditioned CG method to \f$ \bar A\bar x = \bar b\f$ with \f$ \bar A
 * := {E^{-1}}^\mathrm{T} A E^{-1} \f$, \f$ \bar x := Ex\f$ and \f$ \bar b :=
 * {E^{-1}}^\mathrm{T}\f$, where \f$ M^{-1} = {E^{-1}}^\mathrm{T} E^{-1}\f$ is
 * the preconditioner.
 * The maximum Eigenvalue of \f$ M^{-1} A\f$ is the same as the
 * maximum EV of \f$ {E^{-1}}^\mathrm{T} A E^{-1} \f$
 *
* @attention beware the sign: a negative definite matrix does @b not work in Conjugate gradient
* @copydoc hide_ContainerType
 * @ingroup invert
*/
template< class ContainerType>
class EVE
{
  public:
    using container_type = ContainerType;
    using value_type = get_value_type<ContainerType>; //!< value type of the ContainerType class
    ///@brief Allocate nothing, Call \c construct method before usage
    EVE() {}
    ///@copydoc construct()
    EVE( const ContainerType& copyable, unsigned max_iter = 100):r( copyable), p( r), ap( r), m_max_iter( max_iter) {}
    /**
     * @brief Allocate memory for the pcg method
     *
     * @param copyable A ContainerType must be copy-constructible from this
     * @param max_iter Maximum number of iterations to be used
     */
    void construct( const ContainerType& copyable, unsigned max_iter = 100) {
        ap = p = r = copyable;
        m_max_iter = max_iter;
    }
    /// Set maximum number of iterations
    void set_max( unsigned new_max) {
        m_max_iter = new_max;
    }
    /// Get maximum number of iterations
    unsigned get_max() const {   return m_max_iter; }
    /**
     * @brief Unpreconditioned CG to estimate maximum Eigenvalue of \f$ A x = \lambda x\f$
     *
     * This implements the original algorithm suggested in <a href="http://www.iam.fmph.uniba.sk/amuc/ojs/index.php/algoritmy/article/view/421">Tichy, On error estimation in the conjugate gradient method: Normwise backward error, Proceedings of the Conference Algorithmy, 323-332, 2016 </a>
     * @note This is just a regular CG algorithm which updates an estimate for the largest Eigenvalue in each iteration and returns once the change is marginal.
     * This means on output \c x is the same as after the same number of iterations of a regular CG method.
     * @param A A symmetric, positive definit matrix
     * @param x Contains an initial value on input and the solution on output.
     * @param b The right hand side vector. x and b may be the same vector.
     * @param ev_max (output) maximum Eigenvalue on output
     * @param eps_ev The desired relative accuracy of the largest Eigenvalue
     *
     * @return Number of iterations used to achieve desired precision or max_iterations
     * @copydoc hide_matrix
     */
    template< class MatrixType, class ContainerType0, class ContainerType1>
    unsigned solve( MatrixType& A, ContainerType0& x, const ContainerType1& b, value_type& ev_max, value_type eps_ev=1e-16);
    ///@brief DEPRECATED: use solve method instead
    ///@copydetails solve(MatrixType&,ContainerType0&,const ContainerType1&,value_type&,value_type)
    template< class MatrixType, class ContainerType0, class ContainerType1>
    unsigned operator()( MatrixType& A, ContainerType0& x, const ContainerType1& b, value_type& ev_max, value_type eps_ev=1e-16)
    {
        return solve(A,x,b,ev_max,eps_ev);
    }
    /**
     * @brief Preconditioned CG to estimate maximum Eigenvalue of the generalized problem \f$ Ax = \lambda M x\f$
     *
     * where \f$ M^{-1}\f$ is the preconditioner.
     * @note This is just a regular PCG algorithm which updates an estimate for the largest Eigenvalue in each iteration and returns once the change is marginal.
     * This means on output \c x is the same as after the same number of iterations of a regular PCG method.
     * @param A A symmetric, positive definit matrix
     * @param x Contains an initial value on input and the solution on output.
     * @param b The right hand side vector. x and b may be the same vector.
     * @param P The preconditioner (\f$ M^{-1}\f$  in the above notation)
     * @param ev_max (output) maximum Eigenvalue on output
     * @param eps_ev The desired accuracy of the largest Eigenvalue
     *
     * @return Number of iterations used to achieve desired precision or max_iterations
     * @copydoc hide_matrix
     */
    template< class MatrixType, class ContainerType0, class ContainerType1, class Preconditioner>
    unsigned solve( MatrixType& A, ContainerType0& x, const ContainerType1& b, Preconditioner& P, value_type& ev_max, value_type eps_ev = 1e-12);
    ///@brief DEPRECATED: use solve method instead
    ///@copydetails solve(MatrixType&,ContainerType0&,const ContainerType1&,Preconditioner&,value_type&,value_type)
    template< class MatrixType, class ContainerType0, class ContainerType1, class Preconditioner>
    unsigned operator()( MatrixType& A, ContainerType0& x, const ContainerType1& b, Preconditioner& P, value_type& ev_max, value_type eps_ev = 1e-12)
    {
        return solve(A,x,b,P,ev_max,eps_ev);
    }
  private:
    ContainerType r, p, ap;
    unsigned m_max_iter;
};

///@cond
template< class ContainerType>
template< class MatrixType, class ContainerType0, class ContainerType1>
unsigned EVE< ContainerType>::solve( MatrixType& A, ContainerType0& x, const ContainerType1&
b, value_type& ev_max, value_type eps_ev)
{
    blas2::symv( A, x, r);
    blas1::axpby( 1., b, -1., r);
    value_type nrm2r_old = blas1::dot( r,r);
    blas1::copy( r, p);
    value_type nrm2r_new, nrmAp;
    value_type alpha = 1., alpha_inv = 1., delta = 0.;
    value_type evdash, gamma = 0., lambda, omega, beta = 0.;
    value_type ev_est = 0.;
    ev_max = 0.;
    for( unsigned i=1; i<m_max_iter; i++)
    {
        lambda = delta*alpha_inv;       // EVE!
        blas2::symv( A, p, ap);
        nrmAp = blas1::dot( p, ap);
        alpha = nrm2r_old /nrmAp;
        alpha_inv = nrmAp /nrm2r_old;   // EVE!
        lambda += alpha_inv;            // EVE!
        blas1::axpby( alpha, p, 1., x);
        blas1::axpby( -alpha, ap, 1., r);
        nrm2r_new = blas1::dot( r, r);

        delta = nrm2r_new /nrm2r_old;                  // EVE!
        evdash = ev_est -lambda;                       // EVE!
        omega = sqrt( evdash*evdash +4.*beta*gamma);   // EVE!
        gamma = 0.5 *(1. -evdash /omega);              // EVE!
        ev_max += omega*gamma;                         // EVE!
        beta = delta*alpha_inv*alpha_inv;              // EVE!
        if( fabs(ev_est-ev_max) < eps_ev*ev_max) {
            return i;
        }
        blas1::axpby(1., r, delta, p);
        nrm2r_old=nrm2r_new;
        ev_est = ev_max;
    }
    return m_max_iter;
};

template< class ContainerType>
template< class Matrix, class ContainerType0, class ContainerType1, class Preconditioner>
unsigned EVE< ContainerType>::solve( Matrix& A, ContainerType0& x, const ContainerType1& b, Preconditioner& P, value_type& ev_max, value_type eps_ev )
{
    blas2::symv( A,x,r);
    blas1::axpby( 1., b, -1., r);
    blas2::symv( P, r, p );
    value_type nrmzr_old = blas1::dot( p,r);
    value_type nrmzr_new, nrmAp;
    value_type alpha = 1., alpha_inv = 1., delta = 0.;
    value_type evdash, gamma = 0., lambda, omega, beta = 0.;
    value_type ev_est = 0.;
    ev_max = 0.;
    for( unsigned i=1; i<m_max_iter; i++)
    {
        lambda = delta*alpha_inv;    // EVE!
        blas2::symv( A, p, ap);
        nrmAp = blas1::dot( p, ap);
        alpha =  nrmzr_old/nrmAp;
        alpha_inv = nrmAp/nrmzr_old; //EVE!
        lambda+= alpha_inv;          //EVE!
        blas1::axpby( alpha, p, 1.,x);
        blas1::axpby( -alpha, ap, 1., r);
        blas2::symv(P,r,ap);
        nrmzr_new = blas1::dot( ap, r);
        delta = nrmzr_new /nrmzr_old;                  // EVE!
        evdash = ev_est -lambda;                       // EVE!
        omega = sqrt( evdash*evdash +4.*beta*gamma);   // EVE!
        gamma = 0.5 *(1. -evdash /omega);              // EVE!
        ev_max += omega*gamma;                         // EVE!
        beta = delta*alpha_inv*alpha_inv;              // EVE!
        if( fabs( ev_est - ev_max) < eps_ev*ev_max)
            return i;
        blas1::axpby(1.,ap, delta, p );
        nrmzr_old=nrmzr_new;
        ev_est = ev_max;
    }
    return m_max_iter;
}
///@endcond

} //namespace dg
#endif //_DG_EVE_

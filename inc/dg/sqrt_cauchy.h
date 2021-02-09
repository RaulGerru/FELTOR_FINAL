#pragma once
#include <boost/math/special_functions.hpp>

#include "blas.h"
#include "helmholtz.h"
#include "lgmres.h"

//! M_PI is non-standard ... so MSVC complains
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/**
 * @brief Matrix class that represents the operator in the Caucha square root integral formula
 *
 * @ingroup matrixoperators
 *
 * discretization of \f[ (-w^2 I -A) x \f]
 * where \f[ A\f] is matrix and w is scalar and x is a vector.
 */
template< class Matrix, class Container>
struct SqrtCauchyIntOp
{
    public:
    using matrix_type = Matrix;
    using container_type = Container;
    using value_type = dg::get_value_type<Container>;
    ///@brief empty object ( no memory allocation)
    SqrtCauchyIntOp() {}
    /**
     * @brief Construct operator \f[ (-w^2 I -A) \f] in cauchy formula
     *
     * @param A symmetric or non-symmetric Matrix, e.g.: a not_normed Helmholtz operator or a symmetric or non-symmetric tridiagonal matrix
     * @param copyable a copyable container 
     * @param multiply_weights multiply inverse weights in front of matrix A (important if matrix A is not_normed) 
     */
    SqrtCauchyIntOp( const Matrix& A, const Container& copyable, const bool& multiply_weights )
    { 
        construct(A, copyable, multiply_weights);
    }
    /**
     * @brief Construct operator \f[ (-w^2 I -A) \f] in cauchy formula
     *
     * @param A symmetric or non-symmetric Matrix, e.g.: a not_normed Helmholtz operator or a symmetric or non-symmetric tridiagonal matrix
     * @param copyable a copyable container 
     * @param multiply_weights multiply inverse weights in front of matrix A (important if matrix A is not_normed) 
     */
    void construct(const Matrix& A, const Container& copyable, const bool& multiply_weights)
    {
        m_A = A;
        m_precond = copyable;
        m_multiply_weights = multiply_weights;
        m_size = copyable.size();
        dg::blas1::scal(m_precond,0.);
        dg::blas1::plus(m_precond,1.0);
        m_weights = m_inv_weights = m_precond;
        m_w=0.0;
    }
    /**
     * @brief Resize vectors (weights, inverse weights and preconditioner)
     *
     * @param new_max new size of vectors
     * 
     */    
    void new_size( unsigned new_max) { 
        m_weights.resize(new_max);
        m_inv_weights.resize(new_max);
        m_precond.resize(new_max);
        m_size = new_max;
        dg::blas1::scal(m_precond,0.);
        dg::blas1::plus(m_precond,1.0);
        m_weights = m_inv_weights = m_precond;
    }
    /**
     * @brief Return the weights of the Helmholtz operator
     * 
     * @return the  weights of the Helmholtz operator
     */
    const Container& weights()const { return m_weights; }
    /**
     * @brief Return the inverse weights of the Helmholtz operator
     *
     * @return the inverse weights of the Helmholtz operator
     */
    const Container& inv_weights()const { return m_inv_weights; }
    /**
     * @brief Return the default preconditioner to use in conjugate gradient
     * 
     * @return the preconditioner of the Helmholtz operator
     */
    const Container& precond()const { return m_precond; }
    /**
     * @brief Set w in the Cauchy integral
     *
     * @param w w in the Cauchy integral
     */
    void set_w( const value_type& w) { m_w=w; }   
/**
     * @brief Set the Matrix A
     *
     * @param A Matrix
     */
    void set_A( const Matrix& A) {
        m_A = A;
    }  
    /**
     * @brief Set the weights
     *
     * @param weights new weights
     */
    void set_weights( const Container& weights) {
        m_weights = weights;
    }
    /**
     * @brief Set the inverse weights
     *
     * @param inv_weights new inverse weights
     */
    void set_inv_weights( const Container& inv_weights) {
        m_inv_weights = inv_weights;
    }
    /**
     * @brief Set the precond
     *
     * @param precond new precond
     */
    void set_precond( const Container& precond) {
        m_precond = precond;
    }
    /**
     * @brief Compute operator
     *
     * i.e. \f[ y= W (w^2 I +V A)*x \f] if weights are multiplied or  \f[ y= (w^2 I +A)*x \f] otherwise
     * @param x left-hand-side
     * @param y result
     */ 
    void symv( const Container& x, Container& y) 
    {
        dg::blas2::symv(m_A, x, y); // A x
        if (m_multiply_weights == true) dg::blas2::symv(m_w*m_w, m_weights, x, 1., y); //W w^2 x + A x 
        else dg::blas1::axpby(m_w*m_w, x, 1., y); // w^2 x + A x 
    } 
  private:
    Container m_weights, m_inv_weights, m_precond;
    Matrix m_A;
    value_type m_w;
    bool m_multiply_weights;
    unsigned m_size;
};


template< class Matrix, class Container>
struct dg::TensorTraits< SqrtCauchyIntOp< Matrix, Container> >
{
    using value_type  = dg::get_value_type<Container>;
    using tensor_category = dg::SelfMadeMatrixTag;
};


/**
    * @brief Compute the square root matrix - vector product via the Cauchy integral
    * i.e. \f[ \sqrt{A} x=  \frac{- 2 K' \sqrt{m}}{\pi N} A \sum_{j=1}^{N} (w_j^2 I -A)^{-1} cn_j dn_j  x \f]
    * A is the matrix, x is the vector, w is a scalar m is the smallest eigenvalue of A, K' is the conjuated complete elliptic integral and cn dn are the jacobi function
 */
template<class Matrix, class Container>
struct CauchySqrtInt
{
  public:
    using matrix_type = Matrix;
    using container_type = Container;
    using value_type = dg::get_value_type<Container>;
    CauchySqrtInt() { }
    /**
     * @brief Construct Rhs operator
     *
     * @param A Helmholtz operator
     * @param g The grid to use
     * @param eps Accuarcy for CG solve
     * @param multiply_weights multiply inverse weights in front of matrix A 
     * @param symmetric true = symmetric A / false = non-symmetric A
     */
    CauchySqrtInt( const Matrix& A, const Container& copyable, value_type eps, const bool& multiply_weights, const bool& symmetric)
    {
        construct(A, copyable, eps, multiply_weights, symmetric);
    }
    /**
     * @brief Construct Rhs operator
     *
     * @param A Helmholtz operator
     * @param g The grid to use
     * @param eps Accuarcy for CG solve
     * @param multiply_weights multiply inverse weights in front of matrix A 
     * @param symmetric true = symmetric A / false = non-symmetric A
     */
    void construct(const Matrix& A, const Container& copyable, value_type eps, const bool& multiply_weights, const bool& symmetric) 
    {
        m_helper = m_helper2 = m_helper3 = copyable;
        m_A = A;
        m_multiply_weights = multiply_weights;
        m_symmetric = symmetric;
        m_eps = eps;
        m_size = m_helper.size();
        std::cout << " size " << m_size<< "\n";
        m_op.construct(m_A, m_helper, m_multiply_weights);
        if (m_symmetric == true) 
        {
            if (m_multiply_weights==true) m_invert.construct( m_helper, m_size*m_size, eps, 1, true, 1.);
            else m_invert.construct( m_helper, m_size*m_size, eps, 1, false, 1.);
        }
        else m_lgmres.construct( m_helper, 30, 10, 10*m_size*m_size);
    }
    /**
     * @brief Resize matrix and set A and vectors and set new size
     *
     * @param new_max new size
     */
     void new_size( unsigned new_max) { 
        m_helper.resize(new_max);
        m_helper2.resize(new_max);
        m_helper3.resize(new_max);
        if (m_symmetric == true) 
        {
            if (m_multiply_weights==true) m_invert.construct( m_helper, new_max*new_max, m_eps, 1, true, 1.);
            else m_invert.construct( m_helper, new_max*new_max, m_eps, 1, false, 1.);
        }
        else m_lgmres.construct( m_helper, 30, 10, 10*new_max*new_max);
        m_op.new_size(new_max);
        m_size = new_max;
        std::cout << " size " << m_size<< "\n";
    } 
    ///@brief Get the current size of vectors
    ///@return the current vector size
    unsigned get_size() const {return m_size;}
    /**
     * @brief Set the Matrix A
     *
     * @param A Matrix
     */
     void set_A( const Matrix& A) { 
         m_A = A; 
         m_op.set_A(A);
    } 
    /**
     * @brief Set the weights
     *
     * @param weights weights
     */
    void set_weights( const Container& weights) {
        m_op.set_weights(weights);
    }
    /**
     * @brief Set the inverse weights
     *
     * @param inv_weights inverse weights
     */
    void set_inv_weights( const Container& inv_weights) {
        m_op.set_inv_weights(inv_weights);
    }
    /**
     * @brief Set the preconditioner
     *
     * @param precond preconditioner
     */
    void set_precond( const Container& precond) {
        m_op.set_precond(precond);
    }
    /**
     * @brief Compute rhs term (including inversion of lhs) 
     *
     * i.e. \f[ b=  \frac{- 2 K' \sqrt{m}}{\pi N} V A \sum_{j=1}^{N} (w^2 I -V A)^{-1} cn dn  x \f]
     * @param y  is \f[ y\f]
     * @param b is \f[ b\approx \sqrt{V A} x\f]
     * @note The Jacobi elliptic functions are related to the Mathematica functions via jacobi_cn(k,u ) = JacobiCN_(u,k^2), ... and the complete elliptic integral of the first kind via comp_ellint_1(k) = EllipticK(k^2) 
     */
    void operator()(const Container& x, Container& b, const value_type& minEV, const value_type& maxEV, const unsigned& iter)
    {
        dg::blas1::scal(m_helper3, 0.0);
        value_type sn=0.;
        value_type cn=0.;
        value_type dn=0.;
        value_type w = 0.;
        value_type t=0.;
        value_type sqrtminEV = sqrt(minEV);
        const value_type k2 = minEV/maxEV;
        const value_type sqrt1mk2 = sqrt(1.-k2);
        const value_type Ks=boost::math::ellint_1(sqrt1mk2 );

        const value_type fac = -2.* Ks*sqrtminEV/(M_PI*iter);
        for (unsigned j=1; j<iter+1; j++)
        {
            t  = (j-0.5)*Ks/iter; //imaginary part .. 1i missing
            cn = 1./boost::math::jacobi_cn(sqrt1mk2, t); 
            sn = boost::math::jacobi_sn(sqrt1mk2, t)*cn;
            dn = boost::math::jacobi_dn(sqrt1mk2, t)*cn;
            w = sqrtminEV*sn;
            dg::blas1::axpby(cn*dn, x, 0.0 , m_helper); //m_helper = cn dn x
            m_op.set_w(w);
            if (m_symmetric == true) m_invert( m_op, m_helper2, m_helper);      // m_helper2 = (w^2 +V A)^(-1) cn dn x
            else m_lgmres.solve( m_op, m_helper2, m_helper, m_op.inv_weights(), m_op.inv_weights(), m_eps, 1); 
            dg::blas1::axpby(-fac, m_helper2, 1.0, m_helper3); // m_helper3 += -fac  (w^2 +V A)^(-1) cn dn x
        }
        dg::blas2::symv(m_A, m_helper3, b); // - A fac sum (w^2 +V A)^(-1) cn dn x
        if (m_multiply_weights == true) dg::blas1::pointwiseDot(m_op.inv_weights(),  b, b);  // fac V A (-w^2 I -V A)^(-1) cn dn x

    }
  private:
    Container m_helper, m_helper2, m_helper3;
    Matrix m_A;
    SqrtCauchyIntOp< Matrix, Container> m_op;
    unsigned m_size;
    bool m_multiply_weights, m_symmetric;
    value_type m_eps;
    dg::Invert<Container> m_invert;
    dg::LGMRES<Container> m_lgmres;

};

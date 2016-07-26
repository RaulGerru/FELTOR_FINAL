#pragma once

#include "cusp/transpose.h"
#include "dg/backend/grid.h"
#include "dg/backend/interpolation.cuh"


namespace dg
{
namespace refined
{

/**
 * @brief Constant refinement in one direction 
 */
enum direction
{
    XDIR, //!< x-direction
    YDIR, //!< y-direction
    XYDIR //!< both directions
};

namespace detail
{

thrust::host_vector<double> exponential_ref( unsigned add_x, unsigned node, unsigned n, unsigned N, dg::bc bcx)
{
    if( add_x == 0)
    {
        thrust::host_vector<double> w_( n*N, 1);
        return w_;
    }
    assert( node <= N);
    //there are add_x+1 finer cells per refined cell ...
    thrust::host_vector< double> left( n*N+n*add_x, 1), right(left);
    for( unsigned k=0; k<n; k++)//the original cell and the additional ones
        left[k] = pow( 2, add_x);
    for( unsigned i=0; i<add_x; i++) 
        for( unsigned k=0; k<n; k++)
            left[(i+1)*n+k] = pow( 2, add_x-i);
    //mirror left into right
    for( unsigned i=0; i<right.size(); i++)
        right[i] = left[ (left.size()-1)-i];
    thrust::host_vector< double> both( n*N+2*n*add_x, 1);
    for( unsigned i=0; i<left.size(); i++)
        both[i] *= left[i];
    for( unsigned i=0; i<right.size(); i++)
        both[i+n*add_x] *= right[i];
    if(      node == 0     && bcx != dg::PER) { return left; }
    else if( node == N && bcx != dg::PER) { return right; }
    else if((node == N || node == 0) && bcx == dg::PER) { return both; }
    else 
    {
        thrust::host_vector<double> w_ = both;
        //now shift indices so that refinement is around nodes
        for( unsigned i=0; i<both.size(); i++)
            w_[((add_x+node)*n+i)%both.size()] = both[i];
        return w_;
    }
}

/**
 * @brief Normalize the given weights and compute the abscissas of the grid
 *
 * @param g The grid to be refined
 * @param weights the unnormalized weights
 *
 * @return The abscissas of the new grid
 */
thrust::host_vector<double> normalize_weights_and_compute_abscissas( const Grid1d<double>& g, thrust::host_vector<double>& weights)
{
    //normalize weights
    unsigned Nx_new = weights.size()/g.n();
    for( unsigned i=0;i<weights.size(); i++)
        weights[i] *= (double)g.N()/(double)Nx_new;

    thrust::host_vector<double> boundaries(Nx_new+1), abs(g.n()*Nx_new);
    boundaries[0] = g.x0();
    for( unsigned i=0; i<Nx_new; i++)
    {
        boundaries[i+1] = boundaries[i] + g.lx()/(double)Nx_new/weights[g.n()*i];
        for( unsigned j=0; j<g.n(); j++)
        {
            abs[i*g.n()+j] =  (boundaries[i+1]+boundaries[i])/2. + 
                (boundaries[i+1]-boundaries[i])/2.*g.dlt().abscissas()[j];
        }
    }
    return abs;
}

/**
 * @brief Create 1d refinement weights and abscissas for the exponential refinement around a node 
 *
 * There will be two refined cells at the end except if a corner node is 
 * given and the boundary condition is not periodic. We count nodes from
 * 0 (left corner) to N (right corner). 
 * @param add_x number of additional cells in the cells idx-1 and idx
 * @param node The cells node-1 and node will be refined
 * @param g The 1d grid to refine
 *
 * @param weights A 1d vector of size n*(Nx+add_x) for one-sided refinement and n*(Nx+2*add_x)) for two-sided refinement
 * @param abscissas A 1d vector of size n*(Nx+add_x) for one-sided refinement and n*(Nx+2*add_x)) for two-sided refinement
 * @return the new number of cells
 */
int exponential_ref( unsigned add_x, unsigned node, const Grid1d<double>& g, thrust::host_vector<double>& weights, thrust::host_vector<double>& abscissas)
{
    if( add_x == 0)
    {
        thrust::host_vector<double> w_( g.size(), 1);
        thrust::host_vector<double> abs_= dg::create::abscissas(g);
        weights = w_; abscissas = abs_; 
        return g.N();
    }
    weights = exponential_ref( add_x, node, g.n(), g.N(), g.bcx());
    unsigned Nx_new = weights.size()/g.n();
    abscissas = normalize_weights_and_compute_abscissas( g, weights);
    return Nx_new;
}

}//namespace detail

/**
 * @brief Refined grid 
 */
struct Grid2d : public dg::Grid2d<double>
{
    /**
     * @brief Refine a corner of a grid
     *
     * @param c
     * @param add_x Add number of cells to the existing one
     * @param add_y Add number of cells to the existing one
     * @param x0
     * @param x1
     * @param y0
     * @param y1
     * @param n
     * @param Nx
     * @param Ny
     * @param bcx
     * @param bcy
     */
    Grid2d( unsigned node_x, unsigned node_y, unsigned add_x, unsigned add_y, 
            double x0, double x1, double y0, double y1, 
            unsigned n, unsigned Nx, unsigned Ny, bc bcx = dg::PER, bc bcy = dg::PER) : dg::Grid2d<double>( x0, x1, y0, y1, n, n_new(Nx, add_x, bcx), n_new(Ny, add_y, bcy), bcx, bcy), 
        wx_(size()), wy_(size()), absX_(size()), absY_(size()),
        g_assoc_( x0, x1, y0, y1, n, Nx, Ny, bcx, bcy)
    {
        Grid1d<double> gx( x0, x1, n, Nx, bcx);
        Grid1d<double> gy( y0, y1, n, Ny, bcy);
        thrust::host_vector<double> wx, ax, wy, ay;
        detail::exponential_ref( add_x, node_x, gx, wx, ax);
        detail::exponential_ref( add_y, node_y, gy, wy, ay);
        //now make product space
        for( unsigned i=0; i<wy.size(); i++)
            for( unsigned j=0; j<wx.size(); j++)
            {
                wx_[i*wx.size()+j] = wx[j];
                wy_[i*wx.size()+j] = wy[i];
                absX_[i*wx.size()+j] = ax[j];
                absY_[i*wx.size()+j] = ay[i];
            }
    }

    /**
     * @brief The grid that this object refines
     *
     * @return  2d grid
     */
    dg::Grid2d<double> associated()const {return g_assoc_;}
    /**
     * @brief Return the abscissas in X-direction 
     *
     * @return A 2d vector
     */
    const thrust::host_vector<double>& abscissasX() const {return absX_;} 
    /**
     * @brief Return the abscissas in Y-direction 
     *
     * @return A 2d vector
     */
    const thrust::host_vector<double>& abscissasY() const {return absY_;} 
    /**
     * @brief Return the weights in X-direction 
     *
     * @return A 2d vector
     */
    const thrust::host_vector<double>& weightsX() const {return wx_;} 
    /**
     * @brief Return the weights in Y-direction 
     *
     * @return A 2d vector
     */
    const thrust::host_vector<double>& weightsY() const {return wy_;} 

    private:
    unsigned n_new( unsigned N, unsigned factor, dg::bc bc)
    {
        if( bc == dg::PER) return N + 2*factor; 
        return N + factor;
    }
    thrust::host_vector<double> wx_, wy_; //weights
    thrust::host_vector<double> absX_, absY_; //abscissas 
    dg::Grid2d<double> g_assoc_;
};

template< class container>
struct Grid3d : public dg::Grid3d<double>
{

};
}//namespace refined


namespace create{

cusp::coo_matrix<int, double, cusp::host_memory> interpolation( const dg::refined::Grid2d& g_fine)
{
    dg::Grid2d<double> g = g_fine.associated();
    //determine number of refined cells
    thrust::host_vector<double> x = g_fine.abscissasX();
    thrust::host_vector<double> y = g_fine.abscissasY();
    
    return dg::create::interpolation( x,y, g);

}

cusp::coo_matrix<int, double, cusp::host_memory> projection( const dg::refined::Grid2d& g_fine)
{
    cusp::coo_matrix<int, double, cusp::host_memory> temp = interpolation( g_fine), A;
    cusp::transpose( temp, A);
    return A;
}

cusp::coo_matrix<int, double, cusp::host_memory> smoothing( const dg::refined::Grid2d& g)
{
    cusp::coo_matrix<int, double, cusp::host_memory> A = interpolation(g);
    cusp::coo_matrix<int, double, cusp::host_memory> B = projection(g);
    cusp::coo_matrix<int, double, cusp::host_memory> C;
    cusp::multiply( A, B, C);
    C.sort_by_row_and_column();
    return C; 
}
}//namespace create

}//namespace dg

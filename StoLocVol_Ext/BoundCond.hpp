// the MIT License (MIT)

// Copyright (c) 2020 Acanthus Solutions (www.acanthussol.com)

// Permission is hereby granted, free of charge, to any person obtaining a copy 
// of this software and associated documentation files (the "Software"), to deal 
// in the Software without restriction, including without limitation the  rights 
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
// copies of the Software, and to permit persons to whom the Software is furnished 
// to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in all 
// copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE 
// SOFTWARE.

#pragma once
#include < array >

#include "MatrixOperator.hpp"
#include "DomainGeometry.hpp"



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//~ Open Bound Condition with zero diffusive term ~//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

// Open Bound Condition with zero diffusion in vol-axis 
class OpenBoundDiff
{
public:
    template<std::size_t M>
    static short applyP1Bounds(double drift,
        double exp2y,
        double msh_stp_zmin,
        double msh_stp_zmax,
        const std::array<double, M>& leverage_vect,
        TridiagMat<M>& p1_new);

    template<std::size_t P>
    static short applyP2Bounds(double kappa,
        double theta,
        double lambda,
        double msh_stp_ymin,
        double msh_stp_ymax,
        const std::array<double, P>& grid_y,
        TridiagMat<P>& p2_new);

    template<std::size_t M, std::size_t P>
    static short applyP0Bounds(double drift,
        const std::array<double, M>& leverage_vect,
        double kappa,
        double theta,
        double lambda,
        double rho,
        const SpatialGeometry<M, P>& grid_geo,
        P0Coeff<M, P>& p0_new);

};


template<std::size_t M>
short OpenBoundDiff::applyP1Bounds(double drift,
    double exp2y,
    double msh_stp_zmin,
    double msh_stp_zmax,
    const std::array<double, M>& leverage_vect,
    TridiagMat<M>& p1_new)
{
    // lower bound
    p1_new.m_a[0] = 0.;
    p1_new.m_b[0] = (drift - 0.5 * exp2y * leverage_vect[0] * leverage_vect[0]) / msh_stp_zmin;
    p1_new.m_c[0] = (-drift + 0.5 * exp2y * leverage_vect[1] * leverage_vect[1]) / msh_stp_zmin;

    // upper bound
    p1_new.m_a[M - 1] = (drift - 0.5 * exp2y * leverage_vect[M - 2] * leverage_vect[M - 2]) / msh_stp_zmax;
    p1_new.m_b[M - 1] = (-drift + 0.5 * exp2y * leverage_vect[M - 1] * leverage_vect[M - 1]) / msh_stp_zmax;
    p1_new.m_c[M - 1] = 0.;

    return 0;
}

template<std::size_t P>
short OpenBoundDiff::applyP2Bounds(double kappa,
    double theta,
    double lambda,
    double msh_stp_ymin,
    double msh_stp_ymax,
    const std::array<double, P>& grid_y,
    TridiagMat<P>& p2_new)
{

    // lower bound
    p2_new.m_a[0] = 0.;
    p2_new.m_b[0] = kappa * (theta - grid_y[0] + msh_stp_ymin) / msh_stp_ymin;
    p2_new.m_c[0] = -kappa * (theta - grid_y[1]) / msh_stp_ymin;

    // upper bound
    p2_new.m_a[P - 1] = kappa * (theta - grid_y[P - 1]) / msh_stp_ymax;
    p2_new.m_b[P - 1] = -kappa * (theta - grid_y[P - 1] - msh_stp_ymax) / msh_stp_ymax;
    p2_new.m_c[P - 1] = 0.;

    return 0;
}

template<std::size_t M, std::size_t P>
short OpenBoundDiff::applyP0Bounds(double drift,
    const std::array<double, M>& leverage_vect,
    double kappa,
    double theta,
    double lambda,
    double rho,
    const SpatialGeometry<M, P>& grid_geo,
    P0Coeff<M, P>& p0_new)
{

    // mixed derivative is zero on the bound 
    for (std::size_t y_indx = 0; y_indx < P - 2; ++y_indx) {
        // inf z-bound
        p0_new.m_zmin_yinf_bound[y_indx] = 0.;
        p0_new.m_zmin_yinf_inner[y_indx] = 0.;
        p0_new.m_zmin_ysup_bound[y_indx] = 0.;
        p0_new.m_zmin_ysup_inner[y_indx] = 0.;
        // sup z-bound
        p0_new.m_zmax_yinf_bound[y_indx] = 0.;
        p0_new.m_zmax_yinf_inner[y_indx] = 0.;
        p0_new.m_zmax_ysup_bound[y_indx] = 0.;
        p0_new.m_zmax_ysup_inner[y_indx] = 0.;
    }

    for (std::size_t z_indx = 0; z_indx < M - 2; ++z_indx) {
        // inf y-bound
        p0_new.m_ymin_zinf_bound[z_indx] = 0.;
        p0_new.m_ymin_zinf_inner[z_indx] = 0.;
        p0_new.m_ymin_zsup_bound[z_indx] = 0.;
        p0_new.m_ymin_zsup_inner[z_indx] = 0.;
        // sup y-bound
        p0_new.m_ymax_zinf_bound[z_indx] = 0.;
        p0_new.m_ymax_zinf_inner[z_indx] = 0.;
        p0_new.m_ymax_zsup_bound[z_indx] = 0.;
        p0_new.m_ymax_zsup_inner[z_indx] = 0.;
    }

    // corners :
    p0_new.m_zmin_ymin_corner = std::make_tuple(0., 0., 0.);
    p0_new.m_zmin_ymax_corner = std::make_tuple(0., 0., 0.);
    p0_new.m_zmax_ymin_corner = std::make_tuple(0., 0., 0.);
    p0_new.m_zmax_ymax_corner = std::make_tuple(0., 0., 0.);

    return 0;
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//~ Open Bound Condition with zero advective term ~//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//



// Open Bound Condition with zero convection in vol-axis
class OpenBoundConv
{
public:
    template<std::size_t M>
    static short applyP1Bounds(double drift,
        double exp2y,
        double msh_stp_zmin,
        double msh_stp_zmax,
        const std::array<double, M>& leverage_vect,
        TridiagMat<M>& p1_new);

    template<std::size_t P>
    static short applyP2Bounds(double kappa,
        double theta,
        double lambda,
        double msh_stp_ymin,
        double msh_stp_ymax,
        const std::array<double, P>& grid_y,
        TridiagMat<P>& p2_new);

    template<std::size_t M, std::size_t P>
    static short applyP0Bounds(double drift,
        const std::array<double, M>& leverage_vect,
        double kappa,
        double theta,
        double lambda,
        double rho,
        const SpatialGeometry<M, P>& grid_geo,
        P0Coeff<M, P>& p0_new);

};



template<std::size_t M>
short OpenBoundConv::applyP1Bounds(double drift,
    double exp2y,
    double msh_stp_zmin,
    double msh_stp_zmax,
    const std::array<double, M>& leverage_vect,
    TridiagMat<M>& p1_new)
{
    // lower bound
    p1_new.m_a[0] = 0.;
    p1_new.m_b[0] = (drift - 0.5 * exp2y * leverage_vect[0] * leverage_vect[0]) / msh_stp_zmin;
    p1_new.m_c[0] = (-drift + 0.5 * exp2y * leverage_vect[1] * leverage_vect[1]) / msh_stp_zmin;

    // upper bound
    p1_new.m_a[M - 1] = (drift - 0.5 * exp2y * leverage_vect[M - 2] * leverage_vect[M - 2]) / msh_stp_zmax;
    p1_new.m_b[M - 1] = (-drift + 0.5 * exp2y * leverage_vect[M - 1] * leverage_vect[M - 1]) / msh_stp_zmax;
    p1_new.m_c[M - 1] = 0.;

    return 0;
}

template<std::size_t P>
short OpenBoundConv::applyP2Bounds(double kappa,
    double theta,
    double lambda,
    double msh_stp_ymin,
    double msh_stp_ymax,
    const std::array<double, P>& grid_y,
    TridiagMat<P>& p2_new)
{

    // lower bound
    p2_new.m_a[0] = 0.;
    p2_new.m_b[0] = -(lambda * lambda) / (msh_stp_ymin * msh_stp_ymin);
    p2_new.m_c[0] = kappa + (lambda * lambda) / (msh_stp_ymin * msh_stp_ymin);

    // upper bound
    p2_new.m_a[P - 1] = kappa + (lambda * lambda) / (msh_stp_ymax * msh_stp_ymax);
    p2_new.m_b[P - 1] = -(lambda * lambda) / (msh_stp_ymax * msh_stp_ymax);
    p2_new.m_c[P - 1] = 0.;


    return 0;
}

template<std::size_t M, std::size_t P>
short OpenBoundConv::applyP0Bounds(double drift,
    const std::array<double, M>& leverage_vect,
    double kappa,
    double theta,
    double lambda,
    double rho,
    const SpatialGeometry<M, P>& grid_geo,
    P0Coeff<M, P>& p0_new)
{

    // mixed derivative is zero on the bound 
    for (std::size_t y_indx = 0; y_indx < P - 2; ++y_indx) {
        // inf z-bound
        p0_new.m_zmin_yinf_bound[y_indx] = 0.;
        p0_new.m_zmin_yinf_inner[y_indx] = 0.;
        p0_new.m_zmin_ysup_bound[y_indx] = 0.;
        p0_new.m_zmin_ysup_inner[y_indx] = 0.;
        // sup z-bound
        p0_new.m_zmax_yinf_bound[y_indx] = 0.;
        p0_new.m_zmax_yinf_inner[y_indx] = 0.;
        p0_new.m_zmax_ysup_bound[y_indx] = 0.;
        p0_new.m_zmax_ysup_inner[y_indx] = 0.;
    }

    for (std::size_t z_indx = 0; z_indx < M - 2; ++z_indx) {
        // inf y-bound
        p0_new.m_ymin_zinf_bound[z_indx] = 0.;
        p0_new.m_ymin_zinf_inner[z_indx] = 0.;
        p0_new.m_ymin_zsup_bound[z_indx] = 0.;
        p0_new.m_ymin_zsup_inner[z_indx] = 0.;
        // sup y-bound
        p0_new.m_ymax_zinf_bound[z_indx] = 0.;
        p0_new.m_ymax_zinf_inner[z_indx] = 0.;
        p0_new.m_ymax_zsup_bound[z_indx] = 0.;
        p0_new.m_ymax_zsup_inner[z_indx] = 0.;
    }

    // corners :
    p0_new.m_zmin_ymin_corner = std::make_tuple(0., 0., 0.);
    p0_new.m_zmin_ymax_corner = std::make_tuple(0., 0., 0.);
    p0_new.m_zmax_ymin_corner = std::make_tuple(0., 0., 0.);
    p0_new.m_zmax_ymax_corner = std::make_tuple(0., 0., 0.);

    return 0;
}



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//~ No-Flux Boundary Condition  ~//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//


class ZeroFlowCond
{
public:
    template<std::size_t M>
    static short applyP1Bounds(double drift,
        double exp2y,
        double msh_stp_zmin,
        double msh_stp_zmax,
        const std::array<double, M>& leverage_vect,
        TridiagMat<M>& p1_new);

    template<std::size_t P>
    static short applyP2Bounds(double kappa,
        double theta,
        double lambda,
        double msh_stp_ymin,
        double msh_stp_ymax,
        const std::array<double, P>& grid_y,
        TridiagMat<P>& p2_new);

    template<std::size_t M, std::size_t P>
    static short applyP0Bounds(double drift,
        const std::array<double, M>& leverage_vect,
        double kappa,
        double theta,
        double lambda,
        double rho,
        const SpatialGeometry<M, P>& grid_geo,
        P0Coeff<M, P>& p0_new);

};








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
#include < vector >

#include "MatrixOperator.hpp"
#include "DomainGeometry.hpp"
#include "ThomasAlgo.hpp"
#include "BoundCond.hpp"
#include "MQRModel.hpp"



//~~~~~~~~~~~~~~~//
//~ ADI Schemes ~//
//~~~~~~~~~~~~~~~//


// function common to all schemes 
template< std::size_t M, typename T >
short buildNewP1Mat(double drift,
    double exp2y,
    const std::array<double, M>& leverage_vect,
    const std::array<double, M - 1>& msh_stp_z,
    TridiagMat<M>& p1_new)
{
    short status{ 0 };

    T::applyP1Bounds<M>(drift, exp2y, msh_stp_z[0], msh_stp_z[M - 2], leverage_vect, p1_new);

    // for a given j on vol axis,  
    for (std::size_t i = 1; i < M - 1; ++i) {
        p1_new.m_a[i] = msh_stp_z[i] / (msh_stp_z[i - 1] * (msh_stp_z[i] + msh_stp_z[i - 1]));
        p1_new.m_a[i] *= (drift + leverage_vect[i - 1] * leverage_vect[i - 1] * exp2y * (1. / msh_stp_z[i] - 0.5));

        p1_new.m_b[i] = (msh_stp_z[i] - msh_stp_z[i - 1]) / (msh_stp_z[i] * msh_stp_z[i - 1])
            * (-drift + 0.5 * leverage_vect[i] * leverage_vect[i] * exp2y);
        p1_new.m_b[i] -= 1. / (msh_stp_z[i] * msh_stp_z[i - 1]) * leverage_vect[i] * leverage_vect[i] * exp2y;

        p1_new.m_c[i] = msh_stp_z[i - 1] / (msh_stp_z[i] * (msh_stp_z[i] + msh_stp_z[i - 1]));
        p1_new.m_c[i] *= (-drift + leverage_vect[i + 1] * leverage_vect[i + 1] * exp2y * (1. / msh_stp_z[i - 1] + 0.5));
    }

    return status;
}


// function common to all schemes 
template< std::size_t M >
short buildS1Mat(double alpha,
    double dt,
    const TridiagMat<M>& p1_new,
    TridiagMat<M>& s1_mat)
{

    s1_mat.m_a[0] = 0.;
    s1_mat.m_b[0] = 1. - alpha * dt * p1_new.m_b[0];
    s1_mat.m_c[0] = -alpha * dt * p1_new.m_c[0];

    for (std::size_t i = 1; i < M - 1; ++i) {
        s1_mat.m_a[i] = -alpha * dt * p1_new.m_a[i];
        s1_mat.m_b[i] = 1. - alpha * dt * p1_new.m_b[i];
        s1_mat.m_c[i] = -alpha * dt * p1_new.m_c[i];
    }

    s1_mat.m_a[M - 1] = -alpha * dt * p1_new.m_a[M - 1];
    s1_mat.m_b[M - 1] = 1. - alpha * dt * p1_new.m_b[M - 1];
    s1_mat.m_c[M - 1] = 0.;

    return 0;
}


// function common to all schemes 
template< std::size_t P, typename T >
short buildNewP2Mat(double kappa,
    double theta,
    double lambda,
    const std::array<double, P>& grid_y,
    const std::array<double, P - 1>& msh_stp_y,
    TridiagMat<P>& p2_new)
{
    short status{ 0 };

    T::applyP2Bounds<P>(kappa, theta, lambda, msh_stp_y[0], msh_stp_y[P - 2], grid_y, p2_new);

    for (std::size_t j = 1; j < P - 1; ++j) {
        p2_new.m_a[j] = msh_stp_y[j] * kappa * (theta - grid_y[j - 1]) + lambda * lambda;
        p2_new.m_a[j] /= (msh_stp_y[j - 1] * (msh_stp_y[j - 1] + msh_stp_y[j]));

        p2_new.m_b[j] = -kappa * (msh_stp_y[j] - msh_stp_y[j - 1]) * (theta - grid_y[j]) - lambda * lambda;
        p2_new.m_b[j] /= (msh_stp_y[j - 1] * msh_stp_y[j]);

        p2_new.m_c[j] = -msh_stp_y[j - 1] * kappa * (theta - grid_y[j + 1]) + lambda * lambda;
        p2_new.m_c[j] /= (msh_stp_y[j] * (msh_stp_y[j - 1] + msh_stp_y[j]));
    }

    return status;
}


// function common to all schemes 
template< std::size_t P >
short buildS2Mat(double alpha,
    double dt,
    const TridiagMat<P>& p2_new,
    TridiagMat<P>& s2_mat)
{

    s2_mat.m_a[0] = 0.;
    s2_mat.m_b[0] = 1. - alpha * dt * p2_new.m_b[0];
    s2_mat.m_c[0] = -alpha * dt * p2_new.m_c[0];

    for (std::size_t i = 1; i < P - 1; ++i) {
        s2_mat.m_a[i] = -alpha * dt * p2_new.m_a[i];
        s2_mat.m_b[i] = 1. - alpha * dt * p2_new.m_b[i];
        s2_mat.m_c[i] = -alpha * dt * p2_new.m_c[i];
    }

    s2_mat.m_a[P - 1] = -alpha * dt * p2_new.m_a[P - 1];
    s2_mat.m_b[P - 1] = 1. - alpha * dt * p2_new.m_b[P - 1];
    s2_mat.m_c[P - 1] = 0.;

    return 0;
}


// function common to all schemes 
template< std::size_t M, std::size_t P, typename T >
short buildNewP0Mat(double drift,
    double kappa,
    double theta,
    double lambda,
    double rho,
    const std::array<double, M>& leverage_vect,
    const SpatialGeometry<M, P> grid_geo,
    P0Coeff<M, P>& p0_new)
{
    // mixed derivative factor:
    p0_new.m_coef_mult = 0.25 * lambda * rho;

    // boundary conditions
    T::applyP0Bounds(drift, leverage_vect,
        kappa, theta, lambda, rho,
        grid_geo,
        p0_new);

    return 0;
}

// function common to all schemes 
template< std::size_t M, std::size_t P >
short buildR1j(std::size_t y_indx,
    double alpha,
    double dt,
    const TridiagMat< M >& p1_mat,
    const TridiagMat< P >& p2_mat,
    const P0Coeff< M, P >& p0_coeff,
    const std::array<double, M>& leverage_vect,
    const std::array<double, M - 1>& msh_stp_z,
    const std::array<double, P - 1>& msh_stp_y,
    const std::array<double, P>& expy_vect,
    const std::vector<std::array<double, P>>& prv_proba,
    std::array< double, M >& R1_j)
{
    // different matrices contributions
    double id_contr{ 0. };
    double p1_contr{ 0. };
    double p2_contr{ 0. };
    double p0_contr{ 0. };

    // y_min bound
    if (y_indx == 0) {

        // z_indx = 0
        id_contr = prv_proba[0][y_indx];

        p1_contr = p1_mat.m_b[0] * prv_proba[0][y_indx] + p1_mat.m_c[0] * prv_proba[1][y_indx];
        p1_contr *= dt * (1. - alpha);

        p2_contr = p2_mat.m_b[y_indx] * prv_proba[0][y_indx]
            + p2_mat.m_c[y_indx] * prv_proba[0][y_indx + 1];
        p2_contr *= dt;

        // corner: trigonometric direction order
        p0_contr = std::get<0>(p0_coeff.m_zmin_ymin_corner) * prv_proba[1][0]
            + std::get<1>(p0_coeff.m_zmin_ymin_corner) * prv_proba[1][1]
            + std::get<2>(p0_coeff.m_zmin_ymin_corner) * prv_proba[0][1];
        p0_contr *= dt;

        R1_j[0] = id_contr + p1_contr + p2_contr + p0_contr;

        // inner points z_indx 
        for (std::size_t z_indx = 1; z_indx < M - 1; ++z_indx) {
            id_contr = prv_proba[z_indx][y_indx];

            p1_contr = p1_mat.m_a[z_indx] * prv_proba[z_indx - 1][y_indx]
                + p1_mat.m_b[z_indx] * prv_proba[z_indx][y_indx]
                + p1_mat.m_c[z_indx] * prv_proba[z_indx + 1][y_indx];
            p1_contr *= dt * (1. - alpha);

            p2_contr = p2_mat.m_b[y_indx] * prv_proba[z_indx][y_indx]
                + p2_mat.m_c[y_indx] * prv_proba[z_indx][y_indx + 1];
            p2_contr *= dt;

            p0_contr = p0_coeff.m_ymin_zinf_bound[z_indx - 1] * prv_proba[z_indx - 1][0]
                + p0_coeff.m_ymin_zinf_inner[z_indx - 1] * prv_proba[z_indx - 1][1]
                + p0_coeff.m_ymin_zsup_bound[z_indx - 1] * prv_proba[z_indx + 1][0]
                + p0_coeff.m_ymin_zsup_inner[z_indx - 1] * prv_proba[z_indx + 1][1];
            p0_contr *= dt;

            R1_j[z_indx] = id_contr + p1_contr + p2_contr + p0_contr;
        }

        // z_indx = M-1
        id_contr = prv_proba[M - 1][y_indx];

        p1_contr = p1_mat.m_a[M - 1] * prv_proba[M - 2][y_indx] + p1_mat.m_b[M - 1] * prv_proba[M - 1][y_indx];
        p1_contr *= dt * (1. - alpha);

        p2_contr = p2_mat.m_b[y_indx] * prv_proba[M - 1][y_indx]
            + p2_mat.m_c[y_indx] * prv_proba[M - 1][y_indx + 1];
        p2_contr *= dt;

        // corner: trigonometric direction order
        p0_contr = std::get<0>(p0_coeff.m_zmax_ymin_corner) * prv_proba[M - 1][1]
            + std::get<1>(p0_coeff.m_zmax_ymin_corner) * prv_proba[M - 2][1]
            + std::get<2>(p0_coeff.m_zmax_ymin_corner) * prv_proba[M - 2][0];
        p0_contr *= dt;

        R1_j[M - 1] = id_contr + p1_contr + p2_contr + p0_contr;

    } // y_indx == 0

    if ((0 < y_indx) && (y_indx < P - 1)) {

        // z_indx = 0
        id_contr = prv_proba[0][y_indx];

        p1_contr = p1_mat.m_b[0] * prv_proba[0][y_indx] + p1_mat.m_c[0] * prv_proba[1][y_indx];
        p1_contr *= dt * (1. - alpha);

        p2_contr = p2_mat.m_a[y_indx] * prv_proba[0][y_indx - 1]
            + p2_mat.m_b[y_indx] * prv_proba[0][y_indx]
            + p2_mat.m_c[y_indx] * prv_proba[0][y_indx + 1];
        p2_contr *= dt;

        p0_contr = p0_coeff.m_zmin_yinf_bound[y_indx - 1] * prv_proba[0][y_indx - 1]
            + p0_coeff.m_zmin_yinf_inner[y_indx - 1] * prv_proba[1][y_indx - 1]
            + p0_coeff.m_zmin_ysup_bound[y_indx - 1] * prv_proba[0][y_indx + 1]
            + p0_coeff.m_zmin_ysup_inner[y_indx - 1] * prv_proba[1][y_indx + 1];
        p0_contr *= dt;

        R1_j[0] = id_contr + p1_contr + p2_contr + p0_contr;

        // inner points z_indx 
        for (std::size_t z_indx = 1; z_indx < M - 1; ++z_indx) {
            id_contr = prv_proba[z_indx][y_indx];

            p1_contr = p1_mat.m_a[z_indx] * prv_proba[z_indx - 1][y_indx]
                + p1_mat.m_b[z_indx] * prv_proba[z_indx][y_indx]
                + p1_mat.m_c[z_indx] * prv_proba[z_indx + 1][y_indx];
            p1_contr *= dt * (1. - alpha);

            p2_contr = p2_mat.m_a[y_indx] * prv_proba[z_indx][y_indx - 1]
                + p2_mat.m_b[y_indx] * prv_proba[z_indx][y_indx]
                + p2_mat.m_c[y_indx] * prv_proba[z_indx][y_indx + 1];
            p2_contr *= dt;

            p0_contr = leverage_vect[z_indx + 1] * expy_vect[y_indx + 1] * prv_proba[z_indx + 1][y_indx + 1]
                + leverage_vect[z_indx - 1] * expy_vect[y_indx - 1] * prv_proba[z_indx - 1][y_indx - 1]
                - leverage_vect[z_indx + 1] * expy_vect[y_indx - 1] * prv_proba[z_indx + 1][y_indx - 1]
                - leverage_vect[z_indx - 1] * expy_vect[y_indx + 1] * prv_proba[z_indx - 1][y_indx + 1];
            p0_contr /= msh_stp_z[z_indx] * msh_stp_y[y_indx] + msh_stp_z[z_indx - 1] * msh_stp_y[y_indx - 1]
                + msh_stp_z[z_indx - 1] * msh_stp_y[y_indx] + msh_stp_z[z_indx] * msh_stp_y[y_indx - 1];
            p0_contr *= p0_coeff.m_coef_mult * dt;

            R1_j[z_indx] = id_contr + p1_contr + p2_contr + p0_contr;
        }

        // z_indx = M-1
        id_contr = prv_proba[M - 1][y_indx];

        p1_contr = p1_mat.m_a[M - 1] * prv_proba[M - 2][y_indx] + p1_mat.m_b[M - 1] * prv_proba[M - 1][y_indx];
        p1_contr *= dt * (1. - alpha);

        p2_contr = p2_mat.m_a[y_indx] * prv_proba[M - 1][y_indx - 1]
            + p2_mat.m_b[y_indx] * prv_proba[M - 1][y_indx]
            + p2_mat.m_c[y_indx] * prv_proba[M - 1][y_indx + 1];
        p2_contr *= dt;

        p0_contr = p0_coeff.m_zmax_yinf_bound[y_indx - 1] * prv_proba[M - 1][y_indx - 1]
            + p0_coeff.m_zmax_yinf_inner[y_indx - 1] * prv_proba[M - 2][y_indx - 1]
            + p0_coeff.m_zmax_ysup_bound[y_indx - 1] * prv_proba[M - 1][y_indx + 1]
            + p0_coeff.m_zmax_ysup_inner[y_indx - 1] * prv_proba[M - 2][y_indx + 1];
        p0_contr *= dt;

        R1_j[M - 1] = id_contr + p1_contr + p2_contr + p0_contr;

    } // inner y_indx point

    // y_max bound
    if (y_indx == P - 1) {

        // z_indx = 0
        id_contr = prv_proba[0][y_indx];

        p1_contr = p1_mat.m_b[0] * prv_proba[0][y_indx] + p1_mat.m_c[0] * prv_proba[1][y_indx];
        p1_contr *= dt * (1. - alpha);

        p2_contr = p2_mat.m_a[y_indx] * prv_proba[0][y_indx - 1]
            + p2_mat.m_b[y_indx] * prv_proba[0][y_indx];
        p2_contr *= dt;

        // corner: trigonometric direction order
        p0_contr = std::get<0>(p0_coeff.m_zmin_ymax_corner) * prv_proba[0][P - 2]
            + std::get<1>(p0_coeff.m_zmin_ymax_corner) * prv_proba[1][P - 2]
            + std::get<2>(p0_coeff.m_zmin_ymax_corner) * prv_proba[1][P - 1];
        p0_contr *= dt;

        R1_j[0] = id_contr + p1_contr + p2_contr + p0_contr;

        // inner points z_indx 
        for (std::size_t z_indx = 1; z_indx < M - 1; ++z_indx) {
            id_contr = prv_proba[z_indx][y_indx];

            p1_contr = p1_mat.m_a[z_indx] * prv_proba[z_indx - 1][y_indx]
                + p1_mat.m_b[z_indx] * prv_proba[z_indx][y_indx]
                + p1_mat.m_c[z_indx] * prv_proba[z_indx + 1][y_indx];
            p1_contr *= dt * (1. - alpha);

            p2_contr = p2_mat.m_a[y_indx] * prv_proba[z_indx][y_indx - 1]
                + p2_mat.m_b[y_indx] * prv_proba[z_indx][y_indx];
            p2_contr *= dt;

            p0_contr = p0_coeff.m_ymax_zinf_bound[z_indx - 1] * prv_proba[z_indx - 1][P - 1]
                + p0_coeff.m_ymax_zinf_inner[z_indx - 1] * prv_proba[z_indx - 1][P - 2]
                + p0_coeff.m_ymax_zsup_bound[z_indx - 1] * prv_proba[z_indx + 1][P - 1]
                + p0_coeff.m_ymax_zsup_inner[z_indx - 1] * prv_proba[z_indx + 1][P - 2];
            p0_contr *= dt;

            R1_j[z_indx] = id_contr + p1_contr + p2_contr + p0_contr;
        }

        // z_indx = M-1
        id_contr = prv_proba[M - 1][y_indx];

        p1_contr = p1_mat.m_a[M - 1] * prv_proba[M - 2][y_indx] + p1_mat.m_b[M - 1] * prv_proba[M - 1][y_indx];
        p1_contr *= dt * (1. - alpha);

        p2_contr = p2_mat.m_a[y_indx] * prv_proba[M - 1][y_indx - 1]
            + p2_mat.m_b[y_indx] * prv_proba[M - 1][y_indx];
        p2_contr *= dt;

        // corner: trigonometric direction order
        p0_contr = std::get<0>(p0_coeff.m_zmax_ymax_corner) * prv_proba[M - 2][P - 1]
            + std::get<1>(p0_coeff.m_zmax_ymax_corner) * prv_proba[M - 2][P - 2]
            + std::get<2>(p0_coeff.m_zmax_ymax_corner) * prv_proba[M - 1][P - 2];
        p0_contr *= dt;

        R1_j[M - 1] = id_contr + p1_contr + p2_contr + p0_contr;

    } // y_indx == P-1

    return 0;
}


// function common to all schemes 
template< std::size_t M, std::size_t P >
short buildR2(std::size_t z_indx,
    double alpha,
    double dt,
    const TridiagMat<P> p2_mat,
    const std::vector< std::array<double, P>>& uprob_mat,
    const std::vector< std::array<double, M>>& u_tilde,
    std::array<double, P>& R2)
{

    R2[0] = p2_mat.m_b[0] * uprob_mat[z_indx][0] + p2_mat.m_c[0] * uprob_mat[z_indx][1];
    R2[0] *= -alpha * dt;
    R2[0] += u_tilde[0][z_indx];

    for (std::size_t y_indx = 1; y_indx < P - 1; ++y_indx) {
        R2[y_indx] = p2_mat.m_a[y_indx] * uprob_mat[z_indx][y_indx - 1]
            + p2_mat.m_b[y_indx] * uprob_mat[z_indx][y_indx]
            + p2_mat.m_c[y_indx] * uprob_mat[z_indx][y_indx + 1];
        R2[y_indx] *= -alpha * dt;
        R2[y_indx] += u_tilde[y_indx][z_indx];
    }

    R2[P - 1] = p2_mat.m_a[P - 1] * uprob_mat[z_indx][P - 2] + p2_mat.m_b[P - 1] * uprob_mat[z_indx][P - 1];
    R2[P - 1] *= -alpha * dt;
    R2[P - 1] += u_tilde[P - 1][z_indx];

    return 0;
}



//~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//~ Douglas-Rachford Scheme ~//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~//

template< std::size_t M, std::size_t P, typename B >
class DOScheme
{
public:
    typedef std::vector<TridiagMat<M>> p1_mat_vec;
    typedef TridiagMat<P> p2_mat;
    typedef std::vector<std::array<double, P>> proba_mat;


    DOScheme() { initMatrices(); };


    short buildNewMatrices(const TnModelData& model_data_t,
        const std::array<double, M>& leverage_vect,
        const SpatialGeometry<M, P>& grid_geo);

    short solveNext(double dt,
        const TnModelData& model_data_t,
        const std::array<double, M>& leverage_vect,
        const SpatialGeometry<M, P>& grid_geo,
        const proba_mat& prv_proba,
        proba_mat& new_proba);


private:
    // semi implicit
    const double m_alpha = 0.5;

    // P0 coefficients 
    template <std::size_t M, std::size_t P>
    using p0coeff = P0Coeff<M, P>;
    p0coeff<M, P> m_p0_coefs;


    // P1 and P2 matrices
    std::vector<TridiagMat<M>> m_p1_matA;
    std::vector<TridiagMat<M>> m_p1_matB;
    std::vector<std::array<double, M>> m_Utilde;
    TridiagMat<P> m_p2_matA;
    TridiagMat<P> m_p2_matB;

    void initMatrices() {
        std::array<double, M> zero_vect;
        zero_vect.fill(0.);
        m_Utilde.assign(P, zero_vect);
        TridiagMat<M> zero_trid = { zero_vect,zero_vect,zero_vect };
        m_p1_matA.assign(P, zero_trid);
        m_p1_matB.assign(P, zero_trid);
    };


    // P1 and P2 matrices pointers
    const p1_mat_vec* m_rp1_ptr;//r-value
    p1_mat_vec* m_lp1_ptr;      //l-value
    const p2_mat* m_rp2_ptr;    //r-value
    p2_mat* m_lp2_ptr;          //l-value

};



template< std::size_t M, std::size_t P, typename B >
short DOScheme<M, P, B>::buildNewMatrices(const TnModelData& model_data_t,
    const std::array<double, M>& leverage_vect,
    const SpatialGeometry<M, P>& grid_geo)
{

    // pointers initializations:
    m_lp1_ptr = &m_p1_matA;
    m_lp2_ptr = &m_p2_matA;
    m_rp1_ptr = &m_p1_matB;
    m_rp2_ptr = &m_p2_matB;

    // we build the P1 matrices:
    for (std::size_t y_indx = 0; y_indx < P; ++y_indx) {

        buildNewP1Mat<M, B>(model_data_t.m_tn_drift, grid_geo.m_exp2y[y_indx], leverage_vect, grid_geo.m_gridsteps_z,
            (*m_lp1_ptr)[y_indx]);
    }

    // we build the P2 matrix:
    buildNewP2Mat<P, B>(model_data_t.m_tn_kappa, model_data_t.m_tn_theta, model_data_t.m_tn_lambda,
        grid_geo.m_grid_y, grid_geo.m_gridsteps_y, *m_lp2_ptr);

    // we update P0
    buildNewP0Mat< M, P, B >(model_data_t.m_tn_drift, model_data_t.m_tn_kappa, model_data_t.m_tn_theta, model_data_t.m_tn_lambda, model_data_t.m_tn_rho,
        leverage_vect, grid_geo, m_p0_coefs);


    return 0;
}

template< std::size_t M, std::size_t P, typename B >
short DOScheme<M, P, B>::solveNext(double dt,
    const TnModelData& model_data_t,
    const std::array<double, M>& leverage_vect,
    const SpatialGeometry<M, P>& grid_geo,
    const proba_mat& prv_proba,
    proba_mat& new_proba)
{
    short status{ 0 };

    // Matrices for Thomas algorithm
    TridiagMat<M> s1_mat;
    std::array<double, M> r1j_vect;
    std::array<double, M> wrk1_vect;

    TridiagMat<P> s2_mat;
    std::array<double, P> r2i_vect;
    std::array<double, P> wrk2_vect;


    // matrix pointers swaps:
    const p1_mat_vec* tmp1_cptr = m_rp1_ptr;
    m_rp1_ptr = m_lp1_ptr;
    m_lp1_ptr = const_cast<p1_mat_vec*>(tmp1_cptr);

    const p2_mat* tmp2_cptr = m_rp2_ptr;
    m_rp2_ptr = m_lp2_ptr;
    m_lp2_ptr = const_cast<p2_mat*>(tmp2_cptr);

    // step DO1:
    for (std::size_t y_indx = 0; y_indx < P; ++y_indx) {

        // we compute the new P1 and S1 matrices 
        buildNewP1Mat<M, B>(model_data_t.m_tn_drift, grid_geo.m_exp2y[y_indx], leverage_vect, grid_geo.m_gridsteps_z,
            (*m_lp1_ptr)[y_indx]);

        buildS1Mat<M>(m_alpha, dt, (*m_lp1_ptr)[y_indx], s1_mat);

        // we build the rhs term R1j:
        buildR1j<M, P>(y_indx, m_alpha, dt,
            (*m_rp1_ptr)[y_indx], *m_rp2_ptr, m_p0_coefs,
            leverage_vect,
            grid_geo.m_gridsteps_z, grid_geo.m_gridsteps_y, grid_geo.m_expy,
            prv_proba, r1j_vect);

        // we solve for m_Utilde[y_indx]:
        status = thomas_inv(M, s1_mat.m_a, s1_mat.m_b, s1_mat.m_c, r1j_vect, m_Utilde[y_indx], wrk1_vect);
    }

    // step DO2:
    // we compute the new P2 and S2 matrices
    buildNewP2Mat<P, B>(model_data_t.m_tn_kappa, model_data_t.m_tn_theta, model_data_t.m_tn_lambda,
        grid_geo.m_grid_y, grid_geo.m_gridsteps_y, *m_lp2_ptr);
    buildS2Mat<P>(m_alpha, dt, *m_lp2_ptr, s2_mat);

    for (std::size_t z_indx = 0; z_indx < M; ++z_indx) {

        // we build the rhs term R2i:
        buildR2(z_indx, m_alpha, dt, *m_rp2_ptr, prv_proba, m_Utilde, r2i_vect);

        // we solve for new_proba[z_indx]:
        status = thomas_inv(P, s2_mat.m_a, s2_mat.m_b, s2_mat.m_c, r2i_vect, new_proba[z_indx], wrk2_vect);
    }

    // we update P0
    buildNewP0Mat< M, P, B >(model_data_t.m_tn_drift, model_data_t.m_tn_kappa, model_data_t.m_tn_theta, model_data_t.m_tn_lambda, model_data_t.m_tn_rho,
        leverage_vect, grid_geo, m_p0_coefs);

    return status;
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//~ Hundsdorfer-Verwer Scheme ~//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

//template< std::size_t M, std::size_t P, typename B >
//class HVScheme
//{
//#error Hundsdorfer-Verwer scheme not provided, pls contact Acanthus Solutions www.acanthussol.com
//} ;


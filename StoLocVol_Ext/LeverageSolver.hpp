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

#include <array>
#include <memory>
#include <cmath>

// results display
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>  

#include "MQRModel.hpp"
#include "ADISchemes.hpp"
#include "BasicSVISmile.hpp"

//~~~~~~~~~~~~~~//
//~ DumpPolicy ~//
//~~~~~~~~~~~~~~//

struct DumpPolicy
{
    enum dump_type { proba = 'p', leverage = 'l' };
    dump_type m_dumptype;
    std::string m_filename;
    std::vector<double> m_dumpdates;
};


//~~~~~~~~~~~~~~~~~~//
//~ LeverageSolver ~//
//~~~~~~~~~~~~~~~~~~//


// M: spot grid size, P: vol grid size, H: time horizon (in months)
// class B: boundary conditions, class S: ADI Scheme, T nbr of events in time grid
template<   std::size_t M,
    std::size_t P,
    unsigned short H,
    typename B,
    template<std::size_t,
    std::size_t, typename> typename S,
    std::size_t T = timestepNumber<H>() >
    class LeverageSolver
{
    static_assert(P % 2 == 1, "Misalignment: Grid point nbr must be odd");
    static_assert(M % 2 == 1, "Misalignment: Grid point nbr must be odd");

public:
    //typedef std::array<std::array<double, P>, M > p_matrix;
    typedef std::vector<std::array<double, P>> p_matrix;

    LeverageSolver() { initProbaMatrices(); };
    LeverageSolver(const PiecewiseLinearCurve& drift,
        const PiecewiseFlatCurve& lambda,
        const PiecewiseFlatCurve& theta,
        const PiecewiseLinearCurve& correl,
        double kappa) {
        initProbaMatrices(), setMQRModel(drift, lambda, theta, correl, kappa);
    };

    void setMQRModel(const PiecewiseLinearCurve& drift,
        const PiecewiseFlatCurve& lambda,
        const PiecewiseFlatCurve& theta,
        const PiecewiseLinearCurve& correl,
        double kappa);

    const MQRModelParameters<M>& getsetMQRModel() const { return *m_model_ptr; };

    template< std::size_t N >
    short calibrateLeverage(const ImpliedVolSmile<N>& vol_smile, const DumpPolicy* display_proba = nullptr);


private:
    void initProbaMatrices() {
        std::array<double, P> zero_vect;
        zero_vect.fill(0.);
        m_U1.assign(M, zero_vect);
        m_U2.assign(M, zero_vect);
    }

    short firstStep(double dt1,
        const std::array<double, M> locvol_vect,
        const TnModelData& mdl_data_t,
        std::array<double, M>& leverage_vect);


    // ADI Scheme
    S< M, P, B > m_motorADI;

    // MQR model and time mapping:
    std::unique_ptr< MQRModelParameters<M> > m_model_ptr;

    p_matrix m_U1;
    p_matrix m_U2;
    const p_matrix* m_rsol;
    p_matrix* m_lsol;

    // spatial discretization
    const SpatialGeometry<M, P> m_spatial_geo ;

    // temporal discretization
    const TemporalGeometry< H > m_temporal_geo;

    
    //const std::array< double, T > m_grid_t = initTimeGridVect<H>();
    //const std::array< double, T - 1 > m_gridsteps_t = initTimeStepVect(m_grid_t);


    std::array< MatBounds, T > m_drift_tmap;
    std::array< MatBounds, T > m_theta_tmap;
    std::array< MatBounds, T > m_lambda_tmap;
    std::array< MatBounds, T > m_correl_tmap;

};




template< std::size_t M, std::size_t P, unsigned short H, typename B, template<std::size_t, std::size_t, typename> typename S, std::size_t T>
void LeverageSolver< M, P, H, B, S, T>::setMQRModel(const PiecewiseLinearCurve& drift,
    const PiecewiseFlatCurve& lambda,
    const PiecewiseFlatCurve& theta,
    const PiecewiseLinearCurve& correl,
    double kappa)
{
    // unique_ptr initialization 
    m_model_ptr = std::make_unique< MQRModelParameters<M> >(drift, lambda, theta, correl, kappa);

    // event mapping:  
    m_drift_tmap = locateValues(drift, m_temporal_geo.m_grid_t);
    m_lambda_tmap = locateValues(lambda, m_temporal_geo.m_grid_t);
    m_theta_tmap = locateValues(theta, m_temporal_geo.m_grid_t);
    m_correl_tmap = locateValues(correl, m_temporal_geo.m_grid_t);
}


template<std::size_t M, std::size_t P>
short approxDiracMass(double dt_1, double drift, double locvol, double kappa, double theta, double lambda, double corr,
    const std::array<double, M>& grid_z,
    const std::array<double, P>& grid_y,
    std::vector<std::array<double, P>>& init_proba)
{
    const double pi{ 3.14159265359 };
    double z_bar = (drift - 0.5 * locvol * locvol) * dt_1;
    double y_bar = kappa * theta * dt_1;
    double sig_z = locvol * std::sqrt(dt_1);
    double sig_y = lambda * std::sqrt(dt_1);
    double zz{ 0. };
    double yy{ 0. };
    double biv{ 0. };

    constexpr double proba_lim = std::numeric_limits<double>::epsilon();

    for (std::size_t z_indx = 0; z_indx < grid_z.size(); ++z_indx) {
        for (std::size_t y_indx = 0; y_indx < grid_y.size(); ++y_indx) {
            zz = (grid_z[z_indx] - z_bar) / sig_z;
            yy = (grid_y[y_indx] - y_bar) / sig_y;
            biv = std::exp(-0.5 / (1. - corr * corr) * (zz * zz - 2 * corr * zz * yy + yy * yy));
            biv /= (2. * pi * sig_z * sig_y * std::sqrt(1. - corr * corr));

            if (biv < proba_lim) {
                biv = 0.;
            }
            init_proba[z_indx][y_indx] = biv;
        }
    }

    return 0;
}

template<std::size_t M>
short firstLeverageVect(double dt_1, double drift, double kappa, double theta, double lambda, double rho,
    const std::array<double, M>& grid_z,
    const std::array<double, M>& locvol_vect,
    std::array<double, M>& leverage_vect)
{
    constexpr double pi{ 3.14159265359 };
    double locvol = locvol_vect[(M - 1) / 2];

    double z_bar = (drift - 0.5 * locvol * locvol) * dt_1;
    double y_bar = kappa * theta * dt_1;
    double sig_z = locvol * std::sqrt(dt_1);
    double sig_y = lambda * std::sqrt(dt_1);
    double zz{ 0. };

    double inv_psi_zi{ 0. };
    double inv_psi_y = std::exp(-2. * (y_bar + (1. - rho * rho) * sig_y * sig_y));
    std::array<double, M> invpsi_vect;


    for (std::size_t z_indx = 0; z_indx < grid_z.size(); ++z_indx) {
        zz = (grid_z[z_indx] - z_bar) / sig_z;
        inv_psi_zi = std::exp(-2. * rho * sig_y * zz);
        invpsi_vect[z_indx] = std::sqrt(inv_psi_y * inv_psi_zi);
        leverage_vect[z_indx] = locvol_vect[z_indx] * std::sqrt(inv_psi_y * inv_psi_zi);
    }

    return 0;
}


template< std::size_t M, std::size_t P >
std::pair<double, double> computePsiFrac(const std::array<double, P - 1>& mshy_vect,
    const std::array<double, P>& proba_vect,
    const std::array<double, P>& exp2y_vect)
{
    double lev{ 0. };
    double psi_num{ 0. };
    double psi_den{ 0. };

    auto prb_it = proba_vect.cbegin();
    auto e2y_it = exp2y_vect.cbegin();
    auto prb_nit = std::next(proba_vect.cbegin());
    auto e2y_nit = std::next(exp2y_vect.cbegin());

    // we compute Psi numerator and denominator:
    for (auto msh_stp : mshy_vect) {
        psi_num += ((*prb_it * *e2y_it) + (*prb_nit * *e2y_nit)) * msh_stp;
        psi_den += (*prb_it + *prb_nit) * msh_stp;
        ++prb_it, ++e2y_it, ++prb_nit, ++e2y_nit;
    }
    std::pair<double, double> res = { psi_num, psi_den };
    return res;
}


template< std::size_t M, std::size_t P >
short computeLeverageVect(const std::array<double, P - 1>& msh_stp_y,
    const std::vector<std::array<double, P>>& proba_mat,
    const std::array<double, P>& exp2y_vect,
    const std::array<double, M>& locvol_vect,
    std::array<double, M>& leverage_vect)
{
    short status{ 0 };
    static_assert(M % 2 == 1, "Misalignment: Grid point nbr must be odd");
    constexpr double proba_lim = std::numeric_limits<double>::epsilon();
    std::array<double, M> invpsi_vect;

    // we exclude negative proba in the center !
    auto psi_frc = computePsiFrac<M, P>(msh_stp_y, proba_mat[(M - 1) / 2], exp2y_vect);
    status = psi_frc.first < proba_lim ? 1 : status;
    status = psi_frc.second < proba_lim ? 1 : status;
    invpsi_vect[(M - 1) / 2] = std::sqrt(psi_frc.second / psi_frc.first);

    bool neg_proba{ false };

    for (std::size_t psi_indx = 1; psi_indx <= (M - 1) / 2; ++psi_indx) {
        //(M - 1) / 2 + lev_indx
        psi_frc = computePsiFrac<M, P>(msh_stp_y, proba_mat[(M - 1) / 2 + psi_indx], exp2y_vect);
        neg_proba = false;
        if ((psi_frc.first < 0.) || (psi_frc.second < 0.)) {
            neg_proba = true;
            status = 1;
        }

        invpsi_vect[(M - 1) / 2 + psi_indx] = psi_frc.second < proba_lim ?
            invpsi_vect[(M - 1) / 2 + psi_indx - 1] : std::sqrt(psi_frc.second / psi_frc.first);

        //(M - 1) / 2 - lev_indx
        psi_frc = computePsiFrac<M, P>(msh_stp_y, proba_mat[(M - 1) / 2 - psi_indx], exp2y_vect);
        neg_proba = false;
        if ((psi_frc.first < 0.) || (psi_frc.second < 0.)) {
            neg_proba = true;
            status = 1;
        }

        invpsi_vect[(M - 1) / 2 - psi_indx] = psi_frc.second < proba_lim ?
            invpsi_vect[(M - 1) / 2 - psi_indx + 1] : std::sqrt(psi_frc.second / psi_frc.first);
    }

    for (std::size_t psi_indx = 0; psi_indx < M; ++psi_indx) {
        leverage_vect[psi_indx] = locvol_vect[psi_indx] * invpsi_vect[psi_indx];
    }

    return status;
}


template< std::size_t M, std::size_t P >
double computeProbaMass(const std::array<double, M - 1>& msh_stp_z,
    const std::array<double, P - 1>& msh_stp_y,
    const std::vector<std::array<double, P>>& proba_mat)
{

    double mass_proba{ 0. };

    // we compute the probability mass :
    for (std::size_t msh_z_indx = 0; msh_z_indx < msh_stp_z.size(); ++msh_z_indx) {
        for (std::size_t msh_y_indx = 0; msh_y_indx < msh_stp_y.size(); ++msh_y_indx) {

            mass_proba += msh_stp_z[msh_z_indx] * msh_stp_y[msh_y_indx] * 0.25 * (
                proba_mat[msh_z_indx][msh_y_indx] + proba_mat[msh_z_indx + 1][msh_y_indx + 1]
                + proba_mat[msh_z_indx][msh_y_indx + 1] + proba_mat[msh_z_indx + 1][msh_y_indx]);

        }
    }
    return mass_proba;
}


template<std::size_t M, std::size_t P>
void dumpProbability(double dump_date,
    const std::vector<std::array<double, P>> proba,
    const SpatialGeometry<M, P>& grid_geo,
    bool is_first_dump,
    const std::string& dump_filename)
{
    std::ofstream dumpfile;
    if (is_first_dump)
        dumpfile.open(dump_filename.c_str(), std::ofstream::out | std::ofstream::trunc);
    else
        dumpfile.open(dump_filename.c_str(), std::ofstream::out | std::ofstream::app);


    if (dumpfile.is_open()) {

        if (is_first_dump) {
            dumpfile << "date" << ',' << "zvalue" << ',' << "yvalue" << ',' << "proba" << '\n';
        }
        // we save only odd points to reduce dump file size:
        for (std::size_t z_indx = 1; z_indx < M; z_indx += 2) {
            for (std::size_t y_indx = 1; y_indx < P; y_indx += 2) {
                // date 
                dumpfile << std::fixed;
                dumpfile << std::setprecision(3);
                dumpfile << dump_date << ',';

                dumpfile << std::setprecision(6);
                // z-value
                dumpfile << grid_geo.m_grid_z[z_indx] << ',';
                // y-value
                dumpfile << grid_geo.m_grid_y[y_indx] << ',';
                // probability
                dumpfile << proba[z_indx][y_indx] << '\n';
            }
        }

        dumpfile.close();
    }

    return;
}


template< std::size_t N, std::size_t M, std::size_t P >
void dumpLocalVolatilityAndLeverage(const SpatialGeometry<M, P>& grid_geo,
    const ImpliedVolSmile<N>& vol_smile,
    const std::vector<std::pair< double, std::array<double, M>>>& leverage_vect,
    const std::string& dump_filename)
{

    const auto& lev_vect = leverage_vect;
    std::vector<double> store_dates;
    store_dates.reserve(lev_vect.size());
    std::for_each(lev_vect.cbegin(), lev_vect.cend(), [&](const std::pair< double, std::array<double, M>>& p) {
        store_dates.push_back(p.first);
        });
    std::vector< MatBounds > vol_mapping = locateValues(vol_smile, store_dates);
    std::array<double, M> locvol_vect;


    std::ofstream dumpfile;
    dumpfile.open(dump_filename.c_str(), std::ofstream::out | std::ofstream::trunc);


    if (dumpfile.is_open()) {

        dumpfile << "type" << ',' << "maturity" << ',' << "strike" << ',' << "value" << '\n';

        // local volatility
        for (std::size_t mat_indx = 0; mat_indx < lev_vect.size(); ++mat_indx) {
            // local vol on the grid:
            getLocalVolVect(vol_smile, vol_mapping[mat_indx], grid_geo.m_grid_z, locvol_vect);

            // local vol display:
            for (std::size_t strike_indx = 0; strike_indx < M; ++strike_indx) {

                dumpfile << "local vol" << ',';

                // date 
                dumpfile << std::fixed;
                dumpfile << std::setprecision(3);
                dumpfile << store_dates[mat_indx] << ',';

                dumpfile << std::setprecision(6);
                dumpfile << grid_geo.m_grid_z[strike_indx] << ',';
                dumpfile << locvol_vect[strike_indx] << '\n';
            }
        }
        // leverage function
        for (std::size_t mat_indx = 0; mat_indx < lev_vect.size(); ++mat_indx) {
            for (std::size_t strike_indx = 0; strike_indx < M; ++strike_indx) {

                dumpfile << "leverage" << ',';

                // date 
                dumpfile << std::fixed;
                dumpfile << std::setprecision(3);
                dumpfile << store_dates[mat_indx] << ',';

                dumpfile << std::setprecision(6);
                dumpfile << grid_geo.m_grid_z[strike_indx] << ',';
                dumpfile << (lev_vect[mat_indx].second)[strike_indx] << '\n';
            }
        }

        dumpfile.close();
    }

    return;
}



template< std::size_t M, std::size_t P, unsigned short H, typename B, template<std::size_t, std::size_t, typename> typename S, std::size_t T>
short LeverageSolver<M, P, H, B, S, T>::firstStep(double dt1,
    const std::array<double, M> locvol_vect,
    const TnModelData& mdl_data_t,
    std::array<double, M>& leverage_vect)
{
    short status{ 0 };
    // first lvalue, rvalue:
    m_rsol = &m_U1;
    m_lsol = &m_U2;


    double locvol = locvol_vect[(M - 1) / 2];
    status = approxDiracMass<M, P>(dt1, mdl_data_t.m_tn_drift, locvol, mdl_data_t.m_tn_kappa,
        mdl_data_t.m_tn_theta, mdl_data_t.m_tn_lambda, mdl_data_t.m_tn_rho,
        m_spatial_geo.m_grid_z, m_spatial_geo.m_grid_y, *m_lsol);
    status = computeLeverageVect<M, P>(m_spatial_geo.m_gridsteps_y, *m_lsol, m_spatial_geo.m_exp2y, locvol_vect, leverage_vect);

    status = firstLeverageVect<M>(dt1,
        mdl_data_t.m_tn_drift, mdl_data_t.m_tn_kappa, mdl_data_t.m_tn_theta, mdl_data_t.m_tn_lambda, mdl_data_t.m_tn_rho,
        m_spatial_geo.m_grid_z, locvol_vect, leverage_vect);

    m_motorADI.buildNewMatrices(mdl_data_t, leverage_vect, m_spatial_geo);

    return status;
}



template< std::size_t M, std::size_t P, unsigned short H, typename B, template<std::size_t, std::size_t, typename> typename S, std::size_t T>
template< std::size_t N>
short LeverageSolver<M, P, H, B, S, T>::calibrateLeverage(const ImpliedVolSmile<N>& vol_smile, const DumpPolicy* display_pol)
{

    short status{ 0 };
    std::array<double, M> leverage_vect;
    TnModelData mdl_data_t;

    // probability mass check
    double proba_mass{ 0. };

    // local vol time mapping
    auto vol_tmap = locateValues<N, T>(vol_smile, m_temporal_geo.m_grid_t);
    std::array<double, M> locvol_vect;

    // first step index: 
    const std::size_t first_time_indx = 2;

    // leverage function storage
    double store_freq = m_model_ptr->getLevFreq();
    std::size_t lev_slice_nbr = int((M / 12. - m_temporal_geo.m_grid_t[first_time_indx]) / store_freq) + 1;
    m_model_ptr->leverageVect().reserve(std::size_t((M / 12. - m_temporal_geo.m_grid_t[first_time_indx]) / store_freq) + 1);

    // proba display management
    std::size_t curr_diplaydate_indx{ 0 };
    bool first_dump{ true };

    // first time step
    status = getLocalVolVect<std::array<double, M>, N>(vol_smile, vol_tmap[first_time_indx], m_spatial_geo.m_grid_z, locvol_vect);
    mdl_data_t.m_tn_kappa = m_model_ptr->getKappa();
    mdl_data_t.m_tn_drift = getValue(m_model_ptr->getDrift(), m_drift_tmap[first_time_indx]);
    mdl_data_t.m_tn_lambda = getValue(m_model_ptr->getLambda(), m_lambda_tmap[first_time_indx]);
    mdl_data_t.m_tn_theta = getValue(m_model_ptr->getThetaMRev(), m_theta_tmap[first_time_indx]);
    mdl_data_t.m_tn_rho = getValue(m_model_ptr->getCorrelSVol(), m_correl_tmap[first_time_indx]);

    // first leverage storage
    status = firstStep(m_temporal_geo.m_grid_t[first_time_indx], locvol_vect, mdl_data_t, leverage_vect);
    double laststore_date = m_temporal_geo.m_grid_t[first_time_indx];
    m_model_ptr->leverageVect().push_back(std::make_pair(laststore_date, leverage_vect));

    // first proba display ?
    if (display_pol && (display_pol->m_dumptype == DumpPolicy::dump_type::proba)) {
        if ((curr_diplaydate_indx < display_pol->m_dumpdates.size()) &&
            (display_pol->m_dumpdates[curr_diplaydate_indx] <= m_temporal_geo.m_grid_t[first_time_indx])) {

            dumpProbability<M, P>(m_temporal_geo.m_grid_t[first_time_indx], *m_lsol, m_spatial_geo, first_dump, display_pol->m_filename);
            first_dump = false;
            ++curr_diplaydate_indx;
        }
    }


    // temporal grid loop 
    for (std::size_t t_indx = first_time_indx + 1; t_indx < m_temporal_geo.m_grid_t.size(); ++t_indx) {
        status = getLocalVolVect<std::array<double, M>, N>(vol_smile, vol_tmap[t_indx], m_spatial_geo.m_grid_z, locvol_vect);
        mdl_data_t.m_tn_kappa = m_model_ptr->getKappa();
        mdl_data_t.m_tn_drift = getValue(m_model_ptr->getDrift(), m_drift_tmap[t_indx]);
        mdl_data_t.m_tn_lambda = getValue(m_model_ptr->getLambda(), m_lambda_tmap[t_indx]);
        mdl_data_t.m_tn_theta = getValue(m_model_ptr->getThetaMRev(), m_theta_tmap[t_indx]);
        mdl_data_t.m_tn_rho = getValue(m_model_ptr->getCorrelSVol(), m_correl_tmap[t_indx]);

        // probability mass check
        proba_mass = computeProbaMass<M, P>(m_spatial_geo.m_gridsteps_z, m_spatial_geo.m_gridsteps_y, *m_lsol);

        // swap l-value / r-value :
        const p_matrix* tmp_cptr = m_rsol;
        m_rsol = m_lsol;
        m_lsol = const_cast<p_matrix*>(tmp_cptr);

        m_motorADI.solveNext(m_temporal_geo.m_gridsteps_t[t_indx - 1], mdl_data_t, leverage_vect, m_spatial_geo, *m_rsol, *m_lsol);

        // we compute the leverage function and store in   
        status = computeLeverageVect<M, P>(m_spatial_geo.m_gridsteps_y, *m_lsol, m_spatial_geo.m_exp2y, locvol_vect, leverage_vect);

        // leverage storage
        if (m_temporal_geo.m_grid_t[t_indx] > laststore_date + store_freq) {
            laststore_date = m_temporal_geo.m_grid_t[t_indx];
            m_model_ptr->leverageVect().push_back(std::make_pair(laststore_date, leverage_vect));
        }

        // proba display
        if (display_pol && (display_pol->m_dumptype == DumpPolicy::dump_type::proba)) {
            if ((curr_diplaydate_indx < display_pol->m_dumpdates.size()) && (display_pol->m_dumpdates[curr_diplaydate_indx] <= m_temporal_geo.m_grid_t[t_indx])) {
                dumpProbability<M, P>(m_temporal_geo.m_grid_t[t_indx], *m_lsol, m_spatial_geo, first_dump, display_pol->m_filename);
                first_dump = false;
                ++curr_diplaydate_indx;
            }
        }


    }
    // last maturity
    m_model_ptr->leverageVect().push_back(std::make_pair(m_temporal_geo.m_grid_t.back(), leverage_vect));

    if (display_pol && (display_pol->m_dumptype == DumpPolicy::dump_type::leverage)) {
        dumpLocalVolatilityAndLeverage<N, M, P>(m_spatial_geo, vol_smile, m_model_ptr->leverageVect(), display_pol->m_filename);
    }


    return status;
}




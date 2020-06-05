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

#include < vector >
#include < array >


//~~~~~~~~~~~~~~~~~~~~//
//~ Curve base class ~//
//~~~~~~~~~~~~~~~~~~~~//

struct MatBounds
{
    double m_dt_inf;
    double m_dt_sup;
    std::size_t m_infIndx;
    std::size_t m_supIndx;
};

struct Curve
{
    Curve() = default;
    Curve(const Curve& crv) : m_dat(crv.m_dat), m_val(crv.m_val) {}
    Curve(const std::vector<std::pair<double, double>> mat_val_vect) {
        m_dat.reserve(mat_val_vect.size());
        m_val.reserve(mat_val_vect.size());
        for (auto mv : mat_val_vect) {
            m_dat.push_back(mv.first);
            m_val.push_back(mv.second);
        }
    }

    const std::vector< double >& maturities() const { return m_dat; };
    std::vector< double >& maturities() { return m_dat; };
    const std::vector< double >& values() const { return m_val; };
    std::vector< double >& values() { return m_val; };

    std::vector< double > m_dat;
    std::vector< double > m_val;

};

double getValue(const Curve& crv, const MatBounds& mb)
{
    const std::vector< double >& val_vect = crv.values();
    double res = val_vect[mb.m_infIndx] * mb.m_dt_sup + val_vect[mb.m_supIndx] * mb.m_dt_inf;
    res /= (mb.m_dt_sup + mb.m_dt_inf);
    return res;
}



//~~~~~~~~~~~~~~~~~~~~~~~~//
//~ PiecewiseLinearCurve ~//
//~~~~~~~~~~~~~~~~~~~~~~~~//


struct PiecewiseFlatCurve : Curve
{
    PiecewiseFlatCurve() = default;
    PiecewiseFlatCurve(const PiecewiseFlatCurve& pfc) : Curve(pfc) {}
    PiecewiseFlatCurve(const std::vector<std::pair<double, double>> mat_val_vect) : Curve(mat_val_vect) {};
};

template< std::size_t N >
std::array< MatBounds, N > locateValues(const PiecewiseFlatCurve& curv, const std::array< double, N >& disc_mat)
{
    std::array< MatBounds, N > res;
    std::size_t crv_indx = 0;
    const auto& crv_mat = curv.maturities();
    auto res_it = res.begin();

    for (std::size_t dsc_indx = 0; dsc_indx < disc_mat.size(); ++dsc_indx) {
        for (; (crv_indx < crv_mat.size()) && (crv_mat[crv_indx] < disc_mat[dsc_indx]); ++crv_indx);

        if (crv_indx < crv_mat.size()) {
            (*res_it).m_dt_inf = 0.;
            (*res_it).m_dt_sup = 1.;
            (*res_it).m_infIndx = crv_indx;
            (*res_it).m_supIndx = crv_indx;
        }
        else {
            (*res_it).m_dt_inf = 1.;
            (*res_it).m_dt_sup = 0.;
            (*res_it).m_infIndx = crv_mat.size() - 1;
            (*res_it).m_supIndx = crv_mat.size() - 1;
        }
        ++res_it;
    }
    return res;
}



//~~~~~~~~~~~~~~~~~~~~~~~~//
//~ PiecewiseLinearCurve ~//
//~~~~~~~~~~~~~~~~~~~~~~~~//


struct PiecewiseLinearCurve : Curve
{
    PiecewiseLinearCurve() = default;
    PiecewiseLinearCurve(const PiecewiseLinearCurve& plc) : Curve(plc) {};
    PiecewiseLinearCurve(const std::vector<std::pair<double, double>> mat_val_vect) : Curve(mat_val_vect) {};
};

template< std::size_t N >
std::array< MatBounds, N > locateValues(const PiecewiseLinearCurve& curv, const std::array< double, N >& disc_mat)
{
    std::array< MatBounds, N > res;
    std::size_t crv_indx = 0;
    const auto& crv_mat = curv.maturities();

    auto res_it = res.begin();

    MatBounds loc{ 0., 1., 0, 0 };

    for (std::size_t dsc_indx = 0; dsc_indx < disc_mat.size(); ++dsc_indx) {
        for (; (crv_indx < crv_mat.size()) && (crv_mat[crv_indx] < disc_mat[dsc_indx]); ++crv_indx);

        if (crv_indx) {
            if (crv_indx < crv_mat.size()) {
                loc.m_dt_inf = disc_mat[dsc_indx] - crv_mat[crv_indx - 1];
                loc.m_dt_sup = crv_mat[crv_indx] - disc_mat[dsc_indx];
                loc.m_infIndx = crv_indx - 1;
                loc.m_supIndx = crv_indx;
            }
            else {
                // beyond the last curve element
                // the curve is considered flattish
                loc.m_dt_inf = 1.;
                loc.m_dt_sup = 0.;
                loc.m_infIndx = crv_mat.size() - 1;
                loc.m_supIndx = crv_mat.size() - 1;
            }
        }
        (*res_it) = loc;
        ++res_it;
    }
    return res;
}


// theta function calibration
PiecewiseFlatCurve calibrateTheta(const PiecewiseFlatCurve& lambda, double kappa)
{
    // we initialize a vector containing exp(kappa ti)
    const std::vector< double >& mat = lambda.maturities();
    std::vector< double > exp_kdelta_ti;

    exp_kdelta_ti.reserve(mat.size());
    double prev_mat = 0.;

    std::for_each(mat.cbegin(), mat.cend(), [&](double curr_mat) {
        exp_kdelta_ti.push_back(std::exp(kappa * (prev_mat - curr_mat)));
        prev_mat = curr_mat;
        });


    // we construct theta (piecewise flat function)
    const std::vector< double >& lbd_val = lambda.values();
    std::vector< double > theta_val;
    theta_val.reserve(lbd_val.size());
    double sum_theta{ 0. };
    double sum_lambda{ 0. };
    double curr_theta{ 0. };
    double exp_ktj{ 0. };

    for (std::size_t i = 0; i < mat.size(); ++i) {
        sum_theta = 0.;
        sum_lambda = 0.;
        for (std::size_t j = 0; j < i; ++j) {
            exp_ktj = std::exp(kappa * (mat[j] - mat[i]));
            sum_theta -= theta_val[j] * exp_ktj * (1. - exp_kdelta_ti[j]);
            sum_lambda -= 0.5 * lbd_val[j] * lbd_val[j] * exp_ktj * exp_ktj * (1. - exp_kdelta_ti[j] * exp_kdelta_ti[j]);
        }
        sum_theta /= (1. - exp_kdelta_ti[i]);
        sum_lambda /= kappa * (1. - exp_kdelta_ti[i]);

        curr_theta = -0.5 * lbd_val[i] * lbd_val[i] / kappa * (1. + exp_kdelta_ti[i]);
        curr_theta += sum_theta;
        curr_theta += sum_lambda;
        theta_val.push_back(curr_theta);
    }

    std::vector<std::pair<double, double>> mat_val_vect;
    mat_val_vect.reserve(mat.size());
    auto mat_it = mat.cbegin();
    auto theta_it = theta_val.cbegin();

    for (; theta_it != theta_val.cend(); ++theta_it, ++mat_it) {
        mat_val_vect.emplace_back(std::pair<double, double>(*mat_it, *theta_it));
    }

    PiecewiseFlatCurve theta_res(mat_val_vect);
    return theta_res;
}


//~~~~~~~~~~~~~~~~~~~~~~//
//~ MQRModelParameters ~//
//~~~~~~~~~~~~~~~~~~~~~~//


template<std::size_t M>
class MQRModelParameters
{
public:
    MQRModelParameters(const MQRModelParameters& mp) : m_drift_sp(mp.m_drift_sp),
        m_lambda_vol(mp.m_lambda_vol), m_rho_corr(mp.m_rho_corr), m_kappa(mp.m_kappa) {};

    MQRModelParameters(const PiecewiseLinearCurve& drift,
        const PiecewiseFlatCurve& lambda,
        const PiecewiseFlatCurve& theta,
        const PiecewiseLinearCurve& correl,
        double kappa) : m_drift_sp(drift), m_lambda_vol(lambda), m_theta_meanrev(theta), m_rho_corr(correl), m_kappa(kappa) {};

    const PiecewiseLinearCurve& getDrift() const { return m_drift_sp; };
    const PiecewiseFlatCurve& getLambda() const { return m_lambda_vol; };
    const PiecewiseLinearCurve& getCorrelSVol() const { return m_rho_corr; };
    const PiecewiseFlatCurve& getThetaMRev() const { return m_theta_meanrev; };
    double getKappa() const { return m_kappa; };

    PiecewiseFlatCurve& getThetaMeanRev() { return m_theta_meanrev; };

    using leverage_vect = std::vector<std::pair< double, std::array<double, M>>>;
    double getLevFreq() const { return m_lev_freq; };
    leverage_vect& leverageVect() { return m_leverage; };
    const leverage_vect& leverageVect() const { return m_leverage; };

private:

    PiecewiseLinearCurve m_drift_sp;

    double m_kappa;
    PiecewiseFlatCurve m_lambda_vol;
    PiecewiseLinearCurve m_rho_corr;
    PiecewiseFlatCurve m_theta_meanrev;

    const double m_lev_freq = 0.25;
    std::vector<std::pair< double, std::array<double, M>>> m_leverage;

};


//~~~~~~~~~~~~~~~//
//~ TnModelData ~//
//~~~~~~~~~~~~~~~//


struct TnModelData
{
    double m_tn_drift;

    double m_tn_kappa;
    double m_tn_lambda;
    double m_tn_theta;
    double m_tn_rho;
};





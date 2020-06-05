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

#include "MQRModel.hpp"


//~~~~~~~~~~~~~~~~~~~//
//~ ImpliedVolSmile ~//
//~~~~~~~~~~~~~~~~~~~//

struct MatSlice
{
    double m_a;
    double m_b;
    double m_rho;
    double m_m;
    double m_sig;
};

template< std::size_t N >
class ImpliedVolSmile
{
public:
    typedef MatSlice OneSlice;
    typedef std::array< double, N > mat_cont;
    typedef std::array< MatSlice, N > smile_cont;

    ImpliedVolSmile() = default;
    ImpliedVolSmile(const mat_cont& mat, const smile_cont& slices) : m_mat(mat), m_smile(slices) {};

    mat_cont& maturities() { return m_mat; };
    const mat_cont& maturities() const { return m_mat; };
    smile_cont& getVolSlices() { return m_smile; };
    const smile_cont& getVolSlices() const { return m_smile; };
    MatSlice& matSlice(std::size_t i) { return m_smile[i]; };
    const MatSlice& matSlice(std::size_t i) const { return m_smile[i]; };

private:
    std::array< double, N > m_mat;
    std::array< OneSlice, N > m_smile;
};


template< std::size_t M, std::size_t N >
std::array< MatBounds, N > locateValues(const ImpliedVolSmile<M>& smile, const std::array< double, N >& disc_mat)
{
    std::array< MatBounds, N > res;
    std::size_t sml_indx = 0;
    const auto& sml_mat = smile.maturities();

    MatBounds loc{ 0., sml_mat[1] - sml_mat[0], 0, 1 };

    for (std::size_t dsc_indx = 0; dsc_indx < disc_mat.size(); ++dsc_indx) {
        for (; (sml_indx < sml_mat.size()) && (sml_mat[sml_indx] < disc_mat[dsc_indx]); ++sml_indx);
        if (sml_indx) {
            if (sml_indx < sml_mat.size()) {
                loc.m_dt_inf = disc_mat[dsc_indx] - sml_mat[sml_indx - 1];
                loc.m_dt_sup = sml_mat[sml_indx] - disc_mat[dsc_indx];
                loc.m_infIndx = sml_indx - 1;
                loc.m_supIndx = sml_indx;
            }
            else {
                loc.m_dt_inf = sml_mat[sml_mat.size() - 1] - sml_mat[sml_mat.size() - 2];
                loc.m_dt_sup = 0.;
                loc.m_infIndx = sml_mat.size() - 2;
                loc.m_supIndx = sml_mat.size() - 1;
            }
        }
        res[dsc_indx] = loc;
    }
    return res;
}

template< std::size_t M >
std::vector< MatBounds > locateValues(const ImpliedVolSmile<M>& smile, const std::vector< double >& disc_mat)
{
    std::vector< MatBounds > res;
    res.reserve(disc_mat.size());
    std::size_t sml_indx = 0;
    const auto& sml_mat = smile.maturities();

    MatBounds loc{ 0., sml_mat[1] - sml_mat[0], 0, 1 };

    for (std::size_t dsc_indx = 0; dsc_indx < disc_mat.size(); ++dsc_indx) {
        for (; (sml_indx < sml_mat.size()) && (sml_mat[sml_indx] < disc_mat[dsc_indx]); ++sml_indx);
        if (sml_indx) {
            if (sml_indx < sml_mat.size()) {
                loc.m_dt_inf = disc_mat[dsc_indx] - sml_mat[sml_indx - 1];
                loc.m_dt_sup = sml_mat[sml_indx] - disc_mat[dsc_indx];
                loc.m_infIndx = sml_indx - 1;
                loc.m_supIndx = sml_indx;
            }
            else {
                loc.m_dt_inf = sml_mat[sml_mat.size() - 1] - sml_mat[sml_mat.size() - 2];
                loc.m_dt_sup = 0.;
                loc.m_infIndx = sml_mat.size() - 2;
                loc.m_supIndx = sml_mat.size() - 1;
            }
        }
        res.push_back(loc);
    }
    return res;
}


template< std::size_t N >
short getImpliedVol(const ImpliedVolSmile< N >& smile, const MatBounds& loc, double log_mn, double& result)
{
    auto vol_slice_inf = smile.matSlice(loc.m_infIndx);
    auto vol_slice_sup = smile.matSlice(loc.m_supIndx);
    // size check 

    double svi_var_inf, svi_var_sup;

    svi_var_inf = vol_slice_inf.m_a + vol_slice_inf.m_b * (vol_slice_inf.m_rho * (log_mn - vol_slice_inf.m_m)
        + std::sqrt((log_mn - vol_slice_inf.m_m) * (log_mn - vol_slice_inf.m_m) + vol_slice_inf.m_sig * vol_slice_inf.m_sig));

    svi_var_sup = vol_slice_sup.m_a + vol_slice_sup.m_b * (vol_slice_sup.m_rho * (log_mn - vol_slice_sup.m_m)
        + std::sqrt((log_mn - vol_slice_sup.m_m) * (log_mn - vol_slice_sup.m_m) + vol_slice_sup.m_sig * vol_slice_sup.m_sig));

    result = (loc.m_dt_inf * svi_var_sup + loc.m_dt_sup * svi_var_inf) / (loc.m_dt_inf + loc.m_dt_sup);

    // from SVI var smile we compute SVI vol: 
    double mat_vol = (smile.maturities())[loc.m_infIndx]; // inf slice maturity
    mat_vol += loc.m_dt_inf;
    result = std::sqrt(result / mat_vol);

    return 0;
}


template< typename T, std::size_t N >
short getImpliedVolVect(const ImpliedVolSmile< N >& smile, const MatBounds& loc, const T& spot_range, T& results)
{
    auto res_it = results.begin();
    auto vol_slice_inf = smile.matSlice(loc.m_infIndx);
    auto vol_slice_sup = smile.matSlice(loc.m_supIndx);
    // size check 

    double mat_vol = (smile.maturities())[loc.m_infIndx]; // inf slice maturity
    mat_vol += loc.m_dt_inf;

    double svi_var_inf, svi_var_sup;

    std::for_each(spot_range.cbegin(), spot_range.cend(), [&](double log_mn) {

        svi_var_inf = vol_slice_inf.m_a + vol_slice_inf.m_b * (vol_slice_inf.m_rho * (log_mn - vol_slice_inf.m_m)
            + std::sqrt((log_mn - vol_slice_inf.m_m) * (log_mn - vol_slice_inf.m_m) + vol_slice_inf.m_sig * vol_slice_inf.m_sig));

        svi_var_sup = vol_slice_sup.m_a + vol_slice_sup.m_b * (vol_slice_sup.m_rho * (log_mn - vol_slice_sup.m_m)
            + std::sqrt((log_mn - vol_slice_sup.m_m) * (log_mn - vol_slice_sup.m_m) + vol_slice_sup.m_sig * vol_slice_sup.m_sig));

        *res_it = (loc.m_dt_inf * svi_var_sup + loc.m_dt_sup * svi_var_inf) / (loc.m_dt_inf + loc.m_dt_sup);
        *res_it = std::sqrt(*res_it / mat_vol);

        ++res_it; }
    );

    return 0;
}


template< std::size_t N >
short getLocalVol(const ImpliedVolSmile< N >& smile, const MatBounds& loc, double log_mn, double& result)
{

    short status{ 0 };
    const MatSlice& vol_slice_inf = smile.matSlice(loc.m_infIndx);
    const MatSlice& vol_slice_sup = smile.matSlice(loc.m_supIndx);

    double var_inf, var_sup, var;
    double der_y_inf, der_y_sup, der_y;
    double der_yy_inf, der_yy_sup, der_yy;
    double w_inf, w_sup;

    // intermediary variable 
    double y_min_m, y_min_m2, sig2, sqrt_trm;
    double locvol_num, locvol_den;

    y_min_m = log_mn - vol_slice_inf.m_m;
    y_min_m2 = y_min_m * y_min_m;
    sig2 = vol_slice_inf.m_sig * vol_slice_inf.m_sig;
    sqrt_trm = std::sqrt(y_min_m2 + sig2);

    var_inf = vol_slice_inf.m_a + vol_slice_inf.m_b * (vol_slice_inf.m_rho * y_min_m + sqrt_trm);
    der_y_inf = vol_slice_inf.m_b * (vol_slice_inf.m_rho + y_min_m / sqrt_trm);
    der_yy_inf = vol_slice_inf.m_b * ((y_min_m2 + sig2 - 0.5 * y_min_m) / std::pow(sqrt_trm, 3.));

    y_min_m = log_mn - vol_slice_sup.m_m;
    y_min_m2 = y_min_m * y_min_m;
    sig2 = vol_slice_sup.m_sig * vol_slice_sup.m_sig;
    sqrt_trm = std::sqrt(y_min_m2 + sig2);

    var_sup = vol_slice_sup.m_a + vol_slice_sup.m_b * (vol_slice_sup.m_rho * y_min_m + sqrt_trm);
    der_y_sup = vol_slice_sup.m_b * (vol_slice_sup.m_rho + y_min_m / sqrt_trm);
    der_yy_sup = vol_slice_sup.m_b * ((y_min_m2 + sig2 - 0.5 * y_min_m) / std::pow(sqrt_trm, 3.));

    w_inf = loc.m_dt_sup / (loc.m_dt_inf + loc.m_dt_sup);
    w_sup = loc.m_dt_inf / (loc.m_dt_inf + loc.m_dt_sup);

    var = w_inf * var_inf + w_sup * var_sup;
    der_y = w_inf * der_y_inf + w_sup * der_y_sup;
    der_yy = w_inf * der_yy_inf + w_sup * der_yy_sup;

    locvol_num = (var_sup - var_inf) / (loc.m_dt_inf + loc.m_dt_sup);
    locvol_den = (1. - 0.5 * log_mn / var * der_y) * (1. - 0.5 * log_mn / var * der_y)
        - 0.25 * (0.25 + 1. / var) * der_y * der_y + 0.5 * der_yy;

#ifdef AS_DEBUG
    // calendar arbitrage 
    if locvol_num < 0.
        status = 1
    // butterfly arbitrage
    if locvol_den < 0.
        status = 2
#endif
    
    result = std::sqrt(locvol_num / locvol_den);
    return status;
}

template< typename T, std::size_t N >
short getLocalVolVect(const ImpliedVolSmile< N >& smile, const MatBounds& loc, const T& spot_range, T& results)
{

    const MatSlice& vol_slice_inf = smile.matSlice(loc.m_infIndx);
    const MatSlice& vol_slice_sup = smile.matSlice(loc.m_supIndx);

    double var_inf, var_sup, var;
    double der_y_inf, der_y_sup, der_y;
    double der_yy_inf, der_yy_sup, der_yy;
    double w_inf, w_sup;

    // intermediary variable
    double y_min_m, y_min_m2, sig2, sqrt_trm;
    double locvol_num, locvol_den;

    auto log_mn_it = spot_range.cbegin();

    for (auto& res_locvol : results) {
        y_min_m = *log_mn_it - vol_slice_inf.m_m;
        y_min_m2 = y_min_m * y_min_m;
        sig2 = vol_slice_inf.m_sig * vol_slice_inf.m_sig;
        sqrt_trm = std::sqrt(y_min_m2 + sig2);

        var_inf = vol_slice_inf.m_a + vol_slice_inf.m_b * (vol_slice_inf.m_rho * y_min_m + sqrt_trm);
        der_y_inf = vol_slice_inf.m_b * (vol_slice_inf.m_rho + y_min_m / sqrt_trm);
        der_yy_inf = vol_slice_inf.m_b * ((y_min_m2 + sig2 - 0.5 * y_min_m) / std::pow(sqrt_trm, 3.));

        y_min_m = *log_mn_it - vol_slice_sup.m_m;
        y_min_m2 = y_min_m * y_min_m;
        sig2 = vol_slice_sup.m_sig * vol_slice_sup.m_sig;
        sqrt_trm = std::sqrt(y_min_m2 + sig2);

        var_sup = vol_slice_sup.m_a + vol_slice_sup.m_b * (vol_slice_sup.m_rho * y_min_m + sqrt_trm);
        der_y_sup = vol_slice_sup.m_b * (vol_slice_sup.m_rho + y_min_m / sqrt_trm);
        der_yy_sup = vol_slice_sup.m_b * ((y_min_m2 + sig2 - 0.5 * y_min_m) / std::pow(sqrt_trm, 3.));

        w_inf = loc.m_dt_sup / (loc.m_dt_inf + loc.m_dt_sup);
        w_sup = loc.m_dt_inf / (loc.m_dt_inf + loc.m_dt_sup);

        var = w_inf * var_inf + w_sup * var_sup;
        der_y = w_inf * der_y_inf + w_sup * der_y_sup;
        der_yy = w_inf * der_yy_inf + w_sup * der_yy_sup;

        locvol_num = (var_sup - var_inf) / (loc.m_dt_inf + loc.m_dt_sup);
        locvol_den = (1. - 0.5 * *log_mn_it / var * der_y) * (1. - 0.5 * *log_mn_it / var * der_y)
            - 0.25 * (0.25 + 1. / var) * der_y * der_y + 0.5 * der_yy;

#ifdef AS_DEBUG
        // calendar arbitrage 
        if locvol_num < 0.
            status = 1
        // butterfly arbitrage
        if locvol_den < 0.
            status = 2
#endif

        res_locvol = std::sqrt(locvol_num / locvol_den);
        log_mn_it++;

    }
    return 0;
}


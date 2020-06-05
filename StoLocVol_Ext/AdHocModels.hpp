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

#include < cmath >

#include "BasicSVISmile.hpp"
#include "MQRModel.hpp"


ImpliedVolSmile<6> adHocVolatility_new()
{
    ImpliedVolSmile<6> smile;
    // ad hoc implied volatility 
    // 2.5m 6m 1y 2y 3y 5y 
    double a{ 0. };
    double b{ 0. };
    double rho{ 0. };
    double m{ 0. };
    double sig{ 0. };

    MatSlice smile_mat;
    smile.maturities() = { 0.20, 0.50, 1.00, 2.00, 3.00, 5.00 };


    // 2.5m
    a = 0.000132941565217;
    b = 0.02;
    rho = -0.85;
    m = -0.00108241;
    sig = 0.036649321739131;

    smile_mat.m_a = a;
    smile_mat.m_b = b;
    smile_mat.m_rho = rho;
    smile_mat.m_m = m;
    smile_mat.m_sig = sig;
    smile.matSlice(0) = smile_mat;


    // 6m 
    a = 0.000381941;
    b = 0.03;
    rho = -0.85;
    m = -0.00108241;
    sig = 0.0594293;

    smile_mat.m_a = a;
    smile_mat.m_b = b;
    smile_mat.m_rho = rho;
    smile_mat.m_m = m;
    smile_mat.m_sig = sig;
    smile.matSlice(1) = smile_mat;

    // 1y   
    a = 0.001407762086957;
    b = 0.03;
    rho = -0.85;
    m = -0.00108241;
    sig = 0.097395930434783;

    smile_mat.m_a = a;
    smile_mat.m_b = b;
    smile_mat.m_rho = rho;
    smile_mat.m_m = m;
    smile_mat.m_sig = sig;
    smile.matSlice(2) = smile_mat;

    // 2y
    a = 0.00345940426087;
    b = 0.03;
    rho = -0.85;
    m = -0.00108241;
    sig = 0.173329191304348;

    smile_mat.m_a = a;
    smile_mat.m_b = b;
    smile_mat.m_rho = rho;
    smile_mat.m_m = m;
    smile_mat.m_sig = sig;
    smile.matSlice(3) = smile_mat;

    // 3y
    a = 0.008003670956522;
    b = 0.032;
    rho = -0.90;
    m = -0.0014;
    sig = 0.240;

    smile_mat.m_a = a;
    smile_mat.m_b = b;
    smile_mat.m_rho = rho;
    smile_mat.m_m = m;
    smile_mat.m_sig = sig;
    smile.matSlice(4) = smile_mat;

    // 5y
    a = 0.02125;
    b = 0.025;
    rho = -0.90;
    m = -0.00108241;
    sig = 0.320;

    smile_mat.m_a = a;
    smile_mat.m_b = b;
    smile_mat.m_rho = rho;
    smile_mat.m_m = m;
    smile_mat.m_sig = sig;
    smile.matSlice(5) = smile_mat;


    return smile;
}



void initAdHocMQRModel(PiecewiseLinearCurve& drift,
    PiecewiseFlatCurve& lambda,
    PiecewiseFlatCurve& theta,
    PiecewiseLinearCurve& correl,
    double& kappa)
{

    // drift (UsdJpy inspired)
    // 6m 0.5       -0.0005
    // 1y 1.00      -0.0005
    // 2y 2.00      -0.0007
    // 3y 3.00      -0.0010
    // 5y 5.00      -0.0010
    std::vector<std::pair<double, double>> drift_mv_vect;
    drift_mv_vect.reserve(5);
    drift_mv_vect.emplace_back(std::pair<double, double>(0.50, -0.0005));
    drift_mv_vect.emplace_back(std::pair<double, double>(1.00, -0.0005));
    drift_mv_vect.emplace_back(std::pair<double, double>(2.00, -0.0007));
    drift_mv_vect.emplace_back(std::pair<double, double>(3.00, -0.0010));
    drift_mv_vect.emplace_back(std::pair<double, double>(5.00, -0.0010));

    //PiecewiseLinearCurve drift_spot(drift_mv_vect);
    drift = PiecewiseLinearCurve(drift_mv_vect);


    // lambda vol of vol
    // mat: 1w 2w 1m 3m 6m 9m 1y 2y 3y 5y 
    // 1w 0.0192    1.000
    // 2w 0.0384    0.850
    // 1m 0.0833    0.750
    // 3m 0.25      0.450
    // 6m 0.5       0.340
    // 9m 0.75      0.290
    // 1y 1.00      0.260
    // 2y 2.00      0.230
    // 3y 3.00      0.210
    // 5y 5.00      0.180

    std::vector<std::pair<double, double>> lambda_mv_vect;
    lambda_mv_vect.reserve(10);
    lambda_mv_vect.emplace_back(std::pair<double, double>(0.0192, 1.000));
    lambda_mv_vect.emplace_back(std::pair<double, double>(0.0384, 0.850));
    lambda_mv_vect.emplace_back(std::pair<double, double>(0.0833, 0.750));
    lambda_mv_vect.emplace_back(std::pair<double, double>(0.2500, 0.450));
    lambda_mv_vect.emplace_back(std::pair<double, double>(0.5000, 0.340));
    lambda_mv_vect.emplace_back(std::pair<double, double>(0.7500, 0.290));
    lambda_mv_vect.emplace_back(std::pair<double, double>(1.0000, 0.260));
    lambda_mv_vect.emplace_back(std::pair<double, double>(2.0000, 0.230));
    lambda_mv_vect.emplace_back(std::pair<double, double>(3.0000, 0.210));
    lambda_mv_vect.emplace_back(std::pair<double, double>(5.0000, 0.180));

    lambda = PiecewiseFlatCurve(lambda_mv_vect);

    // correl spot/vol (model)
    // 1d 0.02      -0.000 
    // 3m 0.25      -0.150
    // 1y 1.00      -0.200
    // 2y 2.00      -0.250
    // 3y 3.00      -0.250
    // 5y 5.00      -0.300

    std::vector<std::pair<double, double>> correl_mv_vect;

    correl_mv_vect.reserve(6);
    correl_mv_vect.emplace_back(std::pair<double, double>(0.02, -0.000));
    correl_mv_vect.emplace_back(std::pair<double, double>(0.25, -0.200));
    correl_mv_vect.emplace_back(std::pair<double, double>(1.00, -0.300));
    correl_mv_vect.emplace_back(std::pair<double, double>(2.00, -0.350));
    correl_mv_vect.emplace_back(std::pair<double, double>(3.00, -0.350));
    correl_mv_vect.emplace_back(std::pair<double, double>(5.00, -0.400));

    correl = PiecewiseLinearCurve(correl_mv_vect);

    // 3months is the characteristic time 
    kappa = 1. / 0.25;

    theta = PiecewiseFlatCurve(calibrateTheta(lambda, kappa));

}
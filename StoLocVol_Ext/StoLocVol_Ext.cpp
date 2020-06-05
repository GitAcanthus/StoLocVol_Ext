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
// SOFTWARE.// the MIT License (MIT)



#include < cmath >
#include < chrono > 
#include < iostream >
#include < vector >
#include < array >

#include "MQRModel.hpp"
#include "AdHocModels.hpp"
#include "LeverageSolver.hpp"



// scoped based basic timer
struct BasicTimer {
    using Tp = std::chrono::time_point<std::chrono::high_resolution_clock>;

    void stopTime() {
        auto end_point = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double, std::milli> dur = (end_point - m_startpoint);
        double milli = dur.count();
        std::cout << milli << "ms\n";
    }
    BasicTimer() { m_startpoint = std::chrono::high_resolution_clock::now(); }
    ~BasicTimer() { stopTime(); }

    Tp m_startpoint;
};

int main()
{

    PiecewiseLinearCurve drift;
    PiecewiseFlatCurve lambda;
    PiecewiseFlatCurve theta;
    PiecewiseLinearCurve correl;
    double kappa;

    initAdHocMQRModel(drift, lambda, theta, correl, kappa);
    ImpliedVolSmile<6> adhoc_smile = adHocVolatility_new();

    DumpPolicy filedump_pol;
    // local vol and leverage comparison
    filedump_pol.m_dumptype = DumpPolicy::dump_type::leverage;
    filedump_pol.m_filename = "loc_and_lev_file.txt";

    // four probability density dump
    //filedump_pol.m_dumptype = DumpPolicy::dump_type::proba;
    //filedump_pol.m_filename = "proba_file.txt";
    //filedump_pol.m_dumpdates.reserve(4);
    //filedump_pol.m_dumpdates.push_back(0.01);
    //filedump_pol.m_dumpdates.push_back(0.50);
    //filedump_pol.m_dumpdates.push_back(1.50);
    //filedump_pol.m_dumpdates.push_back(4.50);


    {
        // basic timer
        BasicTimer my_timer;
        
        // computational domain geometry:  
        // grid spot 181 pts
        // grid vol 201 pts
        // maturity 60 months (187 time steps)

        // open bounds
        LeverageSolver<181, 201, 60, OpenBoundConv, DOScheme> lev_calibration(drift, lambda, theta, correl, kappa);
        //LeverageSolver<181, 201, 60, OpenBoundDiff, DOScheme> lev_calibration(drift, lambda, theta, correl, kappa);
        
        
        // no-flux 
        // only available on request @ www.acanthussol.com 
        //LeverageSolver<181, 201, 60, ZeroFlowCond, DOScheme> lev_calibration(drift, lambda, theta, correl, kappa);

        lev_calibration.calibrateLeverage<6>(adhoc_smile);
    }


    std::cout << "thats it! \n";

}


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
#include < tuple > 


//~~~~~~~~~~~~~~~~~~~~~//
//~ Operator Matrices ~//
//~~~~~~~~~~~~~~~~~~~~~//

template< std::size_t N >
struct TridiagMat
{
    std::array<double, N> m_a;
    std::array<double, N> m_b;
    std::array<double, N> m_c;
};

template<std::size_t M, std::size_t P>
struct P0Coeff
{
    double m_coef_mult;

    // inf z-bound
    std::array<double, P - 2> m_zmin_yinf_bound;
    std::array<double, P - 2> m_zmin_yinf_inner;
    std::array<double, P - 2> m_zmin_ysup_bound;
    std::array<double, P - 2> m_zmin_ysup_inner;
    // sup z-bound
    std::array<double, P - 2> m_zmax_yinf_bound;
    std::array<double, P - 2> m_zmax_yinf_inner;
    std::array<double, P - 2> m_zmax_ysup_bound;
    std::array<double, P - 2> m_zmax_ysup_inner;
    // inf y-bound
    std::array<double, M - 2> m_ymin_zinf_bound;
    std::array<double, M - 2> m_ymin_zinf_inner;
    std::array<double, M - 2> m_ymin_zsup_bound;
    std::array<double, M - 2> m_ymin_zsup_inner;
    // sup y-bound
    std::array<double, M - 2> m_ymax_zinf_bound;
    std::array<double, M - 2> m_ymax_zinf_inner;
    std::array<double, M - 2> m_ymax_zsup_bound;
    std::array<double, M - 2> m_ymax_zsup_inner;

    // corners
    std::tuple<double, double, double> m_zmin_ymin_corner;
    std::tuple<double, double, double> m_zmin_ymax_corner;
    std::tuple<double, double, double> m_zmax_ymin_corner;
    std::tuple<double, double, double> m_zmax_ymax_corner;
};



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


#pragma once

#include < array >

//~~~~~~~~~~~~~~~~~~~//
//~ SpatialGeometry ~//
//~~~~~~~~~~~~~~~~~~~//

// spatial geometry defined at compile time
template< std::size_t M, std::size_t P >
struct SpatialGeometry
{
    const std::array<double, M - 1> m_gridsteps_z = initGridSteps<M>(0.85, 0.03);
    const std::array<double, M>     m_grid_z = initGridPoints<M>(m_gridsteps_z);
    const std::array<double, P - 1> m_gridsteps_y = initGridSteps<P>(1.200, 0.035);
    const std::array<double, P>     m_grid_y = initGridPoints<P>(m_gridsteps_y);
    const std::array<double, P>     m_expy = initExpyVect(m_grid_y);
    const std::array<double, P>     m_exp2y = initExp2yVect(m_grid_y);
};


//~~~~~~~~~~~~~~~~~~~~//
//~ TemporalGeometry ~//
//~~~~~~~~~~~~~~~~~~~~//

// temporal geometry defined at compile time
template< std::size_t H , std::size_t T = timestepNumber<H>() >
struct TemporalGeometry
{
    const std::array<double, T>     m_grid_t = initTimeGridVect<H>();
    const std::array<double, T-1>   m_gridsteps_t = initTimeStepVect(m_grid_t);
};



template<std::size_t N>
constexpr std::array<double, N - 1> initGridSteps(const double Max, const double Mult)
{
    static_assert(N % 2 == 1, "Misalignment: Grid point nbr must be odd");
    static_assert(N > 15, "Misalignment: Grid point nbr must > 15");

    std::array<double, N - 1> ret{};

    ret[(N - 1) / 2 + 0] = 0.25;
    ret[(N - 1) / 2 - 1] = 0.25;

    ret[(N - 1) / 2 + 1] = 0.25;
    ret[(N - 1) / 2 - 2] = 0.25;

    ret[(N - 1) / 2 + 2] = 0.25;
    ret[(N - 1) / 2 - 3] = 0.25;

    ret[(N - 1) / 2 + 3] = 0.25;
    ret[(N - 1) / 2 - 4] = 0.25;

    double sum{ 1. };

    for (std::size_t i = 4; i < (N - 1) / 2; ++i) {
        ret[(N - 1) / 2 + i] = ret[(N - 1) / 2 + i - 1] * (1. + Mult);
        ret[(N - 1) / 2 - i - 1] = ret[(N - 1) / 2 + i];
        sum += ret[(N - 1) / 2 + i];
    }
    for (std::size_t i = 0; i < N - 1; ++i) {
        ret[i] = ret[i] * Max / sum;
    }

    return ret;
}

template<std::size_t N>
constexpr std::array<double, N> initGridPoints(const std::array<double, N - 1> Steps)
{
    static_assert(N % 2 == 1, "Misalignment: Grid point nbr must be odd");
    static_assert(N > 15, "Misalignment: Grid point nbr must > 15");

    std::array<double, N> ret{};
    ret[(N - 1) / 2] = 0.;

    for (std::size_t i = 0; i < (N - 1) / 2; ++i) {
        ret[(N - 1) / 2 + i + 1] = ret[(N - 1) / 2 + i] + Steps[(N - 1) / 2 + i];
        ret[(N - 1) / 2 - i - 1] = -ret[(N - 1) / 2 + i + 1];
    }

    return ret;
}

template<std::size_t P>
constexpr std::array<double, P> initExpyVect(const std::array<double, P>& grid_y)
{
    std::array<double, P> exp_y;
    auto exp_y_it = exp_y.begin();

    std::for_each(grid_y.cbegin(), grid_y.cend(), [&](double y) {
        *exp_y_it = std::exp(y);
        ++exp_y_it; });

    return exp_y;
}

template<std::size_t P>
constexpr std::array<double, P> initExp2yVect(const std::array<double, P>& grid_y)
{
    std::array<double, P> exp_2y;
    auto exp_2y_it = exp_2y.begin();

    std::for_each(grid_y.cbegin(), grid_y.cend(), [&](double y) {
        *exp_2y_it = std::exp(2 * y);
        ++exp_2y_it; });

    return exp_2y;
}


//initTimeGrid
template< unsigned short M >
constexpr std::size_t timestepNumber()
{
    static_assert(M > 0, "Month number must be >0 !");
    std::size_t ret{ 1 };// time 0
    unsigned short mat_min{ 0 };

    //first month: 5x 1d + 5x 2d + 5x 3d 
    ret += 15;

    // months 2 to 3: 6x 5d    
    mat_min = M < 3 ? M : 3;
    for (unsigned short i = 1; i < mat_min; ++i)
        ret += 6;

    // months 4 to 6: 5x 6d    
    mat_min = M < 6 ? M : 6;
    for (unsigned short i = 3; i < mat_min; ++i)
        ret += 5;

    // months 7 to 12: 4x 7d    
    mat_min = M < 12 ? M : 12;
    for (unsigned short i = 6; i < mat_min; ++i)
        ret += 4;

    // months 13 to 36: 3x 10d    
    mat_min = M < 36 ? M : 36;
    for (unsigned short i = 12; i < mat_min; ++i)
        ret += 3;

    // months 37 and beyond 2x 15d
    for (unsigned short i = 36; i < M; ++i)
        ret += 2;

    return ret;
}


template< unsigned short M>
constexpr std::array<double, timestepNumber<M>()> initTimeGridVect()
{
    static_assert(M > 0, "Month number must be >0 !");

    //std::array<double, Nbr> ret{} ;
    std::array<double, timestepNumber<M>()> ret{};
    unsigned short mat_min;
    std::size_t last_indx;

    ret[0] = 0.;

    //first month: 5x 1d + 5x 2d + 5x 3d 
    ret[1] = 1. * 0.002777777777778;
    ret[2] = 2. * 0.002777777777778;
    ret[3] = 3. * 0.002777777777778;
    ret[4] = 4. * 0.002777777777778;
    ret[5] = 5. * 0.002777777777778;

    ret[6] = ret[5] + 0.005555555555556;
    ret[7] = ret[6] + 0.005555555555556;
    ret[8] = ret[7] + 0.005555555555556;
    ret[9] = ret[8] + 0.005555555555556;
    ret[10] = ret[9] + 0.005555555555556;

    ret[11] = ret[10] + 0.008333333333333;
    ret[12] = ret[11] + 0.008333333333333;
    ret[13] = ret[12] + 0.008333333333333;
    ret[14] = ret[13] + 0.008333333333333;
    ret[15] = ret[14] + 0.008333333333333;
    last_indx = 15;

    // months 2 to 3: 6x 5d    
    mat_min = M < 3 ? M : 3;
    for (unsigned short i = 1; i < mat_min; ++i) {
        ret[last_indx + 1] = ret[last_indx] + 0.013888888888889;
        ret[last_indx + 2] = ret[last_indx + 1] + 0.013888888888889;
        ret[last_indx + 3] = ret[last_indx + 2] + 0.013888888888889;
        ret[last_indx + 4] = ret[last_indx + 3] + 0.013888888888889;
        ret[last_indx + 5] = ret[last_indx + 4] + 0.013888888888889;
        ret[last_indx + 6] = ret[last_indx + 5] + 0.013888888888889;
        last_indx += 6;
    }

    // months 3 to 6: 5x 6d    
    mat_min = M < 6 ? M : 6;
    for (unsigned short i = 3; i < mat_min; ++i) {
        ret[last_indx + 1] = ret[last_indx] + 0.016666666666667;
        ret[last_indx + 2] = ret[last_indx + 1] + 0.016666666666667;
        ret[last_indx + 3] = ret[last_indx + 2] + 0.016666666666667;
        ret[last_indx + 4] = ret[last_indx + 3] + 0.016666666666667;
        ret[last_indx + 5] = ret[last_indx + 4] + 0.016666666666667;
        last_indx += 5;
    }

    // months 7 to 12: 4x 7d    
    mat_min = M < 12 ? M : 12;
    for (unsigned short i = 6; i < mat_min; ++i) {
        ret[last_indx + 1] = ret[last_indx] + 0.020833333333333;
        ret[last_indx + 2] = ret[last_indx + 1] + 0.020833333333333;
        ret[last_indx + 3] = ret[last_indx + 2] + 0.020833333333333;
        ret[last_indx + 4] = ret[last_indx + 3] + 0.020833333333333;
        last_indx += 4;
    }

    // months 13 to 36: 3x 10d    
    mat_min = M < 36 ? M : 36;
    for (unsigned short i = 12; i < mat_min; ++i) {
        ret[last_indx + 1] = ret[last_indx] + 0.027777777777778;
        ret[last_indx + 2] = ret[last_indx + 1] + 0.027777777777778;
        ret[last_indx + 3] = ret[last_indx + 2] + 0.027777777777778;
        last_indx += 3;
    }

    // months 37 and beyond 2x 15d
    for (unsigned short i = 36; i < M; ++i) {
        ret[last_indx + 1] = ret[last_indx] + 0.041666666666667;
        ret[last_indx + 2] = ret[last_indx + 1] + 0.041666666666667;
        last_indx += 2;
    }

    return ret;
}


template<std::size_t tnbr>
constexpr std::array<double, (tnbr - 1)> initTimeStepVect(const std::array<double, tnbr>& tgrid)
{
    std::array<double, (tnbr - 1)> mesh_stp;
    auto grid_it = tgrid.cbegin();
    auto grid_nit = std::next(grid_it);

    for (auto& msh_stp : mesh_stp) {
        msh_stp = *grid_nit - *grid_it;
        ++grid_nit, ++grid_it;
    }
    return mesh_stp;
}


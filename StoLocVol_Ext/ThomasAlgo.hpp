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

// thomas algo for std::vector 
// and std::array containers 

template <typename T>
short thomas_inv(std::size_t size,
    const T& a,  // diag a
    const T& b,  // diag b
    const T& c,  // diag c
    const T& r,  // right term
    T& u,        // solution
    T& w)        // workspace
{
    short status{0};
    
    double bet;
    auto it_a = a.cbegin();
    auto it_b = b.cbegin();
    auto it_c = c.cbegin();
    auto it_r = r.cbegin();
    auto it_u = u.begin();
    auto it_p = it_u;
    auto it_w = w.begin();

    bet = *it_b;
    *it_u = *it_r / bet;


    // forward substitution 
    ++it_w, ++it_a, ++it_b, ++it_r, ++it_u;
    for (; it_a != a.cend();) {

#ifdef AS_DEBUG
        if (bet == 0.)
            status = 1;
#endif

        *it_w = *it_c / bet;
        bet = *it_b - *it_a * *it_w;
        *it_u = (*it_r - *it_a * *(it_p)) / bet;

        ++it_a, ++it_b, ++it_c, ++it_r, ++it_u, ++it_p, ++it_w;
    }
    // backward substitution
    auto it_ru = u.rbegin();
    it_w = std::prev(it_w);

    for (++it_ru; it_ru != u.rend();) {
        *it_ru -= *it_w * *it_p;
        ++it_ru, --it_w, --it_p;
    }

    return status ;
}







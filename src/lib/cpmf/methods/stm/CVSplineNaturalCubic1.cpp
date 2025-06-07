// ===============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -------------------------------------------------------------------------------
//    Copyright (C) 2025 Petr Kulhanek, kulhanek@chemi.muni.cz
//
//     This program is free software; you can redistribute it and/or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation; either version 2 of the License, or
//     (at your option) any later version.
//
//     This program is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
//
//     You should have received a copy of the GNU General Public License along
//     with this program; if not, write to the Free Software Foundation, Inc.,
//     51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
// ===============================================================================

#include <CVSplineNaturalCubic.hpp>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CCVSplineNaturalCubic::CCVSplineNaturalCubic(void)
{
}

//------------------------------------------------------------------------------

CCVSplineNaturalCubic::~CCVSplineNaturalCubic(void)
{
    Clear();
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CCVSplineNaturalCubic::Clear(void)
{
    X.clear();
    Y.clear();
    a.clear();
    b.clear();
    c.clear();
    d.clear();
    NumOfKnots = 0;
}

//------------------------------------------------------------------------------

void CCVSplineNaturalCubic::Allocate(int numofknots)
{
    Clear();

    NumOfKnots = numofknots;
    if( NumOfKnots <= 0 ){
        NumOfKnots = 0;
        return;
    }

    X.resize(NumOfKnots);
    Y.resize(NumOfKnots);

    a.resize(NumOfKnots);
    std::fill(a.begin(), a.end(), 0);

    b.resize(NumOfKnots);
    std::fill(b.begin(), b.end(), 0);

    c.resize(NumOfKnots);
    std::fill(c.begin(), c.end(), 0);

    d.resize(NumOfKnots);
    std::fill(d.begin(), d.end(), 0);
}

//------------------------------------------------------------------------------

bool CCVSplineNaturalCubic::AddPoint(int knotid,double alpha,double cv)
{
    if( (knotid < 0) || (knotid >= NumOfKnots)) {
        return(false);
    }

    X[knotid] = alpha;
    Y[knotid] = cv;

    return(true);
}

//------------------------------------------------------------------------------

bool CCVSplineNaturalCubic::Finalize(void)
{
    if (NumOfKnots <= 1)  return(false);    // not enough data

    if (NumOfKnots == 2) {
        // Two-point case: use linear interpolation
        double dx = X[1] - X[0];
        a[0] = Y[0];
        b[0] = (Y[1] - Y[0]) / dx;
        c[0] = 0.0;
        d[0] = 0.0;
        return(true);
    }

    int n = NumOfKnots;
    std::vector<double> h(n - 1), alpha(n - 1);

    // Step 1: Compute h[i]
    for (int i = 0; i < n - 1; ++i) {
        h[i] = X[i + 1] - X[i];
    }

    // Step 2: Compute alpha[i]
    for (int i = 1; i < n - 1; ++i) {
        alpha[i] = (3.0 / h[i]) * (Y[i + 1] - Y[i]) - (3.0 / h[i - 1]) * (Y[i] - Y[i - 1]);
    }

    // Step 3: Solve tridiagonal system for c[i]
    std::vector<double> l(n), mu(n), z(n);
    l[0] = 1.0;
    mu[0] = z[0] = 0.0;

    for (int i = 1; i < n - 1; ++i) {
        l[i] = 2.0 * (X[i + 1] - X[i - 1]) - h[i - 1] * mu[i - 1];
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }

    l[n - 1] = 1.0;
    z[n - 1] = c[n - 1] = 0.0;

    // Step 4: Back substitution for c[i], b[i], d[i]
    for (int j = n - 2; j >= 0; --j) {
        c[j] = z[j] - mu[j] * c[j + 1];
        b[j] = (Y[j + 1] - Y[j]) / h[j] - h[j] * (2.0 * c[j] + c[j + 1]) / 3.0;
        d[j] = (c[j + 1] - c[j]) / (3.0 * h[j]);
        a[j] = Y[j];
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

double CCVSplineNaturalCubic::GetCV(double alpha)
{
    int n = NumOfKnots;

    // Handle out-of-bounds queries by clamping to the closest interval
    if (alpha <= X[0]) {
        alpha = X[0];
    } else if (alpha >= X[n - 1]) {
        alpha = X[n - 1];
    }

    // Find the right interval [X[i], X[i+1]] using linear search
    int i = 0;
    while (i < n - 2 && alpha > X[i + 1]) {
        ++i;
    }

    double dx = alpha - X[i];
    return a[i] + b[i] * dx + c[i] * dx * dx + d[i] * dx * dx * dx;
}

//------------------------------------------------------------------------------

double CCVSplineNaturalCubic::GetCVFirstDer(double alpha)
{
    int n = NumOfKnots;

    // Handle out-of-bounds queries by clamping to the closest interval
    if (alpha <= X[0]) {
        alpha = X[0];
    } else if (alpha >= X[n - 1]) {
        alpha = X[n - 1];
    }

    // Find the right interval [X[i], X[i+1]] using linear search
    int i = 0;
    while (i < n - 2 && alpha > X[i + 1]) {
        ++i;
    }

    double dx = alpha - X[i];
    return b[i] + 2.0 * c[i] * dx + 3.0 * d[i] * dx * dx;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================


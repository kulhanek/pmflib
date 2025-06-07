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
    x.FreeVector();
    y.FreeVector();
    sa.FreeVector();
    sb.FreeVector();
    sc.FreeVector();
    sd.FreeVector();
    n = -1;
}

//------------------------------------------------------------------------------

void CCVSplineNaturalCubic::Allocate(int numofknots)
{
    Clear();

    n = numofknots-1;
    if( n <= 0 ){
        n = 0;
        return;
    }

    x.CreateVector(n+1); // 0,1,...,n
    x.SetZero();

    y.CreateVector(n+1); // 0,1,...,n
    y.SetZero();

    sa.CreateVector(n+1); // 0,1,...,n
    sa.SetZero();

    sb.CreateVector(n+1); // 0,1,...,n
    sb.SetZero();

    sc.CreateVector(n+1); // 0,1,...,n
    sc.SetZero();

    sd.CreateVector(n+1); // 0,1,...,n
    sd.SetZero();
}

//------------------------------------------------------------------------------

bool CCVSplineNaturalCubic::AddPoint(int knotid,double alpha,double cv)
{
    if( (knotid < 0) || (knotid > n)) {
        return(false);
    }

    x[knotid] = alpha;
    y[knotid] = cv;

    return(true);
}

//------------------------------------------------------------------------------

bool CCVSplineNaturalCubic::Finalize(void)
{
    if (n <= 0)  return(false);    // not enough data

    if (n == 1) {
        // linear interpolation - between two points
        sd[0] = y[0];
        sc[0] = (y[1] - y[0]) / (x[1] - x[0]);
        sb[0] = 0.0;
        sa[0] = 0.0;
        return(true);
    }

// helpers
    CSimpleVector<double>   h, p, q, b;

    h.CreateVector(n+1); // 0,1,...,n
    h.SetZero();

    p.CreateVector(n+1);
    p.SetZero();

    q.CreateVector(n+1);
    q.SetZero();

    b.CreateVector(n+1);
    b.SetZero();

// --------------------
    h[0] = x[1] - x[0];
    for(int i=1; i <= n-1; i++){
        h[i] = x[i+1] - x[i];
        p[i] = 2.0*(x[i+1]-x[i-1]);
        q[i] = 3.0*(y[i+1]-y[i])/h[i] - 3.0*(y[i]-y[i-1])/h[i-1];
    }

// --------------------
    for(int i=2; i <= n-1; i++){
        p[i] = p[i] - h[i-1]*h[i-1]/p[i-1];
        q[i] = q[i] - q[i-1]*h[i-1]/p[i-1];
    }

// --------------------
    b[n-1] = q[n-1]/p[n-1];
    for(int i=2; i <= n-1; i++){
        b[n-i] = (q[n-i]-h[n-i]*b[n-i+1])/p[n-i];
    }

// --------------------
    sa[0] = b[1]/(3.0*h[0]);
    sb[0] = 0.0;
    sc[0] = (y[1]-y[0])/h[0] - b[1]*h[0]/3.0;
    sd[0] = y[0];

    for(int i=1; i <= n-1; i++){
        sa[i] = (b[i+1]-b[i])/(3.0*h[i]);
        sb[i] = b[i];
        sc[i] = (b[i]+b[i-1])*h[i-1] + sc[i-1];
        sd[i] = y[i];
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

double CCVSplineNaturalCubic::GetCV(double alpha)
{
    // Handle out-of-bounds queries by clamping to the closest interval
    if( alpha <= x[0] ) {
        alpha = x[0];
    } else if( alpha >= x[n] ) {
        alpha = x[n];
    }

    // Find the right interval [x[i], x[i+1]] using linear search
    int i = 0;
    while( (i <= n - 2) && (alpha >= x[i + 1]) ) {
        i++;
    }

    double dx = alpha - x[i];
    double sp = sd[i] + sc[i] * dx + sb[i] * dx * dx + sa[i] * dx * dx * dx;

    return(sp);
}

//------------------------------------------------------------------------------

double CCVSplineNaturalCubic::GetCVFirstDer(double alpha)
{
    // Handle out-of-bounds queries by clamping to the closest interval
    if( alpha <= x[0] ) {
        alpha = x[0];
    } else if( alpha >= x[n] ) {
        alpha = x[n];
    }

    // Find the right interval [x[i], x[i+1]] using linear search
    int i = 0;
    while( (i <= n - 2) && (alpha >= x[i + 1]) ) {
        i++;
    }

    double dx = alpha - x[i];
    double de = sc[i] + 2.0 * sb[i] * dx + 3.0 * sa[i] * dx * dx;

    return(de);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================


// ===============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -------------------------------------------------------------------------------
//    Copyright (C) 2012 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2011 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2010 Petr Kulhanek, kulhanek@chemi.muni.cz
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

#include <CVSpline.hpp>
#include <ErrorSystem.hpp>
#include <math.h>
#include <spline.hpp>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CCVSpline::CCVSpline(void)
{
    t = NULL;
    y = NULL;
    ypp = NULL;
}

//------------------------------------------------------------------------------

CCVSpline::~CCVSpline(void)
{
    Clear();
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CCVSpline::Clear(void)
{
    if( t != NULL ) delete[] t;
    if( y != NULL ) delete[] y;
    if( ypp != NULL ) delete[] ypp;
    t = NULL;
    y = NULL;
    ypp = NULL;
    NumOfKnots = 0;
}

//------------------------------------------------------------------------------

void CCVSpline::Allocate(int numofknots)
{
    Clear();
    NumOfKnots = numofknots;
    if( NumOfKnots <= 0 ) return;

    t = new double[NumOfKnots];
    y = new double[NumOfKnots];
    ypp = NULL;
}

//------------------------------------------------------------------------------

void CCVSpline::AddPoint(int knotid,double alpha,double cv)
{
    if( (knotid < 0) || (knotid >= NumOfKnots)) {
        RUNTIME_ERROR("knotid is out-of-legal range");
    }
    t[knotid] = alpha;
    y[knotid] = cv;
}

//------------------------------------------------------------------------------

void CCVSpline::Finalize(void)
{
    ypp = spline_cubic_set(NumOfKnots,t,y,0,0,0,0);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

double CCVSpline::GetCV(double alpha)
{
    if( ypp == NULL ) {
        RUNTIME_ERROR("spline is not set");
    }
    double fder,sder;
    double value = spline_cubic_val(NumOfKnots,t,y,ypp,alpha,&fder,&sder);
    return(value);
}

//------------------------------------------------------------------------------

double CCVSpline::GetCVFirstDer(double alpha)
{
    if( ypp == NULL ) {
        RUNTIME_ERROR("spline is not set");
    }
    double fder,sder;
    spline_cubic_val(NumOfKnots,t,y,ypp,alpha,&fder,&sder);
    return(fder);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================


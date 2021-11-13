// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
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
// =============================================================================

#include <Filter.hpp>
//
//! mass * velocity^2 -> kcal/mol
//! g/mol * A^2 / fs^2 = 0.001 kg/mol * (10^-10)^2 m^2 / ((10^-15)^s^2)
//! 0.001*10^-20/10^-30 * kg/mol*m^2/s^s = 10^-3*10^10 J/mol = 10^7 J/mol
//! 10^7 J/mol = 10^7 / 4184 kcal/mol
//! dt (fs) = vdt * sqrt( 4184 / 10^7) = 0.02045482828087295384
//real(PMFDP), parameter  :: PMF_DT2VDT   = 0.02045482828087295384d0 ! sqrt(pc_cal/1d4)
//real(PMFDP), parameter  :: PMF_VDT2DT   = 1.0d0 / PMF_DT2VDT

const double PMF_DT2VDT   = 0.02045482828087295384;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CFilter::CFilter(void)
{
    TimeStep = 0.0;
    SamplingFreq = 0.0;
    InvFDTX = 0.0;
}

//------------------------------------------------------------------------------

CFilter::~CFilter(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CFilter::SetTimeStep(double timestep)
{
    TimeStep = timestep;

    if( TimeStep <= 0.0 ){
        RUNTIME_ERROR("time step must be grater than zero");
    }

    // set sampling frequencies
    SamplingFreq = 1e6/(29.9792458*TimeStep);

    // FIXME
    InvFDTX = 1.0 / (TimeStep * PMF_DT2VDT);
}

//------------------------------------------------------------------------------

void CFilter::RunFilter(CVectorDataPtr& inout)
{
    size_t len = inout->GetLength();
    if( len == 0 ) return;

    CVectorDataPtr in(new CSimpleVector<double>);
    in->CreateVector(len);
    for(size_t idx=0; idx < len; idx++){
        in->GetRawDataField()[idx] = inout->GetRawDataField()[idx];
    }
    RunFilter(in,inout);
}

//------------------------------------------------------------------------------

void CFilter::RunFilter(const CVectorDataPtr& in, CVectorDataPtr& out)
{
    RUNTIME_ERROR("need to be redefined");
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

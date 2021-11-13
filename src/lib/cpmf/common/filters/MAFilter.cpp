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

#include <MAFilter.hpp>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CMAFilter::CMAFilter(void)
{
    CutoffFreq = 0.0;
    NormFreq = 0.0;
    N = 0;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CMAFilter::SetFilter(double timestep, double fcut)
{
    SetTimeStep(timestep);

    CutoffFreq = fcut;

    if( (CutoffFreq <= 0.0) && (CutoffFreq >= SamplingFreq) ){
        RUNTIME_ERROR("illegal cutoff frequency");
    }

    NormFreq = CutoffFreq / SamplingFreq;

    // setup the length
    // https://dsp.stackexchange.com/questions/9966/what-is-the-cut-off-frequency-of-a-moving-average-filter

    N = sqrt((0.442947/NormFreq)*(0.442947/NormFreq) + 1.0);

    if( N < 2 ){
        RUNTIME_ERROR("illegal value of N");
    }
}

//------------------------------------------------------------------------------

int CMAFilter::GetFrameLen(void) const
{
    return(N);
}

//------------------------------------------------------------------------------

void CMAFilter::RunFilter(const CVectorDataPtr& in, CVectorDataPtr& out)
{
    if( in->GetLength() != out->GetLength() ){
        RUNTIME_ERROR("sizes of in and out arrays differ");
    }
    if( in->GetLength() <= N ){
        RUNTIME_ERROR("size of in array is too small");
    }

    double w = 1.0/N;
    double prev = 0.0;
    for(size_t t=0; t < N; t++){
        prev = prev + in->GetRawDataField()[t]*w;
        out->GetRawDataField()[t] = 0.0;
    }
    out->GetRawDataField()[N-1] = prev;
    for(size_t t=N; t < in->GetLength(); t++){
        prev = prev + (in->GetRawDataField()[t] - in->GetRawDataField()[t-N])*w;
        out->GetRawDataField()[t] = prev;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2008 Petr Kulhanek, kulhanek@enzim.hu
//                       Martin Petrek, petrek@chemi.muni.cz
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

#include <errno.h>
#include <string.h>
#include <ABFAccumulator.hpp>
#include <ErrorSystem.hpp>
#include <XMLElement.hpp>
#include <XMLBinData.hpp>

//------------------------------------------------------------------------------

using namespace std;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CABFAccumulator::CABFAccumulator(void)
{
    NCorr           = 1.0;
}

//------------------------------------------------------------------------------

CABFAccumulator::~CABFAccumulator(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CABFAccumulator::GetNumOfSamples(const CSimpleVector<int>& position) const
{
    int glbindex = 0;
    for(int i=0; i < NumOfCVs; i++) {
        glbindex = glbindex*CVs[i]->GetNumOfBins() + position[i];
    }
    return(GetNumOfSamples(glbindex));
}

//------------------------------------------------------------------------------

int CABFAccumulator::GetNumOfSamples(int ibin) const
{
    return(GetData("NSAMPLES",ibin));
}

//------------------------------------------------------------------------------

double CABFAccumulator::GetValue(int icv,int ibin,EABFAccuValue realm) const
{
    int     nsamples = GetData("NSAMPLES",ibin);
    double  micf     = GetData("MICF",ibin,icv);
    double  m2icf    = GetData("M2ICF",ibin,icv);

    double value = 0.0;
    if( nsamples <= 0 ) return(value);

    switch(realm){
// mean force
        // -------------------
        case(EABF_DG_VALUE):
            return( micf );
        // -------------------
        case(EABF_DG_SIGMA):
            return( sqrt(m2icf / nsamples) );
        // -------------------
        case(EABF_DG_ERROR):
            return( sqrt(NCorr * m2icf) / nsamples );
        // -------------------
        // TDS - FIXME
        default:
            RUNTIME_ERROR("unsupported realm");
    }

    return(value);
}

//------------------------------------------------------------------------------

double CABFAccumulator::GetValue(int ibin,EABFAccuValue realm) const
{
    int     nsamples = GetData("NSAMPLES",ibin);
    double  mepot    = GetData("MEPOT",ibin);
    double  m2epot   = GetData("M2EPOT",ibin);

    double value = 0.0;
    if( nsamples <= 0 ) return(value);

    switch(realm){
        // -------------------
        case(EABF_H_VALUE):
            return( mepot );
        // -------------------
        case(EABF_H_SIGMA):
            return( sqrt(m2epot / nsamples) );
        // -------------------
        case(EABF_H_ERROR):
            return( sqrt(NCorr * m2epot) / nsamples );
        // -------------------
        default:
            RUNTIME_ERROR("unsupported realm");
    }

    return(value);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABFAccumulator::SetNumOfSamples(int ibin,int nsamples)
{
    SetData("NSAMPLES",ibin,nsamples);
}

//------------------------------------------------------------------------------

void CABFAccumulator::SetMaskWeight(int ibin,double weight)
{
    SetData("WEIGTH",ibin,weight);
}

//------------------------------------------------------------------------------

void CABFAccumulator::SetNCorr(double ncorr)
{
    NCorr = ncorr;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================




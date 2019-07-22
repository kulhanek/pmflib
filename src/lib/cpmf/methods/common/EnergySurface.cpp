// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2011 Petr Kulhanek, kulhanek@chemi.muni.cz
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

#include <math.h>
#include <errno.h>
#include <string.h>
#include <EnergySurface.hpp>
#include <ErrorSystem.hpp>
#include <XMLElement.hpp>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CEnergySurface::CEnergySurface(void)
{
    NumOfCVs = 0;
    TotNPoints = 0;
}

//------------------------------------------------------------------------------

CEnergySurface::~CEnergySurface(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CEnergySurface::GetNumberOfCoords(void) const
{
    return(NumOfCVs);
}

//------------------------------------------------------------------------------

int CEnergySurface::GetNumberOfPoints(void) const
{
    return(TotNPoints);
}

//------------------------------------------------------------------------------

const CColVariable* CEnergySurface::GetCoordinate(unsigned int cv) const
{
    return(&Sizes[cv]);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CEnergySurface::Allocate(const CMTDHistory* mtd_hist)
{
    if( NumOfCVs > 0 ) Deallocate();
    if( mtd_hist == NULL ) return;
    if( mtd_hist->GetNumberOfCoords() <= 0 ) return;

// allocate items
    NumOfCVs = mtd_hist->GetNumberOfCoords();
    Sizes.CreateVector(NumOfCVs);

// copy cvs and calculate total number of points
    TotNPoints = 1;
    for(int i=0; i < NumOfCVs; i++) {
        Sizes[i].CopyFrom(mtd_hist->GetCoordinate(i));
        TotNPoints *= Sizes[i].GetNumberOfBins();
    }

    Energy.CreateVector(TotNPoints);
    Error.CreateVector(TotNPoints);
    Samples.CreateVector(TotNPoints);

    Clear();
}

//------------------------------------------------------------------------------

void CEnergySurface::Allocate(const CABFAccumulator* abf_accu)
{
    if( NumOfCVs > 0 ) Deallocate();
    if( abf_accu == NULL ) return;
    if( abf_accu->GetNumberOfCoords() <= 0 ) return;

// allocate items
    NumOfCVs = abf_accu->GetNumberOfCoords();
    Sizes.CreateVector(NumOfCVs);

// copy cvs and calculate total number of points
    TotNPoints = 1;
    for(int i=0; i < NumOfCVs; i++) {
        Sizes[i].CopyFrom(abf_accu->GetCoordinate(i));
        TotNPoints *= Sizes[i].GetNumberOfBins();
    }

    Energy.CreateVector(TotNPoints);
    Error.CreateVector(TotNPoints);
    Samples.CreateVector(TotNPoints);

    Clear();
}

//------------------------------------------------------------------------------

void CEnergySurface::Deallocate(void)
{
    Sizes.FreeVector();
    Energy.FreeVector();
    Error.FreeVector();
    Samples.FreeVector();
    NumOfCVs = 0;
    TotNPoints = 0;
}

//------------------------------------------------------------------------------

void CEnergySurface::Clear(void)
{
    for(int i=0; i < TotNPoints; i++) {
        Energy[i] = 0.0;
        Error[i] = 0.0;
        Samples[i] = 0;
    }
}



//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CEnergySurface::SetEnergy(unsigned int index,const double& value)
{
    Energy[index] = value;
}

//------------------------------------------------------------------------------

const double& CEnergySurface::GetEnergy(unsigned int index) const
{
    return(Energy[index]);
}

//------------------------------------------------------------------------------

void CEnergySurface::SetError(unsigned int index,const double& value)
{
    Error[index] = value;
}

//------------------------------------------------------------------------------

const double& CEnergySurface::GetError(unsigned int index) const
{
    return(Error[index]);
}

//------------------------------------------------------------------------------

double CEnergySurface::GetSigmaF2(void) const
{
    double sigmafsum = 0.0;
    double sigmafsum2 = 0.0;
    double count = 0.0;

    for(int k=0; k < TotNPoints; k++) {
        if( Samples[k] > 0 ){
            double ene = Energy[k];
            sigmafsum += ene;
            sigmafsum2 += ene*ene;
            count++;
        }
    }

    double sigmaf2 = 0.0;

    if( count > 0 ){
        sigmaf2 = count*sigmafsum2 - sigmafsum*sigmafsum;
        if(sigmaf2 > 0) {
            sigmaf2 = sqrt(sigmaf2) / count;
        } else {
            sigmaf2 = 0.0;
        }
    }

    sigmaf2 *= sigmaf2;

    return(sigmaf2);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CEnergySurface::SetNumOfSamples(unsigned int index,const int& value)
{
    Samples[index] = value;
}

//------------------------------------------------------------------------------

const int& CEnergySurface::GetNumOfSamples(unsigned int index) const
{
    return(Samples[index]);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CEnergySurface::GetPoint(unsigned int index,CSimpleVector<double>& point) const
{
    for(int k=NumOfCVs-1; k >= 0; k--) {
        const CColVariable* p_coord = &Sizes[k];
        int ibin = index % p_coord->GetNumberOfBins();
        point[k] = p_coord->GetValue(ibin);
        index = index / p_coord->GetNumberOfBins();
    }
}

//------------------------------------------------------------------------------

void CEnergySurface::GetIPoint(unsigned int index,CSimpleVector<int>& point) const
{
    for(int k=NumOfCVs-1; k >= 0; k--) {
        const CColVariable* p_coord = &Sizes[k];
        int ibin = index % p_coord->GetNumberOfBins();
        point[k] = ibin;
        index = index / p_coord->GetNumberOfBins();
    }
}

//------------------------------------------------------------------------------

double CEnergySurface::GetGlobalMinimumValue(void) const
{
    double minimum = 0.0;

    if(TotNPoints > 0) minimum = Energy[0];

    for(int k=0; k < TotNPoints; k++) {
        if(minimum > Energy[k]) minimum = Energy[k];
    }

    return(minimum);
}

//------------------------------------------------------------------------------

void CEnergySurface::ApplyOffset(double offset)
{
    for(int k=0; k < TotNPoints; k++) {
        Energy[k] += offset;
    }
}

//------------------------------------------------------------------------------

void CEnergySurface::AdaptErrorsToGlobalMinimum(void)
{
    double minimum = 0.0;
    double err_at_minimum = 0.0;
    bool   first = true;

    for(int k=0; k < TotNPoints; k++) {
        if( Samples[k] <= 0 ) continue;
        if( (minimum > Energy[k]) || (first == true) ){
            minimum = Energy[k];
            err_at_minimum = Error[k];
            first = false;
        }
    }

    // adapt errors
    for(int k=0; k < TotNPoints; k++) {
        if( Samples[k] <= 0 ) continue;
        Error[k] = fabs(Error[k] - err_at_minimum);
    }
}

//------------------------------------------------------------------------------

void CEnergySurface::AdaptUnsampledToMaxEnergy(void)
{
    double maxe = 0.0;
    bool   first = true;

    // find maximum in sampled region
    for(int k=0; k < TotNPoints; k++) {
        if( Samples[k] <= 0 ) continue; // skip unsampled
        if( (maxe < Energy[k]) || (first == true) ){
            maxe = Energy[k];
            first = false;
        }
    }

    // adapt unsampled region
    for(int k=0; k < TotNPoints; k++) {
        if( Samples[k] <= 0 ){
            Energy[k] = maxe;
        }
    }
}

//------------------------------------------------------------------------------

void CEnergySurface::AdaptUnsampledToMaxEnergy(double maxene)
{
    // adapt unsampled region
    for(int k=0; k < TotNPoints; k++) {
        if( Samples[k] <= 0 ){
            Energy[k] = maxene;
        }
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CEnergySurface::operator+=(const CEnergySurface& source)
{
    if(TotNPoints != source.TotNPoints) {
        ES_ERROR("surfaceces do not match");
        return;
    }

    for(int k=0; k < TotNPoints; k++) {
        Energy[k] += source.Energy[k];
    }
}

//------------------------------------------------------------------------------

void CEnergySurface::operator/=(const double& number)
{
    for(int k=0; k < TotNPoints; k++) {
        Energy[k] /= number;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================



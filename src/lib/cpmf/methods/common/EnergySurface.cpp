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
    NumOfItems = 0;
    TotNPoints = 0;
}

//------------------------------------------------------------------------------

CEnergySurface::~CEnergySurface(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

unsigned int CEnergySurface::GetNumberOfCoords(void) const
{
    return(NumOfItems);
}

//------------------------------------------------------------------------------

unsigned int CEnergySurface::GetNumberOfPoints(void) const
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
    if( NumOfItems > 0 ) Deallocate();
    if( mtd_hist == NULL ) return;
    if( mtd_hist->GetNumberOfCoords() <= 0 ) return;

// allocate items
    NumOfItems = mtd_hist->GetNumberOfCoords();
    Sizes.CreateVector(NumOfItems);

// copy cvs and calculate total number of points
    TotNPoints = 1;
    for(unsigned int i=0; i < NumOfItems; i++) {
        Sizes[i].CopyFrom(mtd_hist->GetCoordinate(i));
        TotNPoints *= Sizes[i].GetNumberOfBins();
    }

    Energy.CreateVector(TotNPoints);
    Samples.CreateVector(TotNPoints);

    Clear();
}

//------------------------------------------------------------------------------

void CEnergySurface::Allocate(const CABFAccumulator* abf_accu)
{
    if( NumOfItems > 0 ) Deallocate();
    if( abf_accu == NULL ) return;
    if( abf_accu->GetNumberOfCoords() <= 0 ) return;

// allocate items
    NumOfItems = abf_accu->GetNumberOfCoords();
    Sizes.CreateVector(NumOfItems);

// copy cvs and calculate total number of points
    TotNPoints = 1;
    for(unsigned int i=0; i < NumOfItems; i++) {
        Sizes[i].CopyFrom(abf_accu->GetCoordinate(i));
        TotNPoints *= Sizes[i].GetNumberOfBins();
    }

    Energy.CreateVector(TotNPoints);
    Samples.CreateVector(TotNPoints);

    Clear();
}

//------------------------------------------------------------------------------

void CEnergySurface::Deallocate(void)
{
    Sizes.FreeVector();
    Energy.FreeVector();
    Samples.FreeVector();
    NumOfItems = 0;
    TotNPoints = 0;
}

//------------------------------------------------------------------------------

void CEnergySurface::Clear(void)
{
    for(unsigned int i=0; i < TotNPoints; i++) {
        Energy[i] = 0.0;
        Samples[i] = 0;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CEnergySurface::CalculateFES(CMTDHistory& mtd_hist,unsigned int mtd_time)
{
// quick compatibility comparison
    if( NumOfItems != mtd_hist.GetNumberOfCoords() ) {
        RUNTIME_ERROR("numofitems mismatch");
    }

// allocate point
    CSimpleVector<double> point;
    point.CreateVector(NumOfItems);

// calculate surface
    unsigned int loc = 0;
    CalculateFES_Part(mtd_hist,point,mtd_time,loc,0);
}

//------------------------------------------------------------------------------

void CEnergySurface::CalculateFES_Part(CMTDHistory& mtd_hist,
                                       CSimpleVector<double>& point,
                                       unsigned int mtd_time,
                                       unsigned int& loc,
                                       unsigned int cv)
{
    if(cv >= NumOfItems) {
        // calculate value
        double value = - mtd_hist.CalculateValue(point,mtd_time);
        Energy[loc++] = value;
        return;
    }

    const CColVariable* p_coord = &Sizes[cv];

// cycle through variable
    for(unsigned int i = 0; i < p_coord->GetNumberOfBins(); i++) {
        point[cv] = p_coord->GetValue(i);
        CalculateFES_Part(mtd_hist,point,mtd_time,loc,cv+1);
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CEnergySurface::CalculateFES(unsigned int ncoords,CSimpleVector<double>& params)
{
// quick compatibility comparison
    if( NumOfItems != ncoords ) {
        RUNTIME_ERROR("numofitems mismatch");
    }

// allocate point
    CSimpleVector<double> point;
    point.CreateVector(NumOfItems);

// calculate fes
    unsigned int loc = 0;
    CalculateFES_MTDParam_Part(params,point,loc,0);
}

//------------------------------------------------------------------------------

void CEnergySurface::CalculateFES_MTDParam_Part(
    CSimpleVector<double>& params,
    CSimpleVector<double>& point,
    unsigned int& loc,
    unsigned int cv)
{
    if(cv >= NumOfItems) {
        // calculate value
        double value = - CalculateValue(params,point);
        Energy[loc++] = value;
        return;
    }

    const CColVariable* p_coord = &Sizes[cv];

// cycle through variable
    for(unsigned int i = 0; i < p_coord->GetNumberOfBins(); i++) {
        point[cv] = p_coord->GetValue(i);
        CalculateFES_MTDParam_Part(params,point,loc,cv+1);
    }
}

//------------------------------------------------------------------------------

double CEnergySurface::CalculateValue(const CSimpleVector<double>& params,
                                      const CSimpleVector<double>& point)
{
    double     value = 0.0;
    double     fexparg;
    int        num_of_hills = params.GetLength()/(1+2*NumOfItems);
    unsigned int        loc = 0;

    for(int i=0; i < num_of_hills; i++) {
        fexparg = 0.0;
        double height = params[loc++];
        for(unsigned int k=0; k < NumOfItems; k++) {
            double value = params[loc++];
            double width = params[loc++];
            double e = point[k] - value;
            fexparg = fexparg + e*e / (2.0 * width * width);
        }
        value = value + height*exp(-fexparg);
    }

    return(value);
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

void CEnergySurface::GetPoint(unsigned int index,CSimpleVector<double>& point)
{
    for(int k=NumOfItems-1; k >= 0; k--) {
        CColVariable* p_coord = &Sizes[k];
        int ibin = index % p_coord->GetNumberOfBins();
        point[k] = p_coord->GetValue(ibin);
        index = index / p_coord->GetNumberOfBins();
    }
}

//------------------------------------------------------------------------------

double CEnergySurface::GetGlobalMinimumValue(void) const
{
    double minimum = 0.0;

    if(TotNPoints > 0) minimum = Energy[0];

    for(unsigned int k=0; k < TotNPoints; k++) {
        if(minimum > Energy[k]) minimum = Energy[k];
    }

    return(minimum);
}

//------------------------------------------------------------------------------

void CEnergySurface::ApplyOffset(double offset)
{
    for(unsigned int k=0; k < TotNPoints; k++) {
        Energy[k] += offset;
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

    for(unsigned int k=0; k < TotNPoints; k++) {
        Energy[k] += source.Energy[k];
    }
}

//------------------------------------------------------------------------------

void CEnergySurface::operator/=(const double& number)
{
    for(unsigned int k=0; k < TotNPoints; k++) {
        Energy[k] /= number;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================



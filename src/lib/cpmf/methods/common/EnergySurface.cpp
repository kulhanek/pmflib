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

//------------------------------------------------------------------------------

using namespace std;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CEnergySurface::CEnergySurface(void)
{
    NumOfCVs = 0;
    TotNPoints = 0;
    SLevel = 1.0;
}

//------------------------------------------------------------------------------

CEnergySurface::~CEnergySurface(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CEnergySurface::GetNumOfCVs(void) const
{
    return(NumOfCVs);
}

//------------------------------------------------------------------------------

int CEnergySurface::GetNumOfPoints(void) const
{
    return(TotNPoints);
}

//------------------------------------------------------------------------------

const CColVariable* CEnergySurface::GetCV(unsigned int cv) const
{
    return(&CVs[cv]);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CEnergySurface::Allocate(const CMTDHistory* mtd_hist)
{
    if( NumOfCVs > 0 ) Deallocate();
    if( mtd_hist == NULL ) return;
    if( mtd_hist->GetNumOfCVs() <= 0 ) return;

// allocate items
    NumOfCVs = mtd_hist->GetNumOfCVs();
    CVs.CreateVector(NumOfCVs);

// copy cvs and calculate total number of points
    TotNPoints = 1;
    for(int i=0; i < NumOfCVs; i++) {
        CVs[i].CopyFrom(mtd_hist->GetCV(i));
        TotNPoints *= CVs[i].GetNumOfBins();
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
    if( abf_accu->GetNumOfCVs() <= 0 ) return;

// allocate items
    NumOfCVs = abf_accu->GetNumOfCVs();
    CVs.CreateVector(NumOfCVs);

// copy cvs and calculate total number of points
    TotNPoints = 1;
    for(int i=0; i < NumOfCVs; i++) {
        CVs[i].CopyFrom(abf_accu->GetCV(i));
        TotNPoints *= CVs[i].GetNumOfBins();
    }

    Energy.CreateVector(TotNPoints);
    Error.CreateVector(TotNPoints);
    Samples.CreateVector(TotNPoints);

    Clear();
}

//------------------------------------------------------------------------------

void CEnergySurface::Allocate(const CABPAccumulator* abp_accu)
{
    if( NumOfCVs > 0 ) Deallocate();
    if( abp_accu == NULL ) return;
    if( abp_accu->GetNumOfCVs() <= 0 ) return;

// allocate items
    NumOfCVs = abp_accu->GetNumOfCVs();
    CVs.CreateVector(NumOfCVs);

// copy cvs and calculate total number of points
    TotNPoints = 1;
    for(int i=0; i < NumOfCVs; i++) {
        CVs[i].CopyFrom(abp_accu->GetCV(i));
        TotNPoints *= CVs[i].GetNumOfBins();
    }

    Energy.CreateVector(TotNPoints);
    Error.CreateVector(TotNPoints);
    Samples.CreateVector(TotNPoints);

    Clear();
}

//------------------------------------------------------------------------------

void CEnergySurface::Allocate(const CABFAccumulator* abf_accu,const std::vector<bool>& enabled_cvs)
{
    if( NumOfCVs > 0 ) Deallocate();
    if( abf_accu == NULL ) return;
    if( abf_accu->GetNumOfCVs() <= 0 ) return;

// allocate items
    NumOfCVs = 0;
    for(size_t i = 0; i < enabled_cvs.size(); i++){
        if( enabled_cvs[i] ) NumOfCVs++;
    }
    CVs.CreateVector(NumOfCVs);

// copy cvs and calculate total number of points
    TotNPoints = 1;
    size_t i = 0;
    for(size_t k = 0; k < enabled_cvs.size(); k++){
        if( enabled_cvs[k] ){
            CVs[i].CopyFrom(abf_accu->GetCV(k));
            TotNPoints *= CVs[i].GetNumOfBins();
        }
    }

    Energy.CreateVector(TotNPoints);
    Error.CreateVector(TotNPoints);
    Samples.CreateVector(TotNPoints);

    Clear();
}

//------------------------------------------------------------------------------

void CEnergySurface::Allocate(const CEnergySurface* p_surf,const std::vector<bool>& enabled_cvs)
{
    if( NumOfCVs > 0 ) Deallocate();
    if( p_surf == NULL ) return;
    if( p_surf->GetNumOfCVs() <= 0 ) return;

// allocate items
    NumOfCVs = 0;
    for(size_t i = 0; i < enabled_cvs.size(); i++){
        if( enabled_cvs[i] ) NumOfCVs++;
    }
    CVs.CreateVector(NumOfCVs);

// copy cvs and calculate total number of points
    TotNPoints = 1;
    size_t i = 0;
    for(size_t k = 0; k < enabled_cvs.size(); k++){
        if( enabled_cvs[k] ){
            CVs[i].CopyFrom(p_surf->GetCV(k));
            TotNPoints *= CVs[i].GetNumOfBins();
            i++;
        }
    }

    Energy.CreateVector(TotNPoints);
    Error.CreateVector(TotNPoints);
    Samples.CreateVector(TotNPoints);

    Clear();
}

//------------------------------------------------------------------------------

void CEnergySurface::Deallocate(void)
{
    CVs.FreeVector();
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

double CEnergySurface::GetSigmaF2(bool includeglued) const
{
    double sigmafsum = 0.0;
    double sigmafsum2 = 0.0;
    double count = 0.0;

    for(int k=0; k < TotNPoints; k++) {
        if( Samples[k] == 0 ) continue;
        if( includeglued == false ){
            if( Samples[k] < 0 ) continue;
        }
        double ene = Energy[k];
        sigmafsum += ene;
        sigmafsum2 += ene*ene;
        count++;
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

//------------------------------------------------------------------------------

double CEnergySurface::GetSigmaF2p(bool includeglued) const
{
    double sigmafsum = 0.0;
    double sigmafsum2 = 0.0;
    double count = 0.0;

    for(int k=0; k < TotNPoints; k++) {
        if( Samples[k] == 0 ) continue;
        if( includeglued == false ){
            if( Samples[k] < 0 ) continue;
        }
        double ene = Energy[k] + SLevel*Error[k];
        sigmafsum += ene;
        sigmafsum2 += ene*ene;
        count++;
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

//------------------------------------------------------------------------------

double CEnergySurface::GetSigmaF2m(bool includeglued) const
{
    double sigmafsum = 0.0;
    double sigmafsum2 = 0.0;
    double count = 0.0;

    for(int k=0; k < TotNPoints; k++) {
        if( Samples[k] == 0 ) continue;
        if( includeglued == false ){
            if( Samples[k] < 0 ) continue;
        }
        double ene = Energy[k] - SLevel*Error[k];
        sigmafsum += ene;
        sigmafsum2 += ene*ene;
        count++;
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
        const CColVariable* p_coord = &CVs[k];
        int ibin = index % p_coord->GetNumOfBins();
        point[k] = p_coord->GetValue(ibin);
        index = index / p_coord->GetNumOfBins();
    }
}

//------------------------------------------------------------------------------

void CEnergySurface::GetIPoint(unsigned int index,CSimpleVector<int>& point) const
{
    for(int k=NumOfCVs-1; k >= 0; k--) {
        const CColVariable* p_coord = &CVs[k];
        int ibin = index % p_coord->GetNumOfBins();
        point[k] = ibin;
        index = index / p_coord->GetNumOfBins();
    }
}

//------------------------------------------------------------------------------

int CEnergySurface::IPoint2Bin(const CSimpleVector<int>& point)
{
    int idx = 0;

    for(int i=0; i < NumOfCVs; i++) {
        CColVariable cv = CVs[i];
        int idx_local = point[i];
        if( (idx_local < 0) || (idx_local >= cv.GetNumOfBins())){
            return(-1);
        }
        idx = idx*cv.GetNumOfBins() + idx_local;
    }

    return(idx);
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

void CEnergySurface::SetSLevel(double slevel)
{
    SLevel = slevel;
}

//------------------------------------------------------------------------------

double CEnergySurface::GetSLevel(void) const
{
    return(SLevel);
}

//------------------------------------------------------------------------------

void CEnergySurface::AdaptUnsampledToMaxEnergy(void)
{
    double maxe = 0.0;
    bool   first = true;

    // find maximum in sampled region
    for(int k=0; k < TotNPoints; k++) {
        if( Samples[k] == 0 ) continue; // skip unsampled
        if( (maxe < Energy[k]) || (first == true) ){
            maxe = Energy[k];
            first = false;
        }
    }

    // adapt unsampled region
    for(int k=0; k < TotNPoints; k++) {
        if( Samples[k] == 0 ){
            Energy[k] = maxe;
        }
    }
}

//------------------------------------------------------------------------------

void CEnergySurface::AdaptUnsampledToMaxEnergy(double maxene)
{
    // adapt unsampled region
    for(int k=0; k < TotNPoints; k++) {
        if( Samples[k] == 0 ){
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

bool CEnergySurface::ReduceFES(const std::vector<bool>& keepcvs,double temp,CEnergySurface* p_rsurf)
{
    if( p_rsurf == NULL ){
        ES_ERROR("p_rsurf == NULL");
        return(false);
    }
    if( keepcvs.size() != (size_t)NumOfCVs ){
        RUNTIME_ERROR("keepcvs.size() != NumOfCVs");
    }

    p_rsurf->Allocate(this,keepcvs);
    p_rsurf->Clear();

    CSimpleVector<int>    midx;
    midx.CreateVector(NumOfCVs);

    CSimpleVector<int>    ridx;
    ridx.CreateVector(p_rsurf->NumOfCVs);

    const double R = 1.98720425864083e-3;

// calculate weights
    for(size_t mbin = 0; mbin < (size_t)TotNPoints; mbin++){
        double ene = GetEnergy(mbin);
        GetIPoint(mbin,midx);
        ReduceIPoint(keepcvs,midx,ridx);
        int rbin = p_rsurf->IPoint2Bin(ridx);
        if( rbin == -1 ){
            for(size_t i=0; i < ridx.GetLength(); i++){
                cout << ridx[i] << " ";
            }
            cout << endl;
            for(int i=0; i < p_rsurf->GetNumOfCVs(); i++){
                cout << p_rsurf->GetCV(i)->GetNumOfBins() << " ";
            }
            cout << endl;
            RUNTIME_ERROR("rbin == -1");
        }
        double w = exp(-ene/(R*temp));
        p_rsurf->SetEnergy(rbin,p_rsurf->GetEnergy(rbin) + w);
        p_rsurf->SetNumOfSamples(rbin,1);
    }

// transform back to FE
    for(size_t rbin = 0; rbin < (size_t)p_rsurf->TotNPoints; rbin++){
        double w = p_rsurf->GetEnergy(rbin);
        p_rsurf->SetEnergy(rbin,-R*temp*log(w));
    }

// move global minimum
   double gmin = p_rsurf->GetGlobalMinimumValue();
   p_rsurf->ApplyOffset(-gmin);

    return(true);
}

//------------------------------------------------------------------------------

void CEnergySurface::ReduceIPoint(const std::vector<bool>& keepcvs,CSimpleVector<int>& midx,CSimpleVector<int>& ridx)
{
    ridx.SetZero();
    size_t j = 0;
    for(int i = 0; i < NumOfCVs; i++){
        if( keepcvs[i] ){
            ridx[j] = midx[i];
            j++;
        }
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================



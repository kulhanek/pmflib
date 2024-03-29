// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
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
    NumOfBins = 0;
    SLevel = 1.0;

    Temperature = 300;
    TemperatureUnit = "K";
    TemperatureFConv = 1.0;

    EnergyFConv = 1.0;
    EnergyUnit = "kcal/mol";

    GlobalMinSet    = false;
    GPosSet         = false;
    GPosBin         = 0;
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

int CEnergySurface::GetNumOfBins(void) const
{
    return(NumOfBins);
}

//------------------------------------------------------------------------------

const CColVariablePtr CEnergySurface::GetCV(int cv) const
{
    if( (cv < 0) || (cv >= NumOfCVs) ) {
        RUNTIME_ERROR("cv out-of-range");
    }
    return(CVs[cv]);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CEnergySurface::Allocate(CPMFAccumulatorPtr accu)
{
    if( NumOfCVs > 0 ) Deallocate();
    if( accu == NULL ) return;
    if( accu->GetNumOfCVs() <= 0 ) return;

// copy cvs and calculate total number of points
    NumOfBins = 1;
    for(int i=0; i < accu->GetNumOfCVs(); i++) {
        CColVariablePtr cv = accu->GetCV(i);
        NumOfBins *= cv->GetNumOfBins();
        CVs.push_back(cv);
    }
    NumOfCVs = CVs.size();

    Temperature = accu->GetTemperature();
    TemperatureFConv = accu->GetTemperatureFConv();
    TemperatureUnit = accu->GetTemperatureUnit();

    EnergyFConv = accu->GetEnergyFConv();
    EnergyUnit = accu->GetEnergyUnit();

    Energy.CreateVector(NumOfBins);
    Error.CreateVector(NumOfBins);
    Samples.CreateVector(NumOfBins);

    Clear();
}

//------------------------------------------------------------------------------

void CEnergySurface::Allocate(CEnergySurfacePtr surf,const std::vector<bool>& enabled_cvs)
{
    if( NumOfCVs > 0 ) Deallocate();
    if( surf == NULL ) return;
    if( surf->GetNumOfCVs() <= 0 ) return;
    if( (int)enabled_cvs.size() != surf->GetNumOfCVs() ){
        RUNTIME_ERROR("enabled_cvs and surf inconsistent");
    }

// copy cvs and calculate total number of points
    NumOfBins = 1;
    size_t i = 0;
    for(size_t k = 0; k < enabled_cvs.size(); k++){
        if( enabled_cvs[k] ){
            CColVariablePtr cv = surf->GetCV(i);
            NumOfBins *= cv->GetNumOfBins();
            CVs.push_back(cv);
        }
    }
    NumOfCVs = CVs.size();

    Temperature         = surf->GetTemperature();
    TemperatureFConv    = surf->GetTemperatureFConv();
    TemperatureUnit     = surf->GetTemperatureUnit();

    EnergyFConv = surf->GetEnergyFConv();
    EnergyUnit  = surf->GetEnergyUnit();

    Energy.CreateVector(NumOfBins);
    Error.CreateVector(NumOfBins);
    Samples.CreateVector(NumOfBins);

    Clear();
}

//------------------------------------------------------------------------------

void CEnergySurface::Deallocate(void)
{
    CVs.clear();
    Energy.FreeVector();
    Error.FreeVector();
    Samples.FreeVector();
    NumOfCVs = 0;
    NumOfBins = 0;
}

//------------------------------------------------------------------------------

void CEnergySurface::Clear(void)
{
    for(int i=0; i < NumOfBins; i++) {
        Energy[i] = 0.0;
        Error[i] = 0.0;
        Samples[i] = 0;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CEnergySurface::SetEnergy(int ibin,const double value)
{
    if( (ibin < 0) || (ibin >= NumOfBins) ){
        RUNTIME_ERROR("ibin out-of-range");
    }

    Energy[ibin] = value;
}

//------------------------------------------------------------------------------

double CEnergySurface::GetEnergy(int ibin) const
{
    if( (ibin < 0) || (ibin >= NumOfBins) ){
        RUNTIME_ERROR("ibin out-of-range");
    }
    return(Energy[ibin]);
}

//------------------------------------------------------------------------------

double CEnergySurface::GetEnergyRealValue(int ibin) const
{
    return( GetEnergy(ibin) * EnergyFConv);
}

//------------------------------------------------------------------------------

void CEnergySurface::SetError(int ibin,double value)
{
    Error[ibin] = value;
}

//------------------------------------------------------------------------------

double CEnergySurface::GetError(int ibin) const
{
    return(Error[ibin]);
}

//------------------------------------------------------------------------------

double  CEnergySurface::GetErrorRealValue(int ibin) const
{
    return(Error[ibin] * EnergyFConv);
}

//------------------------------------------------------------------------------

double  CEnergySurface::GetErrorRealValueWithSLevel(int ibin) const
{
    return(SLevel * Error[ibin] * EnergyFConv);
}

//------------------------------------------------------------------------------

double CEnergySurface::GetSigmaF2(bool includeglued) const
{
    double mf = 0.0;
    double m2 = 0.0;
    double n  = 0.0;

    // https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Welford's_online_algorithm

    for(int k=0; k < NumOfBins; k++) {
        if( Samples[k] == 0 ) continue;
        if( includeglued == false ){
            if( Samples[k] < 0 ) continue;
        }
        double ene = Energy[k];
        n++;
        double dx1 = (ene - mf);
        mf = mf + dx1/n;
        double dx2 = (ene - mf);
        m2 = m2 + dx1*dx2;
    }

    double sigmaf2 = 0.0;
    if( n > 0 ) {
        sigmaf2 = m2/n;
    }

    return(sigmaf2);
}

//------------------------------------------------------------------------------

double CEnergySurface::GetSigmaF(bool includeglued) const
{
    double mf = 0.0;
    double m2 = 0.0;
    double n  = 0.0;

    // https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Welford's_online_algorithm

    for(int k=0; k < NumOfBins; k++) {
        if( Samples[k] == 0 ) continue;
        if( includeglued == false ){
            if( Samples[k] < 0 ) continue;
        }
        double ene = Energy[k];
        n++;
        double dx1 = (ene - mf);
        mf = mf + dx1/n;
        double dx2 = (ene - mf);
        m2 = m2 + dx1*dx2;
    }

    double sigmaf = 0.0;
    if( n > 0 ) {
        sigmaf = sqrt(m2/n);
    }

    return(sigmaf);
}

//------------------------------------------------------------------------------

double CEnergySurface::GetRMSError(bool includeglued) const
{
    double m2 = 0.0;
    double n  = 0.0;

    for(int k=0; k < NumOfBins; k++) {
        if( Samples[k] == 0 ) continue;
        if( includeglued == false ){
            if( Samples[k] < 0 ) continue;
        }
        double err = SLevel*Error[k];
        m2 = m2 + err*err;
        n++;
    }

    double rmserr = 0.0;
    if( n > 0 ) {
        rmserr = sqrt(m2/n);
    }

    return(rmserr);
}

//------------------------------------------------------------------------------

double CEnergySurface::GetMaxError(bool includeglued) const
{
    double maxerr = 0.0;

    for(int k=0; k < NumOfBins; k++) {
        if( Samples[k] == 0 ) continue;
        if( includeglued == false ){
            if( Samples[k] < 0 ) continue;
        }
        double err = SLevel*Error[k];
        if( maxerr < err ){
            maxerr = err;
        }
    }

    return(maxerr);
}

//------------------------------------------------------------------------------

double CEnergySurface::GetSigmaF2All(void) const
{
    double mf = 0.0;
    double m2 = 0.0;
    double n  = 0.0;

    // https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Welford's_online_algorithm

    for(int k=0; k < NumOfBins; k++) {
        double ene = Energy[k];
        n++;
        double dx1 = (ene - mf);
        mf = mf + dx1/n;
        double dx2 = (ene - mf);
        m2 = m2 + dx1*dx2;
    }

    double sigmaf2 = 0.0;
    if( n > 0 ) {
        sigmaf2 = m2/n;
    }

    return(sigmaf2);
}

//------------------------------------------------------------------------------

double CEnergySurface::GetTemperature(void) const
{
    return(Temperature);
}

//------------------------------------------------------------------------------

double CEnergySurface::GetRealTemperature(void) const
{
    return(Temperature*TemperatureFConv);
}

//------------------------------------------------------------------------------

const CSmallString& CEnergySurface::GetTemperatureUnit(void) const
{
    return(TemperatureUnit);
}

//------------------------------------------------------------------------------

double CEnergySurface::GetTemperatureFConv(void)
{
    return(TemperatureFConv);
}

//------------------------------------------------------------------------------

double CEnergySurface::GetEnergyFConv(void)
{
    return(EnergyFConv);
}

//------------------------------------------------------------------------------

const CSmallString& CEnergySurface::GetEnergyUnit(void) const
{
    return(EnergyUnit);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CEnergySurface::SetNumOfSamples(int ibin,int value)
{
    Samples[ibin] = value;
}

//------------------------------------------------------------------------------

int CEnergySurface::GetNumOfSamples(int ibin) const
{
    return(Samples[ibin]);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CEnergySurface::GetPoint(int ibin,CSimpleVector<double>& point) const
{
    for(int k=NumOfCVs-1; k >= 0; k--) {
        const CColVariablePtr p_coord = CVs[k];
        int cvbin = ibin % p_coord->GetNumOfBins();
        point[k] = p_coord->GetValue(cvbin);
        ibin = ibin / p_coord->GetNumOfBins();
    }
}

//------------------------------------------------------------------------------

void CEnergySurface::GetIPoint(int ibin,CSimpleVector<int>& point) const
{
    for(int k=NumOfCVs-1; k >= 0; k--) {
        const CColVariablePtr p_coord = CVs[k];
        int cvbin = ibin % p_coord->GetNumOfBins();
        point[k] = cvbin;
        ibin = ibin / p_coord->GetNumOfBins();
    }
}

//------------------------------------------------------------------------------

int CEnergySurface::GetGlobalIndex(const CSimpleVector<int>& position) const
{
    int glbindex = 0;
    for(int i=0; i < NumOfCVs; i++) {
        if( position[i] < 0 ) return(-1);
        if( position[i] >= CVs[i]->GetNumOfBins() ) return(-1);
        glbindex = glbindex*CVs[i]->GetNumOfBins() + position[i];
    }
    return(glbindex);
}

//------------------------------------------------------------------------------

int CEnergySurface::IPoint2Bin(const CSimpleVector<int>& point)
{
    int idx = 0;

    for(int i=0; i < NumOfCVs; i++) {
        CColVariablePtr cv = CVs[i];
        int idx_local = point[i];
        if( (idx_local < 0) || (idx_local >= cv->GetNumOfBins())){
            return(-1);
        }
        idx = idx*cv->GetNumOfBins() + idx_local;
    }

    return(idx);
}

//------------------------------------------------------------------------------

double CEnergySurface::GetGlobalMinimumValue(void) const
{
    double minimum = 0.0;

    bool first = true;

    for(int k=0; k < NumOfBins; k++) {
        if( Samples[k] != 0 ) {
            if( (minimum > Energy[k]) || first) {
                minimum = Energy[k];
                first = false;
            }
        }
    }
    return(minimum);
}

//------------------------------------------------------------------------------

void CEnergySurface::ApplyOffset(double offset)
{
    for(int k=0; k < NumOfBins; k++) {
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
    for(int k=0; k < NumOfBins; k++) {
        if( Samples[k] == 0 ) continue; // skip unsampled
        if( (maxe < Energy[k]) || (first == true) ){
            maxe = Energy[k];
            first = false;
        }
    }

    // adapt unsampled region
    for(int k=0; k < NumOfBins; k++) {
        if( Samples[k] == 0 ){
            Energy[k] = maxe;
        }
    }
}

//------------------------------------------------------------------------------

void CEnergySurface::AdaptUnsampledToMaxEnergy(double maxene)
{
    // adapt unsampled region
    for(int k=0; k < NumOfBins; k++) {
        if( Samples[k] == 0 ){
            Energy[k] = maxene;
        }
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CEnergySurface::AddFES(CEnergySurfacePtr source)
{
    if(NumOfBins != source->NumOfBins) {
        ES_ERROR("surfaces do not match");
        return;
    }

    for(int k=0; k < NumOfBins; k++) {
        Energy[k]   += source->Energy[k];
        Samples[k]  += source->Samples[k];
    }
}

//------------------------------------------------------------------------------

void CEnergySurface::DivideFES(double number)
{
    for(int k=0; k < NumOfBins; k++) {
        Energy[k] /= number;
        // do not divide Samples
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CEnergySurfacePtr CEnergySurface::ReduceFES(const std::vector<bool>& keepcvs)
{
    if( keepcvs.size() != (size_t)NumOfCVs ){
        RUNTIME_ERROR("keepcvs.size() != NumOfCVs");
    }

    CEnergySurfacePtr p_rsurf = CEnergySurfacePtr(new CEnergySurface);

// copy cvs and calculate total number of points
    p_rsurf->NumOfBins = 1;
    for(size_t k = 0; k < keepcvs.size(); k++){
        if( keepcvs[k] ){
            CColVariablePtr cv = GetCV(k);
            p_rsurf->CVs.push_back(cv);
            p_rsurf->NumOfBins *= cv->GetNumOfBins();
        }
    }
    p_rsurf->NumOfCVs = p_rsurf->CVs.size();

    p_rsurf->Energy.CreateVector(p_rsurf->NumOfBins);
    p_rsurf->Error.CreateVector(p_rsurf->NumOfBins);
    p_rsurf->Samples.CreateVector(p_rsurf->NumOfBins);

    p_rsurf->Clear();

    CSimpleVector<int>    midx;
    midx.CreateVector(NumOfCVs);

    CSimpleVector<int>    ridx;
    ridx.CreateVector(p_rsurf->NumOfCVs);

    const double R = 1.98720425864083e-3;

// calculate weights
    for(size_t mbin = 0; mbin < (size_t)NumOfBins; mbin++){
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
        double w = exp(-ene/(R*Temperature));
        p_rsurf->SetEnergy(rbin,p_rsurf->GetEnergy(rbin) + w);
        p_rsurf->SetNumOfSamples(rbin,1);
    }

// transform back to FE
    for(size_t rbin = 0; rbin < (size_t)p_rsurf->NumOfBins; rbin++){
        double w = p_rsurf->GetEnergy(rbin);
        p_rsurf->SetEnergy(rbin,-R*Temperature*log(w));
    }

// move global minimum
    double gmin = p_rsurf->GetGlobalMinimumValue();
    p_rsurf->ApplyOffset(-gmin);

    return(p_rsurf);
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

void CEnergySurface::SetGlobalMin(const CSmallString& spec)
{
    if( NumOfCVs == 0 ){
        RUNTIME_ERROR("accumulator is not set for SetGlobalMin");
    }

    GlobalMinSet = true;
    string sspec(spec);

    // remove "x" from the string
    replace (sspec.begin(), sspec.end(), 'x' , ' ');

    // parse values of CVs
    GPos.CreateVector(NumOfCVs);
    stringstream str(sspec);
    for(int i=0; i < NumOfCVs; i++){
        double val;
        str >> val;
        if( ! str ){
            CSmallString error;
            error << "unable to decode CV value for position: " << i+1;
            RUNTIME_ERROR(error);
        }
        GPos[i] = GetCV(i)->GetIntValue(val);
    }

    GPosSet = true;
}

//------------------------------------------------------------------------------

void CEnergySurface::SetGlobalMin(const CSimpleVector<double>& pos)
{
    GlobalMinSet = true;
    GPos = pos;
    GPosSet = true;
}

//------------------------------------------------------------------------------

bool CEnergySurface::IsGlobalMinSet(void)
{
    return(GPosSet);
}

//------------------------------------------------------------------------------

CSimpleVector<double> CEnergySurface::GetGlobalMinPos(void)
{
    if( GPosSet == false ){
        RUNTIME_ERROR("no global min set")
    }
    return(GPos);
}

//------------------------------------------------------------------------------

int CEnergySurface::GetGlobalMinBin(void)
{
    if( GPosSet == false ){
        RUNTIME_ERROR("no global min set")
    }
    return(GPosBin);
}

//------------------------------------------------------------------------------

double CEnergySurface::GetGlobalMinEnergy(void)
{
   return(GetEnergy(GPosBin));
}

//------------------------------------------------------------------------------

void CEnergySurface::FindGlobalMinBin(void)
{
    // find the closest bin
    CSimpleVector<double>   pos;
    pos.CreateVector(NumOfCVs);
    double minv = 0.0;
    GPosBin = 0;
    for(int ibin=0; ibin < NumOfBins; ibin++){
        GetPoint(ibin,pos);
        double dist2 = 0.0;
        for(int cv=0; cv < NumOfCVs; cv++){
            dist2 = dist2 + (pos[cv]-GPos[cv])*(pos[cv]-GPos[cv]);
        }
        if( ibin == 0 ){
            minv = dist2;
            GPosBin = 0;
        }
        if( dist2 < minv ){
            minv = dist2;
            GPosBin = ibin;
        }
    }

    GetPoint(GPosBin,GPos);
    GPosSet = true;
}

//------------------------------------------------------------------------------

void CEnergySurface::FindGlobalMin(void)
{
    GPosSet = false;
    GPos.CreateVector(NumOfCVs);
    bool   first = true;
    double glb_min = 0.0;
    for(int ibin=0; ibin < NumOfBins; ibin++){
        int samples = GetNumOfSamples(ibin);
        if( samples < -1 ) continue;    // include sampled areas and holes but exclude extrapolated areas
        double value = GetEnergy(ibin);
        if( first || (glb_min > value) ){
            glb_min = value;
            first = false;
            GPosBin = ibin;
            GetPoint(ibin,GPos);
            GPosSet = true;
        }
    }

    if( GPosSet == false ){
        // fallback
        GetPoint(0,GPos);
        GPosSet = true;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================



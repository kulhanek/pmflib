// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
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
// =============================================================================

#include <AUSEngineCST_A.hpp>
#include <CSTProxy_dAdx.hpp>
#include <CSTProxy_dU.hpp>
#include <CSTProxy_mTdSdx.hpp>
#include <CSTProxy_Ecorr.hpp>
#include <iomanip>

//------------------------------------------------------------------------------

using namespace std;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CAUSEngineCST_A::CAUSEngineCST_A(void)
{
    RegisterRealm(0, "AUS", "CST", "dA | dU | -TdS");

    NoEnergy = false;
}

//------------------------------------------------------------------------------

CAUSEngineCST_A::~CAUSEngineCST_A(void)
{
}

//------------------------------------------------------------------------------

void CAUSEngineCST_A::SetAccumulator(CPMFAccumulatorPtr accu)
{
    if( accu == NULL ) return;

    // set number of hyprms
    NumOfSigmaF2 = 5;
    NumOfCoVar = 0;
    CGPRKernel::SetAccumulator(accu);
    NumOfSigmaN2 = 2*NumOfCVs + 3;

    A_ES = CEnergySurfacePtr(new CEnergySurface);
    A_ES->Allocate(accu);
    A.SetOutputES(A_ES);

    CCSTProxy_dAdx_Ptr a_proxy = CCSTProxy_dAdx_Ptr(new CCSTProxy_dAdx());
    a_proxy->SetRealm(CST_dLdx);
    a_proxy->Init(accu);
    A.SetInputEnergyDerProxy(a_proxy);

    B_ES = CEnergySurfacePtr(new CEnergySurface);
    B_ES->Allocate(accu);
    B.SetOutputES(B_ES);

    CCSTProxy_mTdSdx_Ptr b_proxy = CCSTProxy_mTdSdx_Ptr(new CCSTProxy_mTdSdx());
    b_proxy->SetRealm(CST_TdS_LT);
    b_proxy->Init(accu);
    B.SetInputEnergyDerProxy(b_proxy);

    C_ES = CEnergySurfacePtr(new CEnergySurface);
    C_ES->Allocate(accu);
    C.SetOutputES(C_ES);

    CCSTProxy_dU_Ptr c_proxy = CCSTProxy_dU_Ptr(new CCSTProxy_dU());
    c_proxy->SetRealm(CST_ETOTFW);
    c_proxy->Init(accu);
    C.SetInputEnergyProxy(c_proxy);

    D_ES = CEnergySurfacePtr(new CEnergySurface);
    D_ES->Allocate(accu);
    D.SetOutputES(D_ES);

    CCSTProxy_dU_Ptr d_proxy = CCSTProxy_dU_Ptr(new CCSTProxy_dU());
    d_proxy->SetRealm(CST_ETOT);
    d_proxy->Init(accu);
    D.SetInputEnergyProxy(d_proxy);

    E_ES = CEnergySurfacePtr(new CEnergySurface);
    E_ES->Allocate(accu);
    E.SetOutputES(E_ES);

    CCSTProxy_Ecorr_Ptr e_proxy = CCSTProxy_Ecorr_Ptr(new CCSTProxy_Ecorr());
    e_proxy->SetRealm(CST_dA_corr);
    e_proxy->Init(accu);
    E.SetInputEnergyProxy(e_proxy);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CAUSEngineCST_A::SetIncludeError(bool set)
{
    A.SetIncludeError(set);
    B.SetIncludeError(set);
    C.SetIncludeError(set);
    D.SetIncludeError(set);
    E.SetIncludeError(set);

    IncludeError = set;
}

//-----------------------------------------------------------------------------

void CAUSEngineCST_A::SetNoEnergy(bool set)
{
    A.SetNoEnergy(set);
    B.SetNoEnergy(set);
    C.SetNoEnergy(set);
    D.SetNoEnergy(set);
    E.SetNoEnergy(set);

    NoEnergy = set;
}

//-----------------------------------------------------------------------------

void CAUSEngineCST_A::IncludeGluedAreas(bool set)
{
    A.IncludeGluedAreas(set);
    B.IncludeGluedAreas(set);
    C.IncludeGluedAreas(set);
    D.IncludeGluedAreas(set);
    E.IncludeGluedAreas(set);
}

//-----------------------------------------------------------------------------

void CAUSEngineCST_A::PrepForHyprmsGrd(bool set)
{
    A.PrepForHyprmsGrd(set);
    B.PrepForHyprmsGrd(set);
    C.PrepForHyprmsGrd(set);
    D.PrepForHyprmsGrd(set);
    E.PrepForHyprmsGrd(set);

   NeedInv |= set;
}

//-----------------------------------------------------------------------------

void CAUSEngineCST_A::SetCalcLogPL(bool set)
{
    A.SetCalcLogPL(set);
    B.SetCalcLogPL(set);
    C.SetCalcLogPL(set);
    D.SetCalcLogPL(set);
    E.SetCalcLogPL(set);

   NeedInv |= set;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CAUSEngineCST_A::RunGPR(CVerboseStr& vout,bool nostat)
{
    vout        << "   AUS proxy ... " << GetFullDescription() << endl;

    bool result = true;

    A.SetLAMethod(Method);
    A.SetRCond(RCond);
    A.SetUseInv(UseInv);
    A.SetKernel(GetKernel());

    B.SetLAMethod(Method);
    B.SetRCond(RCond);
    B.SetUseInv(UseInv);
    B.SetKernel(GetKernel());

    C.SetLAMethod(Method);
    C.SetRCond(RCond);
    C.SetUseInv(UseInv);
    C.SetKernel(GetKernel());
    // C.UseFirstKernelDerivatives(true);

    D.SetLAMethod(Method);
    D.SetRCond(RCond);
    D.SetUseInv(UseInv);
    D.SetKernel(GetKernel());
    // D.UseFirstKernelDerivatives(true);

    E.SetLAMethod(Method);
    E.SetRCond(RCond);
    E.SetUseInv(UseInv);
    E.SetKernel(GetKernel());
    // E.UseFirstKernelDerivatives(true);

    // it must be re-mapped
    int idx = 0;
    A.SetSigmaF2(0,SigmaF2[idx++]);
    B.SetSigmaF2(0,SigmaF2[idx++]);
    C.SetSigmaF2(0,SigmaF2[idx++]);
    D.SetSigmaF2(0,SigmaF2[idx++]);
    E.SetSigmaF2(0,SigmaF2[idx++]);

    idx = 0;
    for(size_t cv=0; cv < NumOfCVs; cv++){
        double wfac = WFac[idx++];
        A.SetWFac(cv,wfac);
        B.SetWFac(cv,wfac);
        C.SetWFac(cv,wfac);
        D.SetWFac(cv,wfac);
        E.SetWFac(cv,wfac);
    }

    idx = 0;
    for(size_t cv=0; cv < NumOfCVs; cv++){
        A.SetSigmaN2(cv,SigmaN2[idx++]);
    }
    for(size_t cv=0; cv < NumOfCVs; cv++){
        B.SetSigmaN2(cv,SigmaN2[idx++]);
    }
    C.SetSigmaN2(0,SigmaN2[idx++]);
    D.SetSigmaN2(0,SigmaN2[idx++]);
    E.SetSigmaN2(0,SigmaN2[idx++]);

    vout << "   ### ------------" << endl;
    result &= A.RunGPR(vout,nostat);  // <lam>
    vout << "   ### ------------" << endl;
    result &= B.RunGPR(vout,nostat);  // Cov(lam,Etot)
    vout << "   ### ------------" << endl;
    result &= C.RunGPR(vout,nostat);  // d<Etot>_FW
    vout << "   ### ------------" << endl;
    result &= D.RunGPR(vout,nostat);  // d<Etot>
    vout << "   ### ------------" << endl;
    result &= E.RunGPR(vout,nostat);  // -TdS{CST}corr
    vout << "   ### ------------" << endl;

    if( ! nostat ){
        vout << "   Total result ..." << endl;
        // and log of marginal likelihood
            vout << "      logML     = " << setprecision(5) << GetLogML() << endl;
        if( NeedInv || UseInv ){
            // and log of pseudo-likelihood
            vout << "      logPL     = " << setprecision(5) << GetLogPL() << endl;
        }
    }

    // finalize ENE surfaces if requested
    if( ! NoEnergy ){
        CalculateEnergy(vout);
    }

    return(result);
}

//------------------------------------------------------------------------------

void CAUSEngineCST_A::CalculateEnergy(CVerboseStr& vout)
{
    vout << "   ### -----------" << endl;
    vout << "   Calculating ENE surfaces ..." << endl;

// FEN
    vout << "   ** FEN" << endl;

    for(size_t ibin=0; ibin < NumOfBins; ibin++){
        int nsamples = A_ES->GetNumOfSamples(ibin);
        double fen1 = A_ES->GetEnergy(ibin);    // <lam>
        double fen2 = E_ES->GetEnergy(ibin);    // -RT*log(FW)
        double fen = fen1+fen2;
        ASurface->SetNumOfSamples(ibin,nsamples);
        ASurface->SetEnergy(ibin,fen);
        if( IncludeError ) {
            double efen1 = A_ES->GetError(ibin);
            double efen2 = E_ES->GetError(ibin);
            double efen = sqrt(efen1*efen1+efen2*efen2);        // RAW estimate - uncorrelated
            ASurface->SetError(ibin,efen);
        }
    }

// update FES
    if( ASurface->IsGlobalMinSet() ){

        CSimpleVector<double> gpos;

        gpos = ASurface->GetGlobalMinPos();
        vout << "      Global minimum provided at: ";
        vout << setprecision(5) << gpos[0];
        for(size_t i=1; i < NumOfCVs; i++){
            vout << "x" << setprecision(5) << gpos[0];
        }
        vout << endl;

        ASurface->FindGlobalMinBin();

        gpos = ASurface->GetGlobalMinPos();
        vout << "      Closest bin found at: ";
        vout << setprecision(5) << gpos[0];
        for(size_t i=1; i < NumOfCVs; i++){
            vout << "x" << setprecision(5) << gpos[0];
        }

        double glb_min = ASurface->GetGlobalMinEnergy();
        vout << " (" << setprecision(5) << glb_min << ")" << endl;

        for(size_t ibin=0; ibin < NumOfBins; ibin++){
            if( ASurface->GetNumOfSamples(ibin) >  0 ){
                ASurface->SetEnergy(ibin,ASurface->GetEnergy(ibin)-glb_min);
            }
        }

        USurface->SetGlobalMin(gpos);
        SSurface->SetGlobalMin(gpos);
    } else {
        // search for global minimum
        ASurface->FindGlobalMin();

        double                glb_min = ASurface->GetGlobalMinEnergy();
        CSimpleVector<double> gpos    = ASurface->GetGlobalMinPos();

        vout << "      Global minimum found at: ";
        vout << setprecision(5) << gpos[0];
        for(size_t i=1; i < NumOfCVs; i++){
            vout << "x" << setprecision(5) << gpos[0];
        }
        vout << " (" << setprecision(5) << glb_min << ")" << endl;
        for(size_t ibin=0; ibin < NumOfBins; ibin++){
            if( ASurface->GetNumOfSamples(ibin) >  0 ){
                ASurface->SetEnergy(ibin,ASurface->GetEnergy(ibin)-glb_min);
            }
        }

        USurface->SetGlobalMin(gpos);
        SSurface->SetGlobalMin(gpos);
    }

        vout << "      dA(x) SigmaF2     = " << setw(10) << setprecision(5) << ASurface->GetSigmaF2() << endl;
        vout << "      dA(x) SigmaF      = " << setw(10) << setprecision(5) << ASurface->GetSigmaF() << endl;

// INT
    vout << "   ** INT" << endl;

    for(size_t ibin=0; ibin < NumOfBins; ibin++){
        int nsamples = C_ES->GetNumOfSamples(ibin);
        double inte = C_ES->GetEnergy(ibin);            // HcFW
        USurface->SetNumOfSamples(ibin,nsamples);
        USurface->SetEnergy(ibin,inte);
        if( IncludeError ) {
            double einte = C_ES->GetError(ibin);
            USurface->SetError(ibin,einte);
        }
    }

    {
        CSimpleVector<double> gpos;

        gpos = USurface->GetGlobalMinPos();
        vout << "      Global minimum provided at: ";
        vout << setprecision(5) << gpos[0];
        for(size_t i=1; i < NumOfCVs; i++){
            vout << "x" << setprecision(5) << gpos[0];
        }
        vout << endl;

        USurface->FindGlobalMinBin();

        gpos = USurface->GetGlobalMinPos();
        vout << "      Closest bin found at: ";
        vout << setprecision(5) << gpos[0];
        for(size_t i=1; i < NumOfCVs; i++){
            vout << "x" << setprecision(5) << gpos[0];
        }

        double glb_min = USurface->GetGlobalMinEnergy();
        vout << " (" << setprecision(5) << glb_min << ")" << endl;

        for(size_t ibin=0; ibin < NumOfBins; ibin++){
            if( USurface->GetNumOfSamples(ibin) >  0 ){
                USurface->SetEnergy(ibin,USurface->GetEnergy(ibin)-glb_min);
            }
        }
    }

        vout << "      dU(x) SigmaF2     = " << setw(10) << setprecision(5) << USurface->GetSigmaF2() << endl;
        vout << "      dU(x) SigmaF      = " << setw(10) << setprecision(5) << USurface->GetSigmaF() << endl;

// TDS
        vout << "   ** TDS" << endl;
    for(size_t ibin=0; ibin < NumOfBins; ibin++){
        int nsamples = B_ES->GetNumOfSamples(ibin);
        double tds1 = B_ES->GetEnergy(ibin);    // (1/RT)Cov(L,Hc)
        double tds2 = C_ES->GetEnergy(ibin);    // HcFW
        double tds3 = D_ES->GetEnergy(ibin);    // Hc
        double tds4 = E_ES->GetEnergy(ibin);    // -RT*log(FW)
        double tds = tds1+tds4-(tds2-tds3);
        SSurface->SetNumOfSamples(ibin,nsamples);
        SSurface->SetEnergy(ibin,tds);

        if( IncludeError ) {
            double etds1 = B_ES->GetError(ibin);
            double etds2 = C_ES->GetError(ibin);
            double etds3 = D_ES->GetError(ibin);
            double etds4 = E_ES->GetError(ibin);
            double etds = sqrt(etds1*etds1+etds2*etds2+etds3*etds3+etds4*etds4);    // RAW estimate - uncorrelated
            SSurface->SetError(ibin,etds);
        }
    }

    {
        CSimpleVector<double> gpos;

        gpos = SSurface->GetGlobalMinPos();
        vout << "      Global minimum provided at: ";
        vout << setprecision(5) << gpos[0];
        for(size_t i=1; i < NumOfCVs; i++){
            vout << "x" << setprecision(5) << gpos[0];
        }
        vout << endl;

        SSurface->FindGlobalMinBin();

        gpos = SSurface->GetGlobalMinPos();
        vout << "      Closest bin found at: ";
        vout << setprecision(5) << gpos[0];
        for(size_t i=1; i < NumOfCVs; i++){
            vout << "x" << setprecision(5) << gpos[0];
        }

        double glb_min = SSurface->GetGlobalMinEnergy();
        vout << " (" << setprecision(5) << glb_min << ")" << endl;

        for(size_t ibin=0; ibin < NumOfBins; ibin++){
            if( SSurface->GetNumOfSamples(ibin) >  0 ){
                SSurface->SetEnergy(ibin,SSurface->GetEnergy(ibin)-glb_min);
            }
        }
    }

        vout << "      -TdS(x) SigmaF2   = " << setw(10) << setprecision(5) << SSurface->GetSigmaF2() << endl;
        vout << "      -TdS(x) SigmaF    = " << setw(10) << setprecision(5) << SSurface->GetSigmaF() << endl;

    vout << "   ### -----------" << endl;
    CalcResiduals(vout,false);
    if( DoBalanceResiduals ) {
        vout << "   ### -----------" << endl;
        BalanceResiduals();
        CalcResiduals(vout,true);
    }
}

//------------------------------------------------------------------------------

int CAUSEngineCST_A::GetNumOfTasks(void)
{
    return(5);
}

//------------------------------------------------------------------------------

bool CAUSEngineCST_A::WriteMFInfo(const CSmallString& name,int task)
{
    if( task == 0 ){
        return(A.WriteMFInfo(name));
    } else if( task == 0 ){
        return(B.WriteMFInfo(name));
    } else if( task == 0 ){
        return(C.WriteMFInfo(name));
    }  else if( task == 0 ){
        return(D.WriteMFInfo(name));
    }  else if( task == 0 ){
        return(E.WriteMFInfo(name));
    }  else {
        RUNTIME_ERROR("task out of legal range");
    }
}

//------------------------------------------------------------------------------

double CAUSEngineCST_A::GetLogML(void)
{
    double ml = 0.0;

    ml += A.GetLogML();  // <lam>
    ml += B.GetLogML();  // Cov(lam,Etot)
    ml += C.GetLogML();  // d<Etot>_FW
    ml += D.GetLogML();  // d<Etot>
    ml += E.GetLogML();  // -TdS{CST}corr

    return(ml);
}

//------------------------------------------------------------------------------

double CAUSEngineCST_A::GetLogPL(void)
{
    double loo = 0.0;

    loo += A.GetLogPL();  // <lam>
    loo += B.GetLogPL();  // Cov(lam,Etot)
    loo += C.GetLogPL();  // d<Etot>_FW
    loo += D.GetLogPL();  // d<Etot>
    loo += E.GetLogPL();  // -TdS{CST}corr

    return(loo);
}

//------------------------------------------------------------------------------

void CAUSEngineCST_A::GetLogMLDerivatives(const std::vector<bool>& flags,CSimpleVector<double>& der)
{
    // it must be re-map here
    std::vector<bool>   flags_a(A.GetNumOfHyprms());
    std::vector<bool>   flags_b(B.GetNumOfHyprms());
    std::vector<bool>   flags_c(C.GetNumOfHyprms());
    std::vector<bool>   flags_d(D.GetNumOfHyprms());
    std::vector<bool>   flags_e(E.GetNumOfHyprms());

    CSimpleVector<double>   der_a;
    der_a.CreateVector(A.GetNumOfHyprms());
    der_a.SetZero();
    CSimpleVector<double>   der_b;
    der_b.CreateVector(B.GetNumOfHyprms());
    der_b.SetZero();
    CSimpleVector<double>   der_c;
    der_c.CreateVector(C.GetNumOfHyprms());
    der_c.SetZero();
    CSimpleVector<double>   der_d;
    der_d.CreateVector(D.GetNumOfHyprms());
    der_d.SetZero();
    CSimpleVector<double>   der_e;
    der_e.CreateVector(E.GetNumOfHyprms());
    der_e.SetZero();

    int gidx = 0;
    int aidx = 0;
    int bidx = 0;
    int cidx = 0;
    int didx = 0;
    int eidx = 0;

    flags_a[aidx++] = flags[gidx++];
    flags_b[bidx++] = flags[gidx++];
    flags_c[cidx++] = flags[gidx++];
    flags_d[didx++] = flags[gidx++];
    flags_e[eidx++] = flags[gidx++];

    for(size_t cv=0; cv < NumOfCVs; cv++){
        bool flag = flags[gidx++];
        flags_a[aidx++] = flag;
        flags_b[bidx++] = flag;
        flags_c[cidx++] = flag;
        flags_d[didx++] = flag;
        flags_e[eidx++] = flag;
    }

    for(size_t cv=0; cv < NumOfCVs; cv++){
        flags_a[aidx++] = flags[gidx++];
    }
    for(size_t cv=0; cv < NumOfCVs; cv++){
        flags_b[bidx++] = flags[gidx++];
    }
    flags_c[cidx++] = flags[gidx++];
    flags_d[didx++] = flags[gidx++];
    flags_e[eidx++] = flags[gidx++];

    A.GetLogMLDerivatives(flags_a,der_a);  // <lam>
    B.GetLogMLDerivatives(flags_b,der_b);  // Cov(lam,Etot)
    C.GetLogMLDerivatives(flags_c,der_c);  // d<Etot>_FW
    D.GetLogMLDerivatives(flags_d,der_d);  // d<Etot>
    E.GetLogMLDerivatives(flags_e,der_e);  // -TdS{CST}corr

    // remap derivatives back
    gidx = 0;
    int lidx = 0;
    aidx = 0;
    bidx = 0;
    cidx = 0;
    didx = 0;
    eidx = 0;

    if( flags[gidx++] ) {
       der[lidx] = der[lidx] + der_a[aidx++];
       lidx++;
    }
    if( flags[gidx++] ) {
       der[lidx] = der[lidx] + der_b[bidx++];
       lidx++;
    }
    if( flags[gidx++] ) {
       der[lidx] = der[lidx] + der_c[cidx++];
       lidx++;
    }
    if( flags[gidx++] ) {
       der[lidx] = der[lidx] + der_d[didx++];
       lidx++;
    }
    if( flags[gidx++] ) {
       der[lidx] = der[lidx] + der_e[eidx++];
       lidx++;
    }

    for(size_t cv=0; cv < NumOfCVs; cv++){
        if( flags[gidx++] ) {
           der[lidx] = der[lidx] + der_a[aidx++];
           der[lidx] = der[lidx] + der_b[bidx++];
           der[lidx] = der[lidx] + der_c[cidx++];
           der[lidx] = der[lidx] + der_d[didx++];
           der[lidx] = der[lidx] + der_e[eidx++];
           lidx++;
        }
    }

    for(size_t cv=0; cv < NumOfCVs; cv++){
        if( flags[gidx++] ) {
           der[lidx] = der[lidx] + der_a[aidx++];
           lidx++;
        }
    }
    for(size_t cv=0; cv < NumOfCVs; cv++){
        if( flags[gidx++] ) {
           der[lidx] = der[lidx] + der_b[bidx++];
           lidx++;
        }
    }
    if( flags[gidx++] ) {
       der[lidx] = der[lidx] + der_c[cidx++];
       lidx++;
    }
    if( flags[gidx++] ) {
       der[lidx] = der[lidx] + der_d[didx++];
       lidx++;
    }
    if( flags[gidx++] ) {
       der[lidx] = der[lidx] + der_e[eidx++];
       lidx++;
    }
}

//------------------------------------------------------------------------------

void CAUSEngineCST_A::GetLogPLDerivatives(const std::vector<bool>& flags,CSimpleVector<double>& der)
{
// it must be re-map here
    std::vector<bool>   flags_a(A.GetNumOfHyprms());
    std::vector<bool>   flags_b(B.GetNumOfHyprms());
    std::vector<bool>   flags_c(C.GetNumOfHyprms());
    std::vector<bool>   flags_d(D.GetNumOfHyprms());
    std::vector<bool>   flags_e(E.GetNumOfHyprms());

    CSimpleVector<double>   der_a;
    der_a.CreateVector(A.GetNumOfHyprms());
    der_a.SetZero();
    CSimpleVector<double>   der_b;
    der_b.CreateVector(B.GetNumOfHyprms());
    der_b.SetZero();
    CSimpleVector<double>   der_c;
    der_c.CreateVector(C.GetNumOfHyprms());
    der_c.SetZero();
    CSimpleVector<double>   der_d;
    der_d.CreateVector(D.GetNumOfHyprms());
    der_d.SetZero();
    CSimpleVector<double>   der_e;
    der_e.CreateVector(E.GetNumOfHyprms());
    der_e.SetZero();

    int gidx = 0;
    int aidx = 0;
    int bidx = 0;
    int cidx = 0;
    int didx = 0;
    int eidx = 0;

    flags_a[aidx++] = flags[gidx++];
    flags_b[bidx++] = flags[gidx++];
    flags_c[cidx++] = flags[gidx++];
    flags_d[didx++] = flags[gidx++];
    flags_e[eidx++] = flags[gidx++];

    for(size_t cv=0; cv < NumOfCVs; cv++){
        bool flag = flags[gidx++];
        flags_a[aidx++] = flag;
        flags_b[bidx++] = flag;
        flags_c[cidx++] = flag;
        flags_d[didx++] = flag;
        flags_e[eidx++] = flag;
    }

    for(size_t cv=0; cv < NumOfCVs; cv++){
        flags_a[aidx++] = flags[gidx++];
    }
    for(size_t cv=0; cv < NumOfCVs; cv++){
        flags_b[bidx++] = flags[gidx++];
    }
    flags_c[cidx++] = flags[gidx++];
    flags_d[didx++] = flags[gidx++];
    flags_e[eidx++] = flags[gidx++];

    A.GetLogPLDerivatives(flags_a,der_a);  // <lam>
    B.GetLogPLDerivatives(flags_b,der_b);  // Cov(lam,Etot)
    C.GetLogPLDerivatives(flags_c,der_c);  // d<Etot>_FW
    D.GetLogPLDerivatives(flags_d,der_d);  // d<Etot>
    E.GetLogPLDerivatives(flags_e,der_e);  // -TdS{CST}corr

    // remap derivatives back
    gidx = 0;
    int lidx = 0;
    aidx = 0;
    bidx = 0;
    cidx = 0;
    didx = 0;
    eidx = 0;

    if( flags[gidx++] ) {
       der[lidx] = der[lidx] + der_a[aidx++];
       lidx++;
    }
    if( flags[gidx++] ) {
       der[lidx] = der[lidx] + der_b[bidx++];
       lidx++;
    }
    if( flags[gidx++] ) {
       der[lidx] = der[lidx] + der_c[cidx++];
       lidx++;
    }
    if( flags[gidx++] ) {
       der[lidx] = der[lidx] + der_d[didx++];
       lidx++;
    }
    if( flags[gidx++] ) {
       der[lidx] = der[lidx] + der_e[eidx++];
       lidx++;
    }

    for(size_t cv=0; cv < NumOfCVs; cv++){
        if( flags[gidx++] ) {
           der[lidx] = der[lidx] + der_a[aidx++];
           der[lidx] = der[lidx] + der_b[bidx++];
           der[lidx] = der[lidx] + der_c[cidx++];
           der[lidx] = der[lidx] + der_d[didx++];
           der[lidx] = der[lidx] + der_e[eidx++];
           lidx++;
        }
    }

    for(size_t cv=0; cv < NumOfCVs; cv++){
        if( flags[gidx++] ) {
           der[lidx] = der[lidx] + der_a[aidx++];
           lidx++;
        }
    }
    for(size_t cv=0; cv < NumOfCVs; cv++){
        if( flags[gidx++] ) {
           der[lidx] = der[lidx] + der_b[bidx++];
           lidx++;
        }
    }
    if( flags[gidx++] ) {
       der[lidx] = der[lidx] + der_c[cidx++];
       lidx++;
    }
    if( flags[gidx++] ) {
       der[lidx] = der[lidx] + der_d[didx++];
       lidx++;
    }
    if( flags[gidx++] ) {
       der[lidx] = der[lidx] + der_e[eidx++];
       lidx++;
    }
}

//------------------------------------------------------------------------------

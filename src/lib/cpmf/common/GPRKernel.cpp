// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2023 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2019 Petr Kulhanek, kulhanek@chemi.muni.cz
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

#include <GPRKernel.hpp>
#include <ErrorSystem.hpp>
#include <string>
#include <vector>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

//------------------------------------------------------------------------------

using namespace std;
using namespace boost;
using namespace boost::algorithm;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CGPRKernel::CGPRKernel(void)
{
    NumOfCVs            = 0;
    Kernel              = EGPRK_NONE;
    Accu                = NULL;
    UseNumDiff          = false;
    Alpha               = 10.0;
}

//------------------------------------------------------------------------------

CGPRKernel::~CGPRKernel(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CGPRKernel::SetAccumulator(CPMFAccumulatorPtr accu)
{
    if( accu == NULL ) return;                 // no-accu

    if( NumOfCVs == 0 ){
        NumOfCVs  = (size_t)accu->GetNumOfCVs();
        NumOfBins = (size_t)accu->GetNumOfBins();
    }

    if( NumOfCVs != (size_t)accu->GetNumOfCVs() ){
        RUNTIME_ERROR("inconsistent NumOfCVs");
    }
    if( NumOfBins != (size_t)accu->GetNumOfBins() ){
        RUNTIME_ERROR("inconsistent NumOfBins");
    }

    Accu = accu;
}

//------------------------------------------------------------------------------

void CGPRKernel::SetUseNumDiff(bool iset)
{
    UseNumDiff = iset;
}

//------------------------------------------------------------------------------

bool CGPRKernel::IsNumDiffEnabled(void)
{
    return(UseNumDiff);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CGPRKernel::SetKernel(const CSmallString& kernel)
{
    if( kernel == "ardse" ){
        Kernel = EGPRK_ARDSE;
    } else if( kernel == "ardmc52" ) {
        Kernel = EGPRK_ARDMC52;
    } else if( kernel == "ardmc32" ) {
        Kernel = EGPRK_ARDMC32;
    } else if( kernel == "ardmc12" ) {
        Kernel = EGPRK_ARDMC12;
    } else if( kernel == "ardrq" ) {
        Kernel = EGPRK_ARDRQ;
    } else if( kernel == "default" ) {
        Kernel = EGPRK_ARDSE;
    } else {
        CSmallString error;
        error << "Specified kernel '" << kernel << "' is not supported. "
                 "Supported kernels are: ardse (ARD squared exponential), ardmc52 (ARD Matern class 5/2), ardmc32 (ARD Matern class 3/2), ardmc12 (ARD Matern class 1/2), default(=ardse)";
        INVALID_ARGUMENT(error);
    }
}

//------------------------------------------------------------------------------

const CSmallString CGPRKernel::GetKernelName(void)
{
    switch(Kernel){
        case(EGPRK_ARDSE):
            return("ARD squared exponential");
        case(EGPRK_ARDMC52):
            return("ARD Matern class 5/2");
        case(EGPRK_ARDMC32):
            return("ARD Matern class 3/2");
        case(EGPRK_ARDMC12):
            return("ARD Matern class 1/2");
        case(EGPRK_ARDRQ):
            return("ARD rational quadratic");
        default:
            RUNTIME_ERROR("not implemented");
    }
}

//------------------------------------------------------------------------------

void CGPRKernel::SetupKernel(void)
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }
    if( NumOfCVs == 0 ){
        RUNTIME_ERROR("NumOfCVs is NULL");
    }

    // GPR setup
    CVLengths2.CreateVector(NumOfCVs);
    for(size_t i=0; i < NumOfCVs; i++){
        double l = WFac[i]*Accu->GetCV(i)->GetRange()/Accu->GetCV(i)->GetNumOfBins();
        CVLengths2[i] = l*l;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CGPRKernel::SetWFac(const CSmallString& spec)
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is not set!");
    }
    if( NumOfCVs == 0 ){
        RUNTIME_ERROR("NumOfCVs is zero");
    }

    string          sspec(spec);
    vector<string>  swfacs;

    split(swfacs,sspec,is_any_of("x"),token_compress_on);

    if( swfacs.size() > NumOfCVs ){
        CSmallString error;
        error << "too many wfacs (" << swfacs.size() << ") than required (" << NumOfCVs << ")";
        RUNTIME_ERROR(error);
    }

    WFac.CreateVector(NumOfCVs);

    // parse values of wfac
    double last_wfac = 3.0;
    for(size_t i=0; i < swfacs.size(); i++){
        stringstream str(swfacs[i]);
        str >> last_wfac;
        if( ! str ){
            CSmallString error;
            error << "unable to decode wfac value for position: " << i+1;
            RUNTIME_ERROR(error);
        }
        WFac[i] = last_wfac;
    }

    // pad the rest with the last value
    for(size_t i=swfacs.size(); i < NumOfCVs; i++){
        WFac[i] = last_wfac;
    }
}

//------------------------------------------------------------------------------

void CGPRKernel::SetWFac(CSimpleVector<double>& wfac)
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is not set!");
    }
    if( NumOfCVs == 0 ){
        RUNTIME_ERROR("NumOfCVs is zero");
    }

    if( wfac.GetLength() != NumOfCVs ){
        RUNTIME_ERROR("ncvs inconsistent in the source and target");
    }

    WFac = wfac;
}

//------------------------------------------------------------------------------

void CGPRKernel::SetWFac(size_t cvind, double value)
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is not set!");
    }
    if( NumOfCVs == 0 ){
        RUNTIME_ERROR("NumOfCVs is zero");
    }

    if( cvind >= NumOfCVs ){
        RUNTIME_ERROR("cvind out-of-range");
    }
    // is wfac initialized?
    if( WFac.GetLength() == 0 ){
        WFac.CreateVector(NumOfCVs);
    }
    WFac[cvind] = value;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

double CGPRKernel::GetKernelValue(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp)
{
    // calculate scaled distance
    double scdist2 = 0.0;
    for(size_t ii=0; ii < NumOfCVs; ii++){
        double du = Accu->GetCV(ii)->GetDifference(ip[ii],jp[ii]);
        double dd = CVLengths2[ii];
        scdist2 += du*du/dd;
    }

    // get kernel value
    switch(Kernel){
        case(EGPRK_ARDSE): {
            return(exp(-0.5*scdist2));
        }
        break;
    // -----------------------
        case(EGPRK_ARDMC52):{
            double scdist = sqrt(scdist2);
            return((1.0+sqrt(5.0)*scdist+(5.0/3.0)*scdist2)*exp(-sqrt(5.0)*scdist));
        }
        break;
    // -----------------------
        case(EGPRK_ARDMC32):{
            double scdist = sqrt(scdist2);
            return((1.0+sqrt(3.0)*scdist)*exp(-sqrt(3.0)*scdist));
        }
        break;
    // -----------------------
        case(EGPRK_ARDMC12):{
            double scdist = sqrt(scdist2);
            return(exp(-scdist));
        }
        break;
    // -----------------------
        case(EGPRK_ARDRQ):{
            double in = 1.0 + 0.5*scdist2/Alpha;
            return(pow(in,-Alpha));
        }
        break;
    // -----------------------
        default:
            RUNTIME_ERROR("not implemented");
        break;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

double CGPRKernel::GetKernelIntI(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp)
{
    if( NumOfCVs != 1 ){
        RUNTIME_ERROR("not implemented");
    }

    double kint = 0.0;

    // get kernel value
    switch(Kernel){
    case(EGPRK_ARDSE):{
        double du = Accu->GetCV(0)->GetDifference(ip[0],jp[0]);
        double dd = 2.0*CVLengths2[0];
        kint = 0.5*sqrt(M_PI)*sqrt(dd)*erf(du/sqrt(dd));
        }
        break;
    default:
        RUNTIME_ERROR("not implemented");
        break;
    }

    return(kint);
}

//------------------------------------------------------------------------------

double CGPRKernel::GetKernelIntIJ(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp)
{
    if( NumOfCVs != 1 ){
        RUNTIME_ERROR("not implemented");
    }

    double kint = 0.0;

    // get kernel value
    switch(Kernel){
    case(EGPRK_ARDSE):{
        double du = Accu->GetCV(0)->GetDifference(ip[0],jp[0]);
        double dd = 2.0*CVLengths2[0];
        kint = -0.5*sqrt(M_PI)*sqrt(dd)*du*erf(du/sqrt(dd))-0.5*dd*exp(-du*du/dd)+1.5;
        }
        break;
    default:
        RUNTIME_ERROR("not implemented");
        break;
    }

    return(kint);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CGPRKernel::GetKernelDerI(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,CSimpleVector<double>& kder)
{
    if( UseNumDiff ){
        GetKernelDerINum(ip,jp,kder);
    } else {
        GetKernelDerIAna(ip,jp,kder);
    }
}

//------------------------------------------------------------------------------

void CGPRKernel::GetKernelDerJ(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,CSimpleVector<double>& kder)
{
    if( UseNumDiff ){
        GetKernelDerJNum(ip,jp,kder);
    } else {
        GetKernelDerJAna(ip,jp,kder);
    }
}

//------------------------------------------------------------------------------

void CGPRKernel::GetKernelDerIJ(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,CFortranMatrix& kblock)
{
    if( UseNumDiff ){
        GetKernelDerIJNum(ip,jp,kblock);
    } else {
        GetKernelDerIJAna(ip,jp,kblock);
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CGPRKernel::GetKernelDerIAna(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,CSimpleVector<double>& kder)
{
    // calculate scaled distance
    double scdist2 = 0.0;
    for(size_t ii=0; ii < NumOfCVs; ii++){
        double du = Accu->GetCV(ii)->GetDifference(ip[ii],jp[ii]);
        double dd = CVLengths2[ii];
        scdist2 += du*du/dd;
    }

    // get kernel value
    switch(Kernel){
    case(EGPRK_ARDSE):{
            double pre = exp(-0.5*scdist2);
            for(size_t ii=0; ii < NumOfCVs; ii++){
                double du = Accu->GetCV(ii)->GetDifference(ip[ii],jp[ii]);
                double dd = CVLengths2[ii];
                kder[ii] = -pre*du/dd;
            }
        }
        break;
    case(EGPRK_ARDMC52):{
            double scdist = sqrt(scdist2);
            double pre = -(5.0/3.0)*exp(-sqrt(5.0)*scdist)*(sqrt(5.0)*scdist+1.0);
            for(size_t ii=0; ii < NumOfCVs; ii++){
                double du = Accu->GetCV(ii)->GetDifference(ip[ii],jp[ii]);
                double dd = CVLengths2[ii];
                kder[ii] = pre*du/dd;
            }
        }
        break;
    case(EGPRK_ARDRQ):{
            double in = 1.0 + 0.5*scdist2/Alpha;
            double pre = pow(in,-(Alpha+1.0));
            for(size_t ii=0; ii < NumOfCVs; ii++){
                double du = Accu->GetCV(ii)->GetDifference(ip[ii],jp[ii]);
                double dd = CVLengths2[ii];
                kder[ii] = -pre*du/dd;
            }
        }
        break;
    default:
        RUNTIME_ERROR("not implemented");
        break;
    }
}

//------------------------------------------------------------------------------

void CGPRKernel::GetKernelDerINum(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,CSimpleVector<double>& kder)
{
    CSimpleVector<double> tip;
    tip.CreateVector(NumOfCVs);

    double  dh = 1e-3;
    double  v1,v2;

    // off diagonal elements
    for(size_t ii=0; ii < NumOfCVs; ii++) {

        tip = ip;
        tip[ii] = ip[ii] - dh;
        v1 = GetKernelValue(tip,jp);

        tip = ip;
        tip[ii] = ip[ii] + dh;
        v2 = GetKernelValue(tip,jp);

        kder[ii] = (v2 - v1)/(2.0*dh);
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CGPRKernel::GetKernelDerJAna(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,CSimpleVector<double>& kder)
{
    // calculate scaled distance
    double scdist2 = 0.0;
    for(size_t ii=0; ii < NumOfCVs; ii++){
        double du = Accu->GetCV(ii)->GetDifference(ip[ii],jp[ii]);
        double dd = CVLengths2[ii];
        scdist2 += du*du/dd;
    }

// get derivative
    switch(Kernel){
    case(EGPRK_ARDSE):{
            double pre = exp(-0.5*scdist2);
            for(size_t ii=0; ii < NumOfCVs; ii++){
                double du = Accu->GetCV(ii)->GetDifference(ip[ii],jp[ii]);
                double dd = CVLengths2[ii];
                kder[ii] = pre*du/dd;
            }
        }
        break;
    case(EGPRK_ARDMC52):{
            double scdist = sqrt(scdist2);
            double pre = -(5.0/3.0)*exp(-sqrt(5.0)*scdist)*(sqrt(5.0)*scdist+1.0);
            for(size_t ii=0; ii < NumOfCVs; ii++){
                double du = Accu->GetCV(ii)->GetDifference(ip[ii],jp[ii]);
                double dd = CVLengths2[ii];
                kder[ii] = -pre*du/dd;
            }
        }
        break;
    case(EGPRK_ARDRQ):{
            double in = 1.0 + 0.5*scdist2/Alpha;
            double pre = pow(in,-(Alpha+1.0));
            for(size_t ii=0; ii < NumOfCVs; ii++){
                double du = Accu->GetCV(ii)->GetDifference(ip[ii],jp[ii]);
                double dd = CVLengths2[ii];
                kder[ii] = pre*du/dd;
            }
        }
        break;
    default:
        RUNTIME_ERROR("not implemented");
        break;
    }
}

//------------------------------------------------------------------------------

void CGPRKernel::GetKernelDerJNum(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,CSimpleVector<double>& kder)
{
    CSimpleVector<double> tjp;
    tjp.CreateVector(NumOfCVs);

    double  dh = 1e-3;
    double  v1,v2;

    // off diagonal elements
    for(size_t ii=0; ii < NumOfCVs; ii++) {

        tjp = jp;
        tjp[ii] = jp[ii] - dh;
        v1 = GetKernelValue(ip,tjp);

        tjp = jp;
        tjp[ii] = jp[ii] + dh;
        v2 = GetKernelValue(ip,tjp);

        kder[ii] = (v2 - v1)/(2.0*dh);
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CGPRKernel::GetKernelDerIJAna(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,CFortranMatrix& kblock)
{
    kblock.SetZero();

// calculate scaled distance
    double scdist2 = 0.0;
    for(size_t ii=0; ii < NumOfCVs; ii++){
        double du = Accu->GetCV(ii)->GetDifference(ip[ii],jp[ii]);
        double dd = CVLengths2[ii];
        scdist2 += du*du/dd;
    }

// get derivative
    switch(Kernel){
    case(EGPRK_ARDSE): {
            double pre = exp(-0.5*scdist2);
            for(size_t ii=0; ii < NumOfCVs; ii++){
                for(size_t jj=0; jj < NumOfCVs; jj++){
                    double du = Accu->GetCV(ii)->GetDifference(ip[ii],jp[ii]) *
                                Accu->GetCV(jj)->GetDifference(ip[jj],jp[jj]);
                    double dd = CVLengths2[ii]*CVLengths2[jj];
                    kblock[ii][jj] -= pre*du/dd;
                    if( ii == jj ){
                        kblock[ii][ii] += pre/CVLengths2[ii];
                    }
                }
            }
        }
        break;
    case(EGPRK_ARDMC52):{
            double scdist = sqrt(scdist2);
            double pr = -(5.0/3.0)*exp(-sqrt(5.0)*scdist);
            double d1 = pr*(sqrt(5.0)*scdist+1.0);
            for(size_t ii=0; ii < NumOfCVs; ii++){
                for(size_t jj=0; jj < NumOfCVs; jj++){
                    double du = Accu->GetCV(ii)->GetDifference(ip[ii],jp[ii]) *
                                Accu->GetCV(jj)->GetDifference(ip[jj],jp[jj]);
                    double dd = CVLengths2[ii]*CVLengths2[jj];
                    kblock[ii][jj] += 5.0*pr*du/dd;
                    if( (ii == jj) ){
                        kblock[ii][jj] -= d1/(CVLengths2[ii]);
                    }
                }
            }
        }
        break;
    case(EGPRK_ARDRQ): {
            double in = 1.0 + 0.5*scdist2/Alpha;
            double pre  = pow(in,-(Alpha+2.0));
            double pre1 = pow(in,-(Alpha+1.0));
            for(size_t ii=0; ii < NumOfCVs; ii++){
                for(size_t jj=0; jj < NumOfCVs; jj++){
                    double du = Accu->GetCV(ii)->GetDifference(ip[ii],jp[ii]) *
                                Accu->GetCV(jj)->GetDifference(ip[jj],jp[jj]);
                    double dd = CVLengths2[ii]*CVLengths2[jj];
                    kblock[ii][jj] -= (Alpha+1.0)*pre*du/(dd*Alpha);
                    if( ii == jj ){
                        kblock[ii][ii] += pre1/CVLengths2[ii];
                    }
                }
            }
        }
        break;
    default:
        RUNTIME_ERROR("not implemented");
        break;
    }
}

//------------------------------------------------------------------------------

void CGPRKernel::GetKernelDerIJNum(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,CFortranMatrix& kblock)
{
    CSimpleVector<double> tip;
    tip.CreateVector(NumOfCVs);
    CSimpleVector<double> tjp;
    tjp.CreateVector(NumOfCVs);

    double  dh = 1e-4;
    double  v1,v2,v3,v4;

    // off diagonal elements
    for(size_t ii=0; ii < NumOfCVs; ii++) {
        for(size_t jj=0; jj < NumOfCVs; jj++) {

            tip = ip;
            tip[ii] = ip[ii] + dh;
            tjp = jp;
            tjp[jj] = jp[jj] + dh;
            v1 = GetKernelValue(tip,tjp);

            tip = ip;
            tip[ii] = ip[ii] + dh;
            tjp = jp;
            tjp[jj] = jp[jj] - dh;
            v2 = GetKernelValue(tip,tjp);

            tip = ip;
            tip[ii] = ip[ii] - dh;
            tjp = jp;
            tjp[jj] = jp[jj] + dh;
            v3 = GetKernelValue(tip,tjp);

            tip = ip;
            tip[ii] = ip[ii] - dh;
            tjp = jp;
            tjp[jj] = jp[jj] - dh;
            v4 = GetKernelValue(tip,tjp);

            kblock[ii][jj] = (v1 - v2 - v3 + v4)/(4.0*dh*dh);
        }
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

double CGPRKernel::GetKernelValueWFacDer(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,size_t cv)
{
    if( UseNumDiff ){
        return( GetKernelValueWFacDerNum(ip,jp,cv) );
    } else {
        return( GetKernelValueWFacDerAna(ip,jp,cv) );
    }
}

//------------------------------------------------------------------------------

void CGPRKernel::GetKernelDerIWFacDer(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,size_t cv,CSimpleVector<double>& kder)
{
    if( UseNumDiff ){
        return( GetKernelDerIWFacDerNum(ip,jp,cv,kder) );
    } else {
        return( GetKernelDerIWFacDerAna(ip,jp,cv,kder) );
    }
}

//------------------------------------------------------------------------------

void CGPRKernel::GetKernelDerJWFacDer(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,size_t cv,CSimpleVector<double>& kder)
{
    if( UseNumDiff ){
        return( GetKernelDerJWFacDerNum(ip,jp,cv,kder) );
    } else {
        return( GetKernelDerJWFacDerAna(ip,jp,cv,kder) );
    }
}

//------------------------------------------------------------------------------

void CGPRKernel::GetKernelDerIJWFacDer(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,size_t cv,CFortranMatrix& kblock)
{
    if( UseNumDiff ){
        return( GetKernelDerIJWFacDerNum(ip,jp,cv,kblock) );
    } else {
        return( GetKernelDerIJWFacDerAna(ip,jp,cv,kblock) );
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

double CGPRKernel::GetKernelValueWFacDerAna(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,size_t cv)
{
// calculate scaled distance
    double scdist2 = 0.0;
    for(size_t ii=0; ii < NumOfCVs; ii++){
        double du = Accu->GetCV(ii)->GetDifference(ip[ii],jp[ii]);
        double dd = CVLengths2[ii];
        scdist2 += du*du/dd;
    }

// get derivative
    switch(Kernel){
        case(EGPRK_ARDSE): {
                double val = exp(-0.5*scdist2);
                double du = Accu->GetCV(cv)->GetDifference(ip[cv],jp[cv]);
                double dd = CVLengths2[cv];
                double wf = WFac[cv];
                double der = val * du*du / (dd * wf);
                return(der);
            }
            break;
    // -----------------------
        case(EGPRK_ARDMC52): {
            // possible division by zero is solved in CalcKderWRTWFac
            double scdist = sqrt(scdist2);
            double vex = exp(-sqrt(5.0)*scdist);
            double val = vex*(1.0+sqrt(5.0)*scdist+(5.0/3.0)*scdist2);
            double du = Accu->GetCV(cv)->GetDifference(ip[cv],jp[cv]);
            double dd = CVLengths2[cv];
            double wf = WFac[cv];
            double dexp =   val * sqrt(5.0) * du*du / (dd * wf) / scdist;       // exponential part
            double dpol = - vex * ( sqrt(5.0) * du*du / (dd * wf) / scdist + 2.0 * (5.0/3.0)* du*du / (dd * wf) );   // polynomial part
            return(dexp+dpol);
            }
            break;
    // -----------------------
        case(EGPRK_ARDMC32): {
            // possible division by zero is solved in CalcKderWRTWFac
            double scdist = sqrt(scdist2);
            double vex = exp(-sqrt(3.0)*scdist);
            double val = vex*(1.0+sqrt(3.0)*scdist);
            double du = Accu->GetCV(cv)->GetDifference(ip[cv],jp[cv]);
            double dd = CVLengths2[cv];
            double wf = WFac[cv];
            double dexp =   val * sqrt(3.0) * du*du / (dd * wf) / scdist;   // exponential part
            double dpol = - vex * sqrt(3.0) * du*du / (dd * wf) / scdist;   // polynomial part
            return(dexp+dpol);
            }
            break;
    // -----------------------
        case(EGPRK_ARDMC12): {
            // possible division by zero is solved in CalcKderWRTWFac
            double scdist = sqrt(scdist2);
            double val = exp(-scdist);
            double du = Accu->GetCV(cv)->GetDifference(ip[cv],jp[cv]);
            double dd = CVLengths2[cv];
            double wf = WFac[cv];
            double der = val * du*du / (dd * wf) / scdist;
            return(der);
            }
            break;
    // -----------------------
        case(EGPRK_ARDRQ): {
            double in  = 1.0 + 0.5*scdist2/Alpha;
            double pre = pow(in,-(Alpha+1.0));
            double du = Accu->GetCV(cv)->GetDifference(ip[cv],jp[cv]);
            double dd = CVLengths2[cv];
            double wf = WFac[cv];
            double der = pre * du*du / (2.0*dd*dd*wf);
            return(der);
            }
            break;
    // -----------------------
        default:
            RUNTIME_ERROR("not implemented");
            break;
    }
}

//------------------------------------------------------------------------------

double CGPRKernel::GetKernelValueWFacDerNum(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,size_t cv)
{
    double kval1,kval2;

    double  dh = 1e-3;

    for(size_t i=0; i < NumOfCVs; i++){
        double l;
        if( i == cv ){
            l = (WFac[i]-dh)*Accu->GetCV(i)->GetRange()/Accu->GetCV(i)->GetNumOfBins();
        } else {
            l = WFac[i]*Accu->GetCV(i)->GetRange()/Accu->GetCV(i)->GetNumOfBins();
        }
        CVLengths2[i] = l*l;
    }
    kval1 = GetKernelValue(ip,jp);

    for(size_t i=0; i < NumOfCVs; i++){
        double l;
        if( i == cv ){
            l = (WFac[i]+dh)*Accu->GetCV(i)->GetRange()/Accu->GetCV(i)->GetNumOfBins();
        } else {
            l = WFac[i]*Accu->GetCV(i)->GetRange()/Accu->GetCV(i)->GetNumOfBins();
        }
        CVLengths2[i] = l*l;
    }
    kval2 = GetKernelValue(ip,jp);

    double kder = (kval2-kval1)/(2.0*dh);

    // restore original CVLengths2
    for(size_t i=0; i < NumOfCVs; i++){
        double l = WFac[i]*Accu->GetCV(i)->GetRange()/Accu->GetCV(i)->GetNumOfBins();
        CVLengths2[i] = l*l;
    }

    return(kder);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CGPRKernel::GetKernelDerIWFacDerAna(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,size_t cv,CSimpleVector<double>& kder)
{
// calculate scaled distance
    double scdist2 = 0.0;
    for(size_t ii=0; ii < NumOfCVs; ii++){
        double du = Accu->GetCV(ii)->GetDifference(ip[ii],jp[ii]);
        double dd = CVLengths2[ii];
        scdist2 += du*du/dd;
    }

// get derivative
    switch(Kernel){
    case(EGPRK_ARDSE):{
            double pre = exp(-0.5*scdist2);
            double du  = Accu->GetCV(cv)->GetDifference(ip[cv],jp[cv]);
            double dd  = CVLengths2[cv];
            double wf  = WFac[cv];
            // -1/2 * du*du * wfac^-2*dbin-^2
            // -1/2 * du*du * (-2)*wfac^-3 * dbin^-2
            //        du*du * wfac^-2 * dbin^-2 * wfac^-1
            //        du*du / (dd * wfac)
            double der = pre * du*du / (dd * wf);

            for(size_t ii=0; ii < NumOfCVs; ii++){
                double du = Accu->GetCV(ii)->GetDifference(ip[ii],jp[ii]);
                double dd = CVLengths2[ii];
                // kder[ii] = -pre*du/dd;
                // [(-pre*du)'*(dd) - (-pre*du)*(dd)']/(dd*dd)
                if( cv == ii ){
                    kder[ii] = (-der*du*dd - (-pre*du*(2.0*dd/wf)) )/(dd*dd);
                } else {
                    kder[ii] = -der*du/dd;
                }
            }
        }
        break;
    // -----------------------
        case(EGPRK_ARDRQ): {
            double in   = 1.0 + 0.5*scdist2/Alpha;
            double pre  = pow(in,-(Alpha+2.0));
            double pre1 = pow(in,-(Alpha+1.0));
            double du = Accu->GetCV(cv)->GetDifference(ip[cv],jp[cv]);
            double dd = CVLengths2[cv];
            double wf = WFac[cv];
            double der = -(Alpha+1.0)*pre*du*du/(Alpha*dd*wf);
            for(size_t ii=0; ii < NumOfCVs; ii++){
                double du = Accu->GetCV(ii)->GetDifference(ip[ii],jp[ii]);
                double dd = CVLengths2[ii];
                if( cv == ii ){
                    kder[ii] = 2.0*du*pre1/(dd*wf) - (Alpha+1.0)*du*du*du*pre/(Alpha*dd*dd*wf);
                } else {
                    kder[ii] = der*du/dd;
                }
            }
        }
        break;
    default:
        RUNTIME_ERROR("not implemented");
        break;
    }
}

//------------------------------------------------------------------------------

void CGPRKernel::GetKernelDerIWFacDerNum(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,size_t cv,CSimpleVector<double>& kder)
{
    CSimpleVector<double> kder1,kder2;

    kder1.CreateVector(NumOfCVs);
    kder2.CreateVector(NumOfCVs);

    double  dh = 1e-3;

    for(size_t i=0; i < NumOfCVs; i++){
        double l;
        if( i == cv ){
            l = (WFac[i]-dh)*Accu->GetCV(i)->GetRange()/Accu->GetCV(i)->GetNumOfBins();
        } else {
            l = WFac[i]*Accu->GetCV(i)->GetRange()/Accu->GetCV(i)->GetNumOfBins();
        }
        CVLengths2[i] = l*l;
    }
    GetKernelDerIAna(ip,jp,kder1);

    for(size_t i=0; i < NumOfCVs; i++){
        double l;
        if( i == cv ){
            l = (WFac[i]+dh)*Accu->GetCV(i)->GetRange()/Accu->GetCV(i)->GetNumOfBins();
        } else {
            l = WFac[i]*Accu->GetCV(i)->GetRange()/Accu->GetCV(i)->GetNumOfBins();
        }
        CVLengths2[i] = l*l;
    }
    GetKernelDerIAna(ip,jp,kder2);

    for(size_t ii=0; ii < NumOfCVs; ii++){
        kder[ii] = (kder2[ii]-kder1[ii])/(2.0*dh);
    }

    // restore original CVLengths2
    for(size_t i=0; i < NumOfCVs; i++){
        double l = WFac[i]*Accu->GetCV(i)->GetRange()/Accu->GetCV(i)->GetNumOfBins();
        CVLengths2[i] = l*l;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CGPRKernel::GetKernelDerJWFacDerAna(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,size_t cv,CSimpleVector<double>& kder)
{
// calculate scaled distance
    double scdist2 = 0.0;
    for(size_t ii=0; ii < NumOfCVs; ii++){
        double du = Accu->GetCV(ii)->GetDifference(ip[ii],jp[ii]);
        double dd = CVLengths2[ii];
        scdist2 += du*du/dd;
    }

// get derivative
    switch(Kernel){
    case(EGPRK_ARDSE):{
            double pre = exp(-0.5*scdist2);
            double du  = Accu->GetCV(cv)->GetDifference(ip[cv],jp[cv]);
            double dd  = CVLengths2[cv];
            double wf  = WFac[cv];
            // -1/2 * du*du * wfac^-2*dbin-^2
            // -1/2 * du*du * (-2)*wfac^-3 * dbin^-2
            //        du*du * wfac^-2 * dbin^-2 * wfac^-1
            //        du*du / (dd * wfac)
            double der = pre * du*du / (dd * wf);

            for(size_t ii=0; ii < NumOfCVs; ii++){
                double du = Accu->GetCV(ii)->GetDifference(ip[ii],jp[ii]);
                double dd = CVLengths2[ii];
                // kder[ii] = pre*du/dd;
                // [(pre*du)'*(dd) - (pre*du)*(dd)']/(dd*dd)
                if( cv == ii ){
                    kder[ii] = (der*du*dd - (pre*du*(2.0*dd/wf)) )/(dd*dd);
                } else {
                    kder[ii] = der*du/dd;
                }
            }
        }
        break;
    // -----------------------
        case(EGPRK_ARDRQ): {
            double in   = 1.0 + 0.5*scdist2/Alpha;
            double pre  = pow(in,-(Alpha+2.0));
            double pre1 = pow(in,-(Alpha+1.0));
            double du = Accu->GetCV(cv)->GetDifference(ip[cv],jp[cv]);
            double dd = CVLengths2[cv];
            double wf = WFac[cv];
            double der = -(Alpha+1.0)*pre*du*du/(Alpha*dd*wf);
            for(size_t ii=0; ii < NumOfCVs; ii++){
                double du = Accu->GetCV(ii)->GetDifference(ip[ii],jp[ii]);
                double dd = CVLengths2[ii];
                if( cv == ii ){
                    kder[ii] = -2.0*du*pre1/(dd*wf) + (Alpha+1.0)*du*du*du*pre/(Alpha*dd*dd*wf);
                } else {
                    kder[ii] = -der*du/dd;
                }
            }
        }
        break;
    default:
        RUNTIME_ERROR("not implemented");
        break;
    }
}
//------------------------------------------------------------------------------

void CGPRKernel::GetKernelDerJWFacDerNum(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,size_t cv,CSimpleVector<double>& kder)
{
CSimpleVector<double> kder1,kder2;

    kder1.CreateVector(NumOfCVs);
    kder2.CreateVector(NumOfCVs);

    double  dh = 1e-3;

    for(size_t i=0; i < NumOfCVs; i++){
        double l;
        if( i == cv ){
            l = (WFac[i]-dh)*Accu->GetCV(i)->GetRange()/Accu->GetCV(i)->GetNumOfBins();
        } else {
            l = WFac[i]*Accu->GetCV(i)->GetRange()/Accu->GetCV(i)->GetNumOfBins();
        }
        CVLengths2[i] = l*l;
    }
    GetKernelDerJAna(ip,jp,kder1);

    for(size_t i=0; i < NumOfCVs; i++){
        double l;
        if( i == cv ){
            l = (WFac[i]+dh)*Accu->GetCV(i)->GetRange()/Accu->GetCV(i)->GetNumOfBins();
        } else {
            l = WFac[i]*Accu->GetCV(i)->GetRange()/Accu->GetCV(i)->GetNumOfBins();
        }
        CVLengths2[i] = l*l;
    }
    GetKernelDerJAna(ip,jp,kder2);

    for(size_t ii=0; ii < NumOfCVs; ii++){
        kder[ii] = (kder2[ii]-kder1[ii])/(2.0*dh);
    }

    // restore original CVLengths2
    for(size_t i=0; i < NumOfCVs; i++){
        double l = WFac[i]*Accu->GetCV(i)->GetRange()/Accu->GetCV(i)->GetNumOfBins();
        CVLengths2[i] = l*l;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CGPRKernel::GetKernelDerIJWFacDerAna(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,size_t cv,CFortranMatrix& kblock)
{
// calculate scaled distance
    double scdist2 = 0.0;
    for(size_t ii=0; ii < NumOfCVs; ii++){
        double du = Accu->GetCV(ii)->GetDifference(ip[ii],jp[ii]);
        double dd = CVLengths2[ii];
        scdist2 += du*du/dd;
    }

    kblock.SetZero();

// get derivative
    double wf = WFac[cv];
    double wd3 = 1.0/(CVLengths2[cv]*wf);
    double wd5 = wd3/CVLengths2[cv];
    double dc = Accu->GetCV(cv)->GetDifference(ip[cv],jp[cv]);

    switch(Kernel){
    case(EGPRK_ARDSE): {
            double arg = exp(-0.5*scdist2);
            double argd = arg*dc*dc*wd3;
            for(size_t ii=0; ii < NumOfCVs; ii++){
                for(size_t jj=0; jj < NumOfCVs; jj++){
                    double du = Accu->GetCV(ii)->GetDifference(ip[ii],jp[ii]) *
                                Accu->GetCV(jj)->GetDifference(ip[jj],jp[jj]);
                    double dd = CVLengths2[ii]*CVLengths2[jj];
                    kblock[ii][jj] -= argd*du/dd;
                    if( (cv == ii) && (cv != jj) ) {
                        kblock[ii][jj] += 2.0*arg*du*wd3/(CVLengths2[jj]);
                    } else if( (cv == jj) && (cv != ii) ) {
                        kblock[ii][jj] += 2.0*arg*du*wd3/(CVLengths2[ii]);
                    } else if( (cv == ii) && (cv == jj) ) {
                        kblock[ii][jj] += 4.0*arg*du*wd5;
                    }
                    if( (ii == jj) ){
                        kblock[ii][ii] += argd/CVLengths2[ii];
                        if( ii == cv ){
                            kblock[ii][ii] -= 2.0*arg*wd3;
                        }
                    }
                }
            }
        }
        break;
    case(EGPRK_ARDMC52):{
            double scdist = sqrt(scdist2);
            double pr = -(5.0/3.0)*exp(-sqrt(5.0)*scdist);
            double prd = 0.0;
            if( scdist > 0 ){
               prd = -(5.0/3.0)*sqrt(5.0)*exp(-sqrt(5.0)*scdist)*(1.0/scdist)*wd3*dc*dc;
            }
            double d1 = pr*(sqrt(5.0)*scdist+1.0);
            double d1d = -(25.0/3.0)*exp(-sqrt(5.0)*scdist)*wd3*dc*dc;
            for(size_t ii=0; ii < NumOfCVs; ii++){
                for(size_t jj=0; jj < NumOfCVs; jj++){
                    double du = Accu->GetCV(ii)->GetDifference(ip[ii],jp[ii]) *
                                Accu->GetCV(jj)->GetDifference(ip[jj],jp[jj]);
                    double dd = CVLengths2[ii]*CVLengths2[jj];
                    kblock[ii][jj] += 5.0*prd*du/dd;
                    if( (cv == ii) && (cv != jj) ) {
                        kblock[ii][jj] -= 10.0*pr*du*wd3/(CVLengths2[jj]);
                    } else if( (cv == jj) && (cv != ii) ) {
                        kblock[ii][jj] -= 10.0*pr*du*wd3/(CVLengths2[ii]);
                    } else if( (cv == ii) && (cv == jj) ) {
                        kblock[ii][jj] -= 20.0*pr*du*wd5;
                    }
                    if( (ii == jj) ){
                        kblock[ii][ii] -= d1d/(CVLengths2[ii]);
                        if( ii == cv ){
                            kblock[ii][ii] += 2.0*d1*wd3;
                        }

                    }
                }
            }
        }
        break;
    case(EGPRK_ARDRQ): {
            double in   = 1.0 + 0.5*scdist2/Alpha;
            double pre1 = pow(in,-(Alpha+1.0));
            double pre2 = pow(in,-(Alpha+2.0));
            double pre3 = pow(in,-(Alpha+3.0));

            double du = Accu->GetCV(cv)->GetDifference(ip[cv],jp[cv]);
            double dd = CVLengths2[cv];

            for(size_t ii=0; ii < NumOfCVs; ii++){
                for(size_t jj=0; jj < NumOfCVs; jj++){
                    double duii = Accu->GetCV(ii)->GetDifference(ip[ii],jp[ii]);
                    double dujj = Accu->GetCV(jj)->GetDifference(ip[jj],jp[jj]);
                    double ddii = CVLengths2[ii];
                    double ddjj = CVLengths2[jj];
                           if( (cv == ii) && (cv != jj) ) {
                        kblock[ii][jj] = -(Alpha+2.0)*(Alpha+1.0)*du*du*du*dujj*pre3/(Alpha*Alpha*dd*dd*wf*ddjj)
                                         +2.0*(Alpha+1.0)*du*dujj*pre2/(Alpha*dd*wf*ddjj);
                    } else if( (cv != ii) && (cv == jj) ) {
                        kblock[ii][jj] = -(Alpha+2.0)*(Alpha+1.0)*du*du*du*duii*pre3/(Alpha*Alpha*dd*dd*wf*ddii)
                                         +2.0*(Alpha+1.0)*du*duii*pre2/(Alpha*dd*wf*ddii);
                    } else if( (cv == ii) && (cv == jj) ) {
                        kblock[ii][jj] = -(Alpha+2.0)*(Alpha+1.0)*du*du*du*du*pre3/(Alpha*Alpha*dd*dd*dd*wf)
                                         +5.0*(Alpha+1.0)*du*du*pre2/(Alpha*dd*dd*wf)-2.0*pre1/(dd*wf);
                    } else {
                        kblock[ii][jj] = -(Alpha+2.0)*(Alpha+1.0)*du*du*dujj*dujj*pre3/(Alpha*Alpha*dd*wf*ddjj*ddjj)
                                         +(Alpha+1.0)*du*du*pre2/(Alpha*dd*wf*ddjj);
                    }
                }
            }
        }
        break;
    default:
        RUNTIME_ERROR("not implemented");
        break;
    }
}

//------------------------------------------------------------------------------

// for internal debugging
void CGPRKernel::GetKernelDerIJWFacDerNum(const CSimpleVector<double>& ip,const CSimpleVector<double>& jp,size_t cv,CFortranMatrix& kblock)
{
    CFortranMatrix kblock1,kblock2;

    kblock1.CreateMatrix(NumOfCVs,NumOfCVs);
    kblock2.CreateMatrix(NumOfCVs,NumOfCVs);

    double  dh = 1e-3;

    for(size_t i=0; i < NumOfCVs; i++){
        double l;
        if( i == cv ){
            l = (WFac[i]-dh)*Accu->GetCV(i)->GetRange()/Accu->GetCV(i)->GetNumOfBins();
        } else {
            l = WFac[i]*Accu->GetCV(i)->GetRange()/Accu->GetCV(i)->GetNumOfBins();
        }
        CVLengths2[i] = l*l;
    }
    GetKernelDerIJAna(ip,jp,kblock1);

    for(size_t i=0; i < NumOfCVs; i++){
        double l;
        if( i == cv ){
            l = (WFac[i]+dh)*Accu->GetCV(i)->GetRange()/Accu->GetCV(i)->GetNumOfBins();
        } else {
            l = WFac[i]*Accu->GetCV(i)->GetRange()/Accu->GetCV(i)->GetNumOfBins();
        }
        CVLengths2[i] = l*l;
    }
    GetKernelDerIJAna(ip,jp,kblock2);

    for(size_t ii=0; ii < NumOfCVs; ii++){
        for(size_t jj=0; jj < NumOfCVs; jj++){
            kblock[ii][jj] = (kblock2[ii][jj]-kblock1[ii][jj])/(2.0*dh);
        }
    }

    // restore original CVLengths2
    for(size_t i=0; i < NumOfCVs; i++){
        double l = WFac[i]*Accu->GetCV(i)->GetRange()/Accu->GetCV(i)->GetNumOfBins();
        CVLengths2[i] = l*l;
    }
}

//------------------------------------------------------------------------------

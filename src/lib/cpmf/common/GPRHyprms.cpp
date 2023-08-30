// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2023 Petr Kulhanek, kulhanek@chemi.muni.cz
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

#include <GPRHyprms.hpp>
#include <string>
#include <vector>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <fstream>
#include <SciLapack.hpp>

//------------------------------------------------------------------------------

using namespace std;
using namespace boost;
using namespace boost::algorithm;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CGPRHyprms::CGPRHyprms(void)
{
    NumOfThreads    = 0;

    Method          = EGPRLA_LU;

    UseInv          = false;
    NeedInv         = false;

    RCond           = 1e-7;

    NumOfSigmaF2    = 0;
    NumOfCoVar      = 0;
    NumOfNCorr      = 0;
    NumOfSigmaN2    = 0;
}

//------------------------------------------------------------------------------

CGPRHyprms::~CGPRHyprms(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CGPRHyprms::SetSigmaF2(const CSmallString& spec)
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is not set!");
    }
    if( NumOfSigmaF2 == 0 ){
        RUNTIME_ERROR("NumOfSigmaF2 is zero");
    }

    string          sspec(spec);
    vector<string>  ssigmaf2;

    split(ssigmaf2,sspec,is_any_of("x"),token_compress_on);

    if( ssigmaf2.size() > NumOfSigmaF2){
        CSmallString error;
        error << "too many sigmaf2 (" << ssigmaf2.size() << ") than required (" << NumOfSigmaF2 << ")";
        RUNTIME_ERROR(error);
    }

    SigmaF2.CreateVector(NumOfSigmaF2);

    // parse values of sigmaf2
    double last_sigmaf2 = 0.1;
    for(size_t i=0; i < ssigmaf2.size(); i++){
        stringstream str(ssigmaf2[i]);
        str >> last_sigmaf2;
        if( ! str ){
            CSmallString error;
            error << "unable to decode sigmaf2 value for position: " << i+1;
            RUNTIME_ERROR(error);
        }
        SigmaF2[i] = last_sigmaf2;
    }

    // pad the rest with the last value
    for(size_t i=ssigmaf2.size(); i < NumOfSigmaF2; i++){
        SigmaF2[i] = last_sigmaf2;
    }
}

//------------------------------------------------------------------------------

void CGPRHyprms::SetSigmaF2(CSimpleVector<double>& sigmaf2)
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is not set!");
    }
    if( NumOfSigmaF2 == 0 ){
        RUNTIME_ERROR("NumOfSigmaF2 is zero");
    }

    if( sigmaf2.GetLength() != NumOfSigmaF2 ){
        RUNTIME_ERROR("dimmension inconsistent in the source and target");
    }
    SigmaF2 = sigmaf2;
}

//------------------------------------------------------------------------------

void CGPRHyprms::SetSigmaF2(size_t cvind, double value)
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is not set!");
    }
    if( NumOfSigmaF2 == 0 ){
        RUNTIME_ERROR("NumOfSigmaF2 is zero");
    }

    if( cvind >= NumOfSigmaF2 ){
        RUNTIME_ERROR("cvind out-of-range");
    }
    // is SigmaF2 initialized?
    if( SigmaF2.GetLength() == 0 ){
        SigmaF2.CreateVector(NumOfSigmaF2);
    }
    SigmaF2[cvind] = value;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CGPRHyprms::SetCoVar(const CSmallString& spec)
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is not set!");
    }
    if( NumOfCoVar == 0 ){
        RUNTIME_ERROR("NumOfCoVar is zero");
    }

    string          sspec(spec);
    vector<string>  scovar;

    split(scovar,sspec,is_any_of("x"),token_compress_on);

    if( scovar.size() > NumOfCoVar){
        CSmallString error;
        error << "too many covariances (" << scovar.size() << ") than required (" << NumOfCoVar << ")";
        RUNTIME_ERROR(error);
    }

    CoVar.CreateVector(NumOfCoVar);

    // parse values of scovar
    double last_covar = 0.0;
    for(size_t i=0; i < scovar.size(); i++){
        stringstream str(scovar[i]);
        str >> last_covar;
        if( ! str ){
            CSmallString error;
            error << "unable to decode scovar value for position: " << i+1;
            RUNTIME_ERROR(error);
        }
        CoVar[i] = last_covar;
    }

    // pad the rest with the last value
    for(size_t i=scovar.size(); i < NumOfCoVar; i++){
        CoVar[i] = last_covar;
    }
}

//------------------------------------------------------------------------------

void CGPRHyprms::SetCoVar(CSimpleVector<double>& covar)
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is not set!");
    }
    if( NumOfCoVar == 0 ){
        RUNTIME_ERROR("NumOfCoVar is zero");
    }

    if( covar.GetLength() != NumOfCoVar ){
        RUNTIME_ERROR("dimmension inconsistent in the source and target");
    }
    CoVar = covar;
}

//------------------------------------------------------------------------------

void CGPRHyprms::SetCoVar(size_t cvind, double value)
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is not set!");
    }
    if( NumOfCoVar == 0 ){
        RUNTIME_ERROR("NumOfCoVar is zero");
    }

    if( cvind >= NumOfCoVar ){
        RUNTIME_ERROR("cvind out-of-range");
    }
    // is CoVar initialized?
    if( CoVar.GetLength() == 0 ){
        CoVar.CreateVector(NumOfCoVar);
    }
    CoVar[cvind] = value;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CGPRHyprms::SetNCorr(const CSmallString& spec)
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is not set!");
    }
    if( NumOfNCorr == 0 ){
        RUNTIME_ERROR("NumOfNCorr is zero");
    }

    string          sspec(spec);
    vector<string>  sncorr;

    split(sncorr,sspec,is_any_of("x"),token_compress_on);

    if( sncorr.size() > NumOfNCorr ){
        CSmallString error;
        error << "too many ncorr (" << sncorr.size() << ") than required (" << NumOfNCorr << ")";
        RUNTIME_ERROR(error);
    }

    NCorr.CreateVector(NumOfNCorr);

    // parse values of ncorr
    double last_ncorr = 1.0;
    for(size_t i=0; i < sncorr.size(); i++){
        stringstream str(sncorr[i]);
        str >> last_ncorr;
        if( ! str ){
            CSmallString error;
            error << "unable to decode ncorr value for position: " << i+1;
            RUNTIME_ERROR(error);
        }
        NCorr[i] = last_ncorr;
    }

    // pad the rest with the last value
    for(size_t i=sncorr.size(); i < NumOfNCorr; i++){
        NCorr[i] = last_ncorr;
    }
}

//------------------------------------------------------------------------------

void CGPRHyprms::SetNCorr(CSimpleVector<double>& ncorr)
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is not set!");
    }
    if( NumOfNCorr == 0 ){
        RUNTIME_ERROR("NumOfNCorr is zero");
    }

    if( ncorr.GetLength() != NumOfNCorr ){
        RUNTIME_ERROR("ncvs inconsistent in the source and target");
    }

    NCorr = ncorr;
}

//------------------------------------------------------------------------------

void CGPRHyprms::SetNCorr(size_t cvind, double value)
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is not set!");
    }
    if( NumOfNCorr == 0 ){
        RUNTIME_ERROR("NumOfNCorr is zero");
    }

    if( cvind >= NumOfNCorr ){
        RUNTIME_ERROR("cvind out-of-range");
    }
    // is NCorr initialized?
    if( NCorr.GetLength() == 0 ){
        NCorr.CreateVector(NumOfNCorr);
    }
    NCorr[cvind] = value;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CGPRHyprms::SetSigmaN2(const CSmallString& spec)
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is not set!");
    }
    if( NumOfSigmaN2 == 0 ){
        RUNTIME_ERROR("NumOfSigmaN2 is zero");
    }

    string          sspec(spec);
    vector<string>  ssigman2;

    split(ssigman2,sspec,is_any_of("x"),token_compress_on);

    if( ssigman2.size() > NumOfSigmaN2){
        CSmallString error;
        error << "too many sigman2 (" << ssigman2.size() << ") than required (" << NumOfSigmaN2 << ")";
        RUNTIME_ERROR(error);
    }

    SigmaN2.CreateVector(NumOfSigmaN2);

    // parse values of sigman2
    double last_sigman2 = 0.0;
    for(size_t i=0; i < ssigman2.size(); i++){
        stringstream str(ssigman2[i]);
        str >> last_sigman2;
        if( ! str ){
            CSmallString error;
            error << "unable to decode sigman2 value for position: " << i+1;
            RUNTIME_ERROR(error);
        }
        SigmaN2[i] = last_sigman2;
    }

    // pad the rest with the last value
    for(size_t i=ssigman2.size(); i < NumOfSigmaN2; i++){
        SigmaN2[i] = last_sigman2;
    }
}

//------------------------------------------------------------------------------

void CGPRHyprms::SetSigmaN2(CSimpleVector<double>& sigman2)
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is not set!");
    }
    if( NumOfSigmaN2 == 0 ){
        RUNTIME_ERROR("NumOfSigmaN2 is zero");
    }

    if( sigman2.GetLength() != NumOfSigmaN2 ){
        RUNTIME_ERROR("ncvs inconsistent in the source and target");
    }

    SigmaN2 = sigman2;
}

//------------------------------------------------------------------------------

void CGPRHyprms::SetSigmaN2(size_t cvind, double value)
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is not set!");
    }
    if( NumOfSigmaN2 == 0 ){
        RUNTIME_ERROR("NumOfSigmaN2 is zero");
    }

    if( cvind >= NumOfSigmaN2 ){
        RUNTIME_ERROR("cvind out-of-range");
    }
    // is SigmaN2 initialized?
    if( SigmaN2.GetLength() == 0 ){
        SigmaN2.CreateVector(NumOfSigmaN2);
    }
    SigmaN2[cvind] = value;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CGPRHyprms::LoadGPRHyprms(const CSmallString& name)
{
    ifstream fin;
    fin.open(name);
    if( ! fin ){
        CSmallString error;
        error << "unable to open file with GPR hyperparameters: " << name;
        RUNTIME_ERROR(error);
    }

    string line;
    while( getline(fin,line) ){
        // is it comment?
        if( (line.size() > 0) && (line[0] == '#') ) continue;

        // parse line
        stringstream str(line);
        string key, buf;
        double value;
        str >> key >> buf >> value;
        if( ! str ){
            CSmallString error;
            error << "GPR hyperparameters file, unable to decode line: " << line.c_str();
            RUNTIME_ERROR(error);
        }
        if( key.find("SigmaF2#") != string::npos ) {
            std::replace( key.begin(), key.end(), '#', ' ');
            stringstream kstr(key);
            string swfac;
            int    cvind;
            kstr >> swfac >> cvind;
            if( ! kstr ){
                CSmallString error;
                error << "GPR hyperparameters file, unable to decode sigmaf2 key: " << key.c_str();
                RUNTIME_ERROR(error);
            }
            cvind--; // transform to 0-based indexing
            SetSigmaF2(cvind,value);
        } else         if( key.find("CoVar#") != string::npos ) {
            std::replace( key.begin(), key.end(), '#', ' ');
            stringstream kstr(key);
            string swfac;
            int    cvind;
            kstr >> swfac >> cvind;
            if( ! kstr ){
                CSmallString error;
                error << "GPR hyperparameters file, unable to decode covar key: " << key.c_str();
                RUNTIME_ERROR(error);
            }
            cvind--; // transform to 0-based indexing
            SetCoVar(cvind,value);
        } else if( key.find("WFac#") != string::npos ) {
            std::replace( key.begin(), key.end(), '#', ' ');
            stringstream kstr(key);
            string swfac;
            int    cvind;
            kstr >> swfac >> cvind;
            if( ! kstr ){
                CSmallString error;
                error << "GPR hyperparameters file, unable to decode wfac key: " << key.c_str();
                RUNTIME_ERROR(error);
            }
            cvind--; // transform to 0-based indexing
            SetWFac(cvind,value);
        } else if( key.find("SigmaN2#") != string::npos ) {
            std::replace( key.begin(), key.end(), '#', ' ');
            stringstream kstr(key);
            string swfac;
            int    cvind;
            kstr >> swfac >> cvind;
            if( ! kstr ){
                CSmallString error;
                error << "GPR hyperparameters file, unable to decode sigman2 key: " << key.c_str();
                RUNTIME_ERROR(error);
            }
            cvind--; // transform to 0-based indexing
            SetSigmaN2(cvind,value);
        } else {
            CSmallString error;
            error << "GPR hyperparameters file, unrecognized key: " << key.c_str();
            RUNTIME_ERROR(error);
        }
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

double CGPRHyprms::GetLogML(void)
{
    double ml = 0.0;
    return(ml);
}

//------------------------------------------------------------------------------

double CGPRHyprms::GetLogPL(void)
{
    double loo = 0.0;
    return(loo);
}

//------------------------------------------------------------------------------

void CGPRHyprms::GetLogMLDerivatives(const std::vector<bool>& flags,CSimpleVector<double>& der)
{
}

//------------------------------------------------------------------------------

void CGPRHyprms::GetLogPLDerivatives(const std::vector<bool>& flags,CSimpleVector<double>& der)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CGPRHyprms::PrintExecInfo(CVerboseStr& vout)
{
    NumOfThreads = 1;

#if defined(_OPENMP)
    {
        NumOfThreads = omp_get_max_threads();
        vout << "   OpenMP - number of threads: " << NumOfThreads << endl;
    }
#else
    vout << "   No OpenMP - sequential mode." << endl;
#endif
    RunBlasLapackPar();
    CSciLapack::PrintExecInfo(vout);
}

//------------------------------------------------------------------------------

void CGPRHyprms::RunBlasLapackSeq(void)
{
    CSciLapack::SetNumThreadsLocal(1);
}

//------------------------------------------------------------------------------

void CGPRHyprms::RunBlasLapackPar(void)
{
    CSciLapack::SetNumThreadsLocal(NumOfThreads);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CGPRHyprms::SetLAMethod(EGPRLAMethod mset)
{
    Method = mset;
}

//------------------------------------------------------------------------------

void CGPRHyprms::SetLAMethod(const CSmallString& method)
{
    if( method == "svd" ){
        SetLAMethod(EGPRLA_SVD);
    } else if( method == "svd2" ){
        SetLAMethod(EGPRLA_SVD2);
    } else if( method == "lu" ) {
        SetLAMethod(EGPRLA_LU);
    } else if( method == "ll" ) {
        SetLAMethod(EGPRLA_LL);
    } else if( method == "default" ) {
        SetLAMethod(EGPRLA_LU);
    } else {
        CSmallString error;
        error << "Specified method '" << method << "' for linear algebra is not supported. "
                 "Supported methods are: svd (simple driver), svd2 (conquer and divide driver), lu, ll (Cholesky decomposition), default (=lu)";
        INVALID_ARGUMENT(error);
    }
}

//------------------------------------------------------------------------------

void CGPRHyprms::SetRCond(double rcond)
{
    RCond = rcond;
}

//------------------------------------------------------------------------------

void CGPRHyprms::SetUseInv(bool iset)
{
    UseInv = iset;
}

//------------------------------------------------------------------------------

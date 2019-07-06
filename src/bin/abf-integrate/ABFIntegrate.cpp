// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2008 Martin Petrek, petrek@chemi.muni.cz
//                       Petr Kulhanek, kulhanek@enzim.hu
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
#include <ErrorSystem.hpp>
#include <ABFIntegratorRFD.hpp>
#include <ABFIntegratorRFD2.hpp>
#include <ABFIntegratorRBF.hpp>
#include <ABFIntegratorGPR.hpp>
#include <EnergySurface.hpp>
#include <ESPrinter.hpp>
#include "ABFIntegrate.hpp"
#include <iomanip>

//------------------------------------------------------------------------------

using namespace std;

//------------------------------------------------------------------------------

MAIN_ENTRY(CABFIntegrate)

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CABFIntegrate::CABFIntegrate(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CABFIntegrate::Init(int argc,char* argv[])
{
// encode program options, all check procedures are done inside of CABFIntOpts
    int result = Options.ParseCmdLine(argc,argv);

// should we exit or was it error?
    if(result != SO_CONTINUE) return(result);

// attach verbose stream to cout and set desired verbosity level
    vout.Attach(Console);
    if( Options.GetOptVerbose() ) {
        vout.Verbosity(CVerboseStr::debug);
    } else {
        vout.Verbosity(CVerboseStr::high);
    }

    StartTime.GetActualTimeAndDate();

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# abf-integrate (PMFLib utility)  started at " << StartTime.GetSDateAndTime() << endl;
    vout << "# Version: " << LibBuildVersion_PMF << endl;
    vout << "# ==============================================================================" << endl;

    if(Options.GetArgABFAccuName() != "-") {
        vout << "# ABF accu file (in)    : " << Options.GetArgABFAccuName() << endl;
    } else {
        vout << "# ABF accu file (in)    : - (standard input)" << endl;
    }
    if(Options.GetArgFEOutputName() != "-") {
        vout << "# Free energy file (out): " << Options.GetArgFEOutputName() << endl;
    } else {
        vout << "# Free energy file (out): - (standard output)" << endl;
    }
    vout << "# ------------------------------------------------" << endl;
        if(Options.GetOptMethod() == "rfd" ) {
            vout << "# Integration method    : RFD (reverse finite differences via csparse)" << endl;
        } else if(Options.GetOptMethod() == "rfd2" ) {
            vout << "# Integration method    : RFD (reverse finite differences via lapack)" << endl;
        } else if( Options.GetOptMethod() == "rbf" ){
            vout << "# Integration method    : RBF (radial basis functions)" << endl;
        } else if( Options.GetOptMethod() == "gpr" ) {
            vout << "# Integration method    : GPR (gaussian process)" << endl;
        } else {
            INVALID_ARGUMENT("method - not implemented");
        }
        if( Options.GetOptWithError() ) {
        vout << "# Integrated domains    : force+error" << endl;
        } else {
        vout << "# Integrated domains    : force only" << endl;
        }
    if( (Options.GetOptMethod() == "rfd") || (Options.GetOptMethod() == "rfd2" ) ){
        vout << "# FD number of points   : " << Options.GetOptFDPoints() << endl;
        vout << "# Periodicity           : " << bool_to_str(Options.GetOptPeriodicity()) << endl;
        if( Options.GetOptMethod() == "rfd2" ){
        vout << "# SVD rcond             : " << setprecision(3) << Options.GetOptRCond() << endl;
        }
    } else if ( Options.GetOptMethod() == "rbf" ){
        if( Options.IsOptRFac2Set() ){
        vout << "# Reduction factor rfac : " << setprecision(3) << Options.GetOptRFac() << " x " << setprecision(3) << Options.GetOptRFac2()  << endl;
        } else {
        vout << "# Reduction factor rfac : " << setprecision(3) << Options.GetOptRFac() << endl;
        }
        if( Options.IsOptWFac2Set() ){
        vout << "# Width factor wfac     : " << setprecision(3) << Options.GetOptWFac() << " x " << setprecision(3) << Options.GetOptWFac2()  << endl;
        } else {
        vout << "# Width factor wfac     : " << setprecision(3) << Options.GetOptWFac() << endl;
        }
        vout << "# SVD rcond             : " << setprecision(3) << Options.GetOptRCond() << endl;
        vout << "# RBF overhang          : " << Options.GetOptOverhang() << endl;
    } else if ( Options.GetOptMethod() == "gpr"  ) {
        vout << "# SigmaF2               : " << setprecision(3) << Options.GetOptSigmaF2() << endl;  
        if( Options.IsOptWFac2Set() ){
        vout << "# Width factor wfac     : " << setprecision(3) << Options.GetOptWFac() << " x " << setprecision(3) << Options.GetOptWFac2()  << endl;
        } else {
        vout << "# Width factor wfac     : " << setprecision(3) << Options.GetOptWFac() << endl;
        }
        vout << "# SVD rcond             : " << setprecision(3) << Options.GetOptRCond() << endl;
    } else {
        ES_ERROR("not implemented method");
        return(SO_USER_ERROR);
    }
    vout << "# ------------------------------------------------" << endl;
    if(Options.GetOptLimit() == 0) {
        vout << "# Sampling limit        : all bins will be taken into account" << endl;
    } else {
        vout << "# Sampling limit        : " << Options.GetOptLimit() << endl;
    }
    if(Options.GetOptEnergyLimit() == -1) {
        vout << "# Energy limit          : not applied" << endl;
    } else {
        vout << "# Energy limit          : " << Options.GetOptEnergyLimit() << endl;
    }
        vout << "# Number of corr. sam.  : " << Options.GetOptNCorr() << endl;
        vout << "# Integration offset    : " << Options.GetOptOffset() << endl;
    vout << "# ------------------------------------------------" << endl;
    vout << "# Output FES format     : " << Options.GetOptOutputFormat() << endl;
    vout << "# No header to output   : " << bool_to_str(Options.GetOptNoHeader()) << endl;
    vout << "# X format              : " << Options.GetOptIXFormat() << endl;
    vout << "# Y format              : " << Options.GetOptOEFormat() << endl;
    vout << "# ------------------------------------------------------------------------------" << endl;

#ifdef HAVE_MKL_PARALLEL
    vout << "# Note: linked with parallel version of MKL" << endl;
#endif

    // open files -----------------------------------
    if( InputFile.Open(Options.GetArgABFAccuName(),"r") == false ){
        ES_ERROR("unable to open input file");
        return(SO_USER_ERROR);
    }
    if( OutputFile.Open(Options.GetArgFEOutputName(),"w") == false ){
        ES_ERROR("unable to open output file");
        return(SO_USER_ERROR);
    }

    return(SO_CONTINUE);
}

//------------------------------------------------------------------------------

bool CABFIntegrate::Run(void)
{
   CABFAccumulator accumulator;

// load accumulator
    vout << endl;
    vout << "1) Loading ABF accumulator: " << Options.GetArgABFAccuName() << endl;
    try {
        accumulator.Load(InputFile);
    } catch(...) {
        ES_ERROR("unable to load the input ABF accumulator file");
        return(false);
    }
    vout << "   Done" << endl;

    // print CVS info
    accumulator.PrintCVSInfo(vout);
    accumulator.SetNCorr(Options.GetOptNCorr());

    // print header
    if((Options.GetOptNoHeader() == false) && (Options.GetOptOutputFormat() != "fes")) {
        fprintf(OutputFile,"# PMFLib version        : %s\n",LibBuildVersion_PMF);
        fprintf(OutputFile,"# data integrated by    : ");
        if(Options.GetOptMethod() == "rfd" ) {
            fprintf(OutputFile,"RFD (reverse finite differences via csparse)\n");
        } else if(Options.GetOptMethod() == "rfd2" ) {
            fprintf(OutputFile,"RFD (reverse finite differences via lapack)\n");
        } else if( Options.GetOptMethod() == "rbf" ){
            fprintf(OutputFile,"RBF (radial basis functions)\n");
        } else if( Options.GetOptMethod() == "gpr" ) {
            fprintf(OutputFile,"GPR (gaussian process)\n");
        } else {
            INVALID_ARGUMENT("method - not implemented");
        }
        if( Options.GetOptWithError() ) {
        fprintf(OutputFile,"# Integrated domains    : force+error\n");
        } else {
        fprintf(OutputFile,"# Integrated domains    : force only\n");
        }
        if( (Options.GetOptMethod() == "rfd") || (Options.GetOptMethod() == "rfd2") ){
            fprintf(OutputFile,"# FD nuber of points    : %d\n", Options.GetOptFDPoints());
            fprintf(OutputFile,"# Periodicity           : %s\n",(const char*)bool_to_str(Options.GetOptPeriodicity()));
            if( Options.GetOptMethod() == "rfd2" ){
            fprintf(OutputFile,"# SVD rcond             : %5.3f\n", Options.GetOptRCond());
            }
        } else if ( Options.GetOptMethod() == "rbf" ){
            if( Options.IsOptRFac2Set() ){
            fprintf(OutputFile,"# Reduction factor rfac : %5.3f x %5.3f\n", Options.GetOptRFac(), Options.GetOptRFac2());
            } else {
            fprintf(OutputFile,"# Reduction factor rfac : %5.3f\n", Options.GetOptRFac());
            }
            if( Options.IsOptWFac2Set() ){
            fprintf(OutputFile,"# Width factor wfac     : %5.3f X %5.3f\n", Options.GetOptWFac(), Options.GetOptWFac2());
            } else {
            fprintf(OutputFile,"# Width factor wfac     : %5.3f\n", Options.GetOptWFac());
            }
            fprintf(OutputFile,"# SVD rcond             : %5.3f\n", Options.GetOptRCond());
            fprintf(OutputFile,"# RBF overhang          : %d\n", Options.GetOptOverhang());
        } else if ( Options.GetOptMethod() == "gpr"  ) {
            fprintf(OutputFile,"# SigmaF2               : %5.3f\n", Options.GetOptSigmaF2());
            if( Options.IsOptWFac2Set() ){
            fprintf(OutputFile,"# Width factor wfac     : %5.3f X %5.3f\n", Options.GetOptWFac(), Options.GetOptWFac2());
            } else {
            fprintf(OutputFile,"# Width factor wfac     : %5.3f\n", Options.GetOptWFac());
            }
            fprintf(OutputFile,"# SVD rcond             : %5.3f\n", Options.GetOptRCond());
        } else {
            ES_ERROR("not implemented method");
        }

        fprintf(OutputFile,"# Sample limit          : %d\n",Options.GetOptLimit());
        fprintf(OutputFile,"# Energy limit          : %f\n",Options.GetOptEnergyLimit());
        fprintf(OutputFile,"# Number of corr. sam.  : %5.3f\n",Options.GetOptNCorr());
        fprintf(OutputFile,"# Number of coordinates : %d\n",accumulator.GetNumberOfCoords());
        fprintf(OutputFile,"# Total number of bins  : %d\n",accumulator.GetNumberOfBins());
    }

// prepare accumulator --------------------------
    vout << endl;
    vout << "2) Preparing ABF accumulator for integration (sampling limit)"<< endl;
    PrepareAccumulatorI(accumulator);
    vout << "   Done" << endl;

    CEnergySurface     fes;
    fes.Allocate(&accumulator);

    if( Options.GetOptEnergyLimit() > 0.0 ){
        vout << endl;
        vout << "3) ABF accumulator integration (" << Options.GetOptEcutMethod() << ")" << endl;
        if( IntegrateForEcut(accumulator,fes) == false ) return(false);

        vout << endl;
        vout << "2) Preparing ABF accumulator for integration (energy limit)"<< endl;
        PrepareAccumulatorII(accumulator,fes);
        vout << "   Done" << endl;

        fes.Clear();
    }

// integrate data ------------------------------
    vout << endl;
    vout << "3) ABF accumulator integration (" << Options.GetOptMethod() << ")" << endl;
    if( Integrate(accumulator,fes) == false ) return(false);
    vout << "   Done" << endl;

 // apply offset
    fes.ApplyOffset(Options.GetOptOffset() - fes.GetGlobalMinimumValue());

    if( (Options.GetOptWithError()) && (Options.GetOptMethod() != "gpr") ){
        vout << endl;
        vout << "3) ABF accumulator integration (errors)"<< endl;

        if( IntegrateErrors(accumulator,fes) == false ) return(false);

        fes.AdaptErrorsToGlobalMinimum();

        vout << "   Done" << endl;
    }

    if( Options.GetOptUnsampledAsMaxE() ){
        fes.AdaptUnsampledToMaxEnergy();
    }

// print result ---------------------------------
    vout << endl;
    vout << "4) Writing results to file: " << Options.GetArgFEOutputName() << endl;
    CESPrinter printer;

    printer.SetXFormat(Options.GetOptIXFormat());
    printer.SetYFormat(Options.GetOptOEFormat());
    if(Options.GetOptOutputFormat() == "plain") {
        printer.SetOutputFormat(EESPF_PLAIN);
    } else if(Options.GetOptOutputFormat() == "gnuplot") {
        printer.SetOutputFormat(EESPF_GNUPLOT);
    } else if(Options.GetOptOutputFormat() == "fes") {
        printer.SetOutputFormat(EESPF_PMF_FES);
    } else {
        INVALID_ARGUMENT("output format - not implemented");
    }

    if(Options.GetOptPrintAll()) {
        printer.SetSampleLimit(0);
    } else {
        printer.SetSampleLimit(Options.GetOptLimit());
    }

    printer.SetIncludeError(Options.GetOptWithError());
    printer.SetPrintedES(&fes);

    try {
        printer.Print(OutputFile);
    } catch(...) {
        ES_ERROR("unable to save the output free energy file");
        return(false);
    }
    vout << "   Done" << endl;

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CABFIntegrate::IntegrateForEcut(CABFAccumulator&accumulator, CEnergySurface& fes)
{
    if(Options.GetOptEcutMethod() == "rfd" ) {
        CABFIntegratorRFD   integrator;

        integrator.SetPeriodicity(Options.GetOptPeriodicity());
        integrator.SetFDPoints(Options.GetOptFDPoints());

        if( Options.GetOptEcutMethod() == Options.GetOptMethod() ){
            if( Options.GetOptLAMethod() == "lu" ){
                // nothing to do - LU is default
            } else if( Options.GetOptLAMethod() == "default" ) {
                // nothing to do - use default method set in constructor of integrator
            } else {
                INVALID_ARGUMENT("algorithm - not implemented");
            }
        }

        integrator.SetInputABFAccumulator(&accumulator);
        integrator.SetOutputFESurface(&fes);

        if(integrator.Integrate(vout,false) == false) {
            ES_ERROR("unable to integrate ABF accumulator");
            return(false);
        }
    } else if(Options.GetOptEcutMethod() == "rfd2" ) {
        CABFIntegratorRFD2   integrator;

        integrator.SetPeriodicity(Options.GetOptPeriodicity());
        integrator.SetFDPoints(Options.GetOptFDPoints());
        integrator.SetRCond(Options.GetOptRCond());

        if( Options.GetOptEcutMethod() == Options.GetOptMethod() ){
            if( Options.GetOptLAMethod() == "svd" ){
                integrator.SetLLSMehod(ERFDLLS_SVD);
            } else if( Options.GetOptLAMethod() == "qr" ) {
                integrator.SetLLSMehod(ERFDLLS_QR);
            } else if( Options.GetOptLAMethod() == "default" ) {
                // nothing to do - use default method set in constructor of integrator
            } else {
                INVALID_ARGUMENT("algorithm - not implemented");
            }
        }

        integrator.SetInputABFAccumulator(&accumulator);
        integrator.SetOutputFESurface(&fes);

        if(integrator.Integrate(vout,false) == false) {
            ES_ERROR("unable to integrate ABF accumulator");
            return(false);
        }
    } else if( Options.GetOptEcutMethod() == "rbf" ){
        CABFIntegratorRBF   integrator;

        integrator.SetWFac1(Options.GetOptWFac());
        integrator.SetWFac2(Options.GetOptWFac2());
        integrator.SetRCond(Options.GetOptRCond());
        integrator.SetRFac1(Options.GetOptRFac());
        integrator.SetRFac2(Options.GetOptRFac2());
        integrator.SetOverhang(Options.GetOptOverhang());

        if( Options.GetOptEcutMethod() == Options.GetOptMethod() ){
            if( Options.GetOptLAMethod() == "svd" ){
                integrator.SetLLSMehod(ERBFLLS_SVD);
            } else if( Options.GetOptLAMethod() == "qr" ) {
                integrator.SetLLSMehod(ERBFLLS_QR);
            } else if( Options.GetOptLAMethod() == "default" ) {
                // nothing to do - use default method set in constructor of integrator
            } else {
                INVALID_ARGUMENT("algorithm - not implemented");
            }
        }

        integrator.SetInputABFAccumulator(&accumulator);
        integrator.SetOutputFESurface(&fes);

        if(integrator.Integrate(vout,false) == false) {
            ES_ERROR("unable to integrate ABF accumulator");
            return(false);
        }
    } else if( Options.GetOptEcutMethod() == "gpr" ){
        CABFIntegratorGPR   integrator;

        integrator.SetWFac1(Options.GetOptWFac());
        integrator.SetWFac2(Options.GetOptWFac2());
        integrator.SetRCond(Options.GetOptRCond());
        integrator.SetSigmaF2(Options.GetOptSigmaF2());
        integrator.SetIncludeError(false);

        if( Options.GetOptEcutMethod() == Options.GetOptMethod() ){
            if( Options.GetOptLAMethod() == "svd" ){
                integrator.SetINVMehod(EGPRINV_SVD);
            } else if( Options.GetOptLAMethod() == "lu" ) {
                integrator.SetINVMehod(EGPRINV_LU);
            } else if( Options.GetOptLAMethod() == "default" ) {
                // nothing to do - use default method set in constructor of integrator
            } else {
                INVALID_ARGUMENT("algorithm - not implemented");
            }
        }

        integrator.SetInputABFAccumulator(&accumulator);
        integrator.SetOutputFESurface(&fes);

        if(integrator.Integrate(vout) == false) {
            ES_ERROR("unable to integrate ABF accumulator");
            return(false);
        }
    } else {
        INVALID_ARGUMENT("method - not implemented");
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CABFIntegrate::Integrate(CABFAccumulator&accumulator, CEnergySurface& fes)
{
    if(Options.GetOptMethod() == "rfd" ) {
        CABFIntegratorRFD   integrator;

        integrator.SetPeriodicity(Options.GetOptPeriodicity());
        integrator.SetFDPoints(Options.GetOptFDPoints());

        if( Options.GetOptLAMethod() == "lu" ){
            // nothing to do - LU is default
        } else if( Options.GetOptLAMethod() == "default" ) {
            // nothing to do - use default method set in constructor of integrator
        } else {
            INVALID_ARGUMENT("algorithm - not implemented");
        }

        integrator.SetInputABFAccumulator(&accumulator);
        integrator.SetOutputFESurface(&fes);

        if(integrator.Integrate(vout,false) == false) {
            ES_ERROR("unable to integrate ABF accumulator");
            return(false);
        }
    } else if(Options.GetOptMethod() == "rfd2" ) {
        CABFIntegratorRFD2   integrator;

        integrator.SetPeriodicity(Options.GetOptPeriodicity());
        integrator.SetFDPoints(Options.GetOptFDPoints());

        integrator.SetRCond(Options.GetOptRCond());

        if( Options.GetOptLAMethod() == "svd" ){
            integrator.SetLLSMehod(ERFDLLS_SVD);
        } else if( Options.GetOptLAMethod() == "qr" ) {
            integrator.SetLLSMehod(ERFDLLS_QR);
        } else if( Options.GetOptLAMethod() == "default" ) {
            // nothing to do - use default method set in constructor of integrator
        } else {
            INVALID_ARGUMENT("algorithm - not implemented");
        }

        integrator.SetInputABFAccumulator(&accumulator);
        integrator.SetOutputFESurface(&fes);

        if(integrator.Integrate(vout,false) == false) {
            ES_ERROR("unable to integrate ABF accumulator");
            return(false);
        }
    } else if( Options.GetOptMethod() == "rbf" ){
        CABFIntegratorRBF   integrator;

        integrator.SetWFac1(Options.GetOptWFac());
        integrator.SetWFac2(Options.GetOptWFac2());
        integrator.SetRCond(Options.GetOptRCond());
        integrator.SetRFac1(Options.GetOptRFac());
        integrator.SetRFac2(Options.GetOptRFac2());
        integrator.SetOverhang(Options.GetOptOverhang());

        if( Options.GetOptLAMethod() == "svd" ){
            integrator.SetLLSMehod(ERBFLLS_SVD);
        } else if( Options.GetOptLAMethod() == "qr" ) {
            integrator.SetLLSMehod(ERBFLLS_QR);
        } else if( Options.GetOptLAMethod() == "default" ) {
            // nothing to do - use default method set in constructor of integrator
        } else {
            INVALID_ARGUMENT("algorithm - not implemented");
        }

        integrator.SetInputABFAccumulator(&accumulator);
        integrator.SetOutputFESurface(&fes);

        if(integrator.Integrate(vout,false) == false) {
            ES_ERROR("unable to integrate ABF accumulator");
            return(false);
        }
    } else if( Options.GetOptMethod() == "gpr" ){
        CABFIntegratorGPR   integrator;

        integrator.SetWFac1(Options.GetOptWFac());
        integrator.SetWFac2(Options.GetOptWFac2());
        integrator.SetRCond(Options.GetOptRCond());
        integrator.SetSigmaF2(Options.GetOptSigmaF2());
        integrator.SetIncludeError(Options.GetOptWithError());
        integrator.SetNoEnergy(Options.GetOptNoEnergy());

        if( Options.GetOptLAMethod() == "svd" ){
            integrator.SetINVMehod(EGPRINV_SVD);
        } else if( Options.GetOptLAMethod() == "lu" ) {
            integrator.SetINVMehod(EGPRINV_LU);
        } else if( Options.GetOptLAMethod() == "default" ) {
            // nothing to do - use default method set in constructor of integrator
        } else {
            INVALID_ARGUMENT("algorithm - not implemented");
        }

        integrator.SetInputABFAccumulator(&accumulator);
        integrator.SetOutputFESurface(&fes);

        if(integrator.Integrate(vout) == false) {
            ES_ERROR("unable to integrate ABF accumulator");
            return(false);
        }
    } else {
        INVALID_ARGUMENT("method - not implemented");
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CABFIntegrate::IntegrateErrors(CABFAccumulator&accumulator, CEnergySurface& fes)
{
    if(Options.GetOptMethod() == "rfd" ) {
        CABFIntegratorRFD   integrator;

        integrator.SetPeriodicity(Options.GetOptPeriodicity());
        integrator.SetFDPoints(Options.GetOptFDPoints());

        if( Options.GetOptLAMethod() == "lu" ){
            // nothing to do - LU is default
        } else if( Options.GetOptLAMethod() == "default" ) {
            // nothing to do - use default method set in constructor of integrator
        } else {
            INVALID_ARGUMENT("algorithm - not implemented");
        }

        integrator.SetInputABFAccumulator(&accumulator);
        integrator.SetOutputFESurface(&fes);

        if(integrator.Integrate(vout,true) == false) {
            ES_ERROR("unable to integrate ABF accumulator");
            return(false);
        }
    } else if(Options.GetOptMethod() == "rfd2" ) {
        CABFIntegratorRFD2   integrator;

        integrator.SetPeriodicity(Options.GetOptPeriodicity());
        integrator.SetFDPoints(Options.GetOptFDPoints());

        integrator.SetRCond(Options.GetOptRCond());

        if( Options.GetOptLAMethod() == "svd" ){
            integrator.SetLLSMehod(ERFDLLS_SVD);
        } else if( Options.GetOptLAMethod() == "qr" ) {
            integrator.SetLLSMehod(ERFDLLS_QR);
        } else if( Options.GetOptLAMethod() == "default" ) {
            // nothing to do - use default method set in constructor of integrator
        }  else {
            INVALID_ARGUMENT("algorithm - not implemented");
        }

        integrator.SetInputABFAccumulator(&accumulator);
        integrator.SetOutputFESurface(&fes);

        if(integrator.Integrate(vout,true) == false) {
            ES_ERROR("unable to integrate ABF accumulator");
            return(false);
        }
    } else if( Options.GetOptMethod() == "rbf" ){
        CABFIntegratorRBF   integrator;

        integrator.SetWFac1(Options.GetOptWFac());
        integrator.SetWFac2(Options.GetOptWFac2());
        integrator.SetRCond(Options.GetOptRCond());
        integrator.SetRFac1(Options.GetOptRFac());
        integrator.SetRFac2(Options.GetOptRFac2());
        integrator.SetOverhang(Options.GetOptOverhang());

        if( Options.GetOptLAMethod() == "svd" ){
            integrator.SetLLSMehod(ERBFLLS_SVD);
        } else if( Options.GetOptLAMethod() == "qr" ) {
            integrator.SetLLSMehod(ERBFLLS_QR);
        } else if( Options.GetOptLAMethod() == "default" ) {
            // nothing to do - use default method set in constructor of integrator
        } else {
            INVALID_ARGUMENT("algorithm - not implemented");
        }

        integrator.SetInputABFAccumulator(&accumulator);
        integrator.SetOutputFESurface(&fes);

        if(integrator.Integrate(vout,true) == false) {
            ES_ERROR("unable to integrate ABF accumulator");
            return(false);
        }
    } else if( Options.GetOptMethod() == "gpr" ){
        vout << "   Already integrated." << endl;
    } else {
        INVALID_ARGUMENT("method - not implemented");
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

// this part performs following tasks:
//    a) bins with number of samples <= limit will be set to zero

void CABFIntegrate::PrepareAccumulatorI(CABFAccumulator& accumulator)
{
    double sampled = 0.0;
    double maxbins = 0.0;
    for(int ibin=0; ibin < accumulator.GetNumberOfBins(); ibin++) {
        maxbins++;
        // erase datapoints not properly sampled
        if( accumulator.GetNumberOfABFSamples(ibin) <= Options.GetOptLimit() ) {
            accumulator.SetNumberOfABFSamples(ibin,0);
        } else {
            sampled++;
        }
    }

    // calculate sampled area
    if( maxbins > 0 ){
        vout << "   Sampled area: " << setw(5) << setprecision(1) << fixed << sampled/maxbins*100 <<"%" << endl;
    }
}

//------------------------------------------------------------------------------

// this part performs following tasks:
//    a) erase data points with large energy

void CABFIntegrate::PrepareAccumulatorII(CABFAccumulator& accumulator,CEnergySurface& fes)
{
    double sampled = 0.0;
    double maxbins = 0.0;
    for(int ibin=0; ibin < accumulator.GetNumberOfBins(); ibin++) {
        maxbins++;

        if( accumulator.GetNumberOfABFSamples(ibin) > Options.GetOptLimit() ) {
            // consider only properly sampled data points
            if( fes.GetEnergy(ibin) > Options.GetOptEnergyLimit() ){
                // erase datapoints with too large energy
                accumulator.SetNumberOfABFSamples(ibin,0);
            } else {
                sampled++;
            }
        }
    }

    // calculate sampled area
    if( maxbins > 0 ){
        vout << "   Sampled area: " << setw(5) << setprecision(1) << fixed << sampled/maxbins*100 <<"%" << endl;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABFIntegrate::Finalize(void)
{
    // close files if they are own by program
    InputFile.Close();
    OutputFile.Close();

    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    CSmallTime dur;
    dur = dt - StartTime;

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# abf-integrate terminated at " << dt.GetSDateAndTime() << ". Total time: " << dur.GetSTimeAndDay() << endl;
    vout << "# ==============================================================================" << endl;

    if( ErrorSystem.IsError() || Options.GetOptVerbose() ){
        ErrorSystem.PrintErrors(vout);
    }

    vout << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================


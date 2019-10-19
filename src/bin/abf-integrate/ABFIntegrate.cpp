// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2019 Petr Kulhanek, kulhanek@chemi.muni.cz
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
#include <algorithm>
#include <boost/format.hpp>

//------------------------------------------------------------------------------

using namespace std;
using namespace boost;

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

    if( (Options.GetNumberOfProgArgs() != 2) && (Options.GetNumberOfProgArgs() != 3) ){
        ES_ERROR("two or three arguments are expected");
        return(SO_OPTS_ERROR);
    }

    ABFAccuName = Options.GetProgArg(0);
    FEOutputName = Options.GetProgArg(1);
    if( Options.GetNumberOfProgArgs() == 3 ) {
        // optional
        FullFEOutputName = Options.GetProgArg(2);
    }

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

    if( ABFAccuName != "-") {
        vout << "# ABF accu file (in)    : " << ABFAccuName << endl;
    } else {
        vout << "# ABF accu file (in)    : - (standard input)" << endl;
    }
    if( FEOutputName != "-") {
        vout << "# Free energy file (out): " << FEOutputName << endl;
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
// ---------------------------------------
        vout << "# Linear algebra        : " << Options.GetOptLAMethod() << endl;
    if( Options.GetOptMethod() == "rfd2" ){
        vout << "# SVD rcond             : " << setprecision(3) << Options.GetOptRCond() << endl;
    } else if ( Options.GetOptMethod() == "rbf" ){
        if( (Options.GetOptLAMethod() == "svd") || (Options.GetOptLAMethod() == "default") ){
        vout << "# SVD rcond             : " << setprecision(3) << Options.GetOptRCond() << endl;
        }
   } else if ( Options.GetOptMethod() == "gpr"  ) {
        if( (Options.GetOptLAMethod() == "svd") || (Options.GetOptLAMethod() == "svd2")  ){
        vout << "# SVD rcond             : " << setprecision(3) << Options.GetOptRCond() << endl;
        }
    }
// ---------------------------------------
        vout << "# ------------------------------------------------" << endl;
    if( (Options.GetOptMethod() == "rfd") || (Options.GetOptMethod() == "rfd2" ) ){
        vout << "# FD number of points   : " << Options.GetOptFDPoints() << endl;
        vout << "# Periodicity           : " << bool_to_str(Options.GetOptPeriodicity()) << endl;
    } else if ( Options.GetOptMethod() == "rbf" ){
        vout << "# Reduction factor rfac : " << Options.GetOptRFac() << endl;
        vout << "# Width factor wfac     : " << Options.GetOptWFac() << endl;
        vout << "# RBF overhang          : " << Options.GetOptOverhang() << endl;
    } else if ( Options.GetOptMethod() == "gpr"  ) {
        if( Options.IsOptLoadHyprmsSet() ){
            vout << "# GPR hyperprms file    : " << Options.GetOptLoadHyprms() << endl;
            // actual values are printed in detailed output from integrator
        } else {
            vout << "# SigmaF2               : " << setprecision(3) << Options.GetOptSigmaF2() << endl;
            vout << "# NCorr                 : " << setprecision(3) << Options.GetOptNCorr() << endl;
            vout << "# Width factor wfac     : " << Options.GetOptWFac() << endl;
        }
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
    if( (Options.GetOptEcutMethod() == "gpr") || (Options.GetOptEcutMethod() == "rbf") ){
    if(Options.GetOptMFMaxZScore() == -1) {
        vout << "# Max MF error Z-score  : not applied" << endl;
    } else {
        vout << "# Max MF error Z-score  : " << Options.GetOptMFMaxZScore() << endl;
        vout << "# Number of MF Z-tests  : " << Options.GetOptMFZTestPasses() << endl;
    }
    }
        vout << "# Glueing FES factor    : " << Options.GetOptGlueingFactor() << endl;
        vout << "# Glue holes on FES     : " << bool_to_str(Options.GetOptGlueHoles()) << endl;
    if(Options.GetOptEnergyLimit() == -1) {
        vout << "# Energy limit          : not applied" << endl;
    } else {
        vout << "# Energy limit          : " << Options.GetOptEnergyLimit() << endl;
    }
        vout << "# Skip flood fill test  : " << bool_to_str(Options.GetOptNoHeader()) << endl;

    vout << "# ------------------------------------------------" << endl;

    if( Options.IsOptGlobalMinSet() ){
    vout << "# Global FES minimum    : " << Options.GetOptGlobalMin() << endl;
    } else {
    vout << "# Global FES minimum    : -auto-" << endl;
    }
    vout << "# Integration offset    : " << Options.GetOptOffset() << endl;
    vout << "# Output FES format     : " << Options.GetOptOutputFormat() << endl;
    vout << "# No header to output   : " << bool_to_str(Options.GetOptNoHeader()) << endl;
    vout << "# Include bin statuses  : " << bool_to_str(Options.GetOptIncludeBinStat()) << endl;
    vout << "# X format              : " << Options.GetOptIXFormat() << endl;
    vout << "# Y format              : " << Options.GetOptOEFormat() << endl;
    vout << "# ------------------------------------------------------------------------------" << endl;

    // open files -----------------------------------
    if( InputFile.Open(ABFAccuName,"r") == false ){
        ES_ERROR("unable to open input file");
        return(SO_USER_ERROR);
    }
    if( OutputFile.Open(FEOutputName,"w") == false ){
        ES_ERROR("unable to open output file");
        return(SO_USER_ERROR);
    }

    return(SO_CONTINUE);
}

//------------------------------------------------------------------------------

bool CABFIntegrate::Run(void)
{
// load accumulator
    int state = 1;

    vout << endl;
    vout << format("%02d|Loading ABF accumulator: %s")%state%string(ABFAccuName) << endl;
    state++;
    try {
        Accumulator.Load(InputFile);
    } catch(...) {
        ES_ERROR("unable to load the input ABF accumulator file");
        return(false);
    }
    vout << "   Done" << endl;

    // print CVS info
    Accumulator.PrintCVSInfo(vout);
    // DO NOT SET IT HERE, Ncorr is now GPR hyperparameter
    // Accumulator.SetNCorr(Options.GetOptNCorr());
    FES.Allocate(&Accumulator);

    if( Options.GetOptMethod() == "gpr" ){
        // test early stage parsing of --globalmin
        CABFIntegratorGPR   integrator;

        integrator.SetInputABFAccumulator(&Accumulator);
        integrator.SetOutputFESurface(&FES);

        if( Options.IsOptGlobalMinSet() ){
            integrator.SetGlobalMin(Options.GetOptGlobalMin());
        }
    }

// sampling limit -------------------------------
    vout << endl;
    vout << format("%02d|Preparing ABF accumulator for integration (sampling limit)")%state << endl;
    state++;
    PrepareAccumulatorI();
    if( ! Options.GetOptSkipFFTest() ){
        FloodFillTest();
    }
    PrintSampledStat();
    vout << "   Done" << endl;

    if( Options.GetOptMFMaxZScore() > 0.0 ){
        for(int i=1; i <= Options.GetOptMFZTestPasses(); i++ ){
            vout << endl;
            vout << format("%02d|ABF accumulator integration (%s) for mean force error Z-score test #%d")%state%string(Options.GetOptEcutMethod())%i << endl;
            state++;
            if( IntegrateForMFZScore(i) == false ) return(false);

            vout << endl;
            vout << format("%02d|Preparing ABF accumulator for integration (mean force error z-score test #%d")%state%i << ")"<< endl;
            state++;
            PrepareAccumulatorI();
            if( ! Options.GetOptSkipFFTest() ){
                FloodFillTest();
            }
            PrintSampledStat();
            vout << "   Done" << endl;

            FES.Clear();
        }
    }

// glue fes ------------------------------------
    if( Options.GetOptGlueHoles() ){
        vout << endl;
        vout << format("%02d|Preparing ABF accumulator for integration (glue holes on FES)")%state << endl;
        state++;
        GlueHoles();
        PrintSampledStat();
        vout << "   Done" << endl;
        FES.Clear();
    }

    if( Options.GetOptGlueingFactor() > 0 ){
        vout << endl;
        vout << format("%02d|Preparing ABF accumulator for integration (glueing FES)")%state << endl;
        state++;
        vout << "   Searching for border regions in close vicinity of sampled areas ..." << endl;
        for(int i=1; i <= Options.GetOptGlueingFactor(); i++ ){
            GlueingFES(i);
        }
        PrintSampledStat();
        vout << "   Done" << endl;

        FES.Clear();
    }

// energy limit --------------------------------

    if( Options.GetOptEnergyLimit() > 0.0 ){
        vout << endl;
        vout << format("%02d|ABF accumulator integration (%s) for energy limit")%state%string(Options.GetOptEcutMethod()) << endl;
        state++;
        if( IntegrateForEcut() == false ) return(false);

        vout << endl;
        vout << format("%02d|Preparing ABF accumulator for integration (energy limit)")%state << endl;
        state++;
        PrepareAccumulatorII();
        if( ! Options.GetOptSkipFFTest() ){
            FloodFillTest();
        }
        PrintSampledStat();
        vout << "   Done" << endl;

        FES.Clear();
    }

// integrate data ------------------------------
    vout << endl;
    vout << format("%02d|ABF accumulator integration (%s)")%state%string(Options.GetOptMethod()) << endl;
    state++;
    if( Integrate() == false ) return(false);
    vout << "   Done" << endl;

 // apply offset
    if( ! Options.IsOptGlobalMinSet() ){
        FES.ApplyOffset(Options.GetOptOffset() - FES.GetGlobalMinimumValue());
    } else {
        FES.ApplyOffset(Options.GetOptOffset());
    }

    if( Options.GetOptUnsampledAsMaxE() ){
        if( Options.IsOptMaxEnergySet()){
            FES.AdaptUnsampledToMaxEnergy(Options.GetOptMaxEnergy());
        } else {
            FES.AdaptUnsampledToMaxEnergy();
        }
    }

// print result ---------------------------------
    vout << endl;
    vout << format("%02d|Writing results to file: %s")%state%string(FEOutputName) << endl;
    state++;
    CESPrinter printer;

    WriteHeader();

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

    printer.IncludeGluedAreas((Options.GetOptGlueingFactor() > 0)||Options.GetOptIncludeGluedRegions());
    printer.SetIncludeError(Options.GetOptWithError());
    printer.SetIncludeBinStat(Options.GetOptIncludeBinStat());
    printer.SetPrintedES(&FES);

    try {
        printer.Print(OutputFile);
    } catch(...) {
        ES_ERROR("unable to save the output free energy file");
        return(false);
    }
    vout << "   Done" << endl;

    if( Options.GetNumberOfProgArgs() == 3 ){
        vout << endl;
        vout << format("%02d|Writing results to file: %s (full version, --printall)")%state%string(FullFEOutputName) << endl;
        state++;

        if( OutputFile.Open(FullFEOutputName,"w") == false ){
            ES_ERROR("unable to open output file");
            return(SO_USER_ERROR);
        }

        WriteHeader();

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

        // print all
        printer.SetSampleLimit(0);
        printer.IncludeGluedAreas((Options.GetOptGlueingFactor() > 0)||Options.GetOptIncludeGluedRegions());
        printer.SetIncludeError(Options.GetOptWithError());
        printer.SetIncludeBinStat(Options.GetOptIncludeBinStat());
        printer.SetPrintedES(&FES);

        try {
            printer.Print(OutputFile);
        } catch(...) {
            ES_ERROR("unable to save the output free energy file");
            return(false);
        }
        vout << "   Done" << endl;
    }

// save accumulator if requested
    if( Options.GetOptSaveABF() != NULL ){
        vout << endl;
        vout << format("%02d|Saving ABF accumulator to : %s")%state%string(Options.GetOptSaveABF()) << endl;
        state++;
        try {
            Accumulator.Save(Options.GetOptSaveABF());
        } catch(...) {
            ES_ERROR("unable to save the ABF accumulator file");
            return(false);
        }
        vout << "   Done" << endl;
    }

    return(true);
}

//------------------------------------------------------------------------------

void CABFIntegrate::WriteHeader()
{
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
        fprintf(OutputFile,"# ------------------------------------------------------------------------------\n");
            fprintf(OutputFile,"# Linear algebra        : %s\n", (const char*)Options.GetOptLAMethod());
        if( Options.GetOptMethod() == "rfd2" ){
            fprintf(OutputFile,"# SVD rcond             : %5.4e\n",Options.GetOptRCond());
        } else if ( Options.GetOptMethod() == "rbf" ){
            if( (Options.GetOptLAMethod() == "svd") || (Options.GetOptLAMethod() == "default") ){
            fprintf(OutputFile,"# SVD rcond             : %5.4e\n",Options.GetOptRCond());
            }
        } else if ( Options.GetOptMethod() == "gpr"  ) {
            if( (Options.GetOptLAMethod() == "svd") || (Options.GetOptLAMethod() == "svd2")  ){
            fprintf(OutputFile,"# SVD rcond             : %5.4e\n",Options.GetOptRCond());
            }
        }
        fprintf(OutputFile,"# ------------------------------------------------------------------------------\n");
        if( (Options.GetOptMethod() == "rfd") || (Options.GetOptMethod() == "rfd2") ){
            fprintf(OutputFile,"# FD nuber of points    : %d\n", Options.GetOptFDPoints());
            fprintf(OutputFile,"# Periodicity           : %s\n",(const char*)bool_to_str(Options.GetOptPeriodicity()));
        } else if ( Options.GetOptMethod() == "rbf" ){
            fprintf(OutputFile,"# Reduction factor rfac : %s\n", (const char*)Options.GetOptRFac());
            fprintf(OutputFile,"# Width factor wfac     : %s\n", (const char*)Options.GetOptWFac());
            fprintf(OutputFile,"# RBF overhang          : %d\n", Options.GetOptOverhang());
        } else if ( Options.GetOptMethod() == "gpr" ) {
            if( Options.IsOptLoadHyprmsSet() ){
                fprintf(OutputFile,"# GPR hyperprms file    : %s\n", (const char*)Options.GetOptLoadHyprms());
                PrintGPRHyprms(OutputFile);
            } else {
                fprintf(OutputFile,"# SigmaF2               : %10.4f\n", Options.GetOptSigmaF2());
                fprintf(OutputFile,"# NCorr                 : %10.4f\n",Options.GetOptNCorr());
                fprintf(OutputFile,"# Width factor wfac     : %s\n", (const char*)Options.GetOptWFac());
            }
        } else {
            ES_ERROR("not implemented method");
        }

        fprintf(OutputFile,"# ------------------------------------------------------------------------------\n");
        fprintf(OutputFile,"# Sample limit          : %d\n",Options.GetOptLimit());
        if ( (Options.GetOptEcutMethod() == "gpr") || (Options.GetOptEcutMethod() == "rbf") ){
            fprintf(OutputFile,"# Max MF error Z-score  : %f\n",Options.GetOptMFMaxZScore());
            fprintf(OutputFile,"# Number of MF Z-tests  : %d\n",Options.GetOptMFZTestPasses());
        }
        fprintf(OutputFile,"# Glueing FES factor    : %d\n",Options.GetOptGlueingFactor());
        fprintf(OutputFile,"# Glue holes on FES     : %s\n", bool_to_str(Options.GetOptGlueHoles()));
        fprintf(OutputFile,"# Energy limit          : %f\n",Options.GetOptEnergyLimit());
        fprintf(OutputFile,"# Skip flood fill test  : %s\n", bool_to_str(Options.GetOptNoHeader()));
        fprintf(OutputFile,"# ------------------------------------------------------------------------------\n");
        if( Options.IsOptGlobalMinSet() ){
        fprintf(OutputFile,"# Global FES minimum    : %s\n",(const char*)Options.GetOptGlobalMin());
        } else {
        fprintf(OutputFile,"# Global FES minimum    : -auto-\n");
        }
        fprintf(OutputFile,"# Integration offset    : %5.3f\n", Options.GetOptOffset());
        fprintf(OutputFile,"# Include bin statuses  : %s\n",bool_to_str(Options.GetOptIncludeBinStat()));
        fprintf(OutputFile,"# Number of coordinates : %d\n",Accumulator.GetNumberOfCoords());
        fprintf(OutputFile,"# Total number of bins  : %d\n",Accumulator.GetNumberOfBins());
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CABFIntegrate::IntegrateForMFZScore(int pass)
{
    if(Options.GetOptEcutMethod() == "rfd" ) {
        ES_ERROR("illegal combination: --emethod=rfd and --maxzscore");
        return(false);
    } else if(Options.GetOptEcutMethod() == "rfd2" ) {
        ES_ERROR("illegal combination: --emethod=rfd2 and --maxzscore");
        return(false);
    } else if( Options.GetOptEcutMethod() == "rbf" ){
        CABFIntegratorRBF   integrator;

        integrator.SetWFac(Options.GetOptWFac());
        integrator.SetRCond(Options.GetOptRCond());
        integrator.SetRFac(Options.GetOptRFac());
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

        integrator.SetInputABFAccumulator(&Accumulator);
        integrator.SetOutputFESurface(&FES);

        if(integrator.Integrate(vout) == false) {
            ES_ERROR("unable to integrate ABF accumulator");
            return(false);
        }

        if( Options.IsOptMFInfoSet() ){
            CSmallString mfinfo = Options.GetOptMFInfo();
            mfinfo << ".mflimit" << pass;
            if( integrator.WriteMFInfo(mfinfo) == false ) return(false);
        }

        // apply mean force limit
        integrator.FilterByMFZScore(Options.GetOptMFMaxZScore(),vout);

    } else if( Options.GetOptEcutMethod() == "gpr" ){
        CABFIntegratorGPR   integrator;

        integrator.SetInputABFAccumulator(&Accumulator);
        integrator.SetOutputFESurface(&FES);

        if( Options.IsOptLoadHyprmsSet() ){
            LoadGPRHyprms(integrator);
        } else {
            integrator.SetSigmaF2(Options.GetOptSigmaF2());
            integrator.SetNCorr(Options.GetOptNCorr());
            integrator.SetWFac(Options.GetOptWFac());
        }

        integrator.SetNumericK(Options.GetOptNumericK());
        integrator.SetIncludeError(false);

        integrator.SetRCond(Options.GetOptRCond());
        if( Options.GetOptEcutMethod() == Options.GetOptMethod() ){
            if( Options.GetOptLAMethod() == "svd" ){
                integrator.SetINVMehod(EGPRINV_SVD);
            } else if( Options.GetOptLAMethod() == "svd2" ){
                integrator.SetINVMehod(EGPRINV_SVD2);
            } else if( Options.GetOptLAMethod() == "lu" ) {
                integrator.SetINVMehod(EGPRINV_LU);
            } else if( Options.GetOptLAMethod() == "ll" ) {
                integrator.SetINVMehod(EGPRINV_LL);
            } else if( Options.GetOptLAMethod() == "default" ) {
                // nothing to do - use default method set in constructor of integrator
            } else {
                INVALID_ARGUMENT("algorithm - not implemented");
            }
        }

        integrator.SetNoEnergy(true);

        if(integrator.Integrate(vout) == false) {
            ES_ERROR("unable to integrate ABF accumulator");
            return(false);
        }

        if( Options.IsOptMFInfoSet() ){
            CSmallString mfinfo = Options.GetOptMFInfo();
            mfinfo << ".mflimit" << pass;
            if( integrator.WriteMFInfo(mfinfo) == false ) return(false);
        }

        // apply mean force limit
        integrator.FilterByMFZScore(Options.GetOptMFMaxZScore(),vout);

    } else {
        INVALID_ARGUMENT("method - not implemented");
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CABFIntegrate::IntegrateForEcut(void)
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

        integrator.SetInputABFAccumulator(&Accumulator);
        integrator.SetOutputFESurface(&FES);

        if(integrator.Integrate(vout) == false) {
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

        integrator.SetInputABFAccumulator(&Accumulator);
        integrator.SetOutputFESurface(&FES);

        if(integrator.Integrate(vout) == false) {
            ES_ERROR("unable to integrate ABF accumulator");
            return(false);
        }
    } else if( Options.GetOptEcutMethod() == "rbf" ){
        CABFIntegratorRBF   integrator;

        integrator.SetWFac(Options.GetOptWFac());
        integrator.SetRCond(Options.GetOptRCond());
        integrator.SetRFac(Options.GetOptRFac());
        integrator.SetOverhang(Options.GetOptOverhang());
        integrator.IncludeGluedAreas((Options.GetOptGlueingFactor() > 0)||Options.GetOptIncludeGluedRegions());

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

        integrator.SetInputABFAccumulator(&Accumulator);
        integrator.SetOutputFESurface(&FES);

        if(integrator.Integrate(vout) == false) {
            ES_ERROR("unable to integrate ABF accumulator");
            return(false);
        }
    } else if( Options.GetOptEcutMethod() == "gpr" ){
        CABFIntegratorGPR   integrator;

        integrator.SetInputABFAccumulator(&Accumulator);
        integrator.SetOutputFESurface(&FES);

        if( Options.IsOptLoadHyprmsSet() ){
            LoadGPRHyprms(integrator);
        } else {
            integrator.SetSigmaF2(Options.GetOptSigmaF2());
            integrator.SetNCorr(Options.GetOptNCorr());
            integrator.SetWFac(Options.GetOptWFac());
        }

        integrator.SetNumericK(Options.GetOptNumericK());
        integrator.IncludeGluedAreas((Options.GetOptGlueingFactor() > 0)||Options.GetOptIncludeGluedRegions());
        integrator.SetIncludeError(false);

        integrator.SetRCond(Options.GetOptRCond());
        if( Options.GetOptEcutMethod() == Options.GetOptMethod() ){
            if( Options.GetOptLAMethod() == "svd" ){
                integrator.SetINVMehod(EGPRINV_SVD);
            } else if( Options.GetOptLAMethod() == "svd2" ){
                integrator.SetINVMehod(EGPRINV_SVD2);
            } else if( Options.GetOptLAMethod() == "lu" ) {
                integrator.SetINVMehod(EGPRINV_LU);
            } else if( Options.GetOptLAMethod() == "ll" ) {
                integrator.SetINVMehod(EGPRINV_LL);
            } else if( Options.GetOptLAMethod() == "default" ) {
                // nothing to do - use default method set in constructor of integrator
            } else {
                INVALID_ARGUMENT("algorithm - not implemented");
            }
        }

        if( Options.IsOptGlobalMinSet() ){
            integrator.SetGlobalMin(Options.GetOptGlobalMin());
        }

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

bool CABFIntegrate::Integrate()
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

        integrator.SetUseOldRFDMode(Options.GetOptUseOldRFD());
        integrator.SetInputABFAccumulator(&Accumulator);
        integrator.SetOutputFESurface(&FES);

        if(integrator.Integrate(vout) == false) {
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

        integrator.SetInputABFAccumulator(&Accumulator);
        integrator.SetOutputFESurface(&FES);

        if(integrator.Integrate(vout) == false) {
            ES_ERROR("unable to integrate ABF accumulator");
            return(false);
        }
    } else if( Options.GetOptMethod() == "rbf" ){
        CABFIntegratorRBF   integrator;

        integrator.SetWFac(Options.GetOptWFac());
        integrator.SetRCond(Options.GetOptRCond());
        integrator.SetRFac(Options.GetOptRFac());
        integrator.SetOverhang(Options.GetOptOverhang());
        integrator.IncludeGluedAreas((Options.GetOptGlueingFactor() > 0)||Options.GetOptIncludeGluedRegions());

        if( Options.GetOptLAMethod() == "svd" ){
            integrator.SetLLSMehod(ERBFLLS_SVD);
        } else if( Options.GetOptLAMethod() == "qr" ) {
            integrator.SetLLSMehod(ERBFLLS_QR);
        } else if( Options.GetOptLAMethod() == "default" ) {
            // nothing to do - use default method set in constructor of integrator
        } else {
            INVALID_ARGUMENT("algorithm - not implemented");
        }

        integrator.SetInputABFAccumulator(&Accumulator);
        integrator.SetOutputFESurface(&FES);

        if(integrator.Integrate(vout) == false) {
            ES_ERROR("unable to integrate ABF accumulator");
            return(false);
        }

        if( Options.IsOptMFInfoSet() ){
            if( integrator.WriteMFInfo(Options.GetOptMFInfo()) == false ) return(false);
        }

    } else if( Options.GetOptMethod() == "gpr" ){
        CABFIntegratorGPR   integrator;

        integrator.SetInputABFAccumulator(&Accumulator);
        integrator.SetOutputFESurface(&FES);

        if( Options.IsOptLoadHyprmsSet() ){
            LoadGPRHyprms(integrator);
        } else {
            integrator.SetSigmaF2(Options.GetOptSigmaF2());
            integrator.SetNCorr(Options.GetOptNCorr());
            integrator.SetWFac(Options.GetOptWFac());
        }

        integrator.SetIncludeError(Options.GetOptWithError());
        integrator.SetNoEnergy(Options.GetOptNoEnergy());
        integrator.SetNumericK(Options.GetOptNumericK());
        integrator.IncludeGluedAreas((Options.GetOptGlueingFactor() > 0)||Options.GetOptIncludeGluedRegions());

        if( Options.IsOptGlobalMinSet() ){
            integrator.SetGlobalMin(Options.GetOptGlobalMin());
        }

        integrator.SetRCond(Options.GetOptRCond());
        if( Options.GetOptLAMethod() == "svd" ){
            integrator.SetINVMehod(EGPRINV_SVD);
        } else if( Options.GetOptLAMethod() == "svd2" ){
                integrator.SetINVMehod(EGPRINV_SVD2);
        } else if( Options.GetOptLAMethod() == "lu" ) {
            integrator.SetINVMehod(EGPRINV_LU);
        } else if( Options.GetOptLAMethod() == "ll" ) {
            integrator.SetINVMehod(EGPRINV_LL);
        } else if( Options.GetOptLAMethod() == "default" ) {
            // nothing to do - use default method set in constructor of integrator
        } else {
            INVALID_ARGUMENT("algorithm - not implemented");
        }

        if(integrator.Integrate(vout) == false) {
            ES_ERROR("unable to integrate ABF accumulator");
            return(false);
        }

        if( Options.IsOptMFInfoSet() ){
            if( integrator.WriteMFInfo(Options.GetOptMFInfo()) == false ) return(false);
        }

    } else {
        INVALID_ARGUMENT("method - not implemented");
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABFIntegrate::PrintGPRHyprms(FILE* p_fout)
{
    ifstream fin;
    fin.open(Options.GetOptLoadHyprms());
    if( ! fin ){
        CSmallString error;
        error << "unable to open file with GPR hyperparameters: " << Options.GetOptLoadHyprms();
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
        if( key == "SigmaF2" ){
            fprintf(p_fout,"# SigmaF2               : %10.4f\n",value);
        } else if( key == "NCorr" ){
            fprintf(p_fout,"# NCorr                 : %10.4f\n",value);
        } else if( key.find("WFac#") != string::npos ) {
            fprintf(p_fout,"# %-7s               : %10.4f\n",key.c_str(),value);
        } else {
            CSmallString error;
            error << "GPR hyperparameters file, unrecognized key: " << key.c_str();
            RUNTIME_ERROR(error);
        }
    }
}

//------------------------------------------------------------------------------

void CABFIntegrate::LoadGPRHyprms(CABFIntegratorGPR& gpr)
{
    ifstream fin;
    fin.open(Options.GetOptLoadHyprms());
    if( ! fin ){
        CSmallString error;
        error << "unable to open file with GPR hyperparameters: " << Options.GetOptLoadHyprms();
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
        if( key == "SigmaF2" ){
            gpr.SetSigmaF2(value);
        } else if( key == "NCorr" ){
            gpr.SetNCorr(value);
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
            gpr.SetWFac(cvind,value);
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

// this part performs following tasks:
//    a) bins with number of samples <= limit will be set to zero

void CABFIntegrate::PrepareAccumulatorI(void)
{
    for(int ibin=0; ibin < Accumulator.GetNumberOfBins(); ibin++) {
        // erase datapoints not properly sampled, preserve glueing
        if( (Accumulator.GetNumberOfABFSamples(ibin) >= 0) && (Accumulator.GetNumberOfABFSamples(ibin) <= Options.GetOptLimit()) ) {
            Accumulator.SetNumberOfABFSamples(ibin,0);
        }
    }
}

//------------------------------------------------------------------------------

// this part performs following tasks:
//    a) erase data points with large energy

void CABFIntegrate::PrepareAccumulatorII(void)
{
    // filter by energy
    for(int ibin=0; ibin < Accumulator.GetNumberOfBins(); ibin++) {
        if( Accumulator.GetNumberOfABFSamples(ibin) != 0 ) {
            // consider only properly sampled data points
            if( FES.GetEnergy(ibin) > Options.GetOptEnergyLimit() ){
                // erase datapoints with too large energy
                Accumulator.SetNumberOfABFSamples(ibin,0);
            }
        }
    }
}

//------------------------------------------------------------------------------

void CABFIntegrate::PrintSampledStat(void)
{
    // calculate sampled area
    double maxbins = Accumulator.GetNumberOfBins();
    int    sampled = 0;
    int    glued = 0;
    for(int ibin=0; ibin < Accumulator.GetNumberOfBins(); ibin++) {
        if( Accumulator.GetNumberOfABFSamples(ibin) > 0 ) {
            sampled++;
        }
        if( Accumulator.GetNumberOfABFSamples(ibin) < 0 ) {
            glued++;
        }
    }
    if( maxbins > 0 ){
        vout << "   Sampled area: "
             << setw(6) << sampled << " / " << (int)maxbins << " | " << setw(5) << setprecision(1) << fixed << sampled/maxbins*100 <<"%" << endl;
    }
    if( glued > 0 ){
        vout << "   Glued area:   "
             << setw(6) << glued << " / " << (int)maxbins << " | " << setw(5) << setprecision(1) << fixed << glued/maxbins*100 <<"%" << endl;
    }
}

//------------------------------------------------------------------------------

void CABFIntegrate::FloodFillTest(void)
{
    vout << "   Searching for discontinuous regions ..." << endl;
    int seedid = 1;

    FFSeeds.CreateVector(Accumulator.GetNumberOfBins());
    FFSeeds.SetZero();
    IPos.CreateVector(Accumulator.GetNumberOfCoords());
    TPos.CreateVector(Accumulator.GetNumberOfCoords());

    double maxbins = Accumulator.GetNumberOfBins();
    int    maxseedid = 0;
    int    maxsampled = 0;
    bool   first = true;

    while( InstallNewSeed(seedid,false) ){
        int sampled = 1;    // for initial seed set by InstallNewSeed
        int newsamples = 0;

        while( (newsamples = FillSeed(seedid,false)) > 0 ){
            sampled += newsamples;
        }

        if( maxbins > 0 ){
            vout << "   Region: " << setw(6) << seedid << " - sampled area: "
                 << setw(6) << sampled << " / " << (int)maxbins << " | " << setw(5) << setprecision(1) << fixed << sampled/maxbins*100 <<"%" << endl;
        }

        if( first || (maxsampled < sampled) ){
            first = false;
            maxsampled = sampled;
            maxseedid = seedid;
        }

        seedid++;
    }
    seedid--;

    // quit if one or none region
    if( seedid <= 1 ){
        vout << "   All is continuous." << endl;
        return;
    }

    vout << "   Clearing all except region: " << setw(3) << maxseedid <<  endl;

    for(int ibin=0; ibin < Accumulator.GetNumberOfBins(); ibin++) {
        if( FFSeeds[ibin] != maxseedid ) {
            Accumulator.SetNumberOfABFSamples(ibin,0);
        }
    }
}

//------------------------------------------------------------------------------

bool CABFIntegrate::InstallNewSeed(int seedid,bool unsampled)
{
    for(int ibin=0; ibin < Accumulator.GetNumberOfBins(); ibin++) {
        if( unsampled ){
            if( (FFSeeds[ibin] == 0) && ( Accumulator.GetNumberOfABFSamples(ibin) == 0 ) ) {
                FFSeeds[ibin] = seedid;
                return(true);
            }
        } else {
            if( (FFSeeds[ibin] == 0) && ( Accumulator.GetNumberOfABFSamples(ibin) > 0 ) ) {
                FFSeeds[ibin] = seedid;
                return(true);
            }
        }
    }

    return(false);
}

//------------------------------------------------------------------------------

int CABFIntegrate::FillSeed(int seedid,bool unsampled)
{
    int newsamples = 0;
    int ndir = 1;
    for(int j=0; j < Accumulator.GetNumberOfCoords(); j++){
        ndir *= 3;
    }

    for(int ibin=0; ibin < Accumulator.GetNumberOfBins(); ibin++) {
        if( unsampled ){
            if( Accumulator.GetNumberOfABFSamples(ibin) > 0 ) continue; // skip sampled regions
        } else {
            if( Accumulator.GetNumberOfABFSamples(ibin) <= 0 ) continue; // skip unsampled regions
        }
        if( FFSeeds[ibin] != seedid ) continue; // skip different regions

        // convert to ipont
        Accumulator.GetIPoint(ibin,IPos);

        // in each direction
        for(int j=0; j < ndir; j++){
            GetTPoint(IPos,j,TPos);
            int tbin = Accumulator.GetGlobalIndex(TPos);
            if( tbin >= 0 ){
                if( FFSeeds[tbin] == 0 ){
                    if( unsampled ){
                        if( Accumulator.GetNumberOfABFSamples(tbin) == 0 ){
                            FFSeeds[tbin] = seedid;
                            newsamples++;
                        }
                    } else {
                        if( Accumulator.GetNumberOfABFSamples(tbin) > 0 ){
                            FFSeeds[tbin] = seedid;
                            newsamples++;
                        }
                    }
                }
            }
        }
    }

    return(newsamples);
}

//------------------------------------------------------------------------------

void CABFIntegrate::GetTPoint(CSimpleVector<int>& ipos,int d,CSimpleVector<int>& tpos)
{
    for(int k=Accumulator.GetNumberOfCoords()-1; k >= 0; k--) {
        int ibin = d % 3 - 1;
        tpos[k] = ibin + ipos[k];
        d = d / 3;
    }
}

//------------------------------------------------------------------------------

void CABFIntegrate::GlueingFES(int factor)
{
    IPos.CreateVector(Accumulator.GetNumberOfCoords());
    TPos.CreateVector(Accumulator.GetNumberOfCoords());

    int ndir = 1;
    for(int j=0; j < Accumulator.GetNumberOfCoords(); j++){
        ndir *= 3;
    }

    vout << "   Gluing FES: factor = " << factor;

    int glued = 0;

    for(int ibin=0; ibin < Accumulator.GetNumberOfBins(); ibin++) {
        if( Accumulator.GetNumberOfABFSamples(ibin) != 0 ) continue; // skip glued or sampled bins

        // convert to ipont
        Accumulator.GetIPoint(ibin,IPos);

        // is sampled or glued region in close vicinty?

        // in each direction
        for(int j=0; j < ndir; j++){
            GetTPoint(IPos,j,TPos);
            int tbin = Accumulator.GetGlobalIndex(TPos);
            if( tbin >= 0 ){
                if( factor == 1 ){
                    if( Accumulator.GetNumberOfABFSamples(tbin) > 0 ){
                        Accumulator.SetNumberOfABFSamples(ibin,-factor);
                        glued++;
                        break;
                    }
                } else {
                    if( - Accumulator.GetNumberOfABFSamples(tbin) == factor - 1 ){
                        Accumulator.SetNumberOfABFSamples(ibin,-factor);
                        glued++;
                        break;
                    }
                }
            }
        }
    }

    vout << ", glued bins = " << glued << endl;
}

//------------------------------------------------------------------------------

void CABFIntegrate::GlueHoles(void)
{
    vout << "   Searching for holes on FES ..." << endl;
    int seedid = 1;

    FFSeeds.CreateVector(Accumulator.GetNumberOfBins());
    FFSeeds.SetZero();
    IPos.CreateVector(Accumulator.GetNumberOfCoords());
    TPos.CreateVector(Accumulator.GetNumberOfCoords());

    double maxbins = Accumulator.GetNumberOfBins();
    int    numofholes = 0;

    SeedSampled(seedid);
    seedid++;

    while( InstallNewSeed(seedid,true) ){
        int sampled = 1;    // for initial seed set by InstallNewSeed
        int newsamples = 0;

        while( (newsamples = FillSeed(seedid,true)) > 0 ){
            sampled += newsamples;
        }

        // detect type of area
        bool hole = false;
        if( IsHole(seedid) ){
            MarkAsHole(seedid);
            numofholes++;
            hole = true;
        }

        if( maxbins > 0 ){
            vout << "   Region: " << setw(6) << seedid << " - sampled area: "
                 << setw(6) << sampled << " / " << (int)maxbins << " | " << setw(5) << setprecision(1) << fixed << sampled/maxbins*100 <<"%";
            if( hole ){
                vout << " hole - glued." << endl;
            } else {
                vout << " area at edge." << endl;
            }
        }
        seedid++;
    }

    // print stat
    vout << "   Number of holes: " <<  numofholes << endl;
}

//------------------------------------------------------------------------------

void CABFIntegrate::SeedSampled(int seedid)
{
    for(int ibin=0; ibin < Accumulator.GetNumberOfBins(); ibin++) {
        if( Accumulator.GetNumberOfABFSamples(ibin) > 0 ) {
            FFSeeds[ibin] = seedid;
        }
    }
}

//------------------------------------------------------------------------------

bool CABFIntegrate::IsHole(int seedid)
{
    IPos.CreateVector(Accumulator.GetNumberOfCoords());
    TPos.CreateVector(Accumulator.GetNumberOfCoords());

    int ndir = 1;
    for(int j=0; j < Accumulator.GetNumberOfCoords(); j++){
        ndir *= 3;
    }

    for(int ibin=0; ibin < Accumulator.GetNumberOfBins(); ibin++) {
        if( FFSeeds[ibin] != seedid ) continue;

        // convert to ipont
        Accumulator.GetIPoint(ibin,IPos);

        // test neighbouring of the point in each direction
        for(int j=0; j < ndir; j++){
            GetTPoint(IPos,j,TPos);
            int tbin = Accumulator.GetGlobalIndex(TPos);
            if( tbin < 0 ){
                // outside of ABF accumulator - it is not a hole
                return(false);
            }
        }
    }

    // it is a hole - all points are confined in sampled regions.
    return(true);
}

//------------------------------------------------------------------------------

void CABFIntegrate::MarkAsHole(int seedid)
{
    for(int ibin=0; ibin < Accumulator.GetNumberOfBins(); ibin++) {
        if( FFSeeds[ibin] == seedid ){
            Accumulator.SetNumberOfABFSamples(ibin,-1);
        }
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


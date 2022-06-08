// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2008 Petr Kulhanek, kulhanek@enzim.hu
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
#include "ACCUCombine.hpp"
#include <ErrorSystem.hpp>
#include <SmallTimeAndDate.hpp>
#include <boost/format.hpp>
#include <iomanip>

//------------------------------------------------------------------------------

using namespace std;
using namespace boost;

//------------------------------------------------------------------------------

MAIN_ENTRY(CACCUCombine)

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CACCUCombine::CACCUCombine(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CACCUCombine::Init(int argc,char* argv[])
{
// encode program options, all check procedures are done inside of CACCUIntOpts
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

    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# accu-combine (PMFLib utility)  started at " << dt.GetSDateAndTime() << endl;
    vout << "# Version: " << LibBuildVersion_PMF << endl;
    vout << "# ==============================================================================" << endl;
    for(int i=0; i < Options.GetNumberOfProgArgs() - 1; i++){
        vout << format("# Input PMF accumulator #%05d : %s")%(i+1)%Options.GetProgArg(i) << endl;
    }
    vout << format("# Output PMF accumulator     : %s")%Options.GetProgArg(Options.GetNumberOfProgArgs()-1) << endl;

    vout << "# ------------------------------------------------------------------------------" << endl;

    return(SO_CONTINUE);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CACCUCombine::Run(void)
{
// load accumulators
    int state = 1;
    vout << endl;
    vout << format("%02d:Loading input PMF accumulators ...")%state << endl;
    state++;
    try {
        for(int i=0; i < Options.GetNumberOfProgArgs() - 1; i++){
            vout << format("  >> %05d %s")%(i+1)%string(Options.GetProgArg(i)) << endl;
            CPMFAccumulatorPtr inaccu = CPMFAccumulatorPtr(new CPMFAccumulator);
            inaccu->Load(Options.GetProgArg(i));
            InAccus.push_back(inaccu);
        }
    } catch(...) {
        ES_ERROR("unable to load the input PMF accumulator");
        return(false);
    }

// -----------------------------------------------------------------------------
    vout << endl;
    vout << format("%02d:Statistics of input accumulators")%state << endl;
    state++;
    PrintSampledStat();
    vout << "   Done." << endl;

// -----------------------------------------------------------------------------

    vout << endl;
    if( Options.GetOptLinear() ){
        vout << format("%02d:Combining PMF accumulators ... (linear mode)")%state << endl;
        CombineLinear();
    } else {
        vout << format("%02d:Combining PMF accumulators ... (tree mode)")%state << endl;
        CombineAsTree();
    }
    state++;

    vout << endl;
    vout << format("%02d:Saving the resulting PMF accumulator: %s")%state%Options.GetProgArg(Options.GetNumberOfProgArgs()-1) << endl;
    state++;

// save final accumulator
    try {
        OutAccu->Save(Options.GetProgArg(Options.GetNumberOfProgArgs()-1));
    } catch(...) {
        ES_ERROR("unable to save final PMF accumulator file");
        return(false);
    }

    return(true);
}

//------------------------------------------------------------------------------

void CACCUCombine::CombineLinear(void)
{
    OutAccu = InAccus.front();  // the first accumulator becomes the output one
    // combine the others into the first one
    for(int i=1; i < Options.GetNumberOfProgArgs() - 1; i++){
        OutAccu->Combine(InAccus[i],Options.GetOptCommonOnly());
    }
}

//------------------------------------------------------------------------------

void CACCUCombine::CombineAsTree(void)
{
    int pass = 1;
    while( InAccus.size() > 1 ){
        vout << format("   Pass: %5d - number of accumulators: %5d")%pass%InAccus.size() << endl;
        std::vector<CPMFAccumulatorPtr> newset;
        std::vector<CPMFAccumulatorPtr>::iterator it = InAccus.begin();
        std::vector<CPMFAccumulatorPtr>::iterator ie = InAccus.end();
        CPMFAccumulatorPtr combined;
        bool first = true;
        while( it != ie ){
            CPMFAccumulatorPtr accu = *it;
            if( first ){
                // first
                combined = accu;
                newset.push_back(combined);
                it++;
                first = false;
            } else {
                // second - combine
                combined->Combine(accu,Options.GetOptCommonOnly());
                it++;
                if(  it != ie ){
                    // third - combine if it is the last item
                    CPMFAccumulatorPtr accu = *it;
                    if( accu == InAccus.back() ){
                        combined->Combine(accu,Options.GetOptCommonOnly());
                        it++;
                    }
                }
                first = true;
            }
        }
        InAccus = newset;
        pass++;
    }

    // final accumulator
    OutAccu = InAccus.front();
}

//------------------------------------------------------------------------------

void CACCUCombine::PrintSampledStat(void)
{
    for(size_t i=0; i < InAccus.size(); i++){
        CPMFAccumulatorPtr accu = InAccus[i];
        vout << format("** PMF Accumulator #%05d ...")%(i+1);
        // calculate sampled area
        double maxbins = accu->GetNumOfBins();
        int    sampled = 0;
        int    limit = 0;
        for(int ibin=0; ibin < accu->GetNumOfBins(); ibin++) {
            if( accu->GetNumOfSamples(ibin) > 0 ) {
                sampled++;
            }
            if( accu->GetNumOfSamples(ibin) > Options.GetOptLimit() ) {
                limit++;
            } else {
                accu->SetNumOfSamples(ibin,0);
            }
        }
        if( maxbins > 0 ){
            vout << " Sampled area: "
                 << setw(6) << sampled << " / " << (int)maxbins << " | " << setw(5) << setprecision(1) << fixed << sampled/maxbins*100 <<"%";
            vout << " Within limit: "
                 << setw(6) << limit << " / " << (int)maxbins << " | " << setw(5) << setprecision(1) << fixed << limit/maxbins*100 <<"%";
        }
        vout << endl;
        if( accu->CheckCVSInfo(InAccus[0]) == false ){
            CSmallString error;
            error << "inconsistent dimensions of two PMF accumulators";
            RUNTIME_ERROR(error);
        }
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CACCUCombine::Finalize(void)
{
    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# accu-combine terminated at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;

    if( ErrorSystem.IsError() || Options.GetOptVerbose() ){
        ErrorSystem.PrintErrors(vout);
    }

    vout << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================


// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2008 Petr Kulhanek, kulhanek@enzim.hu
//    Copyright (C) 2007 Petr Kulhanek, kulhanek@enzim.hu
//    Copyright (C) 2006 Petr Kulhanek, kulhanek@chemi.muni.cz &
//                       Martin Petrek, petrek@chemi.muni.cz
//    Copyright (C) 2005 Petr Kulhanek, kulhanek@chemi.muni.cz
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

#include <stdio.h>
#include <math.h>
#include <ErrorSystem.hpp>
#include "BMIntegrate.hpp"
#include <SmallTimeAndDate.hpp>
#include <errno.h>

//------------------------------------------------------------------------------

using namespace std;

//------------------------------------------------------------------------------

MAIN_ENTRY(CBMIntegrate)

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CBMIntegrate::CBMIntegrate(void)
{
    InputFile = NULL;
    OwnInputFile = false;
    OutputFile = NULL;
    OwnOutputFile = false;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CBMIntegrate::Init(int argc,char* argv[])
{
    if((InputFile != NULL) || (OutputFile != NULL)) return(false);     // files already opened

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

    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# con-integrat (PMFLib utility)  started at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;

    if(Options.GetArgInput() == "-") {
        vout << "# Input file         : - (standard input)" << endl;
    } else {
        vout << "# Input file         : " << Options.GetArgInput() << endl;
    }
    if(Options.GetArgOutput() == "-") {
        vout << "# Output file        : -  (standard output)" << endl;
    } else {
        vout << "# Output file        : " << Options.GetArgOutput() << endl;
    }
    vout << "# ------------------------------------------------------------------------------" << endl;
    vout << "# Skipped lines      : " << Options.GetOptSkipLines() << endl;
    vout << "# Analysed lines     : " << Options.GetOptAnalLines() << endl;
    vout << "# Padding lines      : " << Options.GetOptPadLines() << endl;
    vout << "# Integration offset : " << Options.GetOptOffset() << endl;
    vout << "# ------------------------------------------------------------------------------" << endl;
    vout << "# Sigma values       : " << bool_to_str(!Options.GetOptNoSigma()) << endl;
    vout << "# Print header       : " << bool_to_str(!Options.GetOptNoHeader()) << endl;
    vout << "# Verbose            : " << bool_to_str(Options.GetOptVerbose()) << endl;
    vout << "# ------------------------------------------------------------------------------" << endl;
    vout << "# X format           : " << Options.GetOptIXFormat() << endl;
    vout << "# Y format           : " << Options.GetOptIYFormat() << endl;
    vout << "# sigma(Y) format    : " << Options.GetOptISFormat() << endl;
    vout << "# I format           : " << Options.GetOptOIFormat() << endl;
    vout << "# sigma(I) format    : " << Options.GetOptOEFormat() << endl;
    vout << "# ------------------------------------------------------------------------------" << endl;
    vout << endl;

// open files -----------------------------------
    if(Options.GetArgInput() == "-") {
        InputFile = stdin;
        OwnInputFile = false;
    } else {
        InputFile = fopen(Options.GetArgInput(),"r");
        if(InputFile == NULL) {
            CSmallString error;
            error << "unable to open input file (" << strerror(errno) << ")";
            ES_ERROR(error);
            return(SO_USER_ERROR);
        }
        OwnInputFile = true;
    }

    if(Options.GetArgOutput() == "-") {
        OutputFile = stdout;
        OwnOutputFile = false;
    } else {
        OutputFile = fopen(Options.GetArgOutput(),"w");
        if(OutputFile == NULL) {
            CSmallString error;
            error << "unable to open output file (" << strerror(errno) << ")";
            ES_ERROR(error);
            if((OwnInputFile == true) && (InputFile != NULL)) {
                fclose(InputFile);
                InputFile = NULL;
                OwnInputFile = false;
            }
            return(SO_USER_ERROR);
        }
        OwnOutputFile = true;
    }

// complete output format -----------------------------------
    OutputFormat  = "   ";
    OutputFormat += Options.GetOptIXFormat();
    OutputFormat += " ";
    OutputFormat += Options.GetOptIYFormat();
    OutputFormat += " ";
    OutputFormat += Options.GetOptISFormat();
    OutputFormat += " ";
    OutputFormat += Options.GetOptOIFormat();
    OutputFormat += " ";
    OutputFormat += Options.GetOptOEFormat();
    OutputFormat += "\n";

    return(result);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CBMIntegrate::Run(void)
{
    if(ReadData() == false) {
        fprintf(stderr,"%s: unable to read data from input\n", (const char*)Options.GetProgramName());
        return(false);
    }

    if(IntegrateData() == false) {
        fprintf(stderr,"%s: unable to integrate data\n", (const char*)Options.GetProgramName());
        return(false);
    }

    if(PrintData() == false) {
        fprintf(stderr,"%s: unable to print data\n", (const char*)Options.GetProgramName());
        return(false);
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CBMIntegrate::ReadData(void)
{
// skip lines from input stream
    int line = 1;
    int c;
    int sl =  Options.GetOptSkipLines();
    while((sl != 0) && (feof(InputFile) == false)) {
        c = fgetc(InputFile);
        if(c == '\n') sl--;
    }

// integrate remaining data
    double x;
    double y;
    double s=0.0;

    int  numofdata = 0;

    while(feof(InputFile) == false) {

        // read data line ----------------------------------------------------------

        // read whole line
        CSmallString data_line;
        if(data_line.ReadLineFromFile(InputFile) == false) break;   // no more data

        if(Options.GetOptNoSigma() == true) {
            s = 0.0;
            int nr = sscanf(data_line,"%le %le",&x,&y);
            if(nr <= 0) break;      // no more data
            if(nr != 2) {
                CSmallString error;
                error << "incosistent number of data records on line " << line << ", requested: 2, found: " << nr;
                ES_ERROR(error);
                return(false);
            }
        } else {
            int nr = sscanf(data_line,"%le %le %le",&x,&y,&s);
            if(nr <= 0) break;      // no more data
            if(nr != 3) {
                CSmallString error;
                error << "incosistent number of data records on line " << line << ", requested: 3, found: " << nr;
                ES_ERROR(error);
                return(false);
            }
        }

        // integrate data ----------------------------------------------------------
        bool result = true;
        result &= X.AppendValue(x);
        result &= Y.AppendValue(y);
        result &= S.AppendValue(s);

        if(result == false) {
            CSmallString error;
            error << "unable to add data to memory";
            ES_ERROR(error);
            return(false);
        }

        // did we analyzed requested number of lines?
        if(numofdata == Options.GetOptAnalLines()) break;

        // skip requested number of lines -----------------------------------------
        int sl =  Options.GetOptPadLines();
        while((sl != 0) && (feof(InputFile) == false)) {
            c = fgetc(InputFile);
            if(c == '\n') sl--;
        }
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CBMIntegrate::IntegrateData(void)
{
    if(Options.GetOptVerbose() == true) {
        printf("# Number of points for integration : %d\n",X.GetSize());
    }

// helper variables
    double x,y,s;      // current value
    double px,py,ps;   // previous value
    double dx,ay,as;   // current difference
    double pi=0.0;     // integrated value
    double psi=0.0;    // integrated error

    bool   p_dir = true;   // direction of integration

// init first point
    if(X.GetSize() >= 1) {
        px = X[0];
        py = Y[0];
        ps = S[0];
    }

// append first point
    bool result = true;
    result &= I.AppendValue(pi);
    result &= IE.AppendValue(psi);

    if(result == false) {
        CSmallString error;
        error << "unable to add integrated data to memory";
        ES_ERROR(error);
        return(false);
    }

// integrate the rest of points
    for(int i=1; i < X.GetSize(); i++) {

        x = X[i];
        y = Y[i];
        s = S[i];

        // incremental changes
        dx = x - px;
        ay = y + py;
        as = s + ps;

        // integrate
        pi = 0.5*dx*ay + pi;
        psi = sqrt(0.5*fabs(dx)*as*as + psi*psi);

        // check integration direction
        if(i == 1) {
            p_dir = (dx >= 0.0);
        } else {
            if(p_dir != (dx >= 0.0)) {
                CSmallString error;
                error << "direction of x-values was changed";
                ES_ERROR(error);
                return(false);
            }
        }

        bool result = true;
        result &= I.AppendValue(pi);
        result &= IE.AppendValue(psi);

        if(result == false) {
            CSmallString error;
            error << "unable to add integrated data to memory";
            ES_ERROR(error);
            return(false);
        }

        px = x;
        py = y;
        ps = s;
    }

// find global minima
    double glb = I[0];
    for(int i=0; i < X.GetSize(); i++) {
        if(I[i] < glb) glb = I[i];
    }

// combine with offset
    for(int i=0; i < X.GetSize(); i++) {
        I[i] = I[i] - glb + Options.GetOptOffset();
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CBMIntegrate::PrintData(void)
{
    if(Options.GetOptNoHeader() == false) {
        // shell we update header sizes according to used formats? - maybe in the next version?
        fprintf(OutputFile,"#               X                Y       Sigma(Y)               I       Sigma(I)\n");
        fprintf(OutputFile,"#  --------------- --------------- -------------- --------------- --------------\n");
    }

    double x;
    double y;
    double s;
    double pi;
    double psi;

    for(int i=0; i < X.GetSize(); i++) {
        x = X[i];
        y = Y[i];
        s = S[i];
        pi  = I[i];
        psi = IE[i];

        // print data --------------------------------------------------------------
        if(fprintf(OutputFile,OutputFormat,x,y,s,pi,psi) <= 0) {
            CSmallString error;
            error << "unable to write to output file";
            ES_ERROR(error);
            return(false);
        }
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CBMIntegrate::Finalize(void)
{
// close files if they are own by program
    if((OwnInputFile == true) && (InputFile != NULL)) {
        fclose(InputFile);
        InputFile = NULL;
        OwnInputFile = false;
    }
    if((OwnOutputFile == true) && (OutputFile != NULL)) {
        fclose(OutputFile);
        OutputFile = NULL;
        OwnOutputFile = false;
    }

    CSmallTimeAndDate dt;
    dt.GetActualTimeAndDate();

    vout << endl;
    vout << "# ==============================================================================" << endl;
    vout << "# con-integrate terminated at " << dt.GetSDateAndTime() << endl;
    vout << "# ==============================================================================" << endl;

    if( ErrorSystem.IsError() || Options.GetOptVerbose() ){
        ErrorSystem.PrintErrors(vout);
    }

    vout << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================





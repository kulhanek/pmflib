// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
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

#include <errno.h>
#include <string.h>
#include <ABFAccumulator.hpp>
#include <ErrorSystem.hpp>
#include <XMLElement.hpp>
#include <XMLBinData.hpp>

//------------------------------------------------------------------------------

using namespace std;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CABFAccumulator::CABFAccumulator(void)
{
    NCVs            = 0;
    TotNBins        = 0;
    NCorr           = 1.0;
    EpotAvailable   = false;
    Temperature     = 300.0;
    EnergyFConv     = 1.0;
    EnergyUnit      = "kcal mol^-1";
}

//------------------------------------------------------------------------------

CABFAccumulator::~CABFAccumulator(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABFAccumulator::Load(const CSmallString& name)
{
    FILE* fin = fopen(name, "r");

    if(fin == NULL) {
        CSmallString error;
        error << "unable to open file '" << name << "' (" << strerror(errno) << ")";
        RUNTIME_ERROR(error);
    }

    try {
        Load(fin);
    } catch(...) {
        fclose(fin);
        throw;
    }

    fclose(fin);
}

//------------------------------------------------------------------------------

void CABFAccumulator::Load(FILE* fin)
{
    if( fin == NULL ) {
        INVALID_ARGUMENT("stream is not open");
    }

// read ABF accumulator header ------------------
    char abf_id[4];
    char ver_id[3];

// first line can contain ABF version
    char buffer[80];
    if( fgets(buffer,80,fin) == NULL ) {
        RUNTIME_ERROR("unable to read the first line");
    }

    int nr = sscanf(buffer,"%3s %2s",abf_id,ver_id);
    if( nr != 2 ){
        CSmallString error;
        error << "illegal header - two leading items expected (line: " << buffer << ")";
        RUNTIME_ERROR(error);
    }

    abf_id[3]='\0';
    ver_id[2]='\0';

// check ID string
    if(strcmp(abf_id,"ABF") != 0) {
        CSmallString error;
        error << "'ABF' magic word was not found in the first line (line: " << buffer << ")";
        RUNTIME_ERROR(error);
    }
    if(strcmp(ver_id,"V3") == 0) {
        Load_v3(buffer,fin);
        return;
    }
    if(strcmp(ver_id,"V4") == 0) {
        Load_v4(buffer,fin);
        return;
    }
    if(strcmp(ver_id,"V5") == 0) {
        Load_v5(buffer,fin);
        return;
    }

    CSmallString error;
    error << "Unsupported version of ABF accumulator (line: " << buffer << ")";
    RUNTIME_ERROR(error);
}

//------------------------------------------------------------------------------

void CABFAccumulator::Load_v3(char* fline,FILE* fin)
{
    if( fin == NULL ) {
        INVALID_ARGUMENT("stream is not open");
    }

// read ABF accumulator header ------------------
    char abf_id[4];
    char ver_id[3];
    int  ncvs = 0;
    int  nr;

    nr = sscanf(fline,"%3s %2s %d",abf_id,ver_id,&ncvs);
    if( nr != 3 ){
        CSmallString error;
        error << "illegal header - three items expected (line: " << fline << ")";
        RUNTIME_ERROR(error);
    }

    abf_id[3]='\0';
    ver_id[2]='\0';

// check ID string
    if(strcmp(abf_id,"ABF") != 0) {
        CSmallString error;
        error << "'ABF' magic word was not found in the first line (line: " << fline << ")";
        RUNTIME_ERROR(error);
    }
    if(strcmp(ver_id,"V3") != 0) {
        CSmallString error;
        error << "only ABF V3 version is supported (line: " << fline << ")";
        RUNTIME_ERROR(error);
    }

    if(ncvs <= 0) {
        CSmallString error;
        error << "number of coordinates has to be greater than zero, but " << ncvs << " was found";
        RUNTIME_ERROR(error);
    }

    SetNumOfCVs(ncvs);

// read coordinate specification ----------------
    for(int i=0; i < NCVs; i++) {
        int             id = 0;
        char            type[20];
        char            name[60];
        double          min_value = 0.0;
        double          max_value = 0.0;
        int             nbins = 0;
        int             tr = 0;

        memset(type,0,20);
        memset(name,0,60);

        // read item
        tr = fscanf(fin,"%d %10s %lf %lf %d",&id,type,&min_value,&max_value,&nbins);
        if( tr != 5 ) {
            CSmallString error;
            error << "unable to read coordinate definition, id: " << i+1 << " (" << tr << " != 5)";
            RUNTIME_ERROR(error);
        }

        // some tests
        if(id != i+1) {
            CSmallString error;
            error << "coordinate id does not match, read: " << id << ", expected: " << i+1;
            RUNTIME_ERROR(error);
        }
        if(max_value <= min_value) {
            CSmallString error;
            error << "min value is not smaller than max value, id: " << id;
            RUNTIME_ERROR(error);
        }
        if(nbins <= 0) {
            CSmallString error;
            error << "number of bins has to be grater than zero, id: " << id;
            RUNTIME_ERROR(error);
        }

        // read item
        tr = fscanf(fin,"%d %55s",&id,name);
        if( tr != 2 ) {
            CSmallString error;
            error << "unable to read coordinate definition, id: " << i+1 << " (" << tr << " != 2)";
            RUNTIME_ERROR(error);
        }
        // some tests
        if(id != i+1) {
            CSmallString error;
            error << "coordinate id does not match, read: " << id << ", expected: " << i+1;
            RUNTIME_ERROR(error);
        }

        // init coordinate
        SetCV(i,name,type,min_value,max_value,nbins);
    }

// alloc accumulator data -----------------------
    FinalizeAllocation();

// read accumulator data ------------------------

// ABF forces =================================================================

// samples
    for(int i=0; i < TotNBins; i++) {
        int ns = 0;
        if( fscanf(fin,"%d",&ns) != 1 ) {
            CSmallString error;
            error << "unable to read number of ABF samples for bin " << i;
            RUNTIME_ERROR(error);
        }
        NSamples[i] = ns;
    }

// do i=1,fnitem
//     read(iounit,40,end=100,err=100) (accumulator%ICFSum(i,j),j=1,accumulator%TotNBins)
// end do

// accumulated force
    for(int i = 0; i < NCVs; i++) {
        for(int j = 0; j < TotNBins; j++) {
            double cf = 0.0;
            if(fscanf(fin,"%lf",&cf) != 1) {
                CSmallString error;
                error << "unable to read accumulated ABF force for bin " << j << " and item " << i;
                RUNTIME_ERROR(error);
            }
            ICFSum[map(i,j)] = cf;
        }
    }

// accumulated force square
    for(int i = 0; i < NCVs; i++) {
        for(int j = 0; j < TotNBins; j++) {
            double cf = 0.0;
            if(fscanf(fin,"%lf",&cf) != 1) {
                CSmallString error;
                error << "unable to read accumulated ABF force squares for bin " << j << " and item " << i;
                RUNTIME_ERROR(error);
            }
            ICFSum2[map(i,j)] = cf;
        }
    }
}

//------------------------------------------------------------------------------

void CABFAccumulator::Load_v4(char* fline,FILE* fin)
{
    if( fin == NULL ) {
        INVALID_ARGUMENT("stream is not open");
    }

// read ABF accumulator header ------------------
    char abf_id[4];
    char ver_id[3];
    int  ncvs = 0;
    int  nr;

// first line can contain either two or four records for version 0

    nr = sscanf(fline,"%s %s %d",abf_id,ver_id,&ncvs);
    if( nr != 3 ){
        CSmallString error;
        error << "illegal header - three items expected (line: " << fline << ")";
        RUNTIME_ERROR(error);
    }

    abf_id[3]='\0';
    ver_id[2]='\0';

// check ID string
    if(strcmp(abf_id,"ABF") != 0) {
        CSmallString error;
        error << "'ABF' magic word was not found in the first line (line: " << fline << ")";
        RUNTIME_ERROR(error);
    }
    if(strcmp(ver_id,"V4") != 0) {
        CSmallString error;
        error << "only ABF V4 version is supported (line: " << fline << ")";
        RUNTIME_ERROR(error);
    }

    if(ncvs <= 0) {
        CSmallString error;
        error << "number of coordinates has to be greater than zero, but " << ncvs << " was found";
        RUNTIME_ERROR(error);
    }

    SetNumOfCVs(ncvs);

// read coordinate specification ----------------
    for(int i=0; i < NCVs; i++) {
        int             id = 0;
        char            type[15];
        char            name[60];
        char            unit[40];
        double          min_value = 0.0;
        double          max_value = 0.0;
        double          fconv = 1.0;
        int             nbins = 0;
        int             tr = 0;

        memset(type,0,15);
        memset(name,0,60);
        memset(unit,0,40);

        // read item
        tr = fscanf(fin,"%d %10s %lf %lf %d",&id,type,&min_value,&max_value,&nbins);
        if( tr != 5 ) {
            CSmallString error;
            error << "unable to read coordinate definition, id: " << i+1 << " (" << tr << " != 5)";
            RUNTIME_ERROR(error);
        }

        // some tests
        if(id != i+1) {
            CSmallString error;
            error << "coordinate id does not match, read: " << id << ", expected: " << i+1;
            RUNTIME_ERROR(error);
        }
        if(max_value <= min_value) {
            CSmallString error;
            error << "min value is not smaller than max value, id: " << id;
            RUNTIME_ERROR(error);
        }
        if(nbins <= 0) {
            CSmallString error;
            error << "number of bins has to be grater than zero, id: " << id;
            RUNTIME_ERROR(error);
        }

        // read item
        tr = fscanf(fin,"%d %55s",&id,name);
        if( tr != 2 ) {
            CSmallString error;
            error << "unable to read coordinate definition, id: " << i+1 << " (" << tr << " != 2)";
            RUNTIME_ERROR(error);
        }
        // some tests
        if(id != i+1) {
            CSmallString error;
            error << "coordinate id does not match, read: " << id << ", expected: " << i+1;
            RUNTIME_ERROR(error);
        }

        // read item
        tr = fscanf(fin,"%d %lf %37s",&id,&fconv,unit);
        if( tr != 3 ) {
            CSmallString error;
            error << "unable to read coordinate definition, id: " << i+1 << " (" << tr << " != 3)";
            RUNTIME_ERROR(error);
        }
        // some tests
        if(id != i+1) {
            CSmallString error;
            error << "coordinate id does not match, read: " << id << ", expected: " << i+1;
            RUNTIME_ERROR(error);
        }

        // init coordinate
        SetCV(i,name,type,min_value,max_value,nbins,fconv,unit);
    }

// alloc accumulator data -----------------------
    FinalizeAllocation();

// read accumulator data ------------------------

// ABF forces =================================================================

// samples
    for(int i=0; i < TotNBins; i++) {
        int ns = 0;
        if( fscanf(fin,"%d",&ns) != 1 ) {
            CSmallString error;
            error << "unable to read number of ABF samples for bin " << i;
            RUNTIME_ERROR(error);
        }
        NSamples[i] = ns;
    }

// do i=1,fnitem
//     read(iounit,40,end=100,err=100) (accumulator%ICFSum(i,j),j=1,accumulator%TotNBins)
// end do

// accumulated force
    for(int i = 0; i < NCVs; i++) {
        for(int j = 0; j < TotNBins; j++) {
            double cf = 0.0;
            if(fscanf(fin,"%lf",&cf) != 1) {
                CSmallString error;
                error << "unable to read ICFSum for bin " << j << " and item " << i;
                RUNTIME_ERROR(error);
            }
            ICFSum[map(i,j)] = cf;
        }
    }

// accumulated force square
    for(int i = 0; i < NCVs; i++) {
        for(int j = 0; j < TotNBins; j++) {
            double cf = 0.0;
            if(fscanf(fin,"%lf",&cf) != 1) {
                CSmallString error;
                error << "unable to read ICFSum2 for bin " << j << " and item " << i;
                RUNTIME_ERROR(error);
            }
            ICFSum2[map(i,j)] = cf;
        }
    }

    // EpotSum
    for(int i=0; i < TotNBins; i++) {
        double ns = 0;
        if( fscanf(fin,"%lf",&ns) != 1 ) {
            CSmallString error;
            error << "unable to read EpotSum for bin " << i;
            RUNTIME_ERROR(error);
        }
        EpotSum[i] = ns;
    }
    for(int i=0; i < TotNBins; i++) {
        double ns = 0;
        if( fscanf(fin,"%lf",&ns) != 1 ) {
            CSmallString error;
            error << "unable to read EpotSum2 for bin " << i;
            RUNTIME_ERROR(error);
        }
        EpotSum2[i] = ns;
    }

    EpotAvailable = true;
    ICFEpotSum.SetZero();
    ICFEpotSum2.SetZero();
}

//------------------------------------------------------------------------------

void CABFAccumulator::Load_v5(char* fline,FILE* fin)
{
    if( fin == NULL ) {
        INVALID_ARGUMENT("stream is not open");
    }

// read ABF accumulator header ------------------
    char abf_id[4];
    char ver_id[3];
    int  ncvs = 0;
    int  nr;

// first line can contain either two or four records for version 0

    nr = sscanf(fline,"%s %s %d",abf_id,ver_id,&ncvs);
    if( nr != 3 ){
        CSmallString error;
        error << "illegal header - three items expected (line: " << fline << ")";
        RUNTIME_ERROR(error);
    }

    abf_id[3]='\0';
    ver_id[2]='\0';

// check ID string
    if(strcmp(abf_id,"ABF") != 0) {
        CSmallString error;
        error << "'ABF' magic word was not found in the first line (line: " << fline << ")";
        RUNTIME_ERROR(error);
    }
    if(strcmp(ver_id,"V5") != 0) {
        CSmallString error;
        error << "only ABF V5 version is supported (line: " << fline << ")";
        RUNTIME_ERROR(error);
    }

    if(ncvs <= 0) {
        CSmallString error;
        error << "number of coordinates has to be greater than zero, but " << ncvs << " was found";
        RUNTIME_ERROR(error);
    }

    SetNumOfCVs(ncvs);

    CSmallString key;

    while( key.ReadLineFromFile(fin) ){

        key.Trim();

        if( key == "CVS") {

            // read coordinate specification
            for(int i=0; i < NCVs; i++) {
                int             id = 0;
                char            type[15];
                char            name[60];
                char            unit[40];
                double          min_value = 0.0;
                double          max_value = 0.0;
                double          fconv = 1.0;
                int             nbins = 0;
                int             tr = 0;

                memset(type,0,15);
                memset(name,0,60);
                memset(unit,0,40);

                // read item
                tr = fscanf(fin,"%d %10s %lf %lf %d",&id,type,&min_value,&max_value,&nbins);
                if( tr != 5 ) {
                    CSmallString error;
                    error << "unable to read coordinate definition, id: " << i+1 << " (" << tr << " != 5)";
                    RUNTIME_ERROR(error);
                }

                // some tests
                if(id != i+1) {
                    CSmallString error;
                    error << "coordinate id does not match, read: " << id << ", expected: " << i+1;
                    RUNTIME_ERROR(error);
                }
                if(max_value <= min_value) {
                    CSmallString error;
                    error << "min value is not smaller than max value, id: " << id;
                    RUNTIME_ERROR(error);
                }
                if(nbins <= 0) {
                    CSmallString error;
                    error << "number of bins has to be grater than zero, id: " << id;
                    RUNTIME_ERROR(error);
                }

                // read item
                tr = fscanf(fin,"%d %55s",&id,name);
                if( tr != 2 ) {
                    CSmallString error;
                    error << "unable to read coordinate definition, id: " << i+1 << " (" << tr << " != 2)";
                    RUNTIME_ERROR(error);
                }
                // some tests
                if(id != i+1) {
                    CSmallString error;
                    error << "coordinate id does not match, read: " << id << ", expected: " << i+1;
                    RUNTIME_ERROR(error);
                }

                // read item
                tr = fscanf(fin,"%d %lf %37s",&id,&fconv,unit);
                if( tr != 3 ) {
                    CSmallString error;
                    error << "unable to read coordinate definition, id: " << i+1 << " (" << tr << " != 3)";
                    RUNTIME_ERROR(error);
                }
                // some tests
                if(id != i+1) {
                    CSmallString error;
                    error << "coordinate id does not match, read: " << id << ", expected: " << i+1;
                    RUNTIME_ERROR(error);
                }

                // init coordinate
                SetCV(i,name,type,min_value,max_value,nbins,fconv,unit);
            }

            // alloc accumulator data
            FinalizeAllocation();
    // -----------------------------------------------------
        } else if( key == "TEMPERATURE" ) {
            // read item
            int tr = fscanf(fin,"%lf",&Temperature);
            if( tr != 1 ) {
                CSmallString error;
                error << "unable to read temperature (" << tr << " != 1)";
                RUNTIME_ERROR(error);
            }
    // -----------------------------------------------------
        } else if( key == "ENERGY" ) {
            char            unit[40];
            int             tr = 0;

            memset(unit,0,40);

            // read item
            tr = fscanf(fin,"%lf %37s",&EnergyFConv,unit);
            if( tr != 2 ) {
                CSmallString error;
                error << "unable to read energy unit (" << tr << " != 2)";
                RUNTIME_ERROR(error);
            }
            EnergyUnit = unit;
    // -----------------------------------------------------
        } else if( key == "NSAMPLES" ) {
            // samples
            for(int i=0; i < TotNBins; i++) {
                int ns = 0;
                if( fscanf(fin,"%d",&ns) != 1 ) {
                    CSmallString error;
                    error << "unable to read number of ABF samples for bin " << i;
                    RUNTIME_ERROR(error);
                }
                NSamples[i] = ns;
            }
    // -----------------------------------------------------
        } else if( key == "ICF" ) {
            // accumulated force
            for(int i = 0; i < NCVs; i++) {
                for(int j = 0; j < TotNBins; j++) {
                    double cf = 0.0;
                    if(fscanf(fin,"%lf",&cf) != 1) {
                        CSmallString error;
                        error << "unable to read ICFSum for bin " << j << " and item " << i;
                        RUNTIME_ERROR(error);
                    }
                    ICFSum[map(i,j)] = cf;
                }
            }

            // accumulated force square
            for(int i = 0; i < NCVs; i++) {
                for(int j = 0; j < TotNBins; j++) {
                    double cf = 0.0;
                    if(fscanf(fin,"%lf",&cf) != 1) {
                        CSmallString error;
                        error << "unable to read ICFSum2 for bin " << j << " and item " << i;
                        RUNTIME_ERROR(error);
                    }
                    ICFSum2[map(i,j)] = cf;
                }
            }
    // -----------------------------------------------------
        } else if( key == "EPOT" ) {
            // EpotSum
            for(int i=0; i < TotNBins; i++) {
                double ns = 0;
                if( fscanf(fin,"%lf",&ns) != 1 ) {
                    CSmallString error;
                    error << "unable to read EpotSum for bin " << i;
                    RUNTIME_ERROR(error);
                }
                EpotSum[i] = ns;
            }
            for(int i=0; i < TotNBins; i++) {
                double ns = 0;
                if( fscanf(fin,"%lf",&ns) != 1 ) {
                    CSmallString error;
                    error << "unable to read EpotSum2 for bin " << i;
                    RUNTIME_ERROR(error);
                }
                EpotSum2[i] = ns;
            }
            EpotAvailable = true;
    // -----------------------------------------------------
        } else if( key == "ICF*EPOT" ) {
            // accumulated force
            for(int i = 0; i < NCVs; i++) {
                for(int j = 0; j < TotNBins; j++) {
                    double cf = 0.0;
                    if(fscanf(fin,"%lf",&cf) != 1) {
                        CSmallString error;
                        error << "unable to read ICFEpotSum for bin " << j << " and item " << i;
                        RUNTIME_ERROR(error);
                    }
                    ICFEpotSum[map(i,j)] = cf;
                }
            }

            // accumulated force square
            for(int i = 0; i < NCVs; i++) {
                for(int j = 0; j < TotNBins; j++) {
                    double cf = 0.0;
                    if(fscanf(fin,"%lf",&cf) != 1) {
                        CSmallString error;
                        error << "unable to read ICFEpotSum2 for bin " << j << " and item " << i;
                        RUNTIME_ERROR(error);
                    }
                    ICFEpotSum2[map(i,j)] = cf;
                }
            }
        } else {
            CSmallString error;
            error << "unrecognized ABF accumulator keyword: '" << key << "'";
            RUNTIME_ERROR(error);
        }

        // read the rest of line and reset key
        key.ReadLineFromFile(fin);
        key = NULL;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABFAccumulator::Save(const CSmallString& name)
{
    FILE* fout = fopen(name, "w");

    if(fout == NULL) {
        CSmallString error;
        error << "unable to open file '" << name << "' (" << strerror(errno) << ")";
        RUNTIME_ERROR(error);
    }

    try {
        Save(fout);
    } catch(...) {
        fclose(fout);
        throw;
    }

    fclose(fout);
}

//------------------------------------------------------------------------------

void CABFAccumulator::Save(FILE* fout)
{
    if(fout == NULL) {
        INVALID_ARGUMENT("stream is not open");
    }

    if( TotNBins == 0 ) {
        CSmallString error;
        error << "no data in accumulator";
        RUNTIME_ERROR(error);
    }

// write ABF accumulator header ------------------
    if(fprintf(fout," ABF V5 %2d\n",NCVs) <= 0) {
        CSmallString error;
        error << "unable to write header";
        RUNTIME_ERROR(error);
    }

    if(fprintf(fout,"TEMPERATURE\n%18.11E\n",Temperature) <= 0) {
        CSmallString error;
        error << "unable to write temperature";
        RUNTIME_ERROR(error);
    }

    if(fprintf(fout,"ENERGY\n%18.11E %36s\n",EnergyFConv,(const char*)EnergyUnit) <= 0) {
        CSmallString error;
        error << "unable to write energy unit";
        RUNTIME_ERROR(error);
    }

// write coordinate specification ----------------
    if(fprintf(fout,"CVS\n") <= 0) {
        CSmallString error;
        error << "unable to write CVS header";
        RUNTIME_ERROR(error);
    }
    for(int i=0; i < NCVs; i++) {
        if(fprintf(fout,"%2d %10s %18.11E %18.11E %6d\n",i+1,
                   (const char*)CVs[i].Type,
                   CVs[i].MinValue,CVs[i].MaxValue,CVs[i].NBins) <= 0) {
            CSmallString error;
            error << "unable to write coordinate definition1 id: " << i+1;
            RUNTIME_ERROR(error);
        }
        if(fprintf(fout,"%2d %55s\n",i+1,(const char*)CVs[i].Name) <= 0) {
            CSmallString error;
            error << "unable to write coordinate definition2 id: " << i+1;
            RUNTIME_ERROR(error);
        }
        if(fprintf(fout,"%2d %18.11E %36s\n",i+1,CVs[i].FConv,(const char*)CVs[i].Unit) <= 0) {
            CSmallString error;
            error << "unable to write coordinate definition3 id: " << i+1;
            RUNTIME_ERROR(error);
        }
    }

// ABF forces =================================================================

// samples
    int counter;

    if(fprintf(fout,"NSAMPLES\n") <= 0) {
        CSmallString error;
        error << "unable to write NSAMPLES header";
        RUNTIME_ERROR(error);
    }
    counter = 0;
    for(int i=0; i < TotNBins; i++) {
        if(fprintf(fout,"%9d ",NSamples[i]) <= 0) {
            CSmallString error;
            error << "unable to write number of ABF samples for bin " << i;
            RUNTIME_ERROR(error);
        }
        if(counter % 8 == 7) fprintf(fout,"\n");
        counter++;
    }
    if(counter % 8 != 0) fprintf(fout,"\n");

// ICFSum
    if(fprintf(fout,"ICF\n") <= 0) {
        CSmallString error;
        error << "unable to write ICF header";
        RUNTIME_ERROR(error);
    }
    counter = 0;
    for(int i = 0; i < NCVs; i++) {
        for(int j = 0; j < TotNBins; j++) {
            if(fprintf(fout,"%19.11E ",ICFSum[map(i,j)]) <= 0) {
                CSmallString error;
                error << "unable to write ICFSum for bin " << j << " and item " << i;
                RUNTIME_ERROR(error);
            }
            if(counter % 4 == 3) fprintf(fout,"\n");
            counter++;
        }
    }
    if(counter % 4 != 0) fprintf(fout,"\n");

// ICFSum2
    counter = 0;
    for(int i = 0; i < NCVs; i++) {
        for(int j = 0; j < TotNBins; j++) {
            if(fprintf(fout,"%19.11E ",ICFSum2[map(i,j)]) <= 0) {
                CSmallString error;
                error << "unable to write ICFSum2 for bin " << j << " and item " << i;
                RUNTIME_ERROR(error);
            }
            if(counter % 4 == 3) fprintf(fout,"\n");
            counter++;
        }
    }
    if(counter % 4 != 0) fprintf(fout,"\n");

// EpotSum
    if(fprintf(fout,"EPOT\n") <= 0) {
        CSmallString error;
        error << "unable to write EPOT header";
        RUNTIME_ERROR(error);
    }
    counter = 0;
    for(int i=0; i < TotNBins; i++) {
        if(fprintf(fout,"%19.11E ",EpotSum[i]) <= 0) {
            CSmallString error;
            error << "unable to write EpotSum for bin " << i;
            RUNTIME_ERROR(error);
        }
        if(counter % 4 == 3) fprintf(fout,"\n");
        counter++;
    }
    if(counter % 4 != 0) fprintf(fout,"\n");

// EpotSum2
    counter = 0;
    for(int i=0; i < TotNBins; i++) {
        if(fprintf(fout,"%19.11E ",EpotSum2[i]) <= 0) {
            CSmallString error;
            error << "unable to write EpotSum2 for bin " << i;
            RUNTIME_ERROR(error);
        }
        if(counter % 4 == 3) fprintf(fout,"\n");
        counter++;
    }
    if(counter % 4 != 0) fprintf(fout,"\n");

// ICFEpotSum
    if(fprintf(fout,"ICF*EPOT\n") <= 0) {
        CSmallString error;
        error << "unable to write ICF*EPOT header";
        RUNTIME_ERROR(error);
    }
    counter = 0;
    for(int i = 0; i < NCVs; i++) {
        for(int j = 0; j < TotNBins; j++) {
            if(fprintf(fout,"%19.11E ",ICFEpotSum[map(i,j)]) <= 0) {
                CSmallString error;
                error << "unable to write ICFEpotSum for bin " << j << " and item " << i;
                RUNTIME_ERROR(error);
            }
            if(counter % 4 == 3) fprintf(fout,"\n");
            counter++;
        }
    }
    if(counter % 4 != 0) fprintf(fout,"\n");

// ICFEpotSum2
    counter = 0;
    for(int i = 0; i < NCVs; i++) {
        for(int j = 0; j < TotNBins; j++) {
            if(fprintf(fout,"%19.11E ",ICFEpotSum2[map(i,j)]) <= 0) {
                CSmallString error;
                error << "unable to write ICFEpotSum2 for bin " << j << " and item " << i;
                RUNTIME_ERROR(error);
            }
            if(counter % 4 == 3) fprintf(fout,"\n");
            counter++;
        }
    }
    if(counter % 4 != 0) fprintf(fout,"\n");
}

//------------------------------------------------------------------------------

void CABFAccumulator::SaveMask(FILE* fout)
{
    if(fout == NULL) {
        INVALID_ARGUMENT("stream is not open");
    }

    if( TotNBins == 0 ) {
        CSmallString error;
        error << "no data in accumulator";
        RUNTIME_ERROR(error);
    }

// 10  format(A4,1X,A2,1X,I2)
// 20  format(I2,1X,A10,1X,E18.11,1X,E18.11,1X,I6)
// 25  format(I2,1X,A55)
// 40  format(4(E19.11,1X))

// write ABF accumulator header ------------------
    if(fprintf(fout,"MABF V3 %2d\n",NCVs) <= 0) {
        CSmallString error;
        error << "unable to write header";
        RUNTIME_ERROR(error);
    }

// write coordinate specification ----------------
    for(int i=0; i < NCVs; i++) {
        if(fprintf(fout,"%2d %-10s %18.11e %18.11E %6d\n",i+1,
                   (const char*)CVs[i].Type,
                   CVs[i].MinValue,CVs[i].MaxValue,CVs[i].NBins) <= 0) {
            CSmallString error;
            error << "unable to write coordinate definition id: " << i+1;
            RUNTIME_ERROR(error);
        }
        if(fprintf(fout,"%2d %-55s\n",i+1,(const char*)CVs[i].Name) <= 0) {
            CSmallString error;
            error << "unable to write coordinate definition id: " << i+1;
            RUNTIME_ERROR(error);
        }
    }

// ABF mask =================================================================

// mask
    int counter;

    counter = 0;
    for(int i=0; i < TotNBins; i++) {
        if(fprintf(fout,"%19.11E ",Mask[i]) <= 0) {
            CSmallString error;
            error << "unable to write mask weight for bin " << i;
            RUNTIME_ERROR(error);
        }
        if(counter % 4 == 3) fprintf(fout,"\n");
        counter++;
    }
    if(counter % 4 != 0) fprintf(fout,"\n");
}

//------------------------------------------------------------------------------

void CABFAccumulator::SaveMaskGnuPlot(FILE* fout)
{
    if(fout == NULL) {
        INVALID_ARGUMENT("stream is not open");
    }

    if( TotNBins == 0 ) {
        CSmallString error;
        error << "no data in accumulator";
        RUNTIME_ERROR(error);
    }

    CSimpleVector<double> point;
    point.CreateVector(NCVs);

    for(int i=0; i < TotNBins; i++) {
        GetPoint(i,point);
        for(int j=0; j < NCVs; j++){
            if(fprintf(fout,"%19.11E ",point[j]) <= 0 ){
                CSmallString error;
                error << "unable to write point coordinate for bin " << i;
                RUNTIME_ERROR(error);
            }
        }
        if(fprintf(fout,"%19.11E\n",Mask[i]) <= 0) {
            CSmallString error;
            error << "unable to write mask weight for bin " << i;
            RUNTIME_ERROR(error);
        }
        if( i % CVs[0].GetNumOfBins() == CVs[0].GetNumOfBins() -1 ){
            fprintf(fout,"\n");
        }
    }
}

//------------------------------------------------------------------------------

void CABFAccumulator::GetPoint(unsigned int index,CSimpleVector<double>& point) const
{
    for(int k=NCVs-1; k >= 0; k--) {
        const CColVariable* p_coord = &CVs[k];
        int ibin = index % p_coord->GetNumOfBins();
        point[k] = p_coord->GetValue(ibin);
        index = index / p_coord->GetNumOfBins();
    }
}

void CABFAccumulator::GetPointRValues(unsigned int index,CSimpleVector<double>& point) const
{
    for(int k=NCVs-1; k >= 0; k--) {
        const CColVariable* p_coord = &CVs[k];
        int ibin = index % p_coord->GetNumOfBins();
        point[k] = p_coord->GetRValue(ibin);
        index = index / p_coord->GetNumOfBins();
    }
}

//------------------------------------------------------------------------------

void CABFAccumulator::GetIPoint(unsigned int index,CSimpleVector<int>& point) const
{
    for(int k=NCVs-1; k >= 0; k--) {
        const CColVariable* p_coord = &CVs[k];
        int ibin = index % p_coord->GetNumOfBins();
        point[k] = ibin;
        index = index / p_coord->GetNumOfBins();
    }
}

//------------------------------------------------------------------------------

void CABFAccumulator::SetTemperature(double temp)
{
    Temperature = temp;
}

//------------------------------------------------------------------------------

double CABFAccumulator::GetTemperature(void)
{
    return(Temperature);
}

//------------------------------------------------------------------------------

void CABFAccumulator::SetEnergyUnit(double fconv,const CSmallString& unit)
{
    EnergyFConv = fconv;
    EnergyUnit  = unit;
}

//------------------------------------------------------------------------------

void CABFAccumulator::SetEpotEnabled(bool enabled)
{
    EpotAvailable = enabled;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABFAccumulator::SetNumOfCVs(int ncvs)
{
    if( ncvs < 0 ) {
        INVALID_ARGUMENT("ncvs < 0");
    }

    if(NCVs > 0) {
        // destroy all previous data
        Deallocate();
        CVs.FreeVector();
        NCVs = 0;
    }

// try to allocate CVs array
    if( ncvs > 0 ) {
        CVs.CreateVector(ncvs);
    }

// all seems to be fine - update items
    NCVs = ncvs;
}

//------------------------------------------------------------------------------

void CABFAccumulator::SetCV(int id,
                                    const CSmallString& name,
                                    const CSmallString& type,
                                    double min_value,double max_value,int nbins)
{
    if( CVs == NULL ){
        RUNTIME_ERROR("CVs is NULL");
    }
    if( id < 0 || id >= NCVs ){
        INVALID_ARGUMENT("id out-of-range");
    }

    if( nbins <= 0 ){
        INVALID_ARGUMENT("nbins <= 0");
    }
    if( max_value < min_value ){
        INVALID_ARGUMENT("max_value < min_value");
    }

    if( NSamples != NULL ) {
        // it was already finalized - destroy data
        Deallocate();
    }

    CVs[id].ID = id;
    CVs[id].Name = name;
    CVs[id].Type = type;
    CVs[id].MinValue = min_value;
    CVs[id].MaxValue = max_value;

    CVs[id].NBins = nbins;
    CVs[id].BinWidth = (CVs[id].MaxValue - CVs[id].MinValue)/CVs[id].NBins;
    CVs[id].Width = CVs[id].MaxValue - CVs[id].MinValue;
}

//------------------------------------------------------------------------------

void CABFAccumulator::SetCV(int id,
                                    const CSmallString& name,
                                    const CSmallString& type,
                                    double min_value,double max_value,int nbins,
                                    double fconv, const CSmallString& unit)
{
    if( CVs == NULL ){
        RUNTIME_ERROR("CVs is NULL");
    }
    if( id < 0 || id >= NCVs ){
        INVALID_ARGUMENT("id out-of-range");
    }

    if( nbins <= 0 ){
        INVALID_ARGUMENT("nbins <= 0");
    }
    if( max_value < min_value ){
        INVALID_ARGUMENT("max_value < min_value");
    }

    if( NSamples != NULL ) {
        // it was already finalized - destroy data
        Deallocate();
    }

    CVs[id].ID = id;
    CVs[id].Name = name;
    CVs[id].Type = type;
    CVs[id].Unit = unit;
    CVs[id].MinValue = min_value;
    CVs[id].MaxValue = max_value;
    CVs[id].FConv = fconv;

    CVs[id].NBins = nbins;
    CVs[id].BinWidth = (CVs[id].MaxValue - CVs[id].MinValue)/CVs[id].NBins;
    CVs[id].Width = CVs[id].MaxValue - CVs[id].MinValue;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CABFAccumulator::GetNumOfCVs(void) const
{
    return(NCVs);
}

//------------------------------------------------------------------------------

const CColVariable* CABFAccumulator::GetCV(int cv) const
{
    return(&CVs[cv]);
}

//------------------------------------------------------------------------------

int CABFAccumulator::GetNumOfBins(int limit) const
{
    if( limit == 0 ) return(TotNBins);
    int number = 0;
    for(int i=0; i < TotNBins; i++) {
        if(NSamples[i] >= limit) number++;
    }
    return(number);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CABFAccumulator::GetNumOfSamples(const CSimpleVector<int>& position) const
{
    int glbindex = 0;
    for(int i=0; i < NCVs; i++) {
        glbindex = glbindex*CVs[i].GetNumOfBins() + position[i];
    }
    return(NSamples[glbindex]);
}

//------------------------------------------------------------------------------

int CABFAccumulator::GetNumOfSamples(int ibin) const
{
    return(NSamples[ibin]);
}

//------------------------------------------------------------------------------

double CABFAccumulator::GetMaskWeight(int ibin) const
{
    return(Mask[ibin]);
}

//------------------------------------------------------------------------------

void CABFAccumulator::SetMaskWeight(int ibin,double weight)
{
    Mask[ibin] = weight;
}

//------------------------------------------------------------------------------

const double& CABFAccumulator::GetICFSum(int icv,int ibin) const
{
    int glbindex = map(icv,ibin);
    return(ICFSum[glbindex]);
}

//------------------------------------------------------------------------------

const double& CABFAccumulator::GetICFSum2(int icv,int ibin) const
{
    int glbindex = map(icv,ibin);
    return(ICFSum2[glbindex]);
}

//------------------------------------------------------------------------------

const double& CABFAccumulator::GetEpotSum(int ibin) const
{
    return(EpotSum[ibin]);
}

//------------------------------------------------------------------------------

const double& CABFAccumulator::GetEpotSum2(int ibin) const
{
    return(EpotSum2[ibin]);
}

//------------------------------------------------------------------------------

const double& CABFAccumulator::GetICFEpotSum(int icv,int ibin) const
{
    int glbindex = map(icv,ibin);
    return(ICFEpotSum[glbindex]);
}

//------------------------------------------------------------------------------

const double& CABFAccumulator::GetICFEpotSum2(int icv,int ibin) const
{
    int glbindex = map(icv,ibin);
    return(ICFEpotSum2[glbindex]);
}

//------------------------------------------------------------------------------

void CABFAccumulator::SetICFSum(int icv,int ibin,const double& value)
{
    int glbindex = map(icv,ibin);
    ICFSum[glbindex] = value;
}

//------------------------------------------------------------------------------

void CABFAccumulator::SetICFSum2(int icv,int ibin,const double& value)
{
    int glbindex = map(icv,ibin);
    ICFSum2[glbindex] = value;
}

//------------------------------------------------------------------------------

void CABFAccumulator::SetNumberOfABFSamples(int ibin,const int& value)
{
    NSamples[ibin] = value;
}

//------------------------------------------------------------------------------

int* CABFAccumulator::GetNSamplesArray(void)
{
    return(NSamples);
}

//------------------------------------------------------------------------------

double* CABFAccumulator::GetICFSumArray(void)
{
    return(ICFSum);
}

//------------------------------------------------------------------------------

double* CABFAccumulator::GetICFSum2Array(void)
{
    return(ICFSum2);
}

//------------------------------------------------------------------------------

double* CABFAccumulator::GetEpotSumArray(void)
{
    return(EpotSum);
}

//------------------------------------------------------------------------------

double* CABFAccumulator::GetEpotSum2Array(void)
{
    return(EpotSum2);
}

//------------------------------------------------------------------------------

double* CABFAccumulator::GetICFEpotSumArray(void)
{
    return(ICFEpotSum);
}

//------------------------------------------------------------------------------

double* CABFAccumulator::GetICFEpotSum2Array(void)
{
    return(ICFEpotSum2);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CABFAccumulator::GetGlobalIndex(const CSimpleVector<int>& position) const
{
    int glbindex = 0;
    for(int i=0; i < NCVs; i++) {
        if( position[i] < 0 ) return(-1);
        if( position[i] >= (int)CVs[i].GetNumOfBins() ) return(-1);
        glbindex = glbindex*CVs[i].GetNumOfBins() + position[i];
    }
    return(glbindex);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABFAccumulator::Clear(void)
{
    if(NSamples == NULL) return;   // allocation was not finalized

    for(int i=0; i < TotNBins; i++) NSamples[i] = 0;

    for(int i=0; i < TotNBins*NCVs; i++) ICFSum[i] = 0.0;
    for(int i=0; i < TotNBins*NCVs; i++) ICFSum2[i] = 0.0;

    for(int i=0; i < TotNBins; i++) EpotSum[i] = 0;
    for(int i=0; i < TotNBins; i++) EpotSum2[i] = 0;

    for(int i=0; i < TotNBins*NCVs; i++) ICFEpotSum[i] = 0;
    for(int i=0; i < TotNBins*NCVs; i++) ICFEpotSum2[i] = 0;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABFAccumulator::FinalizeAllocation()
{
    if(NCVs <= 0) return;        // at least one coordinate is required
    if(CVs == NULL) return;       // CVs array is not allocated
    if(NSamples != NULL) return;    // already finalized!

// check if there is non zero number of bins per coordinate
    TotNBins = 1;
    for(int i=0; i < NCVs; i++) {
        TotNBins *= CVs[i].NBins;
    }

// and now allocate data arrays
    NSamples.CreateVector(TotNBins);
    Mask.CreateVector(TotNBins);
    for(int i=0; i < TotNBins; i++){
        Mask[i] = 1.0;
    }
    ICFSum.CreateVector(TotNBins*NCVs);
    ICFSum2.CreateVector(TotNBins*NCVs);

// EpotSum
    EpotSum.CreateVector(TotNBins);
    EpotSum2.CreateVector(TotNBins);

// abfforce*EpotSum
    ICFEpotSum.CreateVector(TotNBins*NCVs);
    ICFEpotSum2.CreateVector(TotNBins*NCVs);

    // reset data
    Clear();
}

//------------------------------------------------------------------------------

void CABFAccumulator::Deallocate(void)
{
// do not destroy CVs array !

// destroy only data arrays
    NSamples.FreeVector();

    ICFSum.FreeVector();
    ICFSum2.FreeVector();

    EpotSum.FreeVector();
    EpotSum2.FreeVector();

    ICFEpotSum.FreeVector();
    ICFEpotSum2.FreeVector();

    TotNBins        = 0;
    EpotAvailable   = false;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABFAccumulator::LoadCVSInfo(CXMLElement* p_iele)
{
    if(p_iele == NULL) {
        INVALID_ARGUMENT("p_iele is NULL");
    }

    Deallocate();

    CXMLElement* p_ele = p_iele->GetFirstChildElement("CVS");
    if(p_ele == NULL) {
        RUNTIME_ERROR("unable to get CVS element");
    }

    int lnitems = 0;
    if( p_ele->GetAttribute("ncvs",lnitems) == false) {
        RUNTIME_ERROR("unable to get header attributes");
    }

    if( lnitems == 0 ) {
        // no data
        return;
    }

    SetNumOfCVs(lnitems);

    CXMLElement*   p_cel = p_ele->GetFirstChildElement("CV");

    int            ccount = 0;

    while(p_cel != NULL) {
        if( ccount >= lnitems ) {
            LOGIC_ERROR("more CV elements than NCVs");
        }
        CVs[ccount].LoadInfo(p_cel);
        ccount++;
        p_cel = p_cel->GetNextSiblingElement("CV");
    }

// alloc accumulator data -----------------------
    FinalizeAllocation();

// clear all data
    Clear();
}

//------------------------------------------------------------------------------

bool CABFAccumulator::CheckCVSInfo(CXMLElement* p_iele) const
{
    if(p_iele == NULL) {
        INVALID_ARGUMENT("p_iele is NULL");
    }

    CXMLElement* p_ele = p_iele->GetFirstChildElement("CVS");
    if(p_ele == NULL) {
        ES_ERROR("unable to get CVS element");
        return(false);
    }

    bool result = true;

    int lnitems;

    result &= p_ele->GetAttribute("ncvs",lnitems);
    if(result == false) {
        ES_ERROR("unable to get header attributes");
        return(false);
    }

    if(lnitems != NCVs) {
        ES_ERROR("mismatch in the number of coordinates");
        return(false);
    }

    CXMLElement*   p_cel = NULL;
    if(p_ele != NULL) p_cel = p_ele->GetFirstChildElement("CV");
    int            ccount = 0;

    while(p_cel != NULL) {
        if(ccount >= lnitems) {
            ES_ERROR("more COORD elements than NCVs");
            return(false);
        }
        if( CVs[ccount].CheckInfo(p_cel) == false ) {
            CSmallString error;
            error << "mismatch in cv: " << ccount+1;
            ES_ERROR(error);
            return(false);
        }
        ccount++;
        p_cel = p_cel->GetNextSiblingElement("CV");
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CABFAccumulator::CheckCVSInfo(const CABFAccumulator* p_accu) const
{
    if(p_accu == NULL) {
        INVALID_ARGUMENT("p_iele is NULL");
    }

    if(p_accu->NCVs != NCVs) {
        ES_ERROR("mismatch in the number of coordinates");
        return(false);
    }

    for(int i=0; i < NCVs; i++) {
        if(CVs[i].CheckInfo(&p_accu->CVs[i]) == false) {
            CSmallString error;
            error << "mismatch in cv: " << i+1;
            ES_ERROR(error);
            return(false);
        }
    }

    return(true);
}

//------------------------------------------------------------------------------

void CABFAccumulator::SaveCVSInfo(CXMLElement* p_tele) const
{
    if(p_tele == NULL) {
        INVALID_ARGUMENT("p_tele is NULL");
    }

    CXMLElement* p_ele = p_tele->CreateChildElement("CVS");

    p_ele->SetAttribute("ncvs",NCVs);

    for(int i=0; i < NCVs; i++) {
        CXMLElement* p_iele = p_ele->CreateChildElement("CV");
        CVs[i].SaveInfo(p_iele);
    }
}

//------------------------------------------------------------------------------

void CABFAccumulator::PrintCVSInfo(std::ostream& vout)
{
    vout  << endl;
    vout << "=== Collective Variables =======================================================" << endl;
    vout << endl;
    vout << "ID P Type       Unit  Name                       Min value   Max value   NBins  " << endl;
    vout << "-- - ---------- ----- -------------------------- ----------- ----------- -------" << endl;

    for(int i=0; i < NCVs; i++) {
        CVs[i].PrintInfo(vout);
    }
}

//------------------------------------------------------------------------------

void CABFAccumulator::ReadABFData(CXMLElement* p_rele)
{
    if(p_rele == NULL) {
        INVALID_ARGUMENT("p_iele is NULL");
    }

    unsigned int inc_nsamples_size      = GetNumOfBins()*sizeof(int);
    unsigned int inc_icfsum_size        = GetNumOfBins()*GetNumOfCVs()*sizeof(double);
    unsigned int inc_epotsum_size       = GetNumOfBins()*sizeof(double);
    unsigned int inc_icfepotsum_size    = GetNumOfBins()*GetNumOfCVs()*sizeof(double);

    if( (inc_nsamples_size == 0) || (inc_icfsum_size == 0) || (inc_epotsum_size == 0) || (inc_icfepotsum_size == 0) ) {
        RUNTIME_ERROR("size(s) is(are) zero");
    }

// --------------------------
    CXMLBinData* p_nisamples = p_rele->GetFirstChildBinData("NSAMPLES");
    if(p_nisamples == NULL) {
        RUNTIME_ERROR("unable to open NSAMPLES element");
    }

    void* p_data = p_nisamples->GetData();
    if((GetNSamplesArray() == NULL) || (p_data == NULL) || (p_nisamples->GetLength() != inc_nsamples_size)) {
        CSmallString error;
        error << "inconsistent NSAMPLES element dat (r: " << p_nisamples->GetLength()
              << ", l: " <<  inc_nsamples_size << ")";
        RUNTIME_ERROR(error);
    }
    memcpy(GetNSamplesArray(),p_data,inc_nsamples_size);

// --------------------------
    CXMLBinData* p_inc_icfsum = p_rele->GetFirstChildBinData("ICFSUM");
    if(p_inc_icfsum == NULL) {
        RUNTIME_ERROR("unable to open ICFSUM element");
    }

    p_data = p_inc_icfsum->GetData();
    if((GetICFSumArray() == NULL) || (p_data == NULL) || (p_inc_icfsum->GetLength() != inc_icfsum_size)) {
        CSmallString error;
        error << "inconsistent ICFSUM element dat (r: " << p_inc_icfsum->GetLength()
              << ", l: " <<  inc_icfsum_size << ")";
        RUNTIME_ERROR(error);
    }
    memcpy(GetICFSumArray(),p_data,inc_icfsum_size);

// --------------------------
    CXMLBinData* p_inc_icfsum2 = p_rele->GetFirstChildBinData("ICFSUM2");
    if(p_inc_icfsum2 == NULL) {
        RUNTIME_ERROR("unable to open ICFSUM2 element");
    }

    p_data = p_inc_icfsum2->GetData();
    if((GetICFSum2Array() == NULL) || (p_data == NULL) || (p_inc_icfsum2->GetLength() != inc_icfsum_size)) {
        CSmallString error;
        error << "inconsistent ICFSUM2 element dat (r: " << p_inc_icfsum2->GetLength()
              << ", l: " <<  inc_icfsum_size << ")";
        RUNTIME_ERROR(error);
    }
    memcpy(GetICFSum2Array(),p_data,inc_icfsum_size);

// --------------------------
    CXMLBinData* p_inc_epotsum = p_rele->GetFirstChildBinData("EPOTSUM");
    if(p_inc_epotsum == NULL) {
        RUNTIME_ERROR("unable to open EPOTSUM element");
    }

    p_data = p_inc_epotsum->GetData();
    if((GetEpotSumArray() == NULL) || (p_data == NULL) || (p_inc_epotsum->GetLength() != inc_epotsum_size)) {
        CSmallString error;
        error << "inconsistent EPOTSUM element dat (r: " << p_inc_epotsum->GetLength()
              << ", l: " <<  inc_epotsum_size << ")";
        RUNTIME_ERROR(error);
    }
    memcpy(GetEpotSumArray(),p_data,inc_epotsum_size);

// --------------------------
    CXMLBinData* p_inc_epotsum2 = p_rele->GetFirstChildBinData("EPOTSUM2");
    if(p_inc_epotsum2 == NULL) {
        RUNTIME_ERROR("unable to open EPOTSUM2 element");
    }

    p_data = p_inc_epotsum2->GetData();
    if((GetEpotSum2Array() == NULL) || (p_data == NULL) || (p_inc_epotsum2->GetLength() != inc_epotsum_size)) {
        CSmallString error;
        error << "inconsistent EPOTSUM2 element dat (r: " << p_inc_epotsum2->GetLength()
              << ", l: " <<  inc_epotsum_size << ")";
        RUNTIME_ERROR(error);
    }
    memcpy(GetEpotSum2Array(),p_data,inc_epotsum_size);

// --------------------------
    CXMLBinData* p_inc_icfepotsum = p_rele->GetFirstChildBinData("ICFEPOTSUM");
    if(p_inc_icfepotsum == NULL) {
        RUNTIME_ERROR("unable to open ICFEPOTSUM element");
    }

    p_data = p_inc_icfepotsum->GetData();
    if((GetICFSumArray() == NULL) || (p_data == NULL) || (p_inc_icfepotsum->GetLength() != inc_icfepotsum_size)) {
        CSmallString error;
        error << "inconsistent ICFEPOTSUM element dat (r: " << p_inc_icfepotsum->GetLength()
              << ", l: " <<  inc_icfepotsum_size << ")";
        RUNTIME_ERROR(error);
    }
    memcpy(GetICFSumArray(),p_data,inc_icfepotsum_size);

// --------------------------
    CXMLBinData* p_inc_icfepotsum2 = p_rele->GetFirstChildBinData("ICFEPOTSUM2");
    if(p_inc_icfepotsum2 == NULL) {
        RUNTIME_ERROR("unable to open ICFEPOTSUM2 element");
    }

    p_data = p_inc_icfepotsum2->GetData();
    if((GetICFSum2Array() == NULL) || (p_data == NULL) || (p_inc_icfepotsum2->GetLength() != inc_icfepotsum_size)) {
        CSmallString error;
        error << "inconsistent ICFEPOTSUM2 element dat (r: " << p_inc_icfepotsum2->GetLength()
              << ", l: " <<  inc_icfepotsum_size << ")";
        RUNTIME_ERROR(error);
    }
    memcpy(GetICFSum2Array(),p_data,inc_icfepotsum_size);
}

//------------------------------------------------------------------------------

void CABFAccumulator::AddABFData(CXMLElement* p_cele)
{
    if(p_cele == NULL) {
        INVALID_ARGUMENT("p_iele is NULL");
    }

    unsigned int inc_nsamples           = GetNumOfBins();
    unsigned int inc_icfsum             = GetNumOfBins()*GetNumOfCVs();
    unsigned int inc_epotsum            = GetNumOfBins();
    unsigned int inc_icfepotsum         = GetNumOfBins()*GetNumOfCVs();

    unsigned int inc_nsamples_size      = GetNumOfBins()*sizeof(int);
    unsigned int inc_icfsum_size        = GetNumOfBins()*GetNumOfCVs()*sizeof(double);
    unsigned int inc_epotsum_size       = GetNumOfBins()*sizeof(double);
    unsigned int inc_icfepotsum_size    = GetNumOfBins()*GetNumOfCVs()*sizeof(double);

// --------------------------
    CXMLBinData* p_isamples = p_cele->GetFirstChildBinData("NSAMPLES");
    if(p_isamples == NULL) {
        RUNTIME_ERROR("unable to open NSAMPLES element");
    }

    void* p_data = p_isamples->GetData();
    if((GetNSamplesArray() == NULL) || (p_data == NULL) || (p_isamples->GetLength() != inc_nsamples_size)) {
        CSmallString error;
        error << "inconsistent NSAMPLES element dat (r: " << p_isamples->GetLength()
              << ", l: " <<  inc_nsamples_size << ")";
        RUNTIME_ERROR(error);
    }

    int* idst = GetNSamplesArray();
    int* isrc = (int*)p_data;
    for(unsigned int i=0; i < inc_nsamples; i++) {
        *idst++ += *isrc++;
    }

// --------------------------
    CXMLBinData* p_inc_icfsum = p_cele->GetFirstChildBinData("ICFSUM");
    if(p_inc_icfsum == NULL) {
        RUNTIME_ERROR("unable to open ICFSUM element");
    }

    p_data = p_inc_icfsum->GetData();
    if((GetICFSumArray() == NULL) || (p_data == NULL) || (p_inc_icfsum->GetLength() != inc_icfsum_size)) {
        CSmallString error;
        error << "inconsistent ICFSUM element dat (r: " << p_inc_icfsum->GetLength()
              << ", l: " <<  inc_icfsum_size << ")";
        RUNTIME_ERROR(error);
    }

    double* jdst = GetICFSumArray();
    double* jsrc = (double*)p_data;
    for(unsigned int i=0; i < inc_icfsum; i++) {
        *jdst++ += *jsrc++;
    }

// --------------------------
    CXMLBinData* p_inc_icfsum2 = p_cele->GetFirstChildBinData("ICFSUM2");
    if(p_inc_icfsum2 == NULL) {
        RUNTIME_ERROR("unable to open ICFSUM2 element");
    }

    p_data = p_inc_icfsum2->GetData();
    if((GetICFSum2Array() == NULL) || (p_data == NULL) || (p_inc_icfsum2->GetLength() != inc_icfsum_size)) {
        CSmallString error;
        error << "inconsistent ICFSUM2 element dat (r: " << p_inc_icfsum2->GetLength()
              << ", l: " <<  inc_icfsum_size << ")";
        RUNTIME_ERROR(error);
    }

    double* kdst = GetICFSum2Array();
    double* ksrc = (double*)p_data;
    for(unsigned int i=0; i < inc_icfsum; i++) {
        *kdst++ += *ksrc++;
    }

// --------------------------
    CXMLBinData* p_inc_EpotSum = p_cele->GetFirstChildBinData("EPOTSUM");
    if(p_inc_EpotSum == NULL) {
        RUNTIME_ERROR("unable to open EPOTSUM element");
    }

    p_data = p_inc_EpotSum->GetData();
    if((GetEpotSumArray() == NULL) || (p_data == NULL) || (p_inc_EpotSum->GetLength() != inc_epotsum_size)) {
        CSmallString error;
        error << "inconsistent EPOTSUM element dat (r: " << p_inc_EpotSum->GetLength()
              << ", l: " <<  inc_epotsum_size << ")";
        RUNTIME_ERROR(error);
    }

    double* ldst = GetEpotSumArray();
    double* lsrc = (double*)p_data;
    for(unsigned int i=0; i < inc_epotsum; i++) {
        *ldst++ += *lsrc++;
    }

// --------------------------
    CXMLBinData* p_inc_epotsum = p_cele->GetFirstChildBinData("EPOTSUM2");
    if(p_inc_epotsum == NULL) {
        RUNTIME_ERROR("unable to open EPOTSUM2 element");
    }

    p_data = p_inc_epotsum->GetData();
    if((GetEpotSum2Array() == NULL) || (p_data == NULL) || (p_inc_epotsum->GetLength() != inc_epotsum_size)) {
        CSmallString error;
        error << "inconsistent EPOTSUM2 element dat (r: " << p_inc_epotsum->GetLength()
              << ", l: " <<  inc_epotsum_size << ")";
        RUNTIME_ERROR(error);
    }

    double* mdst = GetEpotSum2Array();
    double* msrc = (double*)p_data;
    for(unsigned int i=0; i < inc_epotsum; i++) {
        *mdst++ += *msrc++;
    }

// --------------------------
    CXMLBinData* p_inc_icfepotsum = p_cele->GetFirstChildBinData("ICFEPOTSUM");
    if(p_inc_icfepotsum == NULL) {
        RUNTIME_ERROR("unable to open ICFEPOTSUM element");
    }

    p_data = p_inc_icfsum->GetData();
    if((GetICFSumArray() == NULL) || (p_data == NULL) || (p_inc_icfepotsum->GetLength() != inc_icfepotsum_size)) {
        CSmallString error;
        error << "inconsistent ICFEPOTSUM element dat (r: " << p_inc_icfepotsum->GetLength()
              << ", l: " <<  inc_icfepotsum_size << ")";
        RUNTIME_ERROR(error);
    }

    double* ndst = GetICFEpotSumArray();
    double* nsrc = (double*)p_data;
    for(unsigned int i=0; i < inc_icfepotsum; i++) {
        *ndst++ += *nsrc++;
    }

// --------------------------
    CXMLBinData* p_inc_icfepotsum2 = p_cele->GetFirstChildBinData("ICFEPOTSUM2");
    if(p_inc_icfepotsum2 == NULL) {
        RUNTIME_ERROR("unable to open ICFEPOTSUM2 element");
    }

    p_data = p_inc_icfepotsum2->GetData();
    if((GetICFSum2Array() == NULL) || (p_data == NULL) || (p_inc_icfepotsum2->GetLength() != inc_icfepotsum_size)) {
        CSmallString error;
        error << "inconsistent ICFEPOTSUM2 element dat (r: " << p_inc_icfepotsum2->GetLength()
              << ", l: " <<  inc_icfepotsum_size << ")";
        RUNTIME_ERROR(error);
    }

    double* odst = GetICFEpotSum2Array();
    double* osrc = (double*)p_data;
    for(unsigned int i=0; i < inc_icfepotsum; i++) {
        *odst++ += *osrc++;
    }
}

//------------------------------------------------------------------------------

void CABFAccumulator::WriteABFData(CXMLElement* p_rele)
{
    if(p_rele == NULL) {
        INVALID_ARGUMENT("p_rele is NULL");
    }

    unsigned int inc_nsamples_size      = GetNumOfBins()*sizeof(int);
    unsigned int inc_icfsum_size        = GetNumOfBins()*GetNumOfCVs()*sizeof(double);
    unsigned int inc_epotsum_size       = GetNumOfBins()*sizeof(double);
    unsigned int inc_icfepotsum_size    = GetNumOfBins()*GetNumOfCVs()*sizeof(double);

// write collected data to response
    CXMLBinData* p_isamples = p_rele->CreateChildBinData("NSAMPLES");
    p_isamples->CopyData(GetNSamplesArray(),inc_nsamples_size);

    CXMLBinData* p_inc_icfsum = p_rele->CreateChildBinData("ICFSUM");
    p_inc_icfsum->CopyData(GetICFSumArray(),inc_icfsum_size);

    CXMLBinData* p_inc_icfsum2 = p_rele->CreateChildBinData("ICFSUM2");
    p_inc_icfsum2->CopyData(GetICFSum2Array(),inc_icfsum_size);

    CXMLBinData* p_inc_epotsum = p_rele->CreateChildBinData("EPOTSUM");
    p_inc_epotsum->CopyData(GetEpotSumArray(),inc_epotsum_size);

    CXMLBinData* p_inc_epotsum2 = p_rele->CreateChildBinData("EPOTSUM2");
    p_inc_epotsum2->CopyData(GetEpotSum2Array(),inc_epotsum_size);

    CXMLBinData* p_inc_icfepotsum = p_rele->CreateChildBinData("ICFEPOTSUM");
    p_inc_icfepotsum->CopyData(GetICFEpotSumArray(),inc_icfepotsum_size);

    CXMLBinData* p_inc_icfepotsum2 = p_rele->CreateChildBinData("ICFEPOTSUM2");
    p_inc_icfepotsum2->CopyData(GetICFEpotSum2Array(),inc_icfepotsum_size);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CABFAccumulator::map(int item,int bin) const
{
    return(item + bin*NCVs);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABFAccumulator::AddABFAccumulator(const CABFAccumulator* p_accu)
{
    if( p_accu == NULL ) {
        INVALID_ARGUMENT("p_accu is NULL");
    }

    if( CheckCVSInfo(p_accu) == false ) {
        RUNTIME_ERROR("accumulators are not compatible");
    }

    for(int i=0; i < TotNBins; i++) NSamples[i] += p_accu->NSamples[i];

    for(int i=0; i < TotNBins*NCVs; i++) ICFSum[i]  += p_accu->ICFSum[i];
    for(int i=0; i < TotNBins*NCVs; i++) ICFSum2[i] += p_accu->ICFSum2[i];

    for(int i=0; i < TotNBins; i++) EpotSum[i]  += p_accu->EpotSum[i];
    for(int i=0; i < TotNBins; i++) EpotSum2[i] += p_accu->EpotSum2[i];

    for(int i=0; i < TotNBins*NCVs; i++) ICFEpotSum[i]  += p_accu->ICFEpotSum[i];
    for(int i=0; i < TotNBins*NCVs; i++) ICFEpotSum2[i] += p_accu->ICFEpotSum2[i];
}

//------------------------------------------------------------------------------

void CABFAccumulator::SubABFAccumulator(const CABFAccumulator* p_accu)
{
    if( p_accu == NULL ) {
        INVALID_ARGUMENT("p_accu is NULL");
    }

    if( CheckCVSInfo(p_accu) == false ) {
        RUNTIME_ERROR("accumulators are not compatible");
    }

    for(int i=0; i < TotNBins; i++) NSamples[i] -= p_accu->NSamples[i];

    for(int i=0; i < TotNBins*NCVs; i++) ICFSum[i]  -= p_accu->ICFSum[i];
    for(int i=0; i < TotNBins*NCVs; i++) ICFSum2[i] -= p_accu->ICFSum2[i];

    for(int i=0; i < TotNBins; i++) EpotSum[i]  -= p_accu->EpotSum[i];
    for(int i=0; i < TotNBins; i++) EpotSum2[i] -= p_accu->EpotSum2[i];

    for(int i=0; i < TotNBins*NCVs; i++) ICFEpotSum[i]  -= p_accu->ICFEpotSum[i];
    for(int i=0; i < TotNBins*NCVs; i++) ICFEpotSum2[i] -= p_accu->ICFEpotSum2[i];
}

//------------------------------------------------------------------------------

void CABFAccumulator::SetNCorr(double ncorr)
{
    NCorr = ncorr;
}

//------------------------------------------------------------------------------

double CABFAccumulator::GetValue(int icv,int ibin,EABFAccuValue realm) const
{
    double nsamples = GetNumOfSamples(ibin);
    if( nsamples <= 0 ) return(0.0);

    double value = 0.0;

    // abf accumulated data
    double icfsum = GetICFSum(icv,ibin);
    double icfsum_square = GetICFSum2(icv,ibin);

    double epotsum = GetEpotSum(ibin);
//    double epotsum_square = GetEpotSum2(ibin);

    double icfepotsum = GetICFEpotSum(icv,ibin);
    double icfepotsum_square = GetICFEpotSum2(icv,ibin);

    switch(realm){
// mean force
        case(EABF_DG_VALUE): {
                value = icfsum / nsamples;  // mean ABF force
                return(value);
            }
        case(EABF_DG_SIGMA): {
                // sq is variance of ABF force
                double sq = nsamples*icfsum_square - icfsum*icfsum;
                if(sq > 0) {
                    sq = sqrt(sq) / nsamples;
                } else {
                    sq = 0.0;
                }
                return(sq);
            }
        case(EABF_DG_ERROR): {
                // sq is variance of ABF force
                double sq = nsamples*icfsum_square - icfsum*icfsum;
                if(sq > 0) {
                    sq = sqrt(sq) / nsamples;
                } else {
                    sq = 0.0;
                }
                // value is standard error of mean ABF force
                value = sq / sqrt((double)nsamples/(double)NCorr);
                return(value);
            }
// -TdS
        case(EABF_TDS_VALUE): {
                value = icfepotsum/nsamples -(icfsum/nsamples)*(epotsum/nsamples);
                //cout << icfepotsum/nsamples << " " << (icfsum/nsamples)*(epotsum/nsamples) << endl;
                value = value / (Temperature * 1.98720425864083e-3); // R in kcal/mol/K
                return(value);
            }
        case(EABF_TDS_SIGMA): {
                // FIXME
                return(0.0);
            }
        case(EABF_TDS_ERROR): {
                // FIXME - only part
                // sq is variance of ABF force
                double sq = nsamples*icfepotsum_square - icfepotsum*icfepotsum;
                if(sq > 0) {
                    sq = sqrt(sq) / nsamples;
                } else {
                    sq = 0.0;
                }
                sq = sq / (Temperature * 1.98720425864083e-3);
                // value is standard error of ICF*Epot
                value = sq / sqrt((double)nsamples/(double)NCorr);
                return(5.0);
            }
        default:
            RUNTIME_ERROR("unsupported realm");
    }

    return(value);
}

//------------------------------------------------------------------------------

double CABFAccumulator::GetValue(int ibin,EABFAccuValue realm) const
{
    double nsamples = GetNumOfSamples(ibin);
    if( nsamples <= 0 ) return(0.0);

    double value = 0.0;
    double sum = 0.0;
    double sum_square = 0.0;

    // abf accumulated data
    sum = GetEpotSum(ibin);
    sum_square = GetEpotSum2(ibin);

    switch(realm){
        case(EABF_H_VALUE): {
                value = sum / nsamples;  // mean PotEne
                return(value);
            }
        case(EABF_H_SIGMA): {
                // sq is variance of PotEne
                double sq = nsamples*sum_square - sum*sum;
                if(sq > 0) {
                    sq = sqrt(sq) / nsamples;
                } else {
                    sq = 0.0;
                }
                return(sq);
            }
        case(EABF_H_ERROR): {
                // sq is variance of PotEne
                double sq = nsamples*sum_square - sum*sum;
                if(sq > 0) {
                    sq = sqrt(sq) / nsamples;
                } else {
                    sq = 0.0;
                }
                // value is standard error of mean PotEne
                value = sq / sqrt((double)nsamples/(double)NCorr);
                return(value);
            }
        default:
            RUNTIME_ERROR("unsupported realm");
    }

    return(value);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================




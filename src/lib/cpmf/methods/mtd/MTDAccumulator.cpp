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
#include <MTDAccumulator.hpp>
#include <ErrorSystem.hpp>
#include <XMLElement.hpp>
#include <XMLBinData.hpp>

//------------------------------------------------------------------------------

using namespace std;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CMTDAccumulator::CMTDAccumulator(void)
{
    NCVs            = 0;
    TotNBins        = 0;
    Temperature     = 300.0;
    EnergyFConv     = 1.0;
    EnergyUnit      = "kcal mol^-1";
}

//------------------------------------------------------------------------------

CMTDAccumulator::~CMTDAccumulator(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CMTDAccumulator::Load(const CSmallString& name)
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

void CMTDAccumulator::Load(FILE* fin)
{
    if( fin == NULL ) {
        INVALID_ARGUMENT("stream is not open");
    }

// read MTD accumulator header ------------------
    char abf_id[4];
    char ver_id[3];

// first line can contain MTD version
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
    if(strcmp(abf_id,"MTD") != 0) {
        CSmallString error;
        error << "'MTD' magic word was not found in the first line (line: " << buffer << ")";
        RUNTIME_ERROR(error);
    }
    if(strcmp(ver_id,"V6") == 0) {
        Load_v6(buffer,fin);
        return;
    }

    CSmallString error;
    error << "Unsupported version of MTD accumulator (line: " << buffer << ")";
    RUNTIME_ERROR(error);
}

//------------------------------------------------------------------------------

void CMTDAccumulator::Load_v6(char* fline,FILE* fin)
{
    if( fin == NULL ) {
        INVALID_ARGUMENT("stream is not open");
    }

// read MTD accumulator header ------------------
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
    if(strcmp(abf_id,"MTD") != 0) {
        CSmallString error;
        error << "'MTD' magic word was not found in the first line (line: " << fline << ")";
        RUNTIME_ERROR(error);
    }
    if(strcmp(ver_id,"V6") != 0) {
        CSmallString error;
        error << "only MTD V6 version is supported (line: " << fline << ")";
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
        } else if( key == "ENERGY-UNIT" ) {
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
                    error << "unable to read number of MTD samples for bin " << i;
                    RUNTIME_ERROR(error);
                }
                NSamples[i] = ns;
            }
    // -----------------------------------------------------
        } else if( key == "MTDPOT" ) {
            // accumulated potential
            for(int i=0; i < TotNBins; i++) {
                double pot = 0;
                if( fscanf(fin,"%lf",&pot) != 1 ) {
                    CSmallString error;
                    error << "unable to read MTDPOT for bin " << i;
                    RUNTIME_ERROR(error);
                }
                MTDPot[i] = pot;
            }
    // -----------------------------------------------------
        } else if( key == "MTDFORCE" ) {
            // accumulated force
            for(int i = 0; i < NCVs; i++) {
                for(int j = 0; j < TotNBins; j++) {
                    double cf = 0.0;
                    if(fscanf(fin,"%lf",&cf) != 1) {
                        CSmallString error;
                        error << "unable to read MTDFORCE for bin " << j << " and item " << i;
                        RUNTIME_ERROR(error);
                    }
                    MTDForce[map(i,j)] = cf;
                }
            }
        } else {
            CSmallString error;
            error << "unrecognized MTD accumulator keyword: '" << key << "'";
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

void CMTDAccumulator::Save(const CSmallString& name)
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

void CMTDAccumulator::Save(FILE* fout)
{
//    if(fout == NULL) {
//        INVALID_ARGUMENT("stream is not open");
//    }
//
//    if( TotNBins == 0 ) {
//        CSmallString error;
//        error << "no data in accumulator";
//        RUNTIME_ERROR(error);
//    }
//
//// write MTD accumulator header ------------------
//    if(fprintf(fout," MTD V5 %2d\n",NCVs) <= 0) {
//        CSmallString error;
//        error << "unable to write header";
//        RUNTIME_ERROR(error);
//    }
//
//    if(fprintf(fout,"TEMPERATURE\n%18.11E\n",Temperature) <= 0) {
//        CSmallString error;
//        error << "unable to write temperature";
//        RUNTIME_ERROR(error);
//    }
//
//    if(fprintf(fout,"ENERGY\n%18.11E %36s\n",EnergyFConv,(const char*)EnergyUnit) <= 0) {
//        CSmallString error;
//        error << "unable to write energy unit";
//        RUNTIME_ERROR(error);
//    }
//
//// write coordinate specification ----------------
//    if(fprintf(fout,"CVS\n") <= 0) {
//        CSmallString error;
//        error << "unable to write CVS header";
//        RUNTIME_ERROR(error);
//    }
//    for(int i=0; i < NCVs; i++) {
//        if(fprintf(fout,"%2d %10s %18.11E %18.11E %6d\n",i+1,
//                   (const char*)CVs[i].Type,
//                   CVs[i].MinValue,CVs[i].MaxValue,CVs[i].NBins) <= 0) {
//            CSmallString error;
//            error << "unable to write coordinate definition1 id: " << i+1;
//            RUNTIME_ERROR(error);
//        }
//        if(fprintf(fout,"%2d %55s\n",i+1,(const char*)CVs[i].Name) <= 0) {
//            CSmallString error;
//            error << "unable to write coordinate definition2 id: " << i+1;
//            RUNTIME_ERROR(error);
//        }
//        if(fprintf(fout,"%2d %18.11E %36s\n",i+1,CVs[i].FConv,(const char*)CVs[i].Unit) <= 0) {
//            CSmallString error;
//            error << "unable to write coordinate definition3 id: " << i+1;
//            RUNTIME_ERROR(error);
//        }
//    }
//
//// MTD forces =================================================================
//
//// samples
//    int counter;
//
//    if(fprintf(fout,"NSAMPLES\n") <= 0) {
//        CSmallString error;
//        error << "unable to write NSAMPLES header";
//        RUNTIME_ERROR(error);
//    }
//    counter = 0;
//    for(int i=0; i < TotNBins; i++) {
//        if(fprintf(fout,"%9d ",NSamples[i]) <= 0) {
//            CSmallString error;
//            error << "unable to write number of MTD samples for bin " << i;
//            RUNTIME_ERROR(error);
//        }
//        if(counter % 8 == 7) fprintf(fout,"\n");
//        counter++;
//    }
//    if(counter % 8 != 0) fprintf(fout,"\n");
//
//// MICF
//    if(fprintf(fout,"ICF\n") <= 0) {
//        CSmallString error;
//        error << "unable to write ICF header";
//        RUNTIME_ERROR(error);
//    }
//    counter = 0;
//    for(int i = 0; i < NCVs; i++) {
//        for(int j = 0; j < TotNBins; j++) {
//            if(fprintf(fout,"%19.11E ",MICF[map(i,j)]) <= 0) {
//                CSmallString error;
//                error << "unable to write MICF for bin " << j << " and item " << i;
//                RUNTIME_ERROR(error);
//            }
//            if(counter % 4 == 3) fprintf(fout,"\n");
//            counter++;
//        }
//    }
//    if(counter % 4 != 0) fprintf(fout,"\n");
//
//// M2ICF
//    counter = 0;
//    for(int i = 0; i < NCVs; i++) {
//        for(int j = 0; j < TotNBins; j++) {
//            if(fprintf(fout,"%19.11E ",M2ICF[map(i,j)]) <= 0) {
//                CSmallString error;
//                error << "unable to write M2ICF for bin " << j << " and item " << i;
//                RUNTIME_ERROR(error);
//            }
//            if(counter % 4 == 3) fprintf(fout,"\n");
//            counter++;
//        }
//    }
//    if(counter % 4 != 0) fprintf(fout,"\n");
//
//// MEtot
//    if(fprintf(fout,"ETOT\n") <= 0) {
//        CSmallString error;
//        error << "unable to write ETOT header";
//        RUNTIME_ERROR(error);
//    }
//    counter = 0;
//    for(int i=0; i < TotNBins; i++) {
//        if(fprintf(fout,"%19.11E ",MEtot[i]) <= 0) {
//            CSmallString error;
//            error << "unable to write MEtot for bin " << i;
//            RUNTIME_ERROR(error);
//        }
//        if(counter % 4 == 3) fprintf(fout,"\n");
//        counter++;
//    }
//    if(counter % 4 != 0) fprintf(fout,"\n");
//
//// M2Etot
//    counter = 0;
//    for(int i=0; i < TotNBins; i++) {
//        if(fprintf(fout,"%19.11E ",M2Etot[i]) <= 0) {
//            CSmallString error;
//            error << "unable to write M2Etot for bin " << i;
//            RUNTIME_ERROR(error);
//        }
//        if(counter % 4 == 3) fprintf(fout,"\n");
//        counter++;
//    }
//    if(counter % 4 != 0) fprintf(fout,"\n");
//
//// ICFMEtot
//    if(fprintf(fout,"ICF*ETOT\n") <= 0) {
//        CSmallString error;
//        error << "unable to write ICF*ETOT header";
//        RUNTIME_ERROR(error);
//    }
//    counter = 0;
//    for(int i = 0; i < NCVs; i++) {
//        for(int j = 0; j < TotNBins; j++) {
//            if(fprintf(fout,"%19.11E ",ICFMEtot[map(i,j)]) <= 0) {
//                CSmallString error;
//                error << "unable to write ICFMEtot for bin " << j << " and item " << i;
//                RUNTIME_ERROR(error);
//            }
//            if(counter % 4 == 3) fprintf(fout,"\n");
//            counter++;
//        }
//    }
//    if(counter % 4 != 0) fprintf(fout,"\n");
//
//// ICFM2Etot
//    counter = 0;
//    for(int i = 0; i < NCVs; i++) {
//        for(int j = 0; j < TotNBins; j++) {
//            if(fprintf(fout,"%19.11E ",ICFM2Etot[map(i,j)]) <= 0) {
//                CSmallString error;
//                error << "unable to write ICFM2Etot for bin " << j << " and item " << i;
//                RUNTIME_ERROR(error);
//            }
//            if(counter % 4 == 3) fprintf(fout,"\n");
//            counter++;
//        }
//    }
//    if(counter % 4 != 0) fprintf(fout,"\n");
}

//------------------------------------------------------------------------------

void CMTDAccumulator::GetPoint(unsigned int index,CSimpleVector<double>& point) const
{
    for(int k=NCVs-1; k >= 0; k--) {
        const CColVariable* p_coord = &CVs[k];
        int ibin = index % p_coord->GetNumOfBins();
        point[k] = p_coord->GetValue(ibin);
        index = index / p_coord->GetNumOfBins();
    }
}

void CMTDAccumulator::GetPointRValues(unsigned int index,CSimpleVector<double>& point) const
{
    for(int k=NCVs-1; k >= 0; k--) {
        const CColVariable* p_coord = &CVs[k];
        int ibin = index % p_coord->GetNumOfBins();
        point[k] = p_coord->GetRValue(ibin);
        index = index / p_coord->GetNumOfBins();
    }
}

//------------------------------------------------------------------------------

void CMTDAccumulator::GetIPoint(unsigned int index,CSimpleVector<int>& point) const
{
    for(int k=NCVs-1; k >= 0; k--) {
        const CColVariable* p_coord = &CVs[k];
        int ibin = index % p_coord->GetNumOfBins();
        point[k] = ibin;
        index = index / p_coord->GetNumOfBins();
    }
}

//------------------------------------------------------------------------------

void CMTDAccumulator::SetTemperature(double temp)
{
    Temperature = temp;
}

//------------------------------------------------------------------------------

double CMTDAccumulator::GetTemperature(void)
{
    return(Temperature);
}

//------------------------------------------------------------------------------

void CMTDAccumulator::SetEnergyUnit(double fconv,const CSmallString& unit)
{
    EnergyFConv = fconv;
    EnergyUnit  = unit;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CMTDAccumulator::SetNumOfCVs(int ncvs)
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

void CMTDAccumulator::SetCV(int id,
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

void CMTDAccumulator::SetCV(int id,
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

int CMTDAccumulator::GetNumOfCVs(void) const
{
    return(NCVs);
}

//------------------------------------------------------------------------------

const CColVariable* CMTDAccumulator::GetCV(int cv) const
{
    return(&CVs[cv]);
}

//------------------------------------------------------------------------------

int CMTDAccumulator::GetNumOfBins(int limit) const
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

int CMTDAccumulator::GetNumOfSamples(const CSimpleVector<int>& position) const
{
    int glbindex = 0;
    for(int i=0; i < NCVs; i++) {
        glbindex = glbindex*CVs[i].GetNumOfBins() + position[i];
    }
    return(NSamples[glbindex]);
}

//------------------------------------------------------------------------------

int CMTDAccumulator::GetNumOfSamples(int ibin) const
{
    return(NSamples[ibin]);
}

//------------------------------------------------------------------------------

double CMTDAccumulator::GetMTDPot(int ibin) const
{
    return(MTDPot[ibin]);
}

//------------------------------------------------------------------------------

double CMTDAccumulator::GetMTDForce(int icv,int ibin) const
{
    int glbindex = map(icv,ibin);
    return(MTDForce[glbindex]);
}

//------------------------------------------------------------------------------

void CMTDAccumulator::SetMTDPot(int ibin,const double& value)
{
    MTDPot[ibin] = value;
}

//------------------------------------------------------------------------------

void CMTDAccumulator::SetMTDForce(int icv,int ibin,const double& value)
{
    int glbindex = map(icv,ibin);
    MTDForce[glbindex] = value;
}

//------------------------------------------------------------------------------

void CMTDAccumulator::SetNumOfSamples(int ibin,const int& value)
{
    NSamples[ibin] = value;
}

//------------------------------------------------------------------------------

int* CMTDAccumulator::GetNSamplesArray(void)
{
    return(NSamples);
}

//------------------------------------------------------------------------------

double* CMTDAccumulator::GetMTDPotArray(void)
{
    return(MTDPot);
}

//------------------------------------------------------------------------------

double* CMTDAccumulator::GetMTDForceArray(void)
{
    return(MTDForce);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CMTDAccumulator::GetGlobalIndex(const CSimpleVector<int>& position) const
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

void CMTDAccumulator::Clear(void)
{
    if(NSamples == NULL) return;   // allocation was not finalized

    for(int i=0; i < TotNBins; i++) NSamples[i] = 0;
    for(int i=0; i < TotNBins; i++) MTDPot[i] = 0.0;

    for(int i=0; i < TotNBins*NCVs; i++) MTDForce[i] = 0.0;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CMTDAccumulator::FinalizeAllocation()
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
    MTDPot.CreateVector(TotNBins);

    MTDForce.CreateVector(TotNBins*NCVs);

    // reset data
    Clear();
}

//------------------------------------------------------------------------------

void CMTDAccumulator::Deallocate(void)
{
// do not destroy CVs array !

// destroy only data arrays
    NSamples.FreeVector();

    MTDPot.FreeVector();
    MTDForce.FreeVector();

    TotNBins        = 0;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CMTDAccumulator::LoadCVSInfo(CXMLElement* p_iele)
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

bool CMTDAccumulator::CheckCVSInfo(CXMLElement* p_iele) const
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

bool CMTDAccumulator::CheckCVSInfo(const CMTDAccumulator* p_accu) const
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

void CMTDAccumulator::SaveCVSInfo(CXMLElement* p_tele) const
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

void CMTDAccumulator::PrintCVSInfo(std::ostream& vout)
{
    vout  << endl;
    vout  << "Energy unit: " << EnergyUnit << endl;
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

void CMTDAccumulator::PrintCVSInfo(FILE* p_fout)
{
    fprintf(p_fout,"#\n");
    fprintf(p_fout,"# Energy unit: %s\n",(const char*)EnergyUnit);
    fprintf(p_fout,"#\n");
    fprintf(p_fout,"# == Collective Variables =======================================================\n");
    fprintf(p_fout,"#\n");
    fprintf(p_fout,"# ID P Type       Unit  Name                       Min value   Max value   NBins  \n");
    fprintf(p_fout,"# -- - ---------- ----- -------------------------- ----------- ----------- -------\n");

    for(int i=0; i < NCVs; i++) {
        CVs[i].PrintInfo(p_fout);
    }
}

//------------------------------------------------------------------------------

void CMTDAccumulator::ReadMTDData(CXMLElement* p_rele)
{
    if(p_rele == NULL) {
        INVALID_ARGUMENT("p_iele is NULL");
    }

    unsigned int inc_nsamples_size      = GetNumOfBins()*sizeof(int);
    unsigned int inc_mtdpot_size        = GetNumOfBins()*sizeof(double);
    unsigned int inc_mtdforce_size      = GetNumOfBins()*GetNumOfCVs()*sizeof(double);

    if( (inc_nsamples_size == 0) || (inc_mtdpot_size == 0) || (inc_mtdforce_size == 0) ) {
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
    CXMLBinData* p_inc_mtdpot = p_rele->GetFirstChildBinData("MTDPOT");
    if(p_inc_mtdpot == NULL) {
        RUNTIME_ERROR("unable to open MTDPOT element");
    }

    p_data = p_inc_mtdpot->GetData();
    if((GetMTDPotArray() == NULL) || (p_data == NULL) || (p_inc_mtdpot->GetLength() != inc_mtdpot_size)) {
        CSmallString error;
        error << "inconsistent MTDPOT element dat (r: " << p_inc_mtdpot->GetLength()
              << ", l: " <<  inc_mtdpot_size << ")";
        RUNTIME_ERROR(error);
    }
    memcpy(GetMTDPotArray(),p_data,inc_mtdpot_size);

// --------------------------
    CXMLBinData* p_inc_mtdforce = p_rele->GetFirstChildBinData("MTDFORCE");
    if(p_inc_mtdforce == NULL) {
        RUNTIME_ERROR("unable to open MTDFORCE element");
    }

    p_data = p_inc_mtdforce->GetData();
    if((GetMTDForceArray() == NULL) || (p_data == NULL) || (p_inc_mtdforce->GetLength() != inc_mtdforce_size)) {
        CSmallString error;
        error << "inconsistent MTDFORCE element dat (r: " << p_inc_mtdforce->GetLength()
              << ", l: " <<  inc_mtdforce_size << ")";
        RUNTIME_ERROR(error);
    }
    memcpy(GetMTDForceArray(),p_data,inc_mtdforce_size);
}

//------------------------------------------------------------------------------

void CMTDAccumulator::AddMTDData(CXMLElement* p_cele)
{
    if(p_cele == NULL) {
        INVALID_ARGUMENT("p_iele is NULL");
    }

    unsigned int inc_nsamples           = GetNumOfBins();
    unsigned int inc_mtdpot             = GetNumOfBins();
    unsigned int inc_mtdforce           = GetNumOfBins()*GetNumOfCVs();

    unsigned int inc_nsamples_size      = GetNumOfBins()*sizeof(int);
    unsigned int inc_mtdpot_size        = GetNumOfBins()*sizeof(double);
    unsigned int inc_mtdforce_size      = GetNumOfBins()*GetNumOfCVs()*sizeof(double);


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
    CXMLBinData* p_inc_mtdpot = p_cele->GetFirstChildBinData("MTDPOT");
    if(p_inc_mtdpot == NULL) {
        RUNTIME_ERROR("unable to open MTDPOT element");
    }

    p_data = p_inc_mtdpot->GetData();
    if((GetMTDPotArray() == NULL) || (p_data == NULL) || (p_inc_mtdpot->GetLength() != inc_mtdpot_size)) {
        CSmallString error;
        error << "inconsistent MTDPOT element dat (r: " << p_inc_mtdpot->GetLength()
              << ", l: " <<  inc_mtdpot_size << ")";
        RUNTIME_ERROR(error);
    }

    double* jdst = GetMTDPotArray();
    double* jsrc = (double*)p_data;
    for(unsigned int i=0; i < inc_mtdpot; i++) {
        *jdst++ += *jsrc++;
    }

// --------------------------
    CXMLBinData* p_inc_mtdforce = p_cele->GetFirstChildBinData("MTDFORCE");
    if(p_inc_mtdforce == NULL) {
        RUNTIME_ERROR("unable to open MTDFORCE element");
    }

    p_data = p_inc_mtdforce->GetData();
    if((GetMTDForceArray() == NULL) || (p_data == NULL) || (p_inc_mtdforce->GetLength() != inc_mtdforce_size)) {
        CSmallString error;
        error << "inconsistent MTDFORCE element dat (r: " << p_inc_mtdforce->GetLength()
              << ", l: " <<  inc_mtdforce_size << ")";
        RUNTIME_ERROR(error);
    }

    double* kdst = GetMTDForceArray();
    double* ksrc = (double*)p_data;
    for(unsigned int i=0; i < inc_mtdforce; i++) {
        *kdst++ += *ksrc++;
    }
}

//------------------------------------------------------------------------------

void CMTDAccumulator::WriteMTDData(CXMLElement* p_rele)
{
    if(p_rele == NULL) {
        INVALID_ARGUMENT("p_rele is NULL");
    }

    unsigned int inc_nsamples_size      = GetNumOfBins()*sizeof(int);
    unsigned int inc_mtdpot_size        = GetNumOfBins()*sizeof(double);
    unsigned int inc_mtdforce_size      = GetNumOfBins()*GetNumOfCVs()*sizeof(double);

// write collected data to response
    CXMLBinData* p_isamples = p_rele->CreateChildBinData("NSAMPLES");
    p_isamples->CopyData(GetNSamplesArray(),inc_nsamples_size);

    CXMLBinData* p_inc_icfsum = p_rele->CreateChildBinData("MTDPOT");
    p_inc_icfsum->CopyData(GetMTDPotArray(),inc_mtdpot_size);

    CXMLBinData* p_inc_icfsum2 = p_rele->CreateChildBinData("MTDFORCE");
    p_inc_icfsum2->CopyData(GetMTDForceArray(),inc_mtdforce_size);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CMTDAccumulator::map(int item,int bin) const
{
    return(item + bin*NCVs);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CMTDAccumulator::AddMTDAccumulator(const CMTDAccumulator* p_accu)
{
    if( p_accu == NULL ) {
        INVALID_ARGUMENT("p_accu is NULL");
    }

    if( CheckCVSInfo(p_accu) == false ) {
        RUNTIME_ERROR("accumulators are not compatible");
    }

    for(int i=0; i < TotNBins; i++) NSamples[i] += p_accu->NSamples[i];
    for(int i=0; i < TotNBins; i++) MTDPot[i]   += p_accu->MTDPot[i];

    for(int i=0; i < TotNBins*NCVs; i++) MTDForce[i] += p_accu->MTDForce[i];
}

//------------------------------------------------------------------------------

void CMTDAccumulator::SubMTDAccumulator(const CMTDAccumulator* p_accu)
{
    if( p_accu == NULL ) {
        INVALID_ARGUMENT("p_accu is NULL");
    }

    if( CheckCVSInfo(p_accu) == false ) {
        RUNTIME_ERROR("accumulators are not compatible");
    }

    for(int i=0; i < TotNBins; i++) NSamples[i] -= p_accu->NSamples[i];
    for(int i=0; i < TotNBins; i++) MTDPot[i]   -= p_accu->MTDPot[i];

    for(int i=0; i < TotNBins*NCVs; i++) MTDForce[i] -= p_accu->MTDForce[i];
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================




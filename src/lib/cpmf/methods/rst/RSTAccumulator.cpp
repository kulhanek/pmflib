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
#include <RSTAccumulator.hpp>
#include <ErrorSystem.hpp>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CRSTAccumulator::CRSTAccumulator(void)
{
    NCoords     = 0;
    TotNBins  = 0;
}

//------------------------------------------------------------------------------

CRSTAccumulator::~CRSTAccumulator(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CRSTAccumulator::Load(const CSmallString& name)
{
    FILE* fin = fopen(name, "r");

    if(fin == NULL) {
        CSmallString error;
        error << "unable to open file '" << name << "' (" << strerror(errno) << ")";
        RUNTIME_ERROR(error);
    }

    Load(fin);
    fclose(fin);
}

//------------------------------------------------------------------------------

void CRSTAccumulator::Load(FILE* fin)
{
    if(fin == NULL) {
        INVALID_ARGUMENT("stream is not open");
    }

// read RST accumulator header ------------------
    char rst_id[4];
    char rst_ver[3];
    int  numofcoords = 0;

// first line can contain either two or three records
    char buffer[80];
    if( fgets(buffer,80,fin) == NULL ) {
        CSmallString error;
        error << "unable to read the first line";
        RUNTIME_ERROR(error);
    }

    int nr = sscanf(buffer,"%3s %2s %d",rst_id,rst_ver,&numofcoords);
    if( nr != 3) {
        CSmallString error;
        error << "incorrect header record (nr: " << nr << "), line '" << buffer << "'";
        RUNTIME_ERROR(error);
    }

// check ID string
    rst_id[3]='\0';
    if(strcmp(rst_id,"RST") != 0) {
        CSmallString error;
        error << "'RST' magic word was not found in the first line";
        RUNTIME_ERROR(error);
    }

    rst_ver[2]='\0';
    if(strcmp(rst_ver,"V1") != 0) {
        CSmallString error;
        error << "only RST V1 version is supported (line: " << buffer << ")";
        RUNTIME_ERROR(error);
    }

    if(numofcoords <= 0) {
        CSmallString error;
        error << "number of coordinates has to be greater than zero, but " << numofcoords << " was found";
        RUNTIME_ERROR(error);
    }

    SetNumberOfCoords(numofcoords);

// read coordinate specification ----------------
    for(int i=0; i < NCoords; i++) {
        int             id = 0;
        char            type[11];
        char            name[51];
        double          min_value = 0.0;
        double          max_value = 0.0;
        int             nbins = 0;
        int             tr = 0;

        memset(type,0,11);
        memset(name,0,51);

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
        SetCoordinate(i,name,type,min_value,max_value,nbins);
    }

// alloc accumulator data -----------------------
    FinalizeAllocation();

// read accumulator data ------------------------

// samples
    for(int i=0; i < TotNBins; i++) {
        int ns = 0;
        if(fscanf(fin,"%d",&ns) != 1) {
            CSmallString error;
            error << "unable to read number of samples for bin " << i;
            RUNTIME_ERROR(error);
        }
        NSamples[i] = ns;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CRSTAccumulator::Save(const CSmallString& name)
{
    FILE* fout = fopen(name, "w");

    if(fout == NULL) {
        CSmallString error;
        error << "unable to open file '" << name << "' (" << strerror(errno) << ")";
        RUNTIME_ERROR(error);
    }

    Save(fout);
    fclose(fout);
}

//------------------------------------------------------------------------------

void CRSTAccumulator::Save(FILE* fout)
{
    if(fout == NULL) {
        INVALID_ARGUMENT("stream is not open")
    }

    if(TotNBins == 0) {
        RUNTIME_ERROR("no data in accumulator");
    }

// 10  format(A4,1X,I3)
// 20  format(I2,1X,A5,1X,E18.11,1X,E18.11,1X,I6)
// 30  format(8(I9,1X))
// 40  format(4(E19.11,1X))

// write RST accumulator header ------------------
    if(fprintf(fout," RST %3d\n",NCoords) <= 0) {
        CSmallString error;
        error << "unable to write header";
        RUNTIME_ERROR(error);
    }

// write coordinate specification ----------------
    for(int i=0; i < NCoords; i++) {
        if(fprintf(fout,"%2d %5s %18.11e %18.11E %6d\n",i+1,(const char*)Sizes[i].Type,
                   Sizes[i].MinValue,Sizes[i].MaxValue,Sizes[i].NBins) <= 0) {
            CSmallString error;
            error << "unable to write coordinate definition id: " << i+1;
            RUNTIME_ERROR(error);
        }
    }

// write accumulator data ------------------------

// samples
    int counter;

    counter = 0;
    for(int i=0; i < TotNBins; i++) {
        if(fprintf(fout,"%9d ",NSamples[i]) <= 0) {
            CSmallString error;
            error << "unable to write number of samples for bin " << i;
            RUNTIME_ERROR(error);
        }
        if(counter % 8 == 7) fprintf(fout,"\n");
        counter++;
    }
    if(counter % 8 != 0) fprintf(fout,"\n");
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CRSTAccumulator::SetNumberOfCoords(int numofcoords)
{
    if(numofcoords <= 0) {
        INVALID_ARGUMENT("numofcoords <= 0");
    }

    if(NCoords > 0) {
        // destroy all previous data
        Deallocate();
        Sizes.FreeVector();
        NCoords = 0;
    }

// try to allocate Sizes array
    Sizes.CreateVector(numofcoords);

// all seems to be fine - update items
    NCoords = numofcoords;
}

//------------------------------------------------------------------------------

void CRSTAccumulator::SetCoordinate(int id,
                                    const CSmallString& name,
                                    const CSmallString& type,
                                    double min_value,double max_value,int nbins)
{
    if( Sizes == NULL ){
        RUNTIME_ERROR("Sizes is NULL");
    }
    if( id < 0 || id >= NCoords ){
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

    Sizes[id].ID = id;
    Sizes[id].Name = name;
    Sizes[id].Type = type;
    Sizes[id].MinValue = min_value;
    Sizes[id].MaxValue = max_value;

    Sizes[id].NBins = nbins;
    Sizes[id].BinWidth = (Sizes[id].MaxValue - Sizes[id].MinValue)/Sizes[id].NBins;
    Sizes[id].Width = Sizes[id].MaxValue - Sizes[id].MinValue;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CRSTAccumulator::GetNumberOfCoords(void) const
{
    return(NCoords);
}

//------------------------------------------------------------------------------

int CRSTAccumulator::GetNumberOfBins(void) const
{
    return(TotNBins);
}

//------------------------------------------------------------------------------

int CRSTAccumulator::GetNumberOfBinsWithLimit(int limit) const
{
    int number = 0;
    for(int i=0; i < TotNBins; i++) {
        if(NSamples[i] >= limit) number++;
    }
    return(number);
}

//------------------------------------------------------------------------------

const CColVariable* CRSTAccumulator::GetCoordinate(int cv) const
{
    return(&Sizes[cv]);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CRSTAccumulator::GetNumberOfSamples(const CSimpleVector<int>& position) const
{
    int glbindex = 0;
    for(int i=0; i < NCoords; i++) {
        glbindex = glbindex*Sizes[i].GetNumberOfBins() + position[i];
    }
    return(NSamples[glbindex]);
}

//------------------------------------------------------------------------------

int CRSTAccumulator::GetGlobalIndex(const CSimpleVector<int>& position) const
{
    int glbindex = 0;
    for(int i=0; i < NCoords; i++) {
        glbindex = glbindex*Sizes[i].GetNumberOfBins() + position[i];
    }
    return(glbindex);
}

//------------------------------------------------------------------------------

int CRSTAccumulator::GetNumberOfSamples(int ibin) const
{
    return(NSamples[ibin]);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CRSTAccumulator::Clear(void)
{
    if(NSamples == NULL) return;   // allocation was not finalized

    for(int i=0; i < TotNBins; i++) NSamples[i] = 0;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CRSTAccumulator::FinalizeAllocation(void)
{
    if(NCoords <= 0) return;     // at least one coordinate is required
    if(Sizes == NULL) return;   // Sizes array is not allocated
    if(NSamples != NULL) return;    // already finalized!

// check if there is non zero number of bins per coordinate
    TotNBins = 1;
    for(int i=0; i < NCoords; i++) {
        TotNBins *= Sizes[i].NBins;
    }

// and now allocate data arrays
    NSamples.CreateVector(TotNBins);
}

//------------------------------------------------------------------------------

void CRSTAccumulator::Deallocate(void)
{
// do not destroy Sizes array !

// destroy only data arrays
    NSamples.FreeVector();
    TotNBins = 0;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================




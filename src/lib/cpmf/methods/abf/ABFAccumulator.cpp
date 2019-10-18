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
    NCoords     = 0;
    TotNBins    = 0;
    NCorr       = 1.0;
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
    int  numofcoords = 0;
    int  nr;

// first line can contain either two or four records for version 0
    char buffer[80];
    if( fgets(buffer,80,fin) == NULL ) {
        RUNTIME_ERROR("unable to read the first line");
    }

    nr = sscanf(buffer,"%3s %2s %d",abf_id,ver_id,&numofcoords);
    if( nr != 3 ){
        CSmallString error;
        error << "illegal header - three items expected (line: " << buffer << ")";
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
    if(strcmp(ver_id,"V3") != 0) {
        CSmallString error;
        error << "only ABF V3 version is supported (line: " << buffer << ")";
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
//     read(iounit,40,end=100,err=100) (accumulator%ABFForce(i,j),j=1,accumulator%TotNBins)
// end do

// accumulated force
    for(int i = 0; i < NCoords; i++) {
        for(int j = 0; j < TotNBins; j++) {
            double cf = 0.0;
            if(fscanf(fin,"%lf",&cf) != 1) {
                CSmallString error;
                error << "unable to read accumulated ABF force for bin " << j << " and item " << i;
                RUNTIME_ERROR(error);
            }
            ABFForce[map(i,j)] = cf;
        }
    }

// accumulated force square
    for(int i = 0; i < NCoords; i++) {
        for(int j = 0; j < TotNBins; j++) {
            double cf = 0.0;
            if(fscanf(fin,"%lf",&cf) != 1) {
                CSmallString error;
                error << "unable to read accumulated ABF force squares for bin " << j << " and item " << i;
                RUNTIME_ERROR(error);
            }
            ABFForce2[map(i,j)] = cf;
        }
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

// 10  format(A3,1X,A2,1X,I2)
// 20  format(I2,1X,A10,1X,E18.11,1X,E18.11,1X,I6)
// 25  format(I2,1X,A55)
// 30  format(8(I9,1X))
// 40  format(4(E19.11,1X))

// write ABF accumulator header ------------------
    if(fprintf(fout,"ABF V3 %2d\n",NCoords) <= 0) {
        CSmallString error;
        error << "unable to write header";
        RUNTIME_ERROR(error);
    }

// write coordinate specification ----------------
    for(int i=0; i < NCoords; i++) {
        if(fprintf(fout,"%2d %10s %18.11e %18.11E %6d\n",i+1,
                   (const char*)Sizes[i].Type,
                   Sizes[i].MinValue,Sizes[i].MaxValue,Sizes[i].NBins) <= 0) {
            CSmallString error;
            error << "unable to write coordinate definition id: " << i+1;
            RUNTIME_ERROR(error);
        }
        if(fprintf(fout,"%2d %55s\n",i+1,(const char*)Sizes[i].Name) <= 0) {
            CSmallString error;
            error << "unable to write coordinate definition id: " << i+1;
            RUNTIME_ERROR(error);
        }
    }

// ABF forces =================================================================

// samples
    int counter;

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

//  do i=1,fnitem
//     write(iounit,40) (accumulator%ABFForce(i,j),j=1,accumulator%TotNBins)
//  end do

// accumulated force
    counter = 0;
    for(int i = 0; i < NCoords; i++) {
        for(int j = 0; j < TotNBins; j++) {
            if(fprintf(fout,"%19.11E ",ABFForce[map(i,j)]) <= 0) {
                CSmallString error;
                error << "unable to write accumulated ABF force for bin " << j << " and item " << i;
                RUNTIME_ERROR(error);
            }
            if(counter % 4 == 3) fprintf(fout,"\n");
            counter++;
        }
    }
    if(counter % 4 != 0) fprintf(fout,"\n");

// accumulated force square
    counter = 0;
    for(int i = 0; i < NCoords; i++) {
        for(int j = 0; j < TotNBins; j++) {
            if(fprintf(fout,"%19.11E ",ABFForce2[map(i,j)]) <= 0) {
                CSmallString error;
                error << "unable to write accumulated ABF force squares for bin " << j << " and item " << i;
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
    if(fprintf(fout,"MABF V3 %2d\n",NCoords) <= 0) {
        CSmallString error;
        error << "unable to write header";
        RUNTIME_ERROR(error);
    }

// write coordinate specification ----------------
    for(int i=0; i < NCoords; i++) {
        if(fprintf(fout,"%2d %-10s %18.11e %18.11E %6d\n",i+1,
                   (const char*)Sizes[i].Type,
                   Sizes[i].MinValue,Sizes[i].MaxValue,Sizes[i].NBins) <= 0) {
            CSmallString error;
            error << "unable to write coordinate definition id: " << i+1;
            RUNTIME_ERROR(error);
        }
        if(fprintf(fout,"%2d %-55s\n",i+1,(const char*)Sizes[i].Name) <= 0) {
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
    point.CreateVector(NCoords);

    for(int i=0; i < TotNBins; i++) {
        GetPoint(i,point);
        for(int j=0; j < NCoords; j++){
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
        if( i % Sizes[0].GetNumberOfBins() == Sizes[0].GetNumberOfBins() -1 ){
            fprintf(fout,"\n");
        }
    }
}

//------------------------------------------------------------------------------

void CABFAccumulator::GetPoint(unsigned int index,CSimpleVector<double>& point) const
{
    for(int k=NCoords-1; k >= 0; k--) {
        const CColVariable* p_coord = &Sizes[k];
        int ibin = index % p_coord->GetNumberOfBins();
        point[k] = p_coord->GetValue(ibin);
        index = index / p_coord->GetNumberOfBins();
    }
}

//------------------------------------------------------------------------------

void CABFAccumulator::GetIPoint(unsigned int index,CSimpleVector<int>& point) const
{
    for(int k=NCoords-1; k >= 0; k--) {
        const CColVariable* p_coord = &Sizes[k];
        int ibin = index % p_coord->GetNumberOfBins();
        point[k] = ibin;
        index = index / p_coord->GetNumberOfBins();
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABFAccumulator::SetNumberOfCoords(int numofcoords)
{
    if( numofcoords < 0 ) {
        INVALID_ARGUMENT("numofcoords < 0");
    }

    if(NCoords > 0) {
        // destroy all previous data
        Deallocate();
        Sizes.FreeVector();
        NCoords = 0;
    }

// try to allocate Sizes array
    if( numofcoords > 0 ) {
        Sizes.CreateVector(numofcoords);
    }

// all seems to be fine - update items
    NCoords = numofcoords;
}

//------------------------------------------------------------------------------

void CABFAccumulator::SetCoordinate(int id,
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

int CABFAccumulator::GetNumberOfCoords(void) const
{
    return(NCoords);
}

//------------------------------------------------------------------------------

const CColVariable* CABFAccumulator::GetCoordinate(int cv) const
{
    return(&Sizes[cv]);
}

//------------------------------------------------------------------------------

int CABFAccumulator::GetNumberOfBins(void) const
{
    return(TotNBins);
}

//------------------------------------------------------------------------------

int CABFAccumulator::GetNumberOfBinsWithABFLimit(int limit) const
{
    int number = 0;
    for(int i=0; i < TotNBins; i++) {
        if(NSamples[i] >= limit) number++;
    }
    return(number);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CABFAccumulator::GetNumberOfABFSamples(const CSimpleVector<int>& position) const
{
    int glbindex = 0;
    for(int i=0; i < NCoords; i++) {
        glbindex = glbindex*Sizes[i].GetNumberOfBins() + position[i];
    }
    return(NSamples[glbindex]);
}

//------------------------------------------------------------------------------

int CABFAccumulator::GetNumberOfABFSamples(int ibin) const
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

const double& CABFAccumulator::GetABFForceSum(int icoord,int ibin) const
{
    int glbindex = map(icoord,ibin);
    return(ABFForce[glbindex]);
}

//------------------------------------------------------------------------------

const double& CABFAccumulator::GetABFForceSquareSum(int icoord,int ibin) const
{
    int glbindex = map(icoord,ibin);
    return(ABFForce2[glbindex]);
}

//------------------------------------------------------------------------------

void CABFAccumulator::SetABFForceSum(int icoord,int ibin,const double& value)
{
    int glbindex = map(icoord,ibin);
    ABFForce[glbindex] = value;
}

//------------------------------------------------------------------------------

void CABFAccumulator::SetABFForceSquareSum(int icoord,int ibin,const double& value)
{
    int glbindex = map(icoord,ibin);
    ABFForce2[glbindex] = value;
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

double* CABFAccumulator::GetABFForceSumArray(void)
{
    return(ABFForce);
}

//------------------------------------------------------------------------------

double* CABFAccumulator::GetABFForceSquareSumArray(void)
{
    return(ABFForce2);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CABFAccumulator::GetGlobalIndex(const CSimpleVector<int>& position) const
{
    int glbindex = 0;
    for(int i=0; i < NCoords; i++) {
        if( position[i] < 0 ) return(-1);
        if( position[i] >= (int)Sizes[i].GetNumberOfBins() ) return(-1);
        glbindex = glbindex*Sizes[i].GetNumberOfBins() + position[i];
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

    for(int i=0; i < TotNBins*NCoords; i++) ABFForce[i] = 0.0;
    for(int i=0; i < TotNBins*NCoords; i++) ABFForce2[i] = 0.0;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABFAccumulator::FinalizeAllocation(void)
{
    if(NCoords <= 0) return;        // at least one coordinate is required
    if(Sizes == NULL) return;       // Sizes array is not allocated
    if(NSamples != NULL) return;    // already finalized!

// check if there is non zero number of bins per coordinate
    TotNBins = 1;
    for(int i=0; i < NCoords; i++) {
        TotNBins *= Sizes[i].NBins;
    }

// and now allocate data arrays
    NSamples.CreateVector(TotNBins);
    Mask.CreateVector(TotNBins);
    for(int i=0; i < TotNBins; i++){
        Mask[i] = 1.0;
    }
    ABFForce.CreateVector(TotNBins*NCoords);
    ABFForce2.CreateVector(TotNBins*NCoords);
}

//------------------------------------------------------------------------------

void CABFAccumulator::Deallocate(void)
{
// do not destroy Sizes array !

// destroy only data arrays
    NSamples.FreeVector();
    ABFForce.FreeVector();
    ABFForce2.FreeVector();

    TotNBins = 0;
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
    if( p_ele->GetAttribute("NCoords",lnitems) == false) {
        RUNTIME_ERROR("unable to get header attributes");
    }

    if( lnitems == 0 ) {
        // no data
        return;
    }

    SetNumberOfCoords(lnitems);

    CXMLElement*   p_cel = p_ele->GetFirstChildElement("COORD");

    int            ccount = 0;

    while(p_cel != NULL) {
        if( ccount >= lnitems ) {
            LOGIC_ERROR("more COORD elements than NCoords");
        }
        Sizes[ccount].LoadInfo(p_cel);
        ccount++;
        p_cel = p_cel->GetNextSiblingElement("COORD");
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

    result &= p_ele->GetAttribute("NCoords",lnitems);
    if(result == false) {
        ES_ERROR("unable to get header attributes");
        return(false);
    }

    if(lnitems != NCoords) {
        ES_ERROR("mismatch in the number of coordinates");
        return(false);
    }

    CXMLElement*   p_cel = NULL;
    if(p_ele != NULL) p_cel = p_ele->GetFirstChildElement("COORD");
    int            ccount = 0;

    while(p_cel != NULL) {
        if(ccount >= lnitems) {
            ES_ERROR("more COORD elements than NCoords");
            return(false);
        }
        if( Sizes[ccount].CheckInfo(p_cel) == false ) {
            CSmallString error;
            error << "mismatch in cv: " << ccount+1;
            ES_ERROR(error);
            return(false);
        }
        ccount++;
        p_cel = p_cel->GetNextSiblingElement("COORD");
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CABFAccumulator::CheckCVSInfo(const CABFAccumulator* p_accu) const
{
    if(p_accu == NULL) {
        INVALID_ARGUMENT("p_iele is NULL");
    }

    if(p_accu->NCoords != NCoords) {
        ES_ERROR("mismatch in the number of coordinates");
        return(false);
    }

    for(int i=0; i < NCoords; i++) {
        if(Sizes[i].CheckInfo(&p_accu->Sizes[i]) == false) {
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

    p_ele->SetAttribute("NCoords",NCoords);

    for(int i=0; i < NCoords; i++) {
        CXMLElement* p_iele = p_ele->CreateChildElement("COORD");
        Sizes[i].SaveInfo(p_iele);
    }
}

//------------------------------------------------------------------------------

void CABFAccumulator::PrintCVSInfo(std::ostream& vout)
{
    vout  << endl;
    vout << "=== Collective Variables =======================================================" << endl;
    vout << endl;
    vout << "ID  Type      Name                          Min value       Max value     NBins " << endl;
    vout << "-- ---------- -------------------------- --------------- --------------- -------" << endl;

    for(int i=0; i < NCoords; i++) {
        Sizes[i].PrintInfo(vout);
    }
}

//------------------------------------------------------------------------------

void CABFAccumulator::ReadABFData(CXMLElement* p_rele)
{
    if(p_rele == NULL) {
        INVALID_ARGUMENT("p_iele is NULL");
    }

    unsigned int nisamples_size = GetNumberOfBins()*sizeof(int);
    unsigned int iabfforce_size = GetNumberOfBins()*GetNumberOfCoords()*sizeof(double);

    if((nisamples_size == 0) || (iabfforce_size == 0)) {
        RUNTIME_ERROR("size(s) is(are) zero");
    }

    CXMLBinData* p_nisamples = p_rele->GetFirstChildBinData("NISAMPLES");
    if(p_nisamples == NULL) {
        RUNTIME_ERROR("unable to open NISAMPLES element");
    }

    void* p_data = p_nisamples->GetData();
    if((GetNSamplesArray() == NULL) || (p_data == NULL) || (p_nisamples->GetLength() != nisamples_size)) {
        CSmallString error;
        error << "inconsistent NISAMPLES element dat (r: " << p_nisamples->GetLength()
              << ", l: " <<  nisamples_size << ")";
        RUNTIME_ERROR(error);
    }
    memcpy(GetNSamplesArray(),p_data,nisamples_size);

    CXMLBinData* p_iabfforce = p_rele->GetFirstChildBinData("IABFFORCE");
    if(p_iabfforce == NULL) {
        RUNTIME_ERROR("unable to open IABFFORCE element");
    }

    p_data = p_iabfforce->GetData();
    if((GetABFForceSumArray() == NULL) || (p_data == NULL) || (p_iabfforce->GetLength() != iabfforce_size)) {
        CSmallString error;
        error << "inconsistent IABFFORCE element dat (r: " << p_iabfforce->GetLength()
              << ", l: " <<  iabfforce_size << ")";
        RUNTIME_ERROR(error);
    }
    memcpy(GetABFForceSumArray(),p_data,iabfforce_size);


    CXMLBinData* p_iabfforce2 = p_rele->GetFirstChildBinData("IABFFORCE2");
    if(p_iabfforce2 == NULL) {
        RUNTIME_ERROR("unable to open IABFFORCE2 element");
    }

    p_data = p_iabfforce2->GetData();
    if((GetABFForceSquareSumArray() == NULL) || (p_data == NULL) || (p_iabfforce2->GetLength() != iabfforce_size)) {
        CSmallString error;
        error << "inconsistent IABFFORCE element dat (r: " << p_iabfforce2->GetLength()
              << ", l: " <<  iabfforce_size << ")";
        RUNTIME_ERROR(error);
    }
    memcpy(GetABFForceSquareSumArray(),p_data,iabfforce_size);
}

//------------------------------------------------------------------------------

void CABFAccumulator::AddABFData(CXMLElement* p_cele)
{
    if(p_cele == NULL) {
        INVALID_ARGUMENT("p_iele is NULL");
    }

    unsigned int isamples = GetNumberOfBins();
    unsigned int iabfforce = GetNumberOfBins()*GetNumberOfCoords();

    unsigned int isamples_size = GetNumberOfBins()*sizeof(int);
    unsigned int iabfforce_size = GetNumberOfBins()*GetNumberOfCoords()*sizeof(double);

    CXMLBinData* p_isamples = p_cele->GetFirstChildBinData("NISAMPLES");
    if(p_isamples == NULL) {
        RUNTIME_ERROR("unable to open NISAMPLES element");
    }

    void* p_data = p_isamples->GetData();
    if((GetNSamplesArray() == NULL) || (p_data == NULL) || (p_isamples->GetLength() != isamples_size)) {
        CSmallString error;
        error << "inconsistent NISAMPLES element dat (r: " << p_isamples->GetLength()
              << ", l: " <<  isamples_size << ")";
        RUNTIME_ERROR(error);
    }

    int* idst = GetNSamplesArray();
    int* isrc = (int*)p_data;
    for(unsigned int i=0; i < isamples; i++) {
        *idst++ += *isrc++;
    }

    CXMLBinData* p_iabfforce = p_cele->GetFirstChildBinData("IABFFORCE");
    if(p_iabfforce == NULL) {
        RUNTIME_ERROR("unable to open IABFFORCE element");
    }

    p_data = p_iabfforce->GetData();
    if((GetABFForceSumArray() == NULL) || (p_data == NULL) || (p_iabfforce->GetLength() != iabfforce_size)) {
        CSmallString error;
        error << "inconsistent IABFFORCE element dat (r: " << p_iabfforce->GetLength()
              << ", l: " <<  iabfforce_size << ")";
        RUNTIME_ERROR(error);
    }

    double* jdst = GetABFForceSumArray();
    double* jsrc = (double*)p_data;
    for(unsigned int i=0; i < iabfforce; i++) {
        *jdst++ += *jsrc++;
    }

    CXMLBinData* p_iabfforce2 = p_cele->GetFirstChildBinData("IABFFORCE2");
    if(p_iabfforce2 == NULL) {
        RUNTIME_ERROR("unable to open IABFFORCE2 element");
    }

    p_data = p_iabfforce2->GetData();
    if((GetABFForceSquareSumArray() == NULL) || (p_data == NULL) || (p_iabfforce2->GetLength() != iabfforce_size)) {
        CSmallString error;
        error << "inconsistent IABFFORCE element dat (r: " << p_iabfforce2->GetLength()
              << ", l: " <<  iabfforce_size << ")";
        RUNTIME_ERROR(error);
    }

    double* kdst = GetABFForceSquareSumArray();
    double* ksrc = (double*)p_data;
    for(unsigned int i=0; i < iabfforce; i++) {
        *kdst++ += *ksrc++;
    }
}

//------------------------------------------------------------------------------

void CABFAccumulator::WriteABFData(CXMLElement* p_rele)
{
    if(p_rele == NULL) {
        INVALID_ARGUMENT("p_rele is NULL");
    }

    int isamples_size = GetNumberOfBins()*sizeof(int);
    int iabfforce_size = GetNumberOfBins()*GetNumberOfCoords()*sizeof(double);

// write collected data to response
    CXMLBinData* p_isamples = p_rele->CreateChildBinData("NISAMPLES");
    p_isamples->CopyData(GetNSamplesArray(),isamples_size);

    CXMLBinData* p_iabfforce = p_rele->CreateChildBinData("IABFFORCE");
    p_iabfforce->CopyData(GetABFForceSumArray(),iabfforce_size);

    CXMLBinData* p_iabfforce2 = p_rele->CreateChildBinData("IABFFORCE2");
    p_iabfforce2->CopyData(GetABFForceSquareSumArray(),iabfforce_size);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CABFAccumulator::map(int item,int bin) const
{
    return(item + bin*NCoords);
}

//------------------------------------------------------------------------------

int CABFAccumulator::map_g(int group,int item,int bin) const
{
    return(group*NCoords*TotNBins + item + bin*NCoords);
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

    for(int i=0; i < TotNBins*NCoords; i++) ABFForce[i] += p_accu->ABFForce[i];
    for(int i=0; i < TotNBins*NCoords; i++) ABFForce2[i] += p_accu->ABFForce2[i];
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

    for(int i=0; i < TotNBins*NCoords; i++) ABFForce[i] -= p_accu->ABFForce[i];
    for(int i=0; i < TotNBins*NCoords; i++) ABFForce2[i] -= p_accu->ABFForce2[i];
}

//------------------------------------------------------------------------------

void CABFAccumulator::SetNCorr(double ncorr)
{
    NCorr = ncorr;
}

//------------------------------------------------------------------------------

double CABFAccumulator::GetValue(int icoord,int ibin,EABFAccuValue realm) const
{
    int nsamples = GetNumberOfABFSamples(ibin);
    if( nsamples <= 0 ) return(0.0);

    double value = 0.0;
    double sum = 0.0;
    double sum_square = 0.0;

    // abf accumulated data
    sum = GetABFForceSum(icoord,ibin);
    sum_square = GetABFForceSquareSum(icoord,ibin);

    switch(realm){
        case(EABF_MEAN_FORCE_VALUE): {
                value = sum / nsamples;  // mean ABF force
                return(value);
            }
        case(EABF_INST_FORCE_SIGMA): {
                // sq is variance of ABF force
                double sq = nsamples*sum_square - sum*sum;
                if(sq > 0) {
                    sq = sqrt(sq) / nsamples;
                } else {
                    sq = 0.0;
                }
                return(sq);
            }
        case(EABF_MEAN_FORCE_ERROR): {
                // sq is variance of ABF force
                double sq = nsamples*sum_square - sum*sum;
                if(sq > 0) {
                    sq = sqrt(sq) / nsamples;
                } else {
                    sq = 0.0;
                }
                // value is standard error of mean ABF force
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




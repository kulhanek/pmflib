// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2020 Petr Kulhanek, kulhanek@chemi.muni.cz
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
#include <ABPAccumulator.hpp>
#include <ErrorSystem.hpp>
#include <XMLElement.hpp>
#include <XMLBinData.hpp>

//------------------------------------------------------------------------------

using namespace std;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CABPAccumulator::CABPAccumulator(void)
{
    NCoords     = 0;
    TotNBins    = 0;
    Temperature = 0.0;
}

//------------------------------------------------------------------------------

CABPAccumulator::~CABPAccumulator(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABPAccumulator::Load(const CSmallString& name)
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

void CABPAccumulator::Load(FILE* fin)
{
    if( fin == NULL ) {
        INVALID_ARGUMENT("stream is not open");
    }

// read ABP accumulator header ------------------
    char abf_id[4];
    char ver_id[3];
    char kernel[4];
    int  numofcoords = 0;
    int  nr;

// first line can contain either two or four records for version 0
    char buffer[80];
    if( fgets(buffer,80,fin) == NULL ) {
        RUNTIME_ERROR("unable to read the first line");
    }

    nr = sscanf(buffer,"%3s %2s %d %4s %lf",abf_id,ver_id,&numofcoords,kernel,&Temperature);
    if( nr != 5 ){
        CSmallString error;
        error << "illegal header - five items expected (line: " << buffer << ")";
        RUNTIME_ERROR(error);
    }

    abf_id[3]='\0';
    ver_id[2]='\0';
    kernel[3]='\0';

// check ID string
    if(strcmp(abf_id,"ABP") != 0) {
        CSmallString error;
        error << "'ABP' magic word was not found in the first line (line: " << buffer << ")";
        RUNTIME_ERROR(error);
    }
    if(strcmp(ver_id,"V1") != 0) {
        CSmallString error;
        error << "only ABP V1 version is supported (line: " << buffer << ")";
        RUNTIME_ERROR(error);
    }

    if(numofcoords <= 0) {
        CSmallString error;
        error << "number of coordinates has to be greater than zero, but " << numofcoords << " was found";
        RUNTIME_ERROR(error);
    }

    Kernel = kernel;

    SetNumberOfCoords(numofcoords);

// read coordinate specification ----------------
    for(int i=0; i < NCoords; i++) {
        int             id = 0;
        char            type[11];
        char            name[51];
        char            unit[37];
        double          min_value = 0.0;
        double          max_value = 0.0;
        double          fconv = 1.0;
        double          alpha = 0.0;
        int             nbins = 0;
        int             tr = 0;

        memset(type,0,11);
        memset(name,0,51);
        memset(unit,0,37);

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
        tr = fscanf(fin,"%d %lf %lf %37s",&id,&fconv,&alpha,unit);
        if( tr != 4 ) {
            CSmallString error;
            error << "unable to read coordinate definition, id: " << i+1 << " (" << tr << " != 4)";
            RUNTIME_ERROR(error);
        }
        // some tests
        if(id != i+1) {
            CSmallString error;
            error << "coordinate id does not match, read: " << id << ", expected: " << i+1;
            RUNTIME_ERROR(error);
        }

        // init coordinate
        SetCoordinate(i,name,type,min_value,max_value,nbins,alpha,fconv,unit);
    }

// alloc accumulator data -----------------------
    FinalizeAllocation();

// read accumulator data ------------------------

// ABP forces =================================================================

// samples
    for(int i=0; i < TotNBins; i++) {
        int ns = 0;
        if( fscanf(fin,"%d",&ns) != 1 ) {
            CSmallString error;
            error << "unable to read number of ABP samples for bin " << i;
            RUNTIME_ERROR(error);
        }
        NSamples[i] = ns;
    }

// dpop
    for(int i = 0; i < NCoords; i++) {
        for(int j = 0; j < TotNBins; j++) {
            double cf = 0.0;
            if(fscanf(fin,"%lf",&cf) != 1) {
                CSmallString error;
                error << "unable to read accumulated ABP dpop for bin " << j << " and item " << i;
                RUNTIME_ERROR(error);
            }
            DPop[map(i,j)] = cf;
        }
    }

// pop
    for(int i = 0; i < TotNBins; i++) {
        double cf = 0.0;
        if(fscanf(fin,"%lf",&cf) != 1) {
            CSmallString error;
            error << "unable to read accumulated ABP pop for bin " << i;
            RUNTIME_ERROR(error);
        }
        Pop[i] = cf;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABPAccumulator::Save(const CSmallString& name)
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

void CABPAccumulator::Save(FILE* fout)
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
// 26  format(I2,1X,E18.11,1X,A36,1X,E20.11)
// 30  format(8(I9,1X))
// 40  format(4(E19.11,1X))

// write ABP accumulator header ------------------
    if(fprintf(fout," ABP V1 %2d %4s %10.3f\n",NCoords,(const char*)Kernel,Temperature) <= 0) {
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
        if(fprintf(fout,"%2d %18.11E %20.11E %36s\n",i+1,Sizes[i].FConv,Sizes[i].Alpha,(const char*)Sizes[i].Unit) <= 0) {
            CSmallString error;
            error << "unable to write coordinate definition3 id: " << i+1;
            RUNTIME_ERROR(error);
        }
    }

// ABP data =================================================================

// samples
    int counter;

    counter = 0;
    for(int i=0; i < TotNBins; i++) {
        if(fprintf(fout,"%9d ",NSamples[i]) <= 0) {
            CSmallString error;
            error << "unable to write number of ABP samples for bin " << i;
            RUNTIME_ERROR(error);
        }
        if(counter % 8 == 7) fprintf(fout,"\n");
        counter++;
    }
    if(counter % 8 != 0) fprintf(fout,"\n");

// dpop
    counter = 0;
    for(int i = 0; i < NCoords; i++) {
        for(int j = 0; j < TotNBins; j++) {
            if(fprintf(fout,"%19.11E ",DPop[map(i,j)]) <= 0) {
                CSmallString error;
                error << "unable to write accumulated ABP dpop for bin " << j << " and item " << i;
                RUNTIME_ERROR(error);
            }
            if(counter % 4 == 3) fprintf(fout,"\n");
            counter++;
        }
    }
    if(counter % 4 != 0) fprintf(fout,"\n");

// pop
    counter = 0;
    for(int i = 0; i < TotNBins; i++) {
        if(fprintf(fout,"%19.11E ",Pop[i]) <= 0) {
            CSmallString error;
            error << "unable to write accumulated ABP pop for bin " << i;
            RUNTIME_ERROR(error);
        }
        if(counter % 4 == 3) fprintf(fout,"\n");
        counter++;
    }
    if(counter % 4 != 0) fprintf(fout,"\n");
}

//------------------------------------------------------------------------------

void CABPAccumulator::GetPoint(unsigned int index,CSimpleVector<double>& point) const
{
    for(int k=NCoords-1; k >= 0; k--) {
        const CColVariable* p_coord = &Sizes[k];
        int ibin = index % p_coord->GetNumberOfBins();
        point[k] = p_coord->GetValue(ibin);
        index = index / p_coord->GetNumberOfBins();
    }
}

//------------------------------------------------------------------------------

void CABPAccumulator::GetIPoint(unsigned int index,CSimpleVector<int>& point) const
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

void CABPAccumulator::SetNumberOfCoords(int numofcoords)
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

void CABPAccumulator::SetCoordinate(int id,
                                    const CSmallString& name,
                                    const CSmallString& type,
                                    double min_value,double max_value,int nbins,
                                    double alpha,
                                    double fconv, const CSmallString& unit)
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
    Sizes[id].Unit = unit;
    Sizes[id].FConv = fconv;
    Sizes[id].MinValue = min_value;
    Sizes[id].MaxValue = max_value;

    Sizes[id].NBins = nbins;
    Sizes[id].BinWidth = (Sizes[id].MaxValue - Sizes[id].MinValue)/Sizes[id].NBins;
    Sizes[id].Width = Sizes[id].MaxValue - Sizes[id].MinValue;

    Sizes[id].Alpha = alpha;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CABPAccumulator::GetNumberOfCoords(void) const
{
    return(NCoords);
}

//------------------------------------------------------------------------------

const CColVariable* CABPAccumulator::GetCoordinate(int cv) const
{
    return(&Sizes[cv]);
}

//------------------------------------------------------------------------------

int CABPAccumulator::GetNumberOfBins(void) const
{
    return(TotNBins);
}

//------------------------------------------------------------------------------

int CABPAccumulator::GetNumberOfBinsWithABPLimit(int limit) const
{
    int number = 0;
    for(int i=0; i < TotNBins; i++) {
        if(NSamples[i] >= limit) number++;
    }
    return(number);
}

//------------------------------------------------------------------------------

int* CABPAccumulator::GetNSamplesArray(void)
{
   return(NSamples);
}

//-----------------------------------------------------------------------------

double* CABPAccumulator::GetDPopArray(void)
{
   return(DPop);
}

//-----------------------------------------------------------------------------

double* CABPAccumulator::GetPopArray(void)
{
   return(Pop);
}

//------------------------------------------------------------------------------

double CABPAccumulator::GetPop(int ibin) const
{
    return(Pop[ibin]);
}

//------------------------------------------------------------------------------

double CABPAccumulator::GetTemperature(void) const
{
    return(Temperature);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CABPAccumulator::GetNumberOfABPSamples(int ibin) const
{
    return(NSamples[ibin]);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABPAccumulator::Clear(void)
{
    if(NSamples == NULL) return;   // allocation was not finalized

    for(int i=0; i < TotNBins; i++) NSamples[i] = 0;

    for(int i=0; i < TotNBins*NCoords; i++) DPop[i] = 0.0;
    for(int i=0; i < TotNBins; i++) Pop[i] = 1.0;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABPAccumulator::FinalizeAllocation(void)
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
    DPop.CreateVector(TotNBins*NCoords);
    Pop.CreateVector(TotNBins);
}

//------------------------------------------------------------------------------

void CABPAccumulator::Deallocate(void)
{
// do not destroy Sizes array !

// destroy only data arrays
    NSamples.FreeVector();
    DPop.FreeVector();
    Pop.FreeVector();

    TotNBins = 0;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CABPAccumulator::LoadCVSInfo(CXMLElement* p_iele)
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

bool CABPAccumulator::CheckCVSInfo(CXMLElement* p_iele) const
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

bool CABPAccumulator::CheckCVSInfo(const CABPAccumulator* p_accu) const
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

void CABPAccumulator::SaveCVSInfo(CXMLElement* p_tele) const
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

void CABPAccumulator::PrintCVSInfo(std::ostream& vout)
{
    vout  << endl;
    vout << "=== Collective Variables =======================================================" << endl;
    vout << endl;
    vout << "ID P Type       Unit  Name                       Min value   Max value   NBins  " << endl;
    vout << "-- - ---------- ----- -------------------------- ----------- ----------- -------" << endl;

    for(int i=0; i < NCoords; i++) {
        Sizes[i].PrintInfo(vout);
    }
}

//------------------------------------------------------------------------------

void CABPAccumulator::ReadABPData(CXMLElement* p_rele)
{
    if(p_rele == NULL) {
        INVALID_ARGUMENT("p_iele is NULL");
    }

    unsigned int nisamples_size = GetNumberOfBins()*sizeof(int);
    unsigned int idpop_size = GetNumberOfBins()*GetNumberOfCoords()*sizeof(double);
    unsigned int ipop_size = GetNumberOfBins()*sizeof(double);

    if( (nisamples_size == 0) || (idpop_size == 0) || (ipop_size == 0) )  {
        RUNTIME_ERROR("size(s) is(are) zero");
    }

// --------------------------
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

// --------------------------
    CXMLBinData* p_idpop = p_rele->GetFirstChildBinData("IDPOP");
    if(p_idpop == NULL) {
        RUNTIME_ERROR("unable to open IDPOP element");
    }

    p_data = p_idpop->GetData();
    if((GetDPopArray() == NULL) || (p_data == NULL) || (p_idpop->GetLength() != idpop_size)) {
        CSmallString error;
        error << "inconsistent IDPOP element dat (r: " << p_idpop->GetLength()
              << ", l: " <<  idpop_size << ")";
        RUNTIME_ERROR(error);
    }
    memcpy(GetDPopArray(),p_data,idpop_size);

// --------------------------
    CXMLBinData* p_ipop = p_rele->GetFirstChildBinData("IPOP");
    if(p_ipop == NULL) {
        RUNTIME_ERROR("unable to open IPOP element");
    }

    p_data = p_ipop->GetData();
    if((GetPopArray() == NULL) || (p_data == NULL) || (p_ipop->GetLength() != ipop_size)) {
        CSmallString error;
        error << "inconsistent IPOP element dat (r: " << p_ipop->GetLength()
              << ", l: " <<  ipop_size << ")";
        RUNTIME_ERROR(error);
    }
    memcpy(GetPopArray(),p_data,ipop_size);
}

//------------------------------------------------------------------------------

void CABPAccumulator::AddABPData(CXMLElement* p_cele)
{
    if(p_cele == NULL) {
        INVALID_ARGUMENT("p_iele is NULL");
    }

    unsigned int isamples = GetNumberOfBins();
    unsigned int idpop = GetNumberOfBins()*GetNumberOfCoords();
    unsigned int ipop = GetNumberOfBins();

    unsigned int isamples_size = GetNumberOfBins()*sizeof(int);
    unsigned int idpop_size = GetNumberOfBins()*GetNumberOfCoords()*sizeof(double);
    unsigned int ipop_size = GetNumberOfBins()*sizeof(double);

// --------------------------
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

// --------------------------
    CXMLBinData* p_idpop = p_cele->GetFirstChildBinData("IDPOP");
    if(p_idpop == NULL) {
        RUNTIME_ERROR("unable to open IDPOP element");
    }

    p_data = p_idpop->GetData();
    if((GetDPopArray() == NULL) || (p_data == NULL) || (p_idpop->GetLength() != idpop_size)) {
        CSmallString error;
        error << "inconsistent IDPOP element dat (r: " << p_idpop->GetLength()
              << ", l: " <<  idpop_size << ")";
        RUNTIME_ERROR(error);
    }

    double* jdst = GetDPopArray();
    double* jsrc = (double*)p_data;
    for(unsigned int i=0; i < idpop; i++) {
        *jdst++ += *jsrc++;
    }

// --------------------------
    CXMLBinData* p_ipop = p_cele->GetFirstChildBinData("IPOP");
    if(p_ipop == NULL) {
        RUNTIME_ERROR("unable to open IPOP element");
    }

    p_data = p_ipop->GetData();
    if((GetPopArray() == NULL) || (p_data == NULL) || (p_ipop->GetLength() != ipop_size)) {
        CSmallString error;
        error << "inconsistent IPOP element dat (r: " << p_ipop->GetLength()
              << ", l: " <<  ipop_size << ")";
        RUNTIME_ERROR(error);
    }

    double* kdst = GetPopArray();
    double* ksrc = (double*)p_data;
    for(unsigned int i=0; i < ipop; i++) {
        *kdst++ += *ksrc++;
    }
}

//------------------------------------------------------------------------------

void CABPAccumulator::WriteABPData(CXMLElement* p_rele)
{
    if(p_rele == NULL) {
        INVALID_ARGUMENT("p_rele is NULL");
    }

    int isamples_size = GetNumberOfBins()*sizeof(int);
    int idpop_size = GetNumberOfBins()*GetNumberOfCoords()*sizeof(double);
    int ipop_size = GetNumberOfBins()*sizeof(double);

// write collected data to response
    CXMLBinData* p_isamples = p_rele->CreateChildBinData("NISAMPLES");
    p_isamples->CopyData(GetNSamplesArray(),isamples_size);

    CXMLBinData* p_idpop = p_rele->CreateChildBinData("IDPOP");
    p_idpop->CopyData(GetDPopArray(),idpop_size);

    CXMLBinData* p_ipop = p_rele->CreateChildBinData("IPOP");
    p_ipop->CopyData(GetPopArray(),ipop_size);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CABPAccumulator::map(int item,int bin) const
{
    return(item + bin*NCoords);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================




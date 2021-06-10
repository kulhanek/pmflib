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
#include <math.h>
#include <MTDHistory.hpp>
#include <ErrorSystem.hpp>
#include <XMLElement.hpp>
#include <XMLBinData.hpp>
#include <XMLPrinter.hpp>


//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

//void CEnergySurface::CalculateFES(CMTDHistory& mtd_hist,unsigned int mtd_time)
//{
//// quick compatibility comparison
//    if( NumOfCVs != mtd_hist.GetNumOfCVs() ) {
//        RUNTIME_ERROR("numofitems mismatch");
//    }

//// allocate point
//    CSimpleVector<double> point;
//    point.CreateVector(NumOfCVs);

//// calculate surface
//    unsigned int loc = 0;
//    CalculateFES_Part(mtd_hist,point,mtd_time,loc,0);
//}

////------------------------------------------------------------------------------

//void CEnergySurface::CalculateFES_Part(CMTDHistory& mtd_hist,
//                                       CSimpleVector<double>& point,
//                                       unsigned int mtd_time,
//                                       unsigned int& loc,
//                                       unsigned int cv)
//{
//    if(cv >= NumOfCVs) {
//        // calculate value
//        double value = - mtd_hist.CalculateValue(point,mtd_time);
//        Energy[loc++] = value;
//        return;
//    }

//    const CColVariable* p_coord = &CVs[cv];

//// cycle through variable
//    for(unsigned int i = 0; i < p_coord->GetNumOfBins(); i++) {
//        point[cv] = p_coord->GetValue(i);
//        CalculateFES_Part(mtd_hist,point,mtd_time,loc,cv+1);
//    }
//}

////==============================================================================
////------------------------------------------------------------------------------
////==============================================================================

//void CEnergySurface::CalculateFES(int NCVs,CSimpleVector<double>& params)
//{
//// quick compatibility comparison
//    if( NumOfCVs != NCVs ) {
//        RUNTIME_ERROR("NumOfCVs != NCVs");
//    }

//// allocate point
//    CSimpleVector<double> point;
//    point.CreateVector(NumOfCVs);

//// calculate fes
//    unsigned int loc = 0;
//    CalculateFES_MTDParam_Part(params,point,loc,0);
//}

////------------------------------------------------------------------------------

//void CEnergySurface::CalculateFES_MTDParam_Part(
//    CSimpleVector<double>& params,
//    CSimpleVector<double>& point,
//    unsigned int& loc,
//    unsigned int cv)
//{
//    if(cv >= NumOfCVs) {
//        // calculate value
//        double value = - CalculateValue(params,point);
//        Energy[loc++] = value;
//        return;
//    }

//    const CColVariable* p_coord = &CVs[cv];

//// cycle through variable
//    for(unsigned int i = 0; i < p_coord->GetNumOfBins(); i++) {
//        point[cv] = p_coord->GetValue(i);
//        CalculateFES_MTDParam_Part(params,point,loc,cv+1);
//    }
//}

////------------------------------------------------------------------------------

//double CEnergySurface::CalculateValue(const CSimpleVector<double>& params,
//                                      const CSimpleVector<double>& point)
//{
//    double     value = 0.0;
//    double     fexparg;
//    int        num_of_hills = params.GetLength()/(1+2*NumOfCVs);
//    unsigned int        loc = 0;

//    for(int i=0; i < num_of_hills; i++) {
//        fexparg = 0.0;
//        double height = params[loc++];
//        for(int k=0; k < NumOfCVs; k++) {
//            double value = params[loc++];
//            double width = params[loc++];
//            double e = point[k] - value;
//            fexparg = fexparg + e*e / (2.0 * width * width);
//        }
//        value = value + height*exp(-fexparg);
//    }

//    return(value);
//}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CMTDHistory::CMTDHistory(void)
{
    NCVs = 0;
    MaxBufferSize = 1000;
    EnergyUnitFac = 1;
    EnergyUnit = "kcal mol^-1";
}

//------------------------------------------------------------------------------

CMTDHistory::~CMTDHistory(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CMTDHistory::Load(const CSmallString& name)
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

void CMTDHistory::Load(FILE* fin)
{
    if( fin == NULL ) {
        INVALID_ARGUMENT("stream is not open");
    }

// read MTD accumulator header ------------------
    char mtd_id[4];
    char ver_id[3];
    int  numofcoords = 0;
    int  nr;

// first line can contain either two or four records for version 0
    char buffer[80];
    if( fgets(buffer,80,fin) == NULL ) {
        RUNTIME_ERROR("unable to read the first line");
    }

    nr = sscanf(buffer,"%3s %2s %d",mtd_id,ver_id,&numofcoords);
    if( nr != 3 ){
        RUNTIME_ERROR("wrong header");
    }

    mtd_id[3]='\0';
    ver_id[2]='\0';

// check ID string
    if(strcmp(mtd_id,"MTD") != 0) {
        CSmallString error;
        error << "'MTD' magic word was not found in the first line (line: " << buffer << ")";
        RUNTIME_ERROR(error);
    }
    if(strcmp(ver_id,"V3") != 0) {
        CSmallString error;
        error << "'illegal MTD history version (V3 expected) (line: " << buffer << ")";
        RUNTIME_ERROR(error);
    }

    if(numofcoords <= 0) {
        CSmallString error;
        error << "number of coordinates has to be greater than zero, but " << numofcoords << " was found";
        RUNTIME_ERROR(error);
    }

    SetNumOfCVs(numofcoords);

   if( fin == NULL ) {
        INVALID_ARGUMENT("stream is not open");
    }

// read coordinate specification ----------------
    for(int i=0; i < NCVs; i++) {
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
            error << "unable to read coordinate definition, id: " << i+1 << " (number of collums read: " << tr << ")";
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

        tr = fscanf(fin,"%d %55s",&id,name);
        if( tr != 2 ) {
            CSmallString error;
            error << "unable to read coordinate definition, id: " << i+1 << " (number of collums read: " << tr << ")";
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

// now read records -----------------------------

    int             curr_position = 0;
    int             position = 1;
    CMTDBuffer*     p_buff = NULL;

    while(feof(fin) == 0) {
        int nr,level,time;

        if((p_buff == NULL) || (curr_position >= MaxBufferSize)) {
            // allocate new buffer
            p_buff = GetNewBuffer(MaxBufferSize);
            if(p_buff == NULL) {
                CSmallString error;
                error << "unable to allocate new buffer (pos: " << position << ")";
                RUNTIME_ERROR(error);
            }
            curr_position = 0;
        }

        // read level and time
        nr = fscanf(fin,"%d %d",&level,&time);
        if(nr <= 0) break;      // end of file
        if(nr != 2) {
            CSmallString error;
            error << "level and time were not successfully read (pos: " << position << ")";
            RUNTIME_ERROR(error);
        }

        // now read height
        double height;
        nr = fscanf(fin,"%lf",&height);
        if(nr != 1) {
            CSmallString error;
            error << "unable to read height (pos: " << position << ")";
            RUNTIME_ERROR(error);
        }

        p_buff->SetHeight(curr_position,height);

        for(int i=0; i < GetNumOfCVs(); i++) {
            double value,width;
            nr = fscanf(fin,"%lf %lf",&value,&width);
            if(nr != 2) {
                CSmallString error;
                error << "unable to read value and width (pos: " << position << ")";
                RUNTIME_ERROR(error);
            }
            // register value and width
            p_buff->SetValue(curr_position,i,value);
            p_buff->SetWidth(curr_position,i,width);
        }
        p_buff->IncNumberOfHills();
        curr_position++;
        position++;
    }
}

//------------------------------------------------------------------------------

CMTDBuffer* CMTDHistory::GetNewBuffer(int size)
{
    CMTDBuffer* p_buff = new CMTDBuffer();
    p_buff->AllocateData(size,NCVs);
    Buffers.InsertToEnd(p_buff,0,true);
    return(p_buff);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CMTDHistory::Save(const CSmallString& name)
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

void CMTDHistory::Save(FILE* fout)
{
    if(fout == NULL) {
        INVALID_ARGUMENT("stream is not open");
    }

// 10  format(A3,1X,A2,1X,I2)
// 20  format(I2,1X,A10,1X,E18.11,1X,E18.11,1X,I6)
// 30  format(I2,1X,A55)

// write MTD header ------------------
    if( fprintf(fout,"MTD V3 %2d\n",NCVs) <= 0) {
        RUNTIME_ERROR("unable to write header");
    }

// write coordinate specification ----------------
    for(int i=0; i < NCVs; i++) {
        if(fprintf(fout,"%2d %10s %18.11E %18.11E %6d\n",i+1,(const char*)CVs[i].Type,
                   CVs[i].MinValue,CVs[i].MaxValue,CVs[i].NBins) <= 0) {
            CSmallString error;
            error << "unable to write coordinate definition, id: " << i+1;
            RUNTIME_ERROR(error);
        }
        if(fprintf(fout,"%2d %55s\n",i+1,(const char*)CVs[i].Name) <= 0) {
            CSmallString error;
            error << "unable to write coordinate definition, id: " << i+1;
            RUNTIME_ERROR(error);
        }
    }

//  write(MED_RST,'(2X,I5,2X,I10,1X,F10.3,1X)',ADVANCE='NO')
//  do i=1, fnitem
//     write(MED_RST,'(F10.4,1X,F10.4,1X)',ADVANCE='NO')
//  end do

// write records --------------------------------
    CSimpleIterator<CMTDBuffer>  I(Buffers);
    CMTDBuffer*                  p_buf;
    int                          mtd_time;
    int                          mtd_level;

    while((p_buf = I.Current())!= NULL) {
        // write level and MTD time
        mtd_level = p_buf->GetLevel();
        mtd_time  = p_buf->GetStart();

        for(int i=0; i < p_buf->GetNumOfHills(); i++) {
            if(fprintf(fout,"  %5d  %10d ",mtd_level,mtd_time) <= 0) {
                CSmallString error;
                error << "unable to write restart level and MTD time, pos: " << mtd_time;
                RUNTIME_ERROR(error);
            }

            // write height
            if(fprintf(fout,"%10.3f ",p_buf->GetHeight(i)) <= 0) {
                CSmallString error;
                error << "unable to write height, pos: " << mtd_time;
                RUNTIME_ERROR(error);
            }

            // write values and widths
            for(int j=0; j < GetNumOfCVs(); j++) {
                if(fprintf(fout,"%10.4f %10.4f ",p_buf->GetValue(i,j),p_buf->GetWidth(i,j)) <= 0) {
                    CSmallString error;
                    error << "unable to write value and/or width, pos: " << mtd_time << ", item: " << j+1;
                    RUNTIME_ERROR(error);
                }
            }
            fprintf(fout,"\n");
            mtd_time++;
        }
        I++;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CMTDHistory::SetNumOfCVs(int numofcoords)
{
    if(numofcoords <= 0) {
        LOGIC_ERROR("numofcoords <= 0");
    }

    if(NCVs > 0) {
        // destroy all previous data
        Deallocate();
        CVs.FreeVector();
        NCVs = 0;
    }

// try to allocate CVs array
    CVs.CreateVector(numofcoords);

// all seems to be fine - update items
    NCVs = numofcoords;
}

//------------------------------------------------------------------------------

void CMTDHistory::SetEnergyUnit(const CSmallString& unit, double unit_fac)
{
    EnergyUnit = unit;
    EnergyUnitFac = unit_fac;
}

//------------------------------------------------------------------------------

void CMTDHistory::SetCV(int id,
                        const CSmallString& name,
                        const CSmallString& type,
                        double min_value,double max_value,unsigned int nbins)
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

    if(Buffers.NumOfMembers() > 0) {
        // it was already finalized - destroy data
        Deallocate();
    }

    CVs[id].Type = type;
    CVs[id].Name = name;

    CVs[id].MinValue = min_value;
    CVs[id].MaxValue = max_value;

    CVs[id].NBins = nbins;
    CVs[id].BinWidth = (CVs[id].MaxValue - CVs[id].MinValue)/CVs[id].NBins;
    CVs[id].Width = CVs[id].MaxValue - CVs[id].MinValue;
}

//------------------------------------------------------------------------------

void CMTDHistory::Deallocate(void)
{
    Buffers.RemoveAll();   // items are owned by list
}

//------------------------------------------------------------------------------

void CMTDHistory::ReallocateAsSingleBuffer(void)
{
// create new buffer
    CMTDBuffer* p_buffer = new CMTDBuffer;

    int num_of_hills = GetNumOfHills();
    p_buffer->AllocateData(num_of_hills,NCVs);

// copy data
    for(int i=0; i < num_of_hills; i++) {
        p_buffer->SetHeight(i,GetHeight(i));
        for(int j=0; j < NCVs; j++) {
            p_buffer->SetValue(i,j,GetValue(i,j));
            p_buffer->SetWidth(i,j,GetWidth(i,j));
        }
        p_buffer->IncNumberOfHills();
    }

// destroy previous buffers
    Deallocate();

// insert created buffer to the list
    Buffers.InsertToEnd(p_buffer,0,true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CMTDHistory::GetNumOfCVs(void) const
{
    return(NCVs);
}

//------------------------------------------------------------------------------

const CSmallString& CMTDHistory::GetEnergyUnit(void) const
{
    return(EnergyUnit);
}

//------------------------------------------------------------------------------

double CMTDHistory::GetEnergyUnitFac(void) const
{
    return(EnergyUnitFac);
}

//------------------------------------------------------------------------------

int CMTDHistory::GetNumOfHills(void) const
{
    return(GetNumOfHills(Buffers));
}

//------------------------------------------------------------------------------

int CMTDHistory::GetNumOfHills(const CSimpleList<CMTDBuffer>& list) const
{
    unsigned int nsamples = 0;

    CSimpleIteratorC<CMTDBuffer>  I(list);
    const CMTDBuffer*             p_buf;

    while((p_buf = I.Current())!= NULL) {
        nsamples += p_buf->GetNumOfHills();
        I++;
    }

    return(nsamples);
}

//------------------------------------------------------------------------------

const CColVariable* CMTDHistory::GetCV(int cv) const
{
    return(&CVs[cv]);
}

//------------------------------------------------------------------------------

const double& CMTDHistory::GetHeight(int hill_index)
{
    CSimpleIterator<CMTDBuffer>     I(Buffers);
    CMTDBuffer*                     p_buf;
    int                             loc = hill_index;
    static double                   zero = 0.0;

    while((p_buf = I.Current())!= NULL) {
        if(loc < p_buf->GetNumOfHills()) {
            return(p_buf->GetHeight(loc));
        }
        loc = loc - p_buf->GetNumOfHills();
        I++;
    }

    return(zero);
}

//------------------------------------------------------------------------------

const double& CMTDHistory::GetValue(int hill_index, int cv)
{
    CSimpleIterator<CMTDBuffer>     I(Buffers);
    CMTDBuffer*                     p_buf;
    int                             loc = hill_index;
    static double                   zero = 0.0;

    while((p_buf = I.Current())!= NULL) {
        if(loc < p_buf->GetNumOfHills()) {
            return(p_buf->GetValue(loc,cv));
        }
        loc = loc - p_buf->GetNumOfHills();
        I++;
    }

    return(zero);
}

//------------------------------------------------------------------------------

const double& CMTDHistory::GetWidth(int hill_index, int cv)
{
    CSimpleIterator<CMTDBuffer>     I(Buffers);
    CMTDBuffer*                     p_buf;
    int                             loc = hill_index;
    static double                   zero = 0.0;

    while((p_buf = I.Current())!= NULL) {
        if(loc < p_buf->GetNumOfHills()) {
            return(p_buf->GetWidth(loc,cv));
        }
        loc = loc - p_buf->GetNumOfHills();
        I++;
    }

    return(zero);
}

//------------------------------------------------------------------------------

CMTDBuffer* CMTDHistory::GetBuffer(int index)
{
    if(index >= Buffers.NumOfMembers()) return(NULL);
    return(Buffers[index]);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CMTDHistory::GetPVector(CSimpleVector<double>& pvector)
{
    pvector.FreeVector();

    int num_of_params = GetNumOfHills()*(1+2*GetNumOfCVs());
    if(num_of_params <= 0) return;

    pvector.CreateVector(num_of_params);

    CSimpleIterator<CMTDBuffer>    I(Buffers);
    CMTDBuffer*                    p_buf;
    int                            loc = 0;

    while((p_buf = I.Current())!= NULL) {
        for(int i=0; i < p_buf->GetNumOfHills(); i++) {
            pvector[loc++] = p_buf->GetHeight(i);
            for(int k=0; k < GetNumOfCVs(); k++) {
                pvector[loc++] = p_buf->GetValue(i,k);
                pvector[loc++] = p_buf->GetWidth(i,k);
            }
        }
        I++;
    }
}

//------------------------------------------------------------------------------

void CMTDHistory::SetPVector(const CSimpleVector<double>& pvector)
{
    int             curr_position = 0;
    CMTDBuffer*     p_buff = NULL;
    size_t          loc = 0;

    Deallocate();

    while(loc < pvector.GetLength()) {

        if((p_buff == NULL) || (curr_position >= MaxBufferSize)) {
            // allocate new buffer
            p_buff = GetNewBuffer(MaxBufferSize);
            curr_position = 0;
        }

        // now set height
        p_buff->SetHeight(curr_position,pvector[loc++]);

        for(int i=0; i < GetNumOfCVs(); i++) {
            // set value and width
            p_buff->SetValue(curr_position,i,pvector[loc++]);
            p_buff->SetWidth(curr_position,i,pvector[loc++]);
        }
        p_buff->IncNumberOfHills();
        curr_position++;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

double CMTDHistory::CalculateValue(const CSimpleVector<double>& point,int mtdtime)
{
    CSimpleIterator<CMTDBuffer>    I(Buffers);
    CMTDBuffer*                    p_buf;
    double                         value = 0.0;
    int                            processed = 0;
    double                         fexparg;

    while((p_buf = I.Current())!= NULL) {
        for(int i=0; i < p_buf->GetNumOfHills(); i++) {
            fexparg = 0.0;
            for(int k=0; k < GetNumOfCVs(); k++) {
                double e = point[k] - p_buf->GetValue(i,k);
                fexparg = fexparg +
                          e*e / (2.0 * p_buf->GetWidth(i,k) * p_buf->GetWidth(i,k));
            }
            value = value + p_buf->GetHeight(i)*exp(-fexparg);
            processed = processed + 1;
            if((mtdtime > 0) && (processed >= mtdtime)) return(value);
        }
        I++;
    }

    return(value);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CMTDHistory::LoadCVSInfo(CXMLElement* p_iele)
{
    if(p_iele == NULL) {
        INVALID_ARGUMENT("p_iele is NULL");
    }

    Deallocate();

    CXMLElement* p_ele = p_iele->GetFirstChildElement("CVS");

    if(p_ele == NULL) {
        LOGIC_ERROR("unable to get CVS element");
    }

    bool result = true;

    int lnitems;
    result &= p_ele->GetAttribute("nitems",lnitems);

    if(result == false) {
        LOGIC_ERROR("unable to get header attributes");
    }

    SetNumOfCVs(lnitems);

    CXMLElement*   p_cel = NULL;
    if(p_ele != NULL) p_cel = p_ele->GetFirstChildElement("CVS");

    int   ccount = 0;

    while(p_cel != NULL) {
        if(ccount >= NCVs) {
            LOGIC_ERROR("more COORD elements than NCVs");
        }
        CVs[ccount].LoadInfo(p_cel);
        ccount++;
        p_cel = p_cel->GetNextSiblingElement("CV");
    }
}

//------------------------------------------------------------------------------

bool CMTDHistory::CheckCVSInfo(CXMLElement* p_iele)
{
    if(p_iele == NULL) {
        ES_ERROR("p_iele is NULL");
        return(false);
    }

    CXMLElement* p_ele = p_iele->GetFirstChildElement("CVS");

    if(p_ele == NULL) {
        ES_ERROR("unable to get CVS element");
        return(false);
    }

    bool result = true;

    int lnitems;

    result &= p_ele->GetAttribute("nitems",lnitems);

    if(result == false) {
        ES_ERROR("unable to get header attributes");
        return(false);
    }

    if(NCVs != lnitems) {
        ES_ERROR("mismatch in the number of coordinates");
        return(false);
    }

    CXMLElement*   p_cel = NULL;
    if(p_ele != NULL) p_cel = p_ele->GetFirstChildElement("CV");
    int   ccount = 0;

    while(p_cel != NULL) {
        if(ccount >= NCVs) {
            ES_ERROR("more COORD elements than NCVs");
            return(false);
        }
        if(CVs[ccount].CheckInfo(p_cel) == false) {
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

void CMTDHistory::SaveCVSInfo(CXMLElement* p_iele)
{
    if(p_iele == NULL) {
        INVALID_ARGUMENT("p_iele is NULL");
    }

    CXMLElement* p_ele = p_iele->CreateChildElement("CVS");
    p_ele->SetAttribute("nitems",NCVs);

    CColVariable*  p_coord;

    for(int i=0; i < NCVs; i++) {
        CXMLElement* p_iele = p_ele->CreateChildElement("COORD");
        p_coord = &CVs[i];
        p_coord->SaveInfo(p_iele);
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CMTDHistory::ReadMTDData(CXMLElement* p_ele)
{
    if(p_ele == NULL) {
        INVALID_ARGUMENT("p_ele is NULL");
    }

// destroy previous data
    Deallocate();

// read new one
    AddMTDData(p_ele);
}

//------------------------------------------------------------------------------

void CMTDHistory::AddMTDData(CXMLElement* p_ele)
{
    if(p_ele == NULL) {
        INVALID_ARGUMENT("p_ele is NULL");
    }

    CXMLElement* p_history = p_ele->GetFirstChildElement("HISTORY");
    if(p_history == NULL) {
        LOGIC_ERROR("unable to open HISTORY element");
    }

    unsigned int num_of_hills; // this is total number of hills

    if(p_history->GetAttribute("numofhills",num_of_hills) == false) {
        LOGIC_ERROR("unable to get number of hills");
    }

    if(num_of_hills == 0) return;   // there are no data in received history

    CXMLBinData* p_hele = p_history->GetFirstChildBinData("BUFFER");

    unsigned int control_num_of_hills = 0;

    while(p_hele != NULL) {
        unsigned int num_of_values;
        if(p_hele->GetAttribute("numofhills",num_of_values) == false) {
            LOGIC_ERROR("unable to get number of hills in buffer");
        }

        control_num_of_hills += num_of_values;

        // now allocate buffer
        CMTDBuffer* p_buffer = GetNewBuffer(num_of_values);
        p_buffer->ReadBufferData(p_hele);

        p_hele = p_hele->GetNextSiblingBinData("BUFFER");
    }

    if(control_num_of_hills != num_of_hills) {
        LOGIC_ERROR("mismatch in number of hills");
    }
}

//------------------------------------------------------------------------------

CMTDBuffer* CMTDHistory::AddMTDDataAsSingleBuffer(CXMLElement* p_ele)
{
    if(p_ele == NULL) {
        INVALID_ARGUMENT("p_ele is NULL");
    }

    CXMLElement* p_history = p_ele->GetFirstChildElement("HISTORY");
    if(p_history == NULL) {
        LOGIC_ERROR("unable to open HISTORY element");
    }

    int num_of_hills; // this is total number of hills

    if(p_history->GetAttribute("numofhills",num_of_hills) == false) {
        LOGIC_ERROR("unable to get number of hills");
    }

    if(num_of_hills == 0) return(NULL);   // there are no data in received history

    CMTDBuffer* p_buffer = GetNewBuffer(num_of_hills);

    CXMLBinData* p_hele = p_history->GetFirstChildBinData("BUFFER");

    if(p_hele != NULL) {
        // set level and start from the first buffer
        int level_id;
        int start_id;
        if(p_hele->GetAttribute("level",level_id) == false) {
            LOGIC_ERROR("unable to get buffer level id");
        }
        if(p_hele->GetAttribute("start",start_id) == false) {
            LOGIC_ERROR("unable to get buffer start id");
        }
        p_buffer->SetLevel(level_id);
        p_buffer->SetStart(start_id);
    }

    int control_num_of_hills = 0;

    while(p_hele != NULL) {
        int num_of_values;
        if(p_hele->GetAttribute("numofhills",num_of_values) == false) {
            LOGIC_ERROR("unable to get number of hills in buffer");
        }

        int history_size = num_of_values*(1 + 2*NCVs)*sizeof(double);

        if((history_size == 0) || (history_size != (int)p_hele->GetLength())) {
            LOGIC_ERROR("inconsistent history size");
        }

        double* src = (double*)p_hele->GetData();

        if(src == NULL) {
            LOGIC_ERROR("data array is NULL");
        }

        // copy data
        for(int i=0; i < num_of_values; i++) {
            p_buffer->SetHeight(control_num_of_hills+i,*src++);
            for(int j=0; j < NCVs; j++) {
                p_buffer->SetValue(control_num_of_hills+i,j,*src++);
                p_buffer->SetWidth(control_num_of_hills+i,j,*src++);
            }
            p_buffer->IncNumberOfHills();
        }

        control_num_of_hills += num_of_values;

        p_hele = p_hele->GetNextSiblingBinData("BUFFER");
    }

    if(control_num_of_hills != num_of_hills) {
        LOGIC_ERROR("mismatch in number of hills");
    }

    return(p_buffer);
}

//------------------------------------------------------------------------------

void CMTDHistory::WriteMTDData(CXMLElement* p_ele) const
{
    WriteMTDData(p_ele,Buffers);
}

//------------------------------------------------------------------------------

void CMTDHistory::WriteMTDData(CXMLElement* p_ele,const CSimpleList<CMTDBuffer>& list) const
{
    if(p_ele == NULL) {
        INVALID_ARGUMENT("p_ele is NULL");
    }

    CXMLElement* p_history = p_ele->CreateChildElement("HISTORY");

    unsigned int num_of_hills = GetNumOfHills(list);

    p_history->SetAttribute("numofhills",num_of_hills);

    if(num_of_hills == 0) return;   // there are no data to transmit

    CSimpleIteratorC<CMTDBuffer>  I(list);
    const CMTDBuffer*             p_buf;

    while((p_buf = I.Current())!= NULL) {
        CXMLBinData* p_bele = p_history->CreateChildBinData("BUFFER");
        p_buf->WriteBufferData(p_bele);
        I++;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================



// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
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
#include <PMFAccuData.hpp>
#include <ErrorSystem.hpp>
#include <XMLElement.hpp>
#include <XMLBinData.hpp>
#include <iomanip>

//------------------------------------------------------------------------------

using namespace std;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CPMFAccuData::CPMFAccuData(int nbins, int ncvs,int nstlim)
{
    NumOfBins = nbins;
    NumOfCVs = ncvs;
    NSTLimit = nstlim;
}

//------------------------------------------------------------------------------

void CPMFAccuData::InitDataBlock(int len)
{
    // len is used only some cases otherwise it is derived from NumOfBins or NumOfCVs

    Data.clear();

    Size = 0;
    if( Mode == "B" ){
        CVectorDataPtr data = CVectorDataPtr(new CSimpleVector<double>);
        data->CreateVector(NumOfBins);
        Data.push_back(data);
        Size  = NumOfBins;
    } else if ( Mode == "C" ) {
        CVectorDataPtr data = CVectorDataPtr(new CSimpleVector<double>);
        data->CreateVector(NumOfCVs);
        Data.push_back(data);
        Size  = NumOfCVs;
    } else if ( Mode == "D" ) {
        if( len <= 0 ){
            CSmallString error;
            error << "illegal data length for 'D' mode, len = " << len ;
            RUNTIME_ERROR(error);
        }
        CVectorDataPtr data = CVectorDataPtr(new CSimpleVector<double>);
        Size = len;
        data->CreateVector(Size);
        Data.push_back(data);
    } else if ( Mode == "M" ) {
        for(int icv=0; icv < NumOfCVs; icv++){
            CVectorDataPtr data = CVectorDataPtr(new CSimpleVector<double>);
            data->CreateVector(NumOfBins);
            Data.push_back(data);
        }
        Size  = NumOfCVs * NumOfBins;
    } else if ( Mode == "T" ) {
        CVectorDataPtr data = CVectorDataPtr(new CSimpleVector<double>);
        data->CreateVector(NSTLimit);
        Data.push_back(data);
        Size  = NSTLimit;
    } else if ( Mode == "S" ) {
        for(int icv=0; icv < NumOfCVs; icv++){
            CVectorDataPtr data = CVectorDataPtr(new CSimpleVector<double>);
            data->CreateVector(NSTLimit);
            Data.push_back(data);
        }
        Size  = NSTLimit*NumOfCVs;
    } else if ( Mode == "Z" ) {
        for(int icv=0; icv < NumOfCVs; icv++){
            for(int icv=0; icv < NumOfCVs; icv++){
                CVectorDataPtr data = CVectorDataPtr(new CSimpleVector<double>);
                data->CreateVector(NSTLimit);
                Data.push_back(data);
            }
        }
        Size  = NSTLimit*NumOfCVs*NumOfCVs;
    } else {
        CSmallString error;
        error << "unsupported mode '" << Mode << "'";
        RUNTIME_ERROR(error);
    }

    for(size_t idx=0; idx < Data.size(); idx++){
        Data[idx]->SetZero();
    }
}

//------------------------------------------------------------------------------

void CPMFAccuData::Load(FILE* p_fin,const CSmallString& keyline)
{
    CSmallString skey;
    skey.SetLength(21);

    Op.SetLength(2);
    Type.SetLength(1);
    Mode.SetLength(1);

    // 5  format('@',A20,1X,A2,1X,A1,1X,A1,1X,I10[,1X,A20[,1X,A20]])
    int len = 0;
    if( sscanf(keyline,"%21c %2c %1c %1c %10d",skey.GetBuffer(),Op.GetBuffer(),Type.GetBuffer(),Mode.GetBuffer(),&len) != 5 ){
        CSmallString error;
        error << "unable to parse keyline: '" << keyline << "'";
        RUNTIME_ERROR(error);
    }
    if( skey.GetLength() >= 2 ){
        Name = skey.GetSubString(1,skey.GetLength()-1);
    }

    Name.Trim();
    Op.Trim();
    Type.Trim();
    Mode.Trim();

    if( keyline.GetLength() >= 60 ){
        MSName = keyline.GetSubStringFromTo(40,59);
    }
    if( keyline.GetLength() >= 81 ){
        MXName = keyline.GetSubStringFromTo(61,80);
    }
    if( keyline.GetLength() >= 102 ){
        MYName = keyline.GetSubStringFromTo(82,101);
    }
    MSName.Trim();
    MXName.Trim();
    MYName.Trim();

    if( len <= 0 ){
        CSmallString error;
        error << "illegal data size for keyline: '" << keyline << "'";
        RUNTIME_ERROR(error);
    }

    InitDataBlock(len);

    if( len != Size ){
        CSmallString error;
        error << "inconsistent size " << len << " != " << Size;
        RUNTIME_ERROR(error);
    }

    size_t seg = 0;
    size_t idx = 0;
    CVectorDataPtr data = Data[seg];

    // read data
    if( Type == "I" ){
        for(int i=0; i < Size; i++){
            if( idx == data->GetLength() ){
                idx = 0;
                seg++;
                data = Data[seg];
            }
            int value;
            if( fscanf(p_fin,"%d",&value) != 1 ){
                CSmallString error;
                error << "unable to read data record for keyline: '" << keyline << "'";
                RUNTIME_ERROR(error);
            }
            data->GetRawDataField()[idx] = value;
            idx++;
        }
    } else if ( Type == "R" ){
        for(int i=0; i < Size; i++){
            if( idx == data->GetLength() ){
                idx = 0;
                seg++;
                data = Data[seg];
            }
            double value;
            if( fscanf(p_fin,"%lf",&value) != 1 ){
                CSmallString error;
                error << "unable to read data record for keyline: '" << keyline << "'";
                RUNTIME_ERROR(error);
            }
            data->GetRawDataField()[idx] = value;
            idx++;
        }
    } else {
        CSmallString error;
        error << "unsupported data type for keyline: '" << keyline << "'";
        RUNTIME_ERROR(error);
    }

    // finish reading to the end of line
    skey.ReadLineFromFile(p_fin,true,true);
}

//------------------------------------------------------------------------------

void CPMFAccuData::Save(FILE* p_fout)
{
    // write keyline
    // 5  format('@',A20,1X,A2,1X,A1,1X,A1,1X,I10)
    if( fprintf(p_fout,"@%-20s %-2s %1s %1s %10d",(const char*)Name,(const char*)Op,(const char*)Type,(const char*)Mode,Size) <= 0 ){
        CSmallString error;
        error << "unable to write data section keyline";
        RUNTIME_ERROR(error);
    }
    if( MSName != NULL ){
        if( fprintf(p_fout," %-20s",(const char*)MSName) <= 0 ){
            CSmallString error;
            error << "unable to write data section keyline - msname";
            RUNTIME_ERROR(error);
        }
    }
    if( MXName != NULL ){
        if( fprintf(p_fout," %-20s",(const char*)MXName) <= 0 ){
            CSmallString error;
            error << "unable to write data section keyline - mxname";
            RUNTIME_ERROR(error);
        }
    }
    if( MYName != NULL ){
        if( fprintf(p_fout," %-20s",(const char*)MYName) <= 0 ){
            CSmallString error;
            error << "unable to write data section keyline - myname";
            RUNTIME_ERROR(error);
        }
    }
    if( fprintf(p_fout,"\n") <= 0 ){
        CSmallString error;
        error << "unable to write data section keyline - end";
        RUNTIME_ERROR(error);
    }

    size_t seg = 0;
    size_t idx = 0;
    CVectorDataPtr data = Data[seg];

    // write data
    if( Type == "I" ){
        for(int i=0; i < Size; i++){
            if( idx == data->GetLength() ){
                idx = 0;
                seg++;
                data = Data[seg];
                fprintf(p_fout,"\n");
            } else {
                if( (i % 8 == 0) && (i != 0) ) fprintf(p_fout,"\n");
            }
            int value = data->GetRawDataField()[idx];
            idx++;
            // 10  format(8(I9,1X))
            if( fprintf(p_fout,"%9d ",value) <= 0 ){
                CSmallString error;
                error << "unable to write data record";
                RUNTIME_ERROR(error);
            }
        }
    } else if ( Type == "R" ){
        for(int i=0; i < Size; i++){
            if( idx == data->GetLength() ){
                idx = 0;
                seg++;
                data = Data[seg];
                fprintf(p_fout,"\n");
            } else {
                if( (i % 4 == 0) && (i != 0) ) fprintf(p_fout,"\n");
            }

            double value = data->GetRawDataField()[idx];
            idx++;
            // 10  format(4(E19.11,1X))
            if( fprintf(p_fout,"%23.15le ",value) <= 0 ){
                CSmallString error;
                error << "unable to write data record";
                RUNTIME_ERROR(error);
            }
        }
    } else {
        CSmallString error;
        error << "unsupported data type for keyline";
        RUNTIME_ERROR(error);
    }

    fprintf(p_fout,"\n");
}

//------------------------------------------------------------------------------

void CPMFAccuData::Load(CXMLElement* p_ele)
{
    if( p_ele == NULL ) {
        INVALID_ARGUMENT("p_ele is NULL");
    }
    if( p_ele->GetName() != "DATA" ){
        RUNTIME_ERROR("element is not DATA");
    }

// set in constructor
// NumOfCVs
// NumOfBins

    bool result = true;
    int  len = 0;

    result &= p_ele->GetAttribute("name",Name);
    result &= p_ele->GetAttribute("op",Op);
    result &= p_ele->GetAttribute("type",Type);
    result &= p_ele->GetAttribute("mode",Mode);
    result &= p_ele->GetAttribute("size",len);
    result &= p_ele->GetAttribute("msname",MSName);
    result &= p_ele->GetAttribute("mxname",MXName);
    result &= p_ele->GetAttribute("myname",MYName);

    if( result == false ){
        RUNTIME_ERROR("unable to get DATA attributes");
    }

    InitDataBlock(len);

    if( len != Size ){
        CSmallString error;
        error << "inconsistent size " << len << " != " << Size;
        RUNTIME_ERROR(error);
    }

    CXMLBinData* p_bitem = p_ele->GetFirstChildBinData("BLOB");
    if( p_bitem == NULL ){
        RUNTIME_ERROR("unable to get BLOB");
    }
    len = p_bitem->GetLength<double>();
    double* ptr = p_bitem->GetData<double>();

    if( len != Size ){
        CSmallString error;
        error << "inconsistent lengths, blob: " << len << ", section: " << Size;
        RUNTIME_ERROR(error);
    }

    size_t seg = 0;
    size_t idx = 0;
    CVectorDataPtr data = Data[seg];

    for(int i=0; i < len; i++){
        if( idx == data->GetLength() ){
            idx = 0;
            seg++;
            data = Data[seg];
        }
        data->GetRawDataField()[idx] = *ptr;
        ptr++;
        idx++;
    }
}

//------------------------------------------------------------------------------

void CPMFAccuData::Save(CXMLElement* p_ele)
{
    if( p_ele == NULL ) {
        INVALID_ARGUMENT("p_ele is NULL");
    }
    if( p_ele->GetName() != "DATA" ){
        RUNTIME_ERROR("element is not DATA");
    }

// set in constructor
// NumOfCVs
// NumOfBins

    p_ele->SetAttribute("name",Name);
    p_ele->SetAttribute("op",Op);
    p_ele->SetAttribute("type",Type);
    p_ele->SetAttribute("mode",Mode);
    p_ele->SetAttribute("size",Size);
    p_ele->SetAttribute("msname",MSName);
    p_ele->SetAttribute("mxname",MXName);
    p_ele->SetAttribute("myname",MYName);

    CSimpleVector<double>   lcopy;
    lcopy.CreateVector(Size);

    size_t seg = 0;
    size_t idx = 0;
    CVectorDataPtr data = Data[seg];

    for(int i=0; i < Size; i++){
        if( idx == data->GetLength() ){
            idx = 0;
            seg++;
            data = Data[seg];
        }
        lcopy[i] = data->GetRawDataField()[idx];
        idx++;
    }

    CXMLBinData* p_bitem = p_ele->CreateChildBinData("BLOB");
    p_bitem->CopyData(lcopy,lcopy.GetLength()*sizeof(double),EXBDT_DOUBLE);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CPMFAccuData::Reset(void)
{
    for(size_t idx=0; idx < Data.size(); idx++){
        Data[idx]->SetZero();
    }
}

//------------------------------------------------------------------------------

const CSmallString& CPMFAccuData::GetName(void) const
{
    return(Name);
}

//------------------------------------------------------------------------------

const CSmallString& CPMFAccuData::GetMSName(void) const
{
    return(MSName);
}

//------------------------------------------------------------------------------

const CSmallString& CPMFAccuData::GetMXName(void) const
{
    return(MXName);
}

//------------------------------------------------------------------------------

const CSmallString& CPMFAccuData::GetMYName(void) const
{
    return(MYName);
}

//------------------------------------------------------------------------------

int CPMFAccuData::GetLength(void) const
{
    return(Size);
}

//------------------------------------------------------------------------------

int CPMFAccuData::GetNumOfBins(void) const
{
    if( (Mode == 'B') || (Mode == 'M') ) return(NumOfBins);
    return(0);
}

//------------------------------------------------------------------------------

int CPMFAccuData::GetNumOfCVs(void) const
{
    if( (Mode == 'C') || (Mode == 'M') ) return(NumOfCVs);
    return(0);
}

//------------------------------------------------------------------------------

const CSmallString& CPMFAccuData::GetType(void) const
{
    return(Type);
}

//------------------------------------------------------------------------------

const CSmallString& CPMFAccuData::GetMode(void) const
{
    return(Mode);
}

//------------------------------------------------------------------------------

const CSmallString& CPMFAccuData::GetOp(void) const
{
    return(Op);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

double CPMFAccuData::GetData(int indi) const
{
    if( Mode == "M" ) RUNTIME_ERROR("not applicable for the M mode");
    if( Mode == "S" ) RUNTIME_ERROR("not applicable for the S mode");
    if( Mode == "Z" ) RUNTIME_ERROR("not applicable for the Z mode");

    if( Mode == "B" ){
        if( (indi < 0) || (NumOfBins <= indi) ) {
            RUNTIME_ERROR("indi out-of-range for the B mode");
        }
        CVectorDataPtr data = Data[0];
        return( data->GetRawDataField()[indi] );
    }

    if( Mode == "T" ){
        if( (indi < 0) || (NSTLimit <= indi) ) {
            RUNTIME_ERROR("indi out-of-range for the T mode");
        }
        CVectorDataPtr data = Data[0];
        return( data->GetRawDataField()[indi] );
    }

    if( Mode == "C" ){
        if( (indi < 0) || (NumOfCVs <= indi) ) {
            RUNTIME_ERROR("indi out-of-range for the C mode");
        }
        CVectorDataPtr data = Data[0];
        return( data->GetRawDataField()[indi] );
    }

    if( Mode == "D" ){
        if( (indi < 0) || (Size <= indi) ) {
            RUNTIME_ERROR("indi out-of-range for the D mode");
        }
        CVectorDataPtr data = Data[0];
        return( data->GetRawDataField()[indi] );
    }

    CSmallString error;
    error << "unsupported mode '" << Mode << "'";
    RUNTIME_ERROR(error);
}

//------------------------------------------------------------------------------

double CPMFAccuData::GetData(int indi, int icv) const
{
    if( Mode == "C" ) RUNTIME_ERROR("not applicable for the C mode");
    if( Mode == "Z" ) RUNTIME_ERROR("not applicable for the Z mode");

    if( Mode == "B" ){
        if( (indi < 0) || (NumOfBins <= indi) ) {
            RUNTIME_ERROR("indi out-of-range for the B mode");
        }
        // ignore icv
        CVectorDataPtr data = Data[0];
        return( data->GetRawDataField()[indi] );
    }

    if( Mode == "M" ){
        if( (indi < 0) || (NumOfBins <= indi) ) {
            RUNTIME_ERROR("indi out-of-range for the M mode");
        }
        if( (icv < 0) || (NumOfCVs <= icv) ) {
            RUNTIME_ERROR("icv out-of-range for the M mode");
        }
        CVectorDataPtr data = Data[icv];
        return( data->GetRawDataField()[indi] );
    }

    if( Mode == "T" ){
        if( (indi < 0) || (NSTLimit <= indi) ) {
            RUNTIME_ERROR("indi out-of-range for the T mode");
        }
        // ignore icv
        CVectorDataPtr data = Data[0];
        return( data->GetRawDataField()[indi] );
    }

    if( Mode == "S" ){
        if( (indi < 0) || (NSTLimit <= indi) ) {
            RUNTIME_ERROR("indi out-of-range for the S mode");
        }
        if( (icv < 0) || (NumOfCVs <= icv) ) {
            RUNTIME_ERROR("icv out-of-range for the S mode");
        }
        CVectorDataPtr data = Data[icv];
        return( data->GetRawDataField()[indi] );
    }

    CSmallString error;
    error << "unsupported mode '" << Mode << "'";
    RUNTIME_ERROR(error);
}

//------------------------------------------------------------------------------

double CPMFAccuData::GetData(int indi, int icv,int jcv) const
{
    if( Mode == "C" ) RUNTIME_ERROR("not applicable for the C mode");
    if( Mode == "B" ) RUNTIME_ERROR("not applicable for the B mode");
    if( Mode == "M" ) RUNTIME_ERROR("not applicable for the M mode");
    if( Mode == "T" ) RUNTIME_ERROR("not applicable for the T mode");
    if( Mode == "S" ) RUNTIME_ERROR("not applicable for the S mode");

    if( Mode == "Z" ){
        if( (indi < 0) || (NSTLimit <= indi) ) {
            RUNTIME_ERROR("indi out-of-range for the Z mode");
        }
        if( (icv < 0) || (NumOfCVs <= icv) ) {
            RUNTIME_ERROR("icv out-of-range for the Z mode");
        }
        if( (jcv < 0) || (NumOfCVs <= jcv) ) {
            RUNTIME_ERROR("jcv out-of-range for the Z mode");
        }
        size_t seg = icv + NumOfCVs*jcv;
        CVectorDataPtr data = Data[seg];
        return( data->GetRawDataField()[indi] );
    }

    CSmallString error;
    error << "unsupported mode '" << Mode << "'";
    RUNTIME_ERROR(error);
}

//------------------------------------------------------------------------------

void CPMFAccuData::SetData(int indi, double value)
{
    if( Mode == "M" ) RUNTIME_ERROR("not applicable for the M mode");
    if( Mode == "S" ) RUNTIME_ERROR("not applicable for the S mode");
    if( Mode == "Z" ) RUNTIME_ERROR("not applicable for the Z mode");

    if( Mode == "B" ){
        if( (indi < 0) || (NumOfBins <= indi) ) {
            RUNTIME_ERROR("indi out-of-range for the B mode");
        }
        CVectorDataPtr data = Data[0];
        data->GetRawDataField()[indi] = value;
        return;
    }

    if( Mode == "T" ){
        if( (indi < 0) || (NSTLimit <= indi) ) {
            RUNTIME_ERROR("indi out-of-range for the T mode");
        }
        CVectorDataPtr data = Data[0];
        data->GetRawDataField()[indi] = value;
        return;
    }

    if( Mode == "C" ){
        if( (indi < 0) || (NumOfCVs <= indi) ) {
            RUNTIME_ERROR("indi out-of-range for the C mode");
        }
        CVectorDataPtr data = Data[0];
        data->GetRawDataField()[indi] = value;
        return;
    }

    if( Mode == "D" ){
        if( (indi < 0) || (Size <= indi) ) {
            RUNTIME_ERROR("indi out-of-range for the D mode");
        }
        CVectorDataPtr data = Data[0];
        data->GetRawDataField()[indi] = value;
        return;
    }

    CSmallString error;
    error << "unsupported mode '" << Mode << "'";
    RUNTIME_ERROR(error);
}

//------------------------------------------------------------------------------

void CPMFAccuData::SetData(int indi, int icv, double value)
{
    if( Mode == "B" ) RUNTIME_ERROR("not applicable for the B mode");
    if( Mode == "T" ) RUNTIME_ERROR("not applicable for the T mode");
    if( Mode == "C" ) RUNTIME_ERROR("not applicable for the C mode");
    if( Mode == "Z" ) RUNTIME_ERROR("not applicable for the Z mode");

    if( Mode == "M" ){
        if( (indi < 0) || (NumOfBins <= indi) ) {
            RUNTIME_ERROR("indi out-of-range for the M mode");
        }
        if( (icv < 0) || (NumOfCVs <= icv) ) {
            RUNTIME_ERROR("icv out-of-range for the M mode");
        }
        CVectorDataPtr data = Data[icv];
        data->GetRawDataField()[indi] = value;
        return;
    }

    if( Mode == "S" ){
        if( (indi < 0) || (NSTLimit <= indi) ) {
            RUNTIME_ERROR("indi out-of-range for the S mode");
        }
        if( (icv < 0) || (NumOfCVs <= icv) ) {
            RUNTIME_ERROR("icv out-of-range for the S mode");
        }
        CVectorDataPtr data = Data[icv];
        data->GetRawDataField()[indi] = value;
        return;
    }

    CSmallString error;
    error << "unsupported mode '" << Mode << "'";
    RUNTIME_ERROR(error);
}

//------------------------------------------------------------------------------

void CPMFAccuData::GetDataBlob(double* p_blob)
{
    for(size_t seg = 0; seg < Data.size(); seg++){
        CVectorDataPtr data  = Data[seg];
        for(size_t idx = 0; idx < data->GetLength(); idx++){
            *p_blob = data->GetRawDataField()[idx];
            p_blob++;
        }
    }
}

//------------------------------------------------------------------------------

void CPMFAccuData::SetDataBlob(double* p_blob)
{
    for(size_t seg = 0; seg < Data.size(); seg++){
        CVectorDataPtr data  = Data[seg];
        for(size_t idx = 0; idx < data->GetLength(); idx++){
            data->GetRawDataField()[idx] = *p_blob;
            p_blob++;
        }
    }
}

//------------------------------------------------------------------------------

CVectorDataPtr CPMFAccuData::GetDataBlob(int icv,int jcv)
{
    int seg = icv + NumOfCVs*jcv;
    return(Data[seg]);
}

//------------------------------------------------------------------------------

void CPMFAccuData::GetDataBlob(int indi,CSimpleVector<double>& data) const
{
    for(size_t indj=0; indj < data.GetLength(); indj++){
        data[indj] = Data[indj]->GetRawDataField()[indi];
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CPMFAccuData::CheckCompatibility(CPMFAccuDataPtr right) const
{
    if( NumOfBins != right->NumOfBins ){
        CSmallString   error;
        error << "two data sections are not compatible (number of bins): " << NumOfBins << " vs " << right->NumOfBins;
        ES_ERROR(error);
        return(false);
    }
    if( NumOfCVs != right->NumOfCVs ){
        CSmallString   error;
        error << "two data sections are not compatible (number of CVs): " << NumOfCVs << " vs " << right->NumOfCVs;
        ES_ERROR(error);
        return(false);
    }
    if( Name != right->Name ){
        CSmallString   error;
        error << "two data sections are not compatible (name): " << Name << " vs " << right->Name;
        ES_ERROR(error);
        return(false);
    }
    if( Op != right->Op ){
        CSmallString   error;
        error << "two data sections are not compatible (op): " << Op << " vs " << right->Op;
        ES_ERROR(error);
        return(false);
    }
    if( Type != right->Type ){
        CSmallString   error;
        error << "two data sections are not compatible (type): " << Type << " vs " << right->Type;
        ES_ERROR(error);
        return(false);
    }
    if( Mode != right->Mode ){
        CSmallString   error;
        error << "two data sections are not compatible (mode): " << Mode << " vs " << right->Mode;
        ES_ERROR(error);
        return(false);
    }
    if( Size != right->Size ){
        CSmallString   error;
        error << "two data sections are not compatible (mode): " << Size << " vs " << right->Size;
        ES_ERROR(error);
        return(false);
    }

    if( MSName != right->MSName ){
        CSmallString   error;
        error << "two data sections are not compatible (msname): " << MSName << " vs " << right->MSName;
        ES_ERROR(error);
        return(false);
    }
    if( MXName != right->MXName ){
        CSmallString   error;
        error << "two data sections are not compatible (mxname): " << MXName << " vs " << right->MXName;
        ES_ERROR(error);
        return(false);
    }
    if( MYName != right->MYName ){
        CSmallString   error;
        error << "two data sections are not compatible (myname): " << MYName << " vs " << right->MYName;
        ES_ERROR(error);
        return(false);
    }

    return(true);
}

//------------------------------------------------------------------------------

CPMFAccuDataPtr CPMFAccuData::CreateTheSame(void) const
{
    CPMFAccuDataPtr ptr = CPMFAccuDataPtr(new CPMFAccuData(NumOfBins,NumOfCVs,NSTLimit));
    ptr->Name   = Name;
    ptr->Op     = Op;
    ptr->Type   = Type;
    ptr->Mode   = Mode;
    ptr->MSName = MSName;
    ptr->MXName = MXName;
    ptr->MYName = MYName;

    ptr->InitDataBlock(Size);

    if(  ptr->Size != Size ){
        CSmallString error;
        error << "inconsistent size " << ptr->Size << " != " << Size;
        RUNTIME_ERROR(error);
    }

    return(ptr);
}

//------------------------------------------------------------------------------

CPMFAccuDataPtr CPMFAccuData::Duplicate(void) const
{
    CPMFAccuDataPtr dup = CreateTheSame();

    for(size_t seg = 0; seg < Data.size(); seg++){
        CVectorDataPtr indata  = Data[seg];
        CVectorDataPtr outdata = dup->Data[seg];
        for(size_t idx = 0; idx < indata->GetLength(); idx++){
            outdata->GetRawDataField()[idx] = indata->GetRawDataField()[idx];
        }
    }

    return(dup);
}

//------------------------------------------------------------------------------

void CPMFAccuData::CombineAD(CPMFAccuDataPtr left,CPMFAccuDataPtr right)
{
    if( CheckCompatibility(left) == false ){
        RUNTIME_ERROR("incompatible with left");
    }
    if( CheckCompatibility(right) == false ){
        RUNTIME_ERROR("incompatible with right");
    }

    if( (Mode != "B") && (Mode != "M") && (Mode != "C") ){
        CSmallString error;
        error <<  "unsupported mode for data section, name: " << Name << ", mode: " << Mode;
        RUNTIME_ERROR(error);
    }

    for(size_t seg = 0; seg < Data.size(); seg++){
        CVectorDataPtr out    = Data[seg];
        CVectorDataPtr ldata  = left->Data[seg];
        CVectorDataPtr rdata  = right->Data[seg];

        for(size_t idx = 0; idx < out->GetLength(); idx++){
            out->GetRawDataField()[idx] = ldata->GetRawDataField()[idx] + rdata->GetRawDataField()[idx];
        }
    }
}

//------------------------------------------------------------------------------

void CPMFAccuData::CombineSA(CPMFAccuDataPtr left,CPMFAccuDataPtr right)
{
    if( CheckCompatibility(left) == false ){
        RUNTIME_ERROR("incompatible with left");
    }
    if( CheckCompatibility(right) == false ){
        RUNTIME_ERROR("incompatible with right");
    }

    for(size_t seg = 0; seg < Data.size(); seg++){
        CVectorDataPtr out    = Data[seg];
        CVectorDataPtr ldata  = left->Data[seg];
        CVectorDataPtr rdata  = right->Data[seg];

        for(size_t idx = 0; idx < out->GetLength(); idx++){
            if( fabs(ldata->GetRawDataField()[idx] - rdata->GetRawDataField()[idx]) > 1e-7 ){
                stringstream serror;
                serror << "data are not the same for '" << GetName() << "', item: " << (idx+1) << ", values: "
                       << ldata->GetRawDataField()[idx] << ", " << rdata->GetRawDataField()[idx];
                RUNTIME_ERROR(serror.str());
            }
            out->GetRawDataField()[idx] = ldata->GetRawDataField()[idx];
        }
    }
}

//------------------------------------------------------------------------------

void CPMFAccuData::CombineWA(CPMFAccuDataPtr left,CPMFAccuDataPtr left_nsamples,CPMFAccuDataPtr right,CPMFAccuDataPtr right_nsamples)
{
    if( CheckCompatibility(left) == false ){
        RUNTIME_ERROR("incompatible with left");
    }
    if( CheckCompatibility(right) == false ){
        RUNTIME_ERROR("incompatible with right");
    }

    if( Mode == "B" ){
        for(int ibin=0; ibin < NumOfBins; ibin++){
            long double ln = left_nsamples->GetData(ibin);
            long double rn = right_nsamples->GetData(ibin);
            if( (ln + rn) != 0.0 ) {
                long double lw = ln / (ln + rn);
                long double rw = rn / (ln + rn);
                long double ld = left->GetData(ibin);
                long double rd = right->GetData(ibin);
                long double re = lw*ld + rw*rd;
                SetData(ibin,re);
            }
        }
    } else if( Mode == "M" ) {
        for(int icv=0; icv < NumOfCVs; icv++){
            for(int ibin=0; ibin < NumOfBins; ibin++){
                long double ln = left_nsamples->GetData(ibin);
                long double rn = right_nsamples->GetData(ibin);
                if( (ln + rn) != 0.0 ) {
                    long double lw = ln / (ln + rn);
                    long double rw = rn / (ln + rn);
                    long double ld = left->GetData(ibin,icv);
                    long double rd = right->GetData(ibin,icv);
                    long double re = lw*ld + rw*rd;
                    // cout << ln << " " << lw << " " << ld << " " << rn << " " << rw << " " << rd << " " << re << endl;
                    SetData(ibin,icv,re);
                }
            }
        }
    } else if( Mode == "C" ) {
        CSmallString error;
        error <<  "WA operation and C mode are unsupported for data section combine operation, name: " << Name;
        RUNTIME_ERROR(error);
    } else {
        CSmallString error;
        error <<  "unsupported mode for data section, name: " << Name << ", mode: " << Mode;
        RUNTIME_ERROR(error);
    }
}

//------------------------------------------------------------------------------

void CPMFAccuData::CombineM2(CPMFAccuDataPtr left,CPMFAccuDataPtr left_nsamples,CPMFAccuDataPtr left_mean,
                             CPMFAccuDataPtr right,CPMFAccuDataPtr right_nsamples,CPMFAccuDataPtr right_mean)
{
    if( CheckCompatibility(left) == false ){
        RUNTIME_ERROR("incompatible with left");
    }
    if( CheckCompatibility(right) == false ){
        RUNTIME_ERROR("incompatible with right");
    }

    if( Mode == "B" ){
        for(int ibin=0; ibin < NumOfBins; ibin++){
            long double ln = left_nsamples->GetData(ibin);
            long double rn = right_nsamples->GetData(ibin);
            if( (ln + rn) != 0.0 ) {
                long double w  = ln * rn / (ln + rn);
                long double dx = left_mean->GetData(ibin) - right_mean->GetData(ibin);
                long double ld = left->GetData(ibin);
                long double rd = right->GetData(ibin);
                long double re = ld + rd + dx*dx*w;
                SetData(ibin,re);
            }
        }
    } else if( Mode == "M" ) {
        for(int icv=0; icv < NumOfCVs; icv++){
            for(int ibin=0; ibin < NumOfBins; ibin++){
                long double ln = left_nsamples->GetData(ibin);
                long double rn = right_nsamples->GetData(ibin);
                if( (ln + rn) != 0.0 ) {
                    long double w  = ln * rn / (ln + rn);
                    long double dx = left_mean->GetData(ibin,icv) - right_mean->GetData(ibin,icv);
                    long double ld = left->GetData(ibin,icv);
                    long double rd = right->GetData(ibin,icv);
                    long double re = ld + rd + dx*dx*w;
                    SetData(ibin,icv,re);
                }
            }
        }
    } else if( Mode == "C" ) {
        CSmallString error;
        error <<  "M2 operation and C mode are unsupported for data section combine operation, name: " << Name;
        RUNTIME_ERROR(error);
    } else {
        CSmallString error;
        error <<  "unsupported mode for data section, name: " << Name << ", mode: " << Mode;
        RUNTIME_ERROR(error);
    }
}

//------------------------------------------------------------------------------

void CPMFAccuData::CombineCO(CPMFAccuDataPtr left,CPMFAccuDataPtr left_nsamples,CPMFAccuDataPtr left_xmean,CPMFAccuDataPtr left_ymean,
                             CPMFAccuDataPtr right,CPMFAccuDataPtr right_nsamples,CPMFAccuDataPtr right_xmean,CPMFAccuDataPtr right_ymean)
{
    if( CheckCompatibility(left) == false ){
        RUNTIME_ERROR("incompatible with left");
    }
    if( CheckCompatibility(right) == false ){
        RUNTIME_ERROR("incompatible with right");
    }

    if( Mode == "B" ){
        for(int ibin=0; ibin < NumOfBins; ibin++){
            long double ln = left_nsamples->GetData(ibin);
            long double rn = right_nsamples->GetData(ibin);
            if( (ln + rn) != 0.0 ) {
                long double w  = ln * rn / (ln + rn);
                long double dx = left_xmean->GetData(ibin) - right_xmean->GetData(ibin);
                long double dy = left_ymean->GetData(ibin) - right_ymean->GetData(ibin);
                long double ld = left->GetData(ibin);
                long double rd = right->GetData(ibin);
                long double re = ld + rd + dx*dy*w;
                SetData(ibin,re);
            }
        }
    } else if( Mode == "M" ) {
        for(int icv=0; icv < NumOfCVs; icv++){
            for(int ibin=0; ibin < NumOfBins; ibin++){
                long double ln = left_nsamples->GetData(ibin);
                long double rn = right_nsamples->GetData(ibin);
                if( (ln + rn) != 0.0 ) {
                    long double w  = ln * rn / (ln + rn);
                    long double dx = left_xmean->GetData(ibin,icv) - right_xmean->GetData(ibin,icv);
                    long double dy = left_ymean->GetData(ibin,icv) - right_ymean->GetData(ibin,icv);
                    long double ld = left->GetData(ibin,icv);
                    long double rd = right->GetData(ibin,icv);
                    long double re = ld + rd + dx*dy*w;
                    SetData(ibin,icv,re);
                }
            }
        }
    } else if( Mode == "C" ) {
        CSmallString error;
        error <<  "CO operation and C mode are unsupported for data section combine operation, name: " << Name;
        RUNTIME_ERROR(error);
    } else {
        CSmallString error;
        error <<  "unsupported mode for data section, name: " << Name << ", mode: " << Mode;
        RUNTIME_ERROR(error);
    }
}

//------------------------------------------------------------------------------

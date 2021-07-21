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

CPMFAccuData::CPMFAccuData(int nbins, int ncvs)
{
    NumOfBins = nbins;
    NumOfCVs = ncvs;
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
    Size = 0;
    if( sscanf(keyline,"%21c %2c %1c %1c %10d",skey.GetBuffer(),Op.GetBuffer(),Type.GetBuffer(),Mode.GetBuffer(),&Size) != 5 ){
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
        MXName = keyline.GetSubStringFromTo(40,59);
    }
    if( keyline.GetLength() >= 81 ){
        MYName = keyline.GetSubStringFromTo(61,80);
    }
    MXName.Trim();
    MYName.Trim();

    if( Size <= 0 ){
        CSmallString error;
        error << "illegal data size for keyline: '" << keyline << "'";
        RUNTIME_ERROR(error);
    }

    if( ! ( ((Mode == 'B') && (Size == NumOfBins)) ||
            ((Mode == 'M') && (Size == NumOfBins*NumOfCVs)) ||
            ((Mode == 'C') && (Size == NumOfCVs)) ||
            ( Mode == 'D') ) ) {
        CSmallString error;
        error << "data size and mode is not consistent: '" << keyline << "'";
        RUNTIME_ERROR(error);
    }

    Data.CreateVector(Size);

    // read data
    if( Type == "I" ){
        for(int i=0; i < Size; i++){
            int value;
            if( fscanf(p_fin,"%d",&value) != 1 ){
                CSmallString error;
                error << "unable to read data record for keyline: '" << keyline << "'";
                RUNTIME_ERROR(error);
            }
            Data[i] = value;
        }
    } else if ( Type == "R" ){
        for(int i=0; i < Size; i++){
            double value;
            if( fscanf(p_fin,"%lf",&value) != 1 ){
                CSmallString error;
                error << "unable to read data record for keyline: '" << keyline << "'";
                RUNTIME_ERROR(error);
            }
            Data[i] = value;
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

    // write data
    if( Type == "I" ){
        for(int i=0; i < Size; i++){
            if( (i % 8 == 0) && (i != 0) ) fprintf(p_fout,"\n");
            int value = Data[i];
            // 10  format(8(I9,1X))
            if( fprintf(p_fout,"%9d ",value) <= 0 ){
                CSmallString error;
                error << "unable to write data record";
                RUNTIME_ERROR(error);
            }
        }
    } else if ( Type == "R" ){
        for(int i=0; i < Size; i++){
            if( (i % 4 == 0) && (i != 0) ) fprintf(p_fout,"\n");
            double value = Data[i];
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

    result &= p_ele->GetAttribute("name",Name);
    result &= p_ele->GetAttribute("op",Op);
    result &= p_ele->GetAttribute("type",Type);
    result &= p_ele->GetAttribute("mode",Mode);
    result &= p_ele->GetAttribute("size",Size);
    result &= p_ele->GetAttribute("mxname",MXName);
    result &= p_ele->GetAttribute("myname",MYName);

    if( result == false ){
        RUNTIME_ERROR("unable to get DATA attributes");
    }

    CXMLBinData* p_bitem = p_ele->GetFirstChildBinData("BLOB");
    if( p_bitem == NULL ){
        RUNTIME_ERROR("unable to get BLOB");
    }
    int     len = p_bitem->GetLength<double>();
    double* ptr = p_bitem->GetData<double>();

    if( len != Size ){
        CSmallString error;
        error << "inconsistent lengths, blob: " << len << ", section: " << Size;
        RUNTIME_ERROR(error);
    }

    Data.CreateVector(len);
    for(int i=0; i < len; i++){
        Data[i] = *ptr;
        ptr++;
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
    p_ele->SetAttribute("mxname",MXName);
    p_ele->SetAttribute("myname",MYName);

    CXMLBinData* p_bitem = p_ele->CreateChildBinData("BLOB");
    p_bitem->CopyData(Data,Data.GetLength()*sizeof(double),EXBDT_DOUBLE);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CPMFAccuData::Reset(void)
{
    Data.SetZero();
}

//------------------------------------------------------------------------------

const CSmallString& CPMFAccuData::GetName(void) const
{
    return(Name);
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

//------------------------------------------------------------------------------

double CPMFAccuData::GetData(int ibin, int cv) const
{
    if( (ibin < 0) || (NumOfBins <= ibin) ) {
        RUNTIME_ERROR("ibin out-of-range");
    }
    if( (cv < 0) || (NumOfCVs <= cv) ) {
        RUNTIME_ERROR("cv out-of-range");
    }
    int idx = 0;
    if( NumOfBins*NumOfCVs == Size ){
        idx = cv*NumOfBins + ibin;
    } else {
        idx = ibin;
    }
    return( Data[idx] );
}

//------------------------------------------------------------------------------

void CPMFAccuData::SetData(int ibin, double value)
{
    if( (ibin < 0) || (NumOfBins <= ibin) ) {
        RUNTIME_ERROR("ibin out-of-range");
    }
    Data[ibin] = value;
}

//------------------------------------------------------------------------------

void CPMFAccuData::SetData(int ibin, int cv, double value)
{
    if( (ibin < 0) || (NumOfBins <= ibin) ) {
        RUNTIME_ERROR("ibin out-of-range");
    }
    if( (cv < 0) || (NumOfCVs <= cv) ) {
        RUNTIME_ERROR("cv out-of-range");
    }
    int idx = 0;
    if( NumOfBins*NumOfCVs == Size ){
        idx = cv*NumOfBins + ibin;
    } else {
        idx = ibin;
    }
    Data[idx] = value;
}

//------------------------------------------------------------------------------

void CPMFAccuData::GetDataBlob(double* p_blob)
{
    double* p_source = Data;

    for(size_t i=0; i < Data.GetLength(); i++){
        *p_blob = *p_source;
        p_blob++;
        p_source++;
    }
}

//------------------------------------------------------------------------------

void CPMFAccuData::SetDataBlob(double* p_blob)
{
    double* p_dest = Data;

    for(size_t i=0; i < Data.GetLength(); i++){
        *p_dest = *p_blob;
        p_blob++;
        p_dest++;
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
    CPMFAccuDataPtr ptr = CPMFAccuDataPtr(new CPMFAccuData(NumOfBins,NumOfCVs));
    ptr->Name   = Name;
    ptr->Op     = Op;
    ptr->Type   = Type;
    ptr->Mode   = Mode;
    ptr->Size   = Size;
    ptr->MXName = MXName;
    ptr->MYName = MYName;

    ptr->Data.CreateVector(Size);
    ptr->Data.SetZero();

    return(ptr);
}

//------------------------------------------------------------------------------

CPMFAccuDataPtr CPMFAccuData::Duplicate(void) const
{
    CPMFAccuDataPtr data = CreateTheSame();
    for(int i=0; i < Size; i++){
        data->Data[i] = Data[i];
    }
    return(data);
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

    for(int i=0; i < Size; i++){
        Data[i] = left->Data[i] + right->Data[i];
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

    for(int i=0; i < Size; i++){
        // FIXME - should we allow some difference?
        if( left->Data[i] != right->Data[i] ){
            stringstream serror;
            serror << "data are not the same for '" << GetName() << "', item: " << (i+1) << ", values: " << left->Data[i] << ", " << right->Data[i];
            RUNTIME_ERROR(serror.str());
        }
        Data[i] = left->Data[i];
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
            double ln = left_nsamples->GetData(ibin);
            double rn = right_nsamples->GetData(ibin);
            if( (ln + rn) != 0.0 ) {
                double lw = ln / (ln + rn);
                double rw = rn / (ln + rn);
                double ld = left->GetData(ibin);
                double rd = right->GetData(ibin);
                double re = lw*ld + rw*rd;
                SetData(ibin,re);
            }
        }
    } else if( Mode == "M" ) {
        for(int icv=0; icv < NumOfCVs; icv++){
            for(int ibin=0; ibin < NumOfBins; ibin++){
                double ln = left_nsamples->GetData(ibin);
                double rn = right_nsamples->GetData(ibin);
                if( (ln + rn) != 0.0 ) {
                    double lw = ln / (ln + rn);
                    double rw = rn / (ln + rn);
                    double ld = left->GetData(ibin,icv);
                    double rd = right->GetData(ibin,icv);
                    double re = lw*ld + rw*rd;
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
            double ln = left_nsamples->GetData(ibin);
            double rn = right_nsamples->GetData(ibin);
            if( (ln + rn) != 0.0 ) {
                double w  = ln * rn / (ln + rn);
                double dx = left_mean->GetData(ibin) - right_mean->GetData(ibin);
                double ld = left->GetData(ibin);
                double rd = right->GetData(ibin);
                double re = ld + rd + dx*dx*w;
                SetData(ibin,re);
            }
        }
    } else if( Mode == "M" ) {
        for(int icv=0; icv < NumOfCVs; icv++){
            for(int ibin=0; ibin < NumOfBins; ibin++){
                double ln = left_nsamples->GetData(ibin);
                double rn = right_nsamples->GetData(ibin);
                if( (ln + rn) != 0.0 ) {
                    double w  = ln * rn / (ln + rn);
                    double dx = left_mean->GetData(ibin,icv) - right_mean->GetData(ibin,icv);
                    double ld = left->GetData(ibin,icv);
                    double rd = right->GetData(ibin,icv);
                    double re = ld + rd + dx*dx*w;
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
            double ln = left_nsamples->GetData(ibin);
            double rn = right_nsamples->GetData(ibin);
            if( (ln + rn) != 0.0 ) {
                double w  = ln * rn / (ln + rn);
                double dx = left_xmean->GetData(ibin) - right_xmean->GetData(ibin);
                double dy = left_ymean->GetData(ibin) - right_ymean->GetData(ibin);
                double ld = left->GetData(ibin);
                double rd = right->GetData(ibin);
                double re = ld + rd + dx*dy*w;
                SetData(ibin,re);
            }
        }
    } else if( Mode == "M" ) {
        for(int icv=0; icv < NumOfCVs; icv++){
            for(int ibin=0; ibin < NumOfBins; ibin++){
                double ln = left_nsamples->GetData(ibin);
                double rn = right_nsamples->GetData(ibin);
                if( (ln + rn) != 0.0 ) {
                    double w  = ln * rn / (ln + rn);
                    double dx = left_xmean->GetData(ibin,icv) - right_xmean->GetData(ibin,icv);
                    double dy = left_ymean->GetData(ibin,icv) - right_ymean->GetData(ibin,icv);
                    double ld = left->GetData(ibin,icv);
                    double rd = right->GetData(ibin,icv);
                    double re = ld + rd + dx*dy*w;
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

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
#include <PMFAccumulator.hpp>
#include <ErrorSystem.hpp>
#include <XMLElement.hpp>
#include <XMLBinData.hpp>
#include <iomanip>
#include <boost/format.hpp>

//------------------------------------------------------------------------------

using namespace std;
using namespace boost;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CPMFAccumulator::CPMFAccumulator(void)
{
    NumOfCVs            = 0;
    NumOfBins           = 0;
    Temperature         = 300.0;
    TemperatureFConv    = 1.0;
    TemperatureUnit     = "K";
    EnergyFConv         = 1.0;
    EnergyUnit          = "kcal mol^-1";
    Method              = "NONE";
    Version             = LibBuildVersion_PMF;
    NCorr               = 1.0;
}

//------------------------------------------------------------------------------

CPMFAccumulator::~CPMFAccumulator(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CPMFAccumulator::Load(const CSmallString& name)
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

void CPMFAccumulator::Load(FILE* fin)
{
    if( fin == NULL ) {
        INVALID_ARGUMENT("stream is not open");
    }

    CSmallString keyline;

// read sections
    while( keyline.ReadLineFromFile(fin,true,true) ){
        if( IsTrajectorySnapshot(keyline) ) return;
        if( IsHeaderSection(keyline) ){
            ReadHeaderSection(fin,keyline);
        } else {
            ReadDataSection(fin,keyline);
        }
    }
}

//------------------------------------------------------------------------------

std::list<CPMFAccumulatorPtr> CPMFAccumulator::LoadFinalSnapshots(const CSmallString& name, size_t number)
{
    std::list<CPMFAccumulatorPtr> accus;

    FILE* fin = fopen(name, "r");

    if(fin == NULL) {
        CSmallString error;
        error << "unable to open file '" << name << "' (" << strerror(errno) << ")";
        RUNTIME_ERROR(error);
    }

    try {
        accus = LoadFinalSnapshots(fin,number);
    } catch(...) {
        fclose(fin);
        throw;
    }

    fclose(fin);

    return(accus);
}


//------------------------------------------------------------------------------

std::list<CPMFAccumulatorPtr> CPMFAccumulator::LoadFinalSnapshots(FILE* fin, size_t number)
{
    if( fin == NULL ) {
        INVALID_ARGUMENT("stream is not open");
    }

    std::list<CPMFAccumulatorPtr> accus;

    CSmallString keyline;

    CPMFAccumulatorPtr accu;

// read sections
    while( keyline.ReadLineFromFile(fin,true,true) ){
        if( IsTrajectoryHeader(keyline) ) {
            continue;   // skip
        }
        if( IsTrajectorySnapshot(keyline) ){
            accu = CPMFAccumulatorPtr(new CPMFAccumulator);
            accus.push_back(accu);
            if( accus.size() > number ){
                accus.pop_front();
            }
            continue;
        } else {
            if( accu == NULL ){
                accu = CPMFAccumulatorPtr(new CPMFAccumulator);
                accus.push_back(accu);
                if( accus.size() > number ){
                    accus.pop_front();
                }
            }
        }

        if( IsHeaderSection(keyline) ){
            accu->ReadHeaderSection(fin,keyline);
        } else {
            accu->ReadDataSection(fin,keyline);
        }
    }

    return(accus);
}

//------------------------------------------------------------------------------

void CPMFAccumulator::Load(CXMLElement* p_root)
{
    if( p_root == NULL ) {
        INVALID_ARGUMENT("p_root is NULL");
    }
    if( p_root->GetName() != "PMFLIB-V6" ){
        CSmallString error;
        error << "root XML element is not PMFLIB-V6 but: " << p_root->GetName() ;
        RUNTIME_ERROR(error);
    }

// root element
    int ncvs = 0;
    p_root->GetAttribute("numofcvs",ncvs);
    if( ncvs <= 0 ){
        RUNTIME_ERROR("no CVs");
    }
    SetNumOfCVs(ncvs);

    CXMLElement* p_item;

// header items ----------------------------------
    p_item = p_root->GetFirstChildElement("HEADER");
    while( p_item != NULL ){
        CSmallString name;
        bool result = true;
        result &= p_item->GetAttribute("name",name);
    // -------------------------------------------
        if( name == "VERSION" ){
            result &= p_item->GetAttribute("value",Version);
    // -------------------------------------------
        } else if (  name == "METHOD" ) {
            result &= p_item->GetAttribute("value",Method);
    // -------------------------------------------
        } else if (  name == "DRIVER" ) {
            result &= p_item->GetAttribute("value",Driver);
    // -------------------------------------------
        } else if (  name == "TEMPERATURE" ) {
            result &= p_item->GetAttribute("value",Temperature);
    // -------------------------------------------
        } else if (  name == "ENERGY-UNIT" ) {
            result &= p_item->GetAttribute("value",EnergyUnit);
            result &= p_item->GetAttribute("fconv",EnergyFConv);
    // -------------------------------------------
        } else if (  name == "TEMPERATURE-UNIT" ) {
            result &= p_item->GetAttribute("value",TemperatureUnit);
            result &= p_item->GetAttribute("fconv",TemperatureFConv);
    // -------------------------------------------
        } else {
            CSmallString error;
            error << "unrecognized header item: " << name;
            RUNTIME_ERROR(name);
        }
        if( result == false ){
            RUNTIME_ERROR("unable to get some HEADER attributes");
        }
        p_item = p_item->GetNextSiblingElement("HEADER");
    }

// load CVs
    LoadCVSInfo(p_root);

// load data sections
    p_item = p_root->GetFirstChildElement("DATA");
    while( p_item != NULL ){
        CSmallString name;
        p_item->GetAttribute("name",name);
        if( name == NULL ){
            RUNTIME_ERROR("DATA element does not have name");
        }
        CPMFAccuDataPtr data = CPMFAccuDataPtr(new CPMFAccuData(NumOfBins,NumOfCVs));
        data->Load(p_item);
        DataBlocks[name] = data;
        p_item = p_item->GetNextSiblingElement("DATA");
    }
}

//------------------------------------------------------------------------------

void CPMFAccumulator::LoadSnapshot(const CSmallString& name,int index)
{
    FILE* fin = fopen(name, "r");

    if(fin == NULL) {
        CSmallString error;
        error << "unable to open file '" << name << "' (" << strerror(errno) << ")";
        RUNTIME_ERROR(error);
    }

    try {
        LoadSnapshot(fin,index);
    } catch(...) {
        fclose(fin);
        throw;
    }

    fclose(fin);
}

//------------------------------------------------------------------------------

void CPMFAccumulator::LoadSnapshot(FILE* fin,int index)
{
    CSmallString   line;

    if( line.ReadLineFromFile(fin,true,true) == false ) {
        RUNTIME_ERROR("unable to read first line from input stream");
    }

    if( line != "#ACCUTRAJ" ) {
        RUNTIME_ERROR("input file is not PMF accumulator trajectory file");
    }

// move to requested position
    int     position    = 0;
    bool    found       = false;

    while( line.ReadLineFromFile(fin,true,true) ) {
        if( line.FindSubString("#ACCUSNAP") == 0 ) {
            position++;
            if(position == index) {
                found = true;
                break;
            }
        }
    }

    if( found == false ){
        CSmallString error;
        error << "PMF accumulator trajectory contains less snapshots than " << index;
        RUNTIME_ERROR(error);
    }

// read accumulator
    Load(fin);
}

//------------------------------------------------------------------------------

void CPMFAccumulator::Combine(CPMFAccumulatorPtr right)
{
// methods must be the same
    if( GetMethod() != right->GetMethod() ){
        CSmallString   error;
        error << "two PMF accumulators are not compatible (methods): " << GetMethod() << " vs " << right->GetMethod();
        RUNTIME_ERROR(error);
    }

// temperature must be the same
    if( fabs(GetTemperature() - right->GetTemperature()) > 0.1 ){
        stringstream serror;
        serror << format("two PMF accumulators are not compatible (temperature): %6.2lf %s vs %6.2lf %s")
                  %GetRealTemperature()%GetTemperatureUnit()
                  %right->GetRealTemperature()%right->GetTemperatureUnit() << endl;
        RUNTIME_ERROR(serror.str());
    }

// do not check for
// DRIVER
// VERSION

// size must be consistent
    if( CheckCVSInfo(right) == false ){
        RUNTIME_ERROR("two PMF accumulators are not size consistent");
    }

// combine data blocks
    if( DataBlocks.size() != right->DataBlocks.size() ){
        CSmallString   error;
        error << "two PMF accumulators are not compatible (number of data sections): " << DataBlocks.size()  << " vs " << right->DataBlocks.size();
        RUNTIME_ERROR(error);
    }

    std::map<CSmallString,CPMFAccuDataPtr>::iterator  it = DataBlocks.begin();
    std::map<CSmallString,CPMFAccuDataPtr>::iterator  ie = DataBlocks.end();

    std::map<CSmallString,CPMFAccuDataPtr>  combined;

    while( it != ie ){
        CPMFAccuDataPtr ldb = it->second;
        CPMFAccuDataPtr rdb = right->GetSectionData(it->first);
        CPMFAccuDataPtr result = Combine(right,ldb,rdb);
        combined[it->first] = result;
        it++;
    }

    DataBlocks.clear();
    DataBlocks = combined;
}

//------------------------------------------------------------------------------

CPMFAccuDataPtr CPMFAccumulator::Combine(CPMFAccumulatorPtr right,CPMFAccuDataPtr ldb,CPMFAccuDataPtr rdb)
{
    CPMFAccuDataPtr result = ldb->CreateTheSame();

// add
    if( result->GetOp() == "AD" ){
        result->CombineAD(ldb,rdb);
// keep the same value
    } else if( result->GetOp() == "SA" ){
        result->CombineSA(ldb,rdb);
// weighted average
    } else if( result->GetOp() == "WA" ){
        // we need two NSAMPLES blocks
        CPMFAccuDataPtr lns = GetSectionData("NSAMPLES");
        if( lns == NULL ){
            CSmallString error;
            error <<  "left NSAMPLES is required for WA operation, data section name: " << result->GetName();
            RUNTIME_ERROR(error);
        }
        CPMFAccuDataPtr rns = right->GetSectionData("NSAMPLES");
        if( rns == NULL ){
            CSmallString error;
            error <<  "right NSAMPLES is required for WA operation, data section name: " << result->GetName();
            RUNTIME_ERROR(error);
        }
        result->CombineWA(ldb,lns,rdb,rns);
// M2 - combine
    } else if( result->GetOp() == "M2" ){
        // we need two NSAMPLES blocks
        CPMFAccuDataPtr lns = GetSectionData("NSAMPLES");
        if( lns == NULL ){
            CSmallString error;
            error <<  "left NSAMPLES is required for M2 operation, data section name: " << result->GetName();
            RUNTIME_ERROR(error);
        }
        CPMFAccuDataPtr rns = right->GetSectionData("NSAMPLES");
        if( rns == NULL ){
            CSmallString error;
            error <<  "right NSAMPLES is required for M2 operation, data section name: " << result->GetName();
            RUNTIME_ERROR(error);
        }
        // in addition, we also need two MEAN blocks
        CPMFAccuDataPtr lmv = GetSectionData(result->GetMXName());
        if( lmv == NULL ){
            CSmallString error;
            error <<  "left mean data section (" << result->GetMXName() << ") is required for M2 operation, data section name: " << result->GetName();
            RUNTIME_ERROR(error);
        }
        CPMFAccuDataPtr rmv = right->GetSectionData(result->GetMXName());
        if( rmv == NULL ){
            CSmallString error;
            error <<  "right mean data section (" << result->GetMXName() << ") is required for M2 operation, data section name: " << result->GetName();
            RUNTIME_ERROR(error);
        }
        result->CombineM2(ldb,lns,lmv,rdb,rns,rmv);
// CO - combine
    } else if( result->GetOp() == "CO" ){
        // we need two NSAMPLES blocks
        CPMFAccuDataPtr lns = GetSectionData("NSAMPLES");
        if( lns == NULL ){
            CSmallString error;
            error <<  "left NSAMPLES is required for CO operation, data section name: " << result->GetName();
            RUNTIME_ERROR(error);
        }
        CPMFAccuDataPtr rns = right->GetSectionData("NSAMPLES");
        if( rns == NULL ){
            CSmallString error;
            error <<  "right NSAMPLES is required for CO operation, data section name: " << result->GetName();
            RUNTIME_ERROR(error);
        }
        // in addition, we also need two MEAN blocks for X and Y
        CPMFAccuDataPtr lmvx = GetSectionData(result->GetMXName());
        if( lmvx == NULL ){
            CSmallString error;
            error <<  "left X mean data section (" << result->GetMXName() << ") is required for CO operation, data section name: " << result->GetName();
            RUNTIME_ERROR(error);
        }
        CPMFAccuDataPtr rmvx = right->GetSectionData(result->GetMXName());
        if( rmvx == NULL ){
            CSmallString error;
            error <<  "right X mean data section (" << result->GetMXName() << ") is required for CO operation, data section name: " << result->GetName();
            RUNTIME_ERROR(error);
        }
        CPMFAccuDataPtr lmvy = GetSectionData(result->GetMYName());
        if( lmvy == NULL ){
            CSmallString error;
            error <<  "left Y mean data section (" << result->GetMYName() << ") is required for CO operation, data section name: " << result->GetName();
            RUNTIME_ERROR(error);
        }
        CPMFAccuDataPtr rmvy = right->GetSectionData(result->GetMYName());
        if( rmvy == NULL ){
            CSmallString error;
            error <<  "right Y mean data section (" << result->GetMYName() << ") is required for CO operation, data section name: " << result->GetName();
            RUNTIME_ERROR(error);
        }
        result->CombineCO(ldb,lns,lmvx,lmvy,rdb,rns,rmvx,rmvy);
// unsupported operation
    } else {
        CSmallString error;
        error <<  "unsupported operation for data section, name: " << result->GetName() << ", op: " << result->GetOp();
        RUNTIME_ERROR(error);
    }

    return(result);
}

//------------------------------------------------------------------------------

void CPMFAccumulator::Combine(CXMLElement* p_ele)
{
    CPMFAccumulatorPtr right = CPMFAccumulatorPtr(new CPMFAccumulator);
    right->Load(p_ele);
    Combine(right);
}

//------------------------------------------------------------------------------

bool CPMFAccumulator::IsTrajectoryHeader(const CSmallString& keyline)
{
    if( keyline == "#ACCUTRAJ" ) return(true);
    return(false);
}

//------------------------------------------------------------------------------

bool CPMFAccumulator::IsTrajectorySnapshot(const CSmallString& keyline)
{
    if( keyline.FindSubString("#ACCUSNAP") == 0 ) return(true);
    return(false);
}

//------------------------------------------------------------------------------

bool CPMFAccumulator::IsHeaderSection(const CSmallString& keyline)
{
    if( keyline.GetLength() < 1 ) return(false);
    if( keyline[0] == '%' ) return(true);
    return(false);
}

//------------------------------------------------------------------------------

const CSmallString CPMFAccumulator::GetSectionName(const CSmallString& keyline) const
{
    if( keyline.GetLength() <= 1 ) return("");
    int first = 1;
    int last = keyline.Scan(" \t\n");
    if( (last - 1) >= 1 ){
        last = last - 1;
    } else {
        last = keyline.GetLength() - 1;
    }
    return( keyline.GetSubStringFromTo(first,last) );
}

//------------------------------------------------------------------------------

void CPMFAccumulator::ReadHeaderSection(FILE* fin,const CSmallString& keyline)
{
    if( fin == NULL ) {
        INVALID_ARGUMENT("stream is not open");
    }

    CSmallString    key = GetSectionName(keyline);
    CSmallString    sbuff;

// -----------------------------------------------
    if( key == "PMFLIB-V6") {
        CSmallString ncvs;
        ncvs.ReadLineFromFile(fin,true,true);
        ncvs.Trim();
        SetNumOfCVs(ncvs.ToInt());
// -----------------------------------------------
    } else if( key == "VERSION") {
        Version.ReadLineFromFile(fin,true,true);
        Version.Trim();
// -----------------------------------------------
    } else if( key == "METHOD") {
        Method.ReadLineFromFile(fin,true,true);
        Method.Trim();
// -----------------------------------------------
    } else if( key == "DRIVER") {
        Driver.ReadLineFromFile(fin,true,true);
        Driver.Trim();
// -----------------------------------------------
    } else if( key == "CVS") {
        // read coordinate specification
        for(int i=0; i < NumOfCVs; i++) {
            CSmallString    type;
            CSmallString    name;
            CSmallString    unit;

            int             id = 0;
            double          min_value = 0.0;
            double          max_value = 0.0;
            double          fconv = 1.0;
            int             nbins = 0;
            int             tr = 0;

//            20  format(I2,1X,E18.11,1X,E18.11,1X,I6,1X,A10)
            // read item
            sbuff.ReadLineFromFile(fin,true,true);
            if( sbuff.GetLength() != 58 ){
                CSmallString error;
                error << "unable to read coordinate definition: '" << sbuff << "'";
                RUNTIME_ERROR(error);
            }
            tr = sscanf(sbuff,"%d %lf %lf %d",&id,&min_value,&max_value,&nbins);
            if( tr != 4 ) {
                CSmallString error;
                error << "unable to read coordinate definition, id: " << i+1 << " (" << tr << " != 4)";
                RUNTIME_ERROR(error);
            }
            type = sbuff.GetSubStringFromTo(48,57);
            type.Trim();
            // some tests
            if(id != i+1) {
                CSmallString error;
                error << "coordinate id does not match, read: " << id << ", expected: " << i+1;
                RUNTIME_ERROR(error);
            }
            if(max_value < min_value) {
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
//            25  format(I2,1X,A55)
            sbuff.ReadLineFromFile(fin,true,true);
            if( sbuff.GetLength() != 58 ){
                CSmallString error;
                error << "unable to read coordinate definition: '" << sbuff << "'";
                RUNTIME_ERROR(error);
            }
            tr = sscanf(sbuff,"%d",&id);
            if( tr != 1 ) {
                CSmallString error;
                error << "unable to read coordinate definition, id: " << i+1 << " (" << tr << " != 1)";
                RUNTIME_ERROR(error);
            }
            name = sbuff.GetSubStringFromTo(3,57);
            name.Trim();
            // some tests
            if(id != i+1) {
                CSmallString error;
                error << "coordinate id does not match, read: " << id << ", expected: " << i+1;
                RUNTIME_ERROR(error);
            }

            // read item
//            26  format(I2,1X,E18.11,1X,A36)
            sbuff.ReadLineFromFile(fin,true,true);
            if( sbuff.GetLength() != 58 ){
                CSmallString error;
                error << "unable to read coordinate definition: '" << sbuff << "'";
                RUNTIME_ERROR(error);
            }
            tr = sscanf(sbuff,"%d %lf",&id,&fconv);
            if( tr != 2 ) {
                CSmallString error;
                error << "unable to read coordinate definition, id: " << i+1 << " (" << tr << " != 2)";
                RUNTIME_ERROR(error);
            }
            unit = sbuff.GetSubStringFromTo(22,57);
            unit.Trim();
            // some tests
            if(id != i+1) {
                CSmallString error;
                error << "coordinate id does not match, read: " << id << ", expected: " << i+1;
                RUNTIME_ERROR(error);
            }

            // init CV
            SetCV(i,name,type,min_value,max_value,nbins,fconv,unit);
        }

    // update NumOfBins
    UpdateNumOfBins();

// -----------------------------------------------------
    } else if( key == "TEMPERATURE" ) {
        // read item
        sbuff.ReadLineFromFile(fin,true,true);
        int tr = sscanf(sbuff,"%lf",&Temperature);
        if( tr != 1 ) {
            CSmallString error;
            error << "unable to read temperature (" << tr << " != 1)";
            RUNTIME_ERROR(error);
        }
// -----------------------------------------------------
    } else if( key == "ENERGY-UNIT" ) {
        int             tr = 0;

// 40  format(3X,E18.11,1X,A36)
        // read item
        sbuff.ReadLineFromFile(fin,true,true);
        tr = sscanf(sbuff,"%lf",&EnergyFConv);
        if( tr != 1 ) {
            CSmallString error;
            error << "unable to read energy unit (" << tr << " != 1)";
            RUNTIME_ERROR(error);
        }
        if( sbuff.GetLength() != 58 ){
            CSmallString error;
            error << "unable to read energy unit: '" << sbuff << "'";
            RUNTIME_ERROR(error);
        }
        EnergyUnit = sbuff.GetSubStringFromTo(22,57);
        EnergyUnit.Trim();

    // -----------------------------------------------------
    } else if( key == "TEMPERATURE-UNIT" ) {
        int             tr = 0;

// 40  format(3X,E18.11,1X,A36)
        // read item
        sbuff.ReadLineFromFile(fin,true,true);
        tr = sscanf(sbuff,"%lf",&TemperatureFConv);
        if( tr != 1 ) {
            CSmallString error;
            error << "unable to read temperature unit (" << tr << " != 1)";
            RUNTIME_ERROR(error);
        }
        if( sbuff.GetLength() != 58 ){
            CSmallString error;
            error << "unable to read temperature unit: '" << sbuff << "'";
            RUNTIME_ERROR(error);
        }
        TemperatureUnit = sbuff.GetSubStringFromTo(22,57);
        TemperatureUnit.Trim();

    }else {
        CSmallString error;
        error << "unrecognized ABF accumulator header keyword: '" << key << "'";
        RUNTIME_ERROR(error);
    }
}


//------------------------------------------------------------------------------

void CPMFAccumulator::UpdateNumOfBins(void)
{
// update NumOfBins
    NumOfBins = 1;
    for(int i=0; i < NumOfCVs; i++) {
        NumOfBins *= CVs[i]->NumOfBins;
    }
}

//------------------------------------------------------------------------------

void CPMFAccumulator::ReadDataSection(FILE* fin,const CSmallString& keyline)
{
    CPMFAccuDataPtr data = CPMFAccuDataPtr(new CPMFAccuData(NumOfBins,NumOfCVs));
    data->Load(fin,keyline);
    DataBlocks[data->GetName()] = data;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CPMFAccumulator::Save(const CSmallString& name)
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

void CPMFAccumulator::Save(CXMLElement* p_ele)
{
    if(p_ele == NULL) {
        INVALID_ARGUMENT("p_ele is NULL");
    }
    if( NumOfBins == 0 ) {
        CSmallString error;
        error << "no data in accumulator";
        RUNTIME_ERROR(error);
    }

// root element
    CXMLElement* p_root = p_ele->CreateChildElement("PMFLIB-V6");
    p_root->SetAttribute("numofcvs",NumOfCVs);

    CXMLElement* p_item;

// header item -----------------------------------
    p_item = p_root->CreateChildElement("HEADER");
    p_item->SetAttribute("name","VERSION");
    p_item->SetAttribute("value",Version);

// header item -----------------------------------
    p_item = p_root->CreateChildElement("HEADER");
    p_item->SetAttribute("name","METHOD");
    p_item->SetAttribute("value",Method);

// header item -----------------------------------
    p_item = p_root->CreateChildElement("HEADER");
    p_item->SetAttribute("name","DRIVER");
    p_item->SetAttribute("value",Driver);

// header item -----------------------------------
    p_item = p_root->CreateChildElement("HEADER");
    p_item->SetAttribute("name","TEMPERATURE");
    p_item->SetAttribute("value",Temperature);

// header item -----------------------------------
    p_item = p_root->CreateChildElement("HEADER");
    p_item->SetAttribute("name","ENERGY-UNIT");
    p_item->SetAttribute("value",EnergyUnit);
    p_item->SetAttribute("fconv",EnergyFConv);

// header item -----------------------------------
    p_item = p_root->CreateChildElement("HEADER");
    p_item->SetAttribute("name","TEMPERATURE-UNIT");
    p_item->SetAttribute("value",TemperatureUnit);
    p_item->SetAttribute("fconv",TemperatureFConv);

// save CVs
    SaveCVSInfo(p_root);

// write data sections
    std::map<CSmallString,CPMFAccuDataPtr>::iterator  dit = DataBlocks.begin();
    std::map<CSmallString,CPMFAccuDataPtr>::iterator  die = DataBlocks.end();

    while( dit != die ){
        CXMLElement*    p_block = p_root->CreateChildElement("DATA");
        CPMFAccuDataPtr ds = dit->second;
        ds->Save(p_block);
        dit++;
    }
}

//------------------------------------------------------------------------------

void CPMFAccumulator::Save(FILE* fout)
{
    if(fout == NULL) {
        INVALID_ARGUMENT("stream is not open");
    }

    if( NumOfBins == 0 ) {
        CSmallString error;
        error << "no data in accumulator";
        RUNTIME_ERROR(error);
    }

// write ABF accumulator header ------------------
    if(fprintf(fout,"%%PMFLIB-V6\n%3d\n",NumOfCVs) <= 0) {
        CSmallString error;
        error << "unable to write header";
        RUNTIME_ERROR(error);
    }

    if(fprintf(fout,"%%VERSION\n%s\n",(const char*)Version) <= 0) {
        CSmallString error;
        error << "unable to write version";
        RUNTIME_ERROR(error);
    }

    if(fprintf(fout,"%%METHOD\n%s\n",(const char*)Method) <= 0) {
        CSmallString error;
        error << "unable to write method";
        RUNTIME_ERROR(error);
    }

    if(fprintf(fout,"%%DRIVER\n%s\n",(const char*)Driver) <= 0) {
        CSmallString error;
        error << "unable to write driver";
        RUNTIME_ERROR(error);
    }

    if(fprintf(fout,"%%TEMPERATURE\n%18.11E\n",Temperature) <= 0) {
        CSmallString error;
        error << "unable to write temperature";
        RUNTIME_ERROR(error);
    }

    // 40  format(3X,E18.11,1X,A36)
    if(fprintf(fout,"%%ENERGY-UNIT\n   %18.11E %36s\n",EnergyFConv,(const char*)EnergyUnit) <= 0) {
        CSmallString error;
        error << "unable to write energy unit";
        RUNTIME_ERROR(error);
    }

    // 40  format(3X,E18.11,1X,A36)
    if(fprintf(fout,"%%TEMPERATURE-UNIT\n   %18.11E %36s\n",TemperatureFConv,(const char*)TemperatureUnit) <= 0) {
        CSmallString error;
        error << "unable to write temperature unit";
        RUNTIME_ERROR(error);
    }

// write coordinate specification ----------------
    if(fprintf(fout,"%%CVS\n") <= 0) {
        CSmallString error;
        error << "unable to write CVS header";
        RUNTIME_ERROR(error);
    }

    //30  format(I2,1X,E18.11,1X,E18.11,1X,I6,1X,A10)
    //31  format(I2,1X,A55)
    //32  format(I2,1X,E18.11,1X,A36)

    for(int i=0; i < NumOfCVs; i++) {
        if(fprintf(fout,"%2d %18.11E %18.11E %6d %10s\n",i+1,
                   CVs[i]->MinValue,CVs[i]->MaxValue,CVs[i]->NumOfBins,
                   (const char*)CVs[i]->Type) <= 0) {
            CSmallString error;
            error << "unable to write coordinate definition I id: " << i+1;
            RUNTIME_ERROR(error);
        }
        if(fprintf(fout,"%2d %55s\n",i+1,(const char*)CVs[i]->Name) <= 0) {
            CSmallString error;
            error << "unable to write coordinate definition II id: " << i+1;
            RUNTIME_ERROR(error);
        }
        if(fprintf(fout,"%2d %18.11E %36s\n",i+1,CVs[i]->FConv,(const char*)CVs[i]->Unit) <= 0) {
            CSmallString error;
            error << "unable to write coordinate definition III id: " << i+1;
            RUNTIME_ERROR(error);
        }
    }

// write data section
    std::map<CSmallString,CPMFAccuDataPtr>::iterator  it = DataBlocks.begin();
    std::map<CSmallString,CPMFAccuDataPtr>::iterator  ie = DataBlocks.end();

    while( it != ie ){
        CPMFAccuDataPtr ds = it->second;
        ds->Save(fout);
        it++;
    }
}

//------------------------------------------------------------------------------

void CPMFAccumulator::SetHeaders(const CSmallString& method, const CSmallString& version, const CSmallString& driver,
                double temp, const CSmallString& temp_unit, double temp_fconv,
                const CSmallString& ene_unit, double ene_fconv)
{
    Method = method;
    Version = version;
    Driver = driver;
    Temperature = temp;
    TemperatureUnit = temp_unit;
    TemperatureFConv = temp_fconv;
    EnergyUnit  = ene_unit;
    EnergyFConv = ene_fconv;
}

//------------------------------------------------------------------------------

void CPMFAccumulator::GetPoint(unsigned int index,CSimpleVector<double>& point) const
{
    for(int k=NumOfCVs-1; k >= 0; k--) {
        const CColVariablePtr p_coord = CVs[k];
        int ibin = index % p_coord->GetNumOfBins();
        point[k] = p_coord->GetValue(ibin);
        index = index / p_coord->GetNumOfBins();
    }
}

void CPMFAccumulator::GetPointRValues(unsigned int index,CSimpleVector<double>& point) const
{
    for(int k=NumOfCVs-1; k >= 0; k--) {
        const CColVariablePtr p_coord = CVs[k];
        int ibin = index % p_coord->GetNumOfBins();
        point[k] = p_coord->GetRValue(ibin);
        index = index / p_coord->GetNumOfBins();
    }
}

//------------------------------------------------------------------------------

void CPMFAccumulator::GetIPoint(unsigned int index,CSimpleVector<int>& point) const
{
    for(int k=NumOfCVs-1; k >= 0; k--) {
        const CColVariablePtr p_coord = CVs[k];
        int ibin = index % p_coord->GetNumOfBins();
        point[k] = ibin;
        index = index / p_coord->GetNumOfBins();
    }
}

//------------------------------------------------------------------------------

double CPMFAccumulator::GetTemperature(void) const
{
    return(Temperature);
}

//------------------------------------------------------------------------------

double CPMFAccumulator::GetRealTemperature(void) const
{
    return(Temperature*TemperatureFConv);
}

//------------------------------------------------------------------------------

const CSmallString& CPMFAccumulator::GetTemperatureUnit(void) const
{
    return(TemperatureUnit);
}

//------------------------------------------------------------------------------

double CPMFAccumulator::GetEnergyFConv(void)
{
    return(EnergyFConv);
}

//------------------------------------------------------------------------------

int CPMFAccumulator::GetNumOfSamples(int ibin) const
{
    return(GetData("NSAMPLES",ibin));
}

//------------------------------------------------------------------------------

int CPMFAccumulator::GetTotalNumOfSamples(void) const
{
    int nsamples = 0;
    for(int ibin=0; ibin < NumOfBins; ibin++){
        nsamples += GetData("NSAMPLES",ibin);
    }
    return(nsamples);
}

//------------------------------------------------------------------------------

double CPMFAccumulator::GetEnergyRealValue(double value) const
{
    return(value * EnergyFConv);
}

//------------------------------------------------------------------------------

void CPMFAccumulator::SetNCorr(double ncorr)
{
    NCorr = ncorr;
}

//------------------------------------------------------------------------------

double CPMFAccumulator::GetNCorr(void) const
{
    return(NCorr);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CPMFAccumulator::SetNumOfCVs(int ncvs)
{
    if( ncvs < 0 ) {
        INVALID_ARGUMENT("ncvs < 0");
    }

    if(NumOfCVs > 0) {
        // destroy all previous data
        Clear();
        CVs.clear();
        NumOfCVs = 0;
    }

// try to allocate CVs array
    if( ncvs > 0 ) {
        for(int i=0; i < ncvs; i++){
            CVs.push_back(CColVariablePtr(new CColVariable));
        }
    }

// all seems to be fine - update items
    NumOfCVs = ncvs;
}

//------------------------------------------------------------------------------

void CPMFAccumulator::SetCV(int id,
                            const CSmallString& name,
                            const CSmallString& type,
                            double min_value,double max_value,int nbins)
{
    if( CVs.size() == 0 ){
        RUNTIME_ERROR("no CVs");
    }
    if( id < 0 || id >= NumOfCVs ){
        INVALID_ARGUMENT("id out-of-range");
    }

    if( nbins <= 0 ){
        INVALID_ARGUMENT("nbins <= 0");
    }
    if( max_value < min_value ){
        INVALID_ARGUMENT("max_value < min_value");
    }

    if( DataBlocks.size() != 0  ) {
        // it was already finalized - destroy data
        Clear();
    }

    CVs[id]->ID = id;
    CVs[id]->Name = name;
    CVs[id]->Type = type;
    CVs[id]->MinValue = min_value;
    CVs[id]->MaxValue = max_value;

    CVs[id]->NumOfBins = nbins;
    CVs[id]->BinWidth = (max_value - min_value)/nbins;
    CVs[id]->Width = max_value - min_value;
}

//------------------------------------------------------------------------------

void CPMFAccumulator::SetCV(int id,
                            const CSmallString& name,
                            const CSmallString& type,
                            double min_value,double max_value,int nbins,
                            double fconv, const CSmallString& unit)
{
    if( CVs.size() == 0 ){
        RUNTIME_ERROR("no CVs");
    }
    if( id < 0 || id >= NumOfCVs ){
        INVALID_ARGUMENT("id out-of-range");
    }

    if( nbins <= 0 ){
        INVALID_ARGUMENT("nbins <= 0");
    }
    if( max_value < min_value ){
        INVALID_ARGUMENT("max_value < min_value");
    }

    if( DataBlocks.size() != 0  ) {
        // it was already finalized - destroy data
        Clear();
    }

    CVs[id]->ID = id;
    CVs[id]->Name = name;
    CVs[id]->Type = type;
    CVs[id]->Unit = unit;
    CVs[id]->MinValue = min_value;
    CVs[id]->MaxValue = max_value;
    CVs[id]->FConv = fconv;

    CVs[id]->NumOfBins = nbins;
    CVs[id]->BinWidth = (max_value - min_value)/nbins;
    CVs[id]->Width = max_value - min_value;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

const CSmallString& CPMFAccumulator::GetMethod(void) const
{
    return(Method);
}

//------------------------------------------------------------------------------

const CSmallString& CPMFAccumulator::GetDriver(void) const
{
    return(Driver);
}

//------------------------------------------------------------------------------

int CPMFAccumulator::GetNumOfBins(void) const
{
    return(NumOfBins);
}

//------------------------------------------------------------------------------

int CPMFAccumulator::GetNumOfCVs(void) const
{
    return(NumOfCVs);
}

//------------------------------------------------------------------------------

const CColVariablePtr CPMFAccumulator::GetCV(int cv) const
{
    return(CVs[cv]);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CPMFAccumulator::GetGlobalIndex(const CSimpleVector<int>& position) const
{
    int glbindex = 0;
    for(int i=0; i < NumOfCVs; i++) {
        if( position[i] < 0 ) return(-1);
        if( position[i] >= CVs[i]->GetNumOfBins() ) return(-1);
        glbindex = glbindex*CVs[i]->GetNumOfBins() + position[i];
    }
    return(glbindex);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CPMFAccumulator::Clear(void)
{
// do not destroy CVs array !

// destroy only data arrays
    DataBlocks.clear();
    NumOfBins        = 0;
}

//------------------------------------------------------------------------------

void CPMFAccumulator::Reset(void)
{
    if( DataBlocks.size() != 0  ) {
        return;
    }

    std::map<CSmallString,CPMFAccuDataPtr>::iterator    it =  DataBlocks.begin();
    std::map<CSmallString,CPMFAccuDataPtr>::iterator    ie =  DataBlocks.end();

    while( it != ie ){
        it->second->Reset();
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CPMFAccumulator::LoadCVSInfo(CXMLElement* p_iele)
{
    if(p_iele == NULL) {
        INVALID_ARGUMENT("p_iele is NULL");
    }
    if( NumOfCVs == 0 ){
        RUNTIME_ERROR("no CVs allocated");
    }

    CXMLElement* p_ele = p_iele->GetFirstChildElement("CVS");
    if(p_ele == NULL) {
        RUNTIME_ERROR("unable to get CVS element");
    }

    CXMLElement*   p_cel = p_ele->GetFirstChildElement("CV");
    int            ccount = 0;

    while(p_cel != NULL) {
        if( ccount >= NumOfCVs ) {
            LOGIC_ERROR("more CV elements than NumOfCVs");
        }
        CVs[ccount]->LoadInfo(p_cel);
        ccount++;
        p_cel = p_cel->GetNextSiblingElement("CV");
    }

    // update NumOfBins
    UpdateNumOfBins();
}

//------------------------------------------------------------------------------

bool CPMFAccumulator::CheckCVSInfo(CXMLElement* p_iele) const
{
    if(p_iele == NULL) {
        INVALID_ARGUMENT("p_iele is NULL");
    }
    if( NumOfCVs == 0 ){
        RUNTIME_ERROR("no CVs allocated");
    }

    CXMLElement* p_ele = p_iele->GetFirstChildElement("CVS");
    if(p_ele == NULL) {
        ES_ERROR("unable to get CVS element");
        return(false);
    }

    CXMLElement*   p_cel = p_ele->GetFirstChildElement("CV");
    int            ccount = 0;

    while(p_cel != NULL) {
        if( ccount >= NumOfCVs ) {
            ES_ERROR("more CV elements than NumOfCVs");
            return(false);
        }
        if( CVs[ccount]->CheckInfo(p_cel) == false ) {
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

bool CPMFAccumulator::CheckCVSInfo(CPMFAccumulatorPtr p_accu) const
{
    if(p_accu == NULL) {
        INVALID_ARGUMENT("p_accu is NULL");
    }

    if(p_accu->NumOfCVs != NumOfCVs) {
        ES_ERROR("mismatch in the number of coordinates");
        return(false);
    }

    for(int i=0; i < NumOfCVs; i++) {
        if(CVs[i]->CheckInfo(p_accu->CVs[i]) == false) {
            CSmallString error;
            error << "mismatch in cv: " << i+1;
            ES_ERROR(error);
            return(false);
        }
    }

    return(true);
}

//------------------------------------------------------------------------------

void CPMFAccumulator::SaveCVSInfo(CXMLElement* p_tele) const
{
    if(p_tele == NULL) {
        INVALID_ARGUMENT("p_tele is NULL");
    }

    CXMLElement* p_ele = p_tele->CreateChildElement("CVS");

    for(int i=0; i < NumOfCVs; i++) {
        CXMLElement* p_iele = p_ele->CreateChildElement("CV");
        CVs[i]->SaveInfo(p_iele);
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CPMFAccumulator::PrintInfo(std::ostream& vout)
{
    PrintAccuInfo(vout);
    PrintCVSInfo(vout);
}

//------------------------------------------------------------------------------

void CPMFAccumulator::PrintInfo(FILE* p_fout)
{
    PrintAccuInfo(p_fout);
    PrintCVSInfo(p_fout);
}

//------------------------------------------------------------------------------

void CPMFAccumulator::PrintAccuInfo(std::ostream& vout)
{
    vout  << endl;
    vout << "=== General Info ===============================================================" << endl;
    vout << endl;
    vout << "Method:        " << Method << endl;
    vout << "Driver:        " << Driver << endl;
    vout << "Version:       " << Version << endl;
    vout << "Temperature:   " << fixed << setprecision(2) << Temperature*TemperatureFConv << " " << TemperatureUnit << endl;
    vout << "Energy unit:   " << EnergyUnit << endl;
    vout << "Data sections: " << DataBlocks.size() << endl;
}

//------------------------------------------------------------------------------

void CPMFAccumulator::PrintCVSInfo(std::ostream& vout)
{
    vout << endl;
    vout << "=== Collective Variables =======================================================" << endl;
    vout << endl;
    vout << "ID P Type       Unit  Name                       Min value   Max value   NBins  " << endl;
    vout << "-- - ---------- ----- -------------------------- ----------- ----------- -------" << endl;

    for(int i=0; i < NumOfCVs; i++) {
        CVs[i]->PrintInfo(vout);
    }
}

//------------------------------------------------------------------------------

void CPMFAccumulator::PrintAccuInfo(FILE* p_fout)
{
    fprintf(p_fout,"#\n");
    fprintf(p_fout,"# === General Info ===============================================================\n");
    fprintf(p_fout,"#\n");
    fprintf(p_fout,"# Method:        %s\n", (const char*)Method);
    fprintf(p_fout,"# Driver:        %s\n", (const char*)Driver);
    fprintf(p_fout,"# Version:       %s\n", (const char*)Version);
    fprintf(p_fout,"# Temperature:   %5.2lf %s\n", Temperature*TemperatureFConv, (const char*)TemperatureUnit);
    fprintf(p_fout,"# Energy unit:   %s\n", (const char*)EnergyUnit);
    fprintf(p_fout,"# Data sections: %d\n", (int)DataBlocks.size());
}

//------------------------------------------------------------------------------

void CPMFAccumulator::PrintCVSInfo(FILE* p_fout)
{
    fprintf(p_fout,"#\n");
    fprintf(p_fout,"# == Collective Variables =======================================================\n");
    fprintf(p_fout,"#\n");
    fprintf(p_fout,"# ID P Type       Unit  Name                       Min value   Max value   NBins  \n");
    fprintf(p_fout,"# -- - ---------- ----- -------------------------- ----------- ----------- -------\n");

    for(int i=0; i < NumOfCVs; i++) {
        CVs[i]->PrintInfo(p_fout);
    }
}

//------------------------------------------------------------------------------

 void CPMFAccumulator::ListSections(std::ostream& vout)
 {
    vout << endl;
    vout << "Number of sections : " << right << setw(10) << DataBlocks.size() << endl;
    vout << "Number of CVs      : " << right << setw(10) << NumOfCVs << endl;
    vout << "Number of bins     : " << right << setw(10) << NumOfBins << endl;
    vout << endl;
    vout << "=== Data Sections ==============================================================" << endl;
    vout << endl;
    vout << "Name                 Size         T M Op NBins      NCVs      " << endl;
    vout << "-------------------- ------------ - - -- ---------- ----------" << endl;

    std::map<CSmallString,CPMFAccuDataPtr>::iterator it = DataBlocks.begin();
    std::map<CSmallString,CPMFAccuDataPtr>::iterator ie = DataBlocks.end();

    while( it != ie ){
        CPMFAccuDataPtr db = it->second;
        vout << setw(20) << left << db->GetName() << " " << right;
        vout << setw(12) << db->GetLength() << " ";
        vout << setw(1)  << db->GetType() << " ";
        vout << setw(1)  << db->GetMode() << " ";
        vout << setw(2)  << db->GetOp() << " ";
        vout << setw(10)  << db->GetNumOfBins() << " ";
        vout << setw(10)  << db->GetNumOfCVs();
        vout << endl;
        it++;
    }
 }

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CPMFAccumulator::HasSectionData(const CSmallString& name) const
{
    std::map<CSmallString,CPMFAccuDataPtr>::const_iterator isec = DataBlocks.find(name);
    if( isec == DataBlocks.end() ) {
        return(false);
    }
    return(true);
}

//------------------------------------------------------------------------------

void CPMFAccumulator::DeleteSectionData(const CSmallString& name)
{
    DataBlocks.erase(name);
}

//------------------------------------------------------------------------------

void CPMFAccumulator::CreateSectionData(const CSmallString& name,const CSmallString& op,const CSmallString& type,const CSmallString& mode,int len)
{
    CPMFAccuDataPtr data = CPMFAccuDataPtr(new CPMFAccuData(NumOfBins,NumOfCVs));
    data->Name  = name;
    data->Op    = op;
    data->Type  = type;
    data->Mode  = mode;

    data->Size  = 0;
    if( data->Mode == "B" ){
        data->Size  = NumOfBins;
    } else if ( data->Mode == "C" ) {
        data->Size  = NumOfCVs;
    } else if ( data->Mode == "M" ) {
        data->Size  = NumOfCVs * NumOfBins;
    } else {
        data->Size = len;
    }
    if( data->Size > 0 ) {
        data->Data.CreateVector(data->Size);
        data->Data.SetZero();
    }

    DataBlocks[name] = data;
}

//------------------------------------------------------------------------------

CPMFAccuDataPtr CPMFAccumulator::GetSectionData(const CSmallString& name) const
{
    std::map<CSmallString,CPMFAccuDataPtr>::const_iterator isec = DataBlocks.find(name);
    if( isec == DataBlocks.end() ) {
        CSmallString error;
        error << "unable to find '" << name << "' data section";
        RUNTIME_ERROR(error);
    }
    const CPMFAccuDataPtr sec = isec->second;
    return(sec);
}

//------------------------------------------------------------------------------

double CPMFAccumulator::GetData(const CSmallString& name, int ibin, int cv) const
{
    std::map<CSmallString,CPMFAccuDataPtr>::const_iterator isec = DataBlocks.find(name);
    if( isec == DataBlocks.end() ) {
        CSmallString error;
        error << "unable to find '" << name << "' data section";
        RUNTIME_ERROR(error);
    }
    const CPMFAccuDataPtr sec = isec->second;
    return( sec->GetData(ibin,cv) );
}

//------------------------------------------------------------------------------

void CPMFAccumulator::SetData(const CSmallString& name, int ibin, double value)
{
    std::map<CSmallString,CPMFAccuDataPtr>::iterator isec = DataBlocks.find(name);
    if( isec == DataBlocks.end() ) {
        CSmallString error;
        error << "unable to find '" << name << "' data section";
        RUNTIME_ERROR(error);
    }
    CPMFAccuDataPtr sec = isec->second;
    return( sec->SetData(ibin,value) );
}

//------------------------------------------------------------------------------

void CPMFAccumulator::SetData(const CSmallString& name, int ibin, int cv, double value)
{
    std::map<CSmallString,CPMFAccuDataPtr>::iterator isec = DataBlocks.find(name);
    if( isec == DataBlocks.end() ) {
        CSmallString error;
        error << "unable to find '" << name << "' data section";
        RUNTIME_ERROR(error);
    }
    CPMFAccuDataPtr sec = isec->second;
    return( sec->SetData(ibin,cv,value) );
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CPMFAccumulator::map(int item,int bin) const
{
    return(item + bin*NumOfCVs);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================




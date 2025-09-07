// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2025 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
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

#include <BaseProxy.hpp>
#include <algorithm>
#include <boost/algorithm/string/join.hpp>

//------------------------------------------------------------------------------

using namespace std;
using namespace boost;


//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CProxyRealmDescr::CProxyRealmDescr(void)
{
    RealmID = -1;
}

//------------------------------------------------------------------------------

bool CProxyRealmDescr::Compare(const CProxyRealmDescr& left, const CProxyRealmDescr& right)
{
    if( left.Method == right.Method ) return( left.Description < right.Description);
    return(left.Method < right.Method);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CBaseProxy::CBaseProxy(void)
{
    RealmID = -1;
}

//------------------------------------------------------------------------------

CBaseProxy::~CBaseProxy(void)
{

}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CBaseProxy::RegisterRealm(int realmid,const CSmallString& realm,const CSmallString& method,const CSmallString& descr)
{
    CProxyRealmDescr rdesc;
    rdesc.RealmID = realmid;
    rdesc.Realm = realm;
    rdesc.Method = method;
    rdesc.Description = descr;

    if( SupportedRealms.count(realm) != 0 ) {
        CSmallString error;
        error << "realm: '" << realm << "/" << method << "/" << descr << "' already registered";
        RUNTIME_ERROR("error");
    }

    SupportedRealms[realm] = rdesc;
    Requires.insert(method);
}

//------------------------------------------------------------------------------

void CBaseProxy::Init(CPMFAccumulatorPtr accu)
{
    if( std::find(Requires.begin(),Requires.end(),string(accu->GetMethod())) == Requires.end() ) {
        CSmallString error;
   //     error << "PMF accumulator '" << accu->GetMethod() << "' is inconsistent with EnergyDerProxy requirements '" << join(Requires,",") << "'";
        RUNTIME_ERROR(error)
    }
    Accu = accu;
}

//------------------------------------------------------------------------------

bool CBaseProxy::IsCompatible(CPMFAccumulatorPtr accu)
{
    return( std::find(Requires.begin(), Requires.end(), std::string(accu->GetMethod())) != Requires.end());
}

//------------------------------------------------------------------------------

bool CBaseProxy::SetRealm(const CSmallString& realm)
{
    if( SupportedRealms.count(realm) == 0 ) return(false);
    RealmID = SupportedRealms[realm].RealmID;
    return(true);
}

//------------------------------------------------------------------------------

void CBaseProxy::SetRealm(int realmid)
{
    std::map<CSmallString,CProxyRealmDescr>::iterator rit = SupportedRealms.begin();
    std::map<CSmallString,CProxyRealmDescr>::iterator rie = SupportedRealms.end();

    CSmallString realm;

    while(rit != rie){
        if( rit->second.RealmID == realmid ){
            RealmID = rit->second.RealmID;
            return;
        }
        rit++;
    }

    RUNTIME_ERROR("unsupported realmID");
}

//------------------------------------------------------------------------------

void CBaseProxy::EnumerateRealms(std::list<CProxyRealmDescr>& dlist)
{
    std::set<CSmallString>::iterator  mit = Requires.begin();
    std::set<CSmallString>::iterator  mie = Requires.end();

    while(mit != mie){
        CSmallString method = *mit;

        std::map<CSmallString,CProxyRealmDescr>::iterator rit = SupportedRealms.begin();
        std::map<CSmallString,CProxyRealmDescr>::iterator rie = SupportedRealms.end();

        while(rit != rie){
            dlist.push_back(rit->second);
            rit++;
        }
        mit++;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CPMFAccumulatorPtr CBaseProxy::GetAccu(void)
{
    return(Accu);
}

//------------------------------------------------------------------------------

CSmallString CBaseProxy::GetRealm(void)
{
    std::map<CSmallString,CProxyRealmDescr>::iterator rit = SupportedRealms.begin();
    std::map<CSmallString,CProxyRealmDescr>::iterator rie = SupportedRealms.end();

    CSmallString realm;

    while(rit != rie){
        if( rit->second.RealmID == RealmID ){
            realm = rit->second.Realm;
        }
        rit++;
    }

    return(realm);
}

//------------------------------------------------------------------------------

CSmallString CBaseProxy::GetMethods(void)
{
    std::set<CSmallString>::iterator it = Requires.begin();
    std::set<CSmallString>::iterator ie = Requires.end();

    CSmallString sm;
    while( it != ie ){
        if( it != Requires.begin() ) sm << ",";
        sm << *it;
        it++;
    }
    return(sm);
}

//------------------------------------------------------------------------------

CSmallString CBaseProxy::GetDescription(void)
{
    std::map<CSmallString,CProxyRealmDescr>::iterator rit = SupportedRealms.begin();
    std::map<CSmallString,CProxyRealmDescr>::iterator rie = SupportedRealms.end();

    CSmallString descr;

    while(rit != rie){
        if( rit->second.RealmID == RealmID ){
            descr = rit->second.Description;
        }
        rit++;
    }

    return(descr);
}

//------------------------------------------------------------------------------

int CBaseProxy::GetNumOfCVs(void) const
{
    if( Accu == NULL ) return(0);
    return(Accu->GetNumOfCVs());
}

//------------------------------------------------------------------------------

int CBaseProxy::GetNumOfBins(void) const
{
    if( Accu == NULL ) return(0);
    return(Accu->GetNumOfBins());
}

//------------------------------------------------------------------------------

int CBaseProxy::GetNumOfSamples(int ibin) const
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }
    return(Accu->GetData("NSAMPLES",ibin));
}

//------------------------------------------------------------------------------

void CBaseProxy::SetNumOfSamples(int ibin,int nsamples)
{
    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }
    Accu->SetData("NSAMPLES",ibin,nsamples);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

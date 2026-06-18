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

CSmallString CBaseProxy::GetFullDescription(void)
{
    CSmallString    desc;
    desc << GetRealm() << " / [" << GetMethods() << "] / " << GetDescription();
    return(desc);
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

void CBaseProxy::GetMeanValue(const CSmallString& section_name,
                              double& mean, double& sd, double& sem, bool value_only, int ibin, int icv) const
{
    mean = 0.0;
    sd = 0.0;
    sem = 0.0;

    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }

    CSmallString m_sname;
    m_sname << "M" << section_name;
    CPMFAccuDataPtr msec = Accu->GetSectionData(m_sname);
    if( ! msec ){
        CSmallString error;
        error << "unable to get (" << m_sname << ") section for mean";
        RUNTIME_ERROR(error);
    }
    if( msec->GetOp() != "WA" ){
        CSmallString error;
        error << "section (" << m_sname << ") does not have the correct (" << msec->GetOp() << ") operation assigned, expecting (WA)";
        RUNTIME_ERROR(error); 
    }

    // https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
    // online algorithm for the mean

    // mean from online summation algorithm 

    if( icv >= 0 ) {
        if( msec->GetMode() == "M" ){
            mean = msec->GetData(ibin,icv);
        } else {
            CSmallString error;
            error << "unsupported M section (" << msec->GetName() << ") data mode (" << msec->GetMode()
                << "), expected (M)";
            RUNTIME_ERROR(error); 
        }
    } else {
        if( msec->GetMode() == "B" ){
            mean = msec->GetData(ibin);
        } else {
            CSmallString error;
            error << "unsupported M section (" << msec->GetName() << ") data mode (" << msec->GetMode()
                << "), expected (B)";
            RUNTIME_ERROR(error); 
        }
    }
    
    if( value_only ) return;

    CSmallString m2_sname;
    m2_sname << "M2" << section_name;

    CPMFAccuDataPtr m2sec = Accu->GetSectionData(m2_sname);
    if( ! m2sec ){
        // no error data available
        return;
    }
    if( m2sec->GetOp() != "M2" ){
        CSmallString error;
        error << "section (" << m2_sname << ") does not have the correct (" << m2sec->GetOp() << ") operation assigned, expecting (M2)";
        RUNTIME_ERROR(error); 
    }

    CSmallString nsam_sname = msec->GetMSName();
    CPMFAccuDataPtr nsamsec = Accu->GetSectionData(nsam_sname);
    if( ! nsamsec ){
        CSmallString error;
        error << "unable to get (" << nsam_sname << ") section for sample standard deviation (sd)";
        RUNTIME_ERROR(error);
    }

    if( msec->GetMSName() != m2sec->GetMSName() ){
        CSmallString error;
        error << "sections (" << m_sname << "," << m2_sname << ") do not have the same MS record (" << msec->GetMSName() << "," << m2sec->GetMSName() << ")";
        RUNTIME_ERROR(error); 
    }

    if( m2sec->GetMXName() != m_sname ){
        CSmallString error;
        error << "section (" << m2_sname << ") does not have the correct MX record (" << m2sec->GetMXName() << "), expected (" << m_sname << ")";
        RUNTIME_ERROR(error); 
    }

    // https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
    // online algorithm for the variance

    double m2 = 0.0;

    if( icv >= 0 ){
        if( m2sec->GetMode() == "M" ){
            m2 = m2sec->GetData(ibin,icv);
        } else {
            CSmallString error;
            error << "unsupported M2 section (" << m2sec->GetName() << ") data mode (" << m2sec->GetMode()
                << "), expected (M)";
            RUNTIME_ERROR(error); 
        }
    } else {
        if( m2sec->GetMode() == "B" ){
            m2 = m2sec->GetData(ibin);
        } else {
            CSmallString error;
            error << "unsupported M2 section (" << m2sec->GetName() << ") data mode (" << m2sec->GetMode()
                << "), expected (B)";
            RUNTIME_ERROR(error); 
        }
    }

    double nsamples = nsamsec->GetData(ibin);

    if( nsamples <= 1 ){
        // no error data available
        return;
    }

    // unbiased sample standard deviation
    sd = sqrt( m2 / (nsamples - 1.0) );

    // standard error of the mean
    sem = sd / sqrt(nsamples);
}

//------------------------------------------------------------------------------

void CBaseProxy::GetWMeanValue(const CSmallString& section_name,const CSmallString& weight_section_name,
                            double& mean, double& sd,double& sem, bool value_only, int ibin, int icv) const
{
    mean = 0.0;
    sd = 0.0;
    sem = 0.0;

    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }

    CSmallString m_sname;
    m_sname << "M" << section_name;

    CPMFAccuDataPtr msec = Accu->GetSectionData(m_sname);
    if( ! msec ){
        CSmallString error;
        error << "unable to get (" << m_sname << ") section for mean";
        RUNTIME_ERROR(error);
    }
    if( msec->GetOp() != "WA" ){
        CSmallString error;
        error << "section (" << m_sname << ") does not have the correct (" << msec->GetOp() << ") operation assigned, expecting (WA)";
        RUNTIME_ERROR(error); 
    }

    CSmallString wsum_sname;
    wsum_sname << weight_section_name << "SUM";

    if( msec->GetMSName() != wsum_sname ){
        CSmallString error;
        error << "section (" << m_sname << ") does not have the correct MS record (" << msec->GetMSName() << "), expected (" << wsum_sname << ")";
        RUNTIME_ERROR(error); 
    }

    // https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
    // online algorithm for the mean

    // mean from online summation algorithm 

    if( icv >= 0 ) {
        if( msec->GetMode() == "M" ){
            mean = msec->GetData(ibin,icv);
        } else {
            CSmallString error;
            error << "unsupported M section (" << msec->GetName() << ") data mode (" << msec->GetMode()
                << "), expected (M)";
            RUNTIME_ERROR(error); 
        }
    } else {
        if( msec->GetMode() == "B" ){
            mean = msec->GetData(ibin);
        } else {
            CSmallString error;
            error << "unsupported M section (" << msec->GetName() << ") data mode (" << msec->GetMode()
                << "), expected (B)";
            RUNTIME_ERROR(error); 
        }
    }

    if( value_only ) return;

    CSmallString m2_sname;
    m2_sname << "M2" << section_name;

    CPMFAccuDataPtr m2sec = Accu->GetSectionData(m2_sname);
    if( ! m2sec ){
        // no error data available
        return;
    }
    if( m2sec->GetOp() != "M2" ){
        CSmallString error;
        error << "section (" << m2_sname << ") does not have the correct (" << m2sec->GetOp() << ") operation assigned, expecting (M2)";
        RUNTIME_ERROR(error); 
    }

    if( m2sec->GetMSName() != wsum_sname ){
        CSmallString error;
        error << "section (" << m2_sname << ") does not have the correct MS record (" << m2sec->GetMSName() << "), expected (" << wsum_sname << ")";
        RUNTIME_ERROR(error); 
    }

    if( m2sec->GetMXName() != m_sname ){
        CSmallString error;
        error << "section (" << m2_sname << ") does not have the correct MX record (" << m2sec->GetMXName() << "), expected (" << m_sname << ")";
        RUNTIME_ERROR(error); 
    }

    CPMFAccuDataPtr wsum_sec = Accu->GetSectionData(wsum_sname);
    if( ! wsum_sec ){
        CSmallString error;
        error << "unable to get (" << wsum_sname << ") section for sample standard deviation (sd)";
        RUNTIME_ERROR(error);
    }

    CSmallString w2sum_sname;
    w2sum_sname << weight_section_name << "SUM2";

    CPMFAccuDataPtr w2sum_sec = Accu->GetSectionData(w2sum_sname);
    if( ! w2sum_sec ){
        CSmallString error;
        error << "unable to get (" << w2sum_sname << ") section for sample standard deviation (sd)";
        RUNTIME_ERROR(error);
    }

    // https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
    // online algorithm for the variance

    // https://seismo.berkeley.edu/~kirchner/Toolkits/Toolkit_12.pdf  

    double m2 = 0.0;

    if( icv >= 0 ){
        if( m2sec->GetMode() == "M" ){
            m2 = m2sec->GetData(ibin,icv);
        } else {
            CSmallString error;
            error << "unsupported M2 section (" << m2sec->GetName() << ") data mode (" << m2sec->GetMode()
                << "), expected (M)";
            RUNTIME_ERROR(error); 
        }
    } else {
        if( m2sec->GetMode() == "B" ){
            m2 = m2sec->GetData(ibin);
        } else {
            CSmallString error;
            error << "unsupported M2 section (" << m2sec->GetName() << ") data mode (" << m2sec->GetMode()
                << "), expected (B)";
            RUNTIME_ERROR(error); 
        }
    }

    double wsum     = wsum_sec->GetData(ibin);
    double w2sum    = w2sum_sec->GetData(ibin);

    if( w2sum <= 0.0 ){
        // no error data available
        return;
    }

    // number of effective measurements
    // Kish effective sample size
    double neff = (wsum * wsum) / w2sum;

    if( neff <= 1.0 ){
        // no error data available
        return;
    }

    // unbiased weighted sample standard deviation
    sd  = sqrt( (m2 / wsum) * neff / (neff - 1.0) );

    // standard error of the mean
    sem = sd / sqrt(neff);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CBaseProxy::GetCovarianceValue(const CSmallString& section_name,
                    double& cval, double& sd, double& sem, bool value_only, int ibin, int icv) const
{
    cval = 0.0;
    sd = 0.0;
    sem = 0.0;

    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }

    CPMFAccuDataPtr csec = Accu->GetSectionData(section_name);
    if( ! csec ){
        CSmallString error;
        error << "unable to get (" << section_name << ") section for mean";
        RUNTIME_ERROR(error);
    }
    if( csec->GetOp() != "CO" ){
        CSmallString error;
        error << "section (" << section_name << ") does not have the correct (" << csec->GetOp() << ") operation assigned, expecting (CO)";
        RUNTIME_ERROR(error); 
    }

    CSmallString nsam_sname = csec->GetMSName();
    CPMFAccuDataPtr nsamsec = Accu->GetSectionData(nsam_sname);
    if( ! nsamsec ){
        CSmallString error;
        error << "unable to get MS (" << nsam_sname << ") section for sample standard deviation (sd)";
        RUNTIME_ERROR(error);
    }

    double nsamples = nsamsec->GetData(ibin);

    if( nsamples <= 1 ){
        // no error data available
        return;
    }

// covariance - for PMF, the covariance will be biased, 1/nsamples

    if( icv >= 0 ) {
        if( csec->GetMode() == "M" ){
            cval = csec->GetData(ibin,icv) / nsamples;
        } else {
            CSmallString error;
            error << "unsupported section (" << csec->GetName() << ") data mode (" << csec->GetMode()
                << "), expected (M)";
            RUNTIME_ERROR(error); 
        }
    } else {
        if( csec->GetMode() == "B" ){
            cval = csec->GetData(ibin) / nsamples;
        } else {
            CSmallString error;
            error << "unsupported section (" << csec->GetName() << ") data mode (" << csec->GetMode()
                << "), expected (B)";
            RUNTIME_ERROR(error); 
        }
    }

    if( value_only ) return;

// std. deviation
    CSmallString mx_sname = csec->GetMXName();
    CSmallString m2x_sname;
    m2x_sname << "M2" << mx_sname.GetSubString(1,-1);

    CPMFAccuDataPtr m2x = Accu->GetSectionData(m2x_sname);
    if( ! m2x ){
        CSmallString error;
        error << "unable to get M2X (" << m2x_sname << ") section for sample standard deviation (sd)";
        RUNTIME_ERROR(error);
    }
    if( m2x->GetMSName() != csec->GetMSName() ){
        CSmallString error;
        error << "sections (" << section_name << "," << m2x_sname << ") do not have the same MS record (" << csec->GetMSName() << "," << m2x->GetMSName() << ")";
        RUNTIME_ERROR(error); 
    }


    bool m_mode = false;

    double xvar = 0.0;

    // m2x and m2y cannot be both in M mode, valid modes are: 1) B,B for (icv < 0) and 2) B,M or M,B for icv >= 0
    if( (m2x->GetMode() == "M") && (icv >= 0) ){  // if icv >= 0, it could be both M or B
        xvar = m2x->GetData(ibin,icv);
        m_mode = true;
    } else if( m2x->GetMode() == "B" ){
        xvar = m2x->GetData(ibin);
    } else {
        CSmallString error;
        error << "unsupported M2X data mode (" << m2x->GetMode()
              << "), expected (M) or (B)";
        RUNTIME_ERROR(error);
    }

    xvar = xvar / nsamples;

    CSmallString my_sname = csec->GetMYName();
    CSmallString m2y_sname;
    m2y_sname << "M2" << my_sname.GetSubString(1,-1);

    CPMFAccuDataPtr m2y = Accu->GetSectionData(m2y_sname);
    if( ! m2y ){
        CSmallString error;
        error << "unable to get M2Y (" << m2y_sname << ") section for sample standard deviation (sd)";
        RUNTIME_ERROR(error);
    }
    if( m2y->GetMSName() != csec->GetMSName() ){
        CSmallString error;
        error << "sections (" << section_name << "," << m2y_sname << ") do not have the same MS record (" << csec->GetMSName() << "," << m2y->GetMSName() << ")";
        RUNTIME_ERROR(error); 
    }

    double yvar = 0.0;

    // m2x and m2y cannot be both in M mode, valid modes are: 1) B,B for (icv < 0) and 2) B,M or M,B for icv >= 0
    if( (m2y->GetMode() == "M") && (icv >= 0) && (m_mode == false)  ){ // if icv >= 0 and m_mode == false, it could be both M or B
        yvar = m2y->GetData(ibin,icv);
        m_mode = true;
    } else if( m2y->GetMode() == "B" ){
        yvar = m2y->GetData(ibin);
    } else {
        CSmallString error;
        error << "unsupported M2Y data mode (" << m2y->GetMode()
              << "), expected (M) or (B)";
        RUNTIME_ERROR(error);
    }

    if( (icv >= 0) && (m_mode == false) ){
        CSmallString error;
        error << "at least one M2 section (" << m2x->GetMode() << ", " << m2y->GetMode()
              << ") must be in the (M) mode";
        RUNTIME_ERROR(error);
    }

    yvar = yvar / nsamples;

    sd = sqrt(xvar * yvar + cval * cval);

    // stadart error of the mean
    sem = sd / sqrt(nsamples);
}

//------------------------------------------------------------------------------

void CBaseProxy::GetWCovarianceValue(const CSmallString& section_name,const CSmallString& weight_section_name,
                    double& cval, double& sd, double& sem, bool value_only, int ibin, int icv) const
{
    cval = 0.0;
    sd   = 0.0;
    sem  = 0.0;

    if( Accu == NULL ){
        RUNTIME_ERROR("Accu is NULL");
    }

    CPMFAccuDataPtr csec = Accu->GetSectionData(section_name);
    if( ! csec ){
        CSmallString error;
        error << "unable to get (" << section_name << ") section for covariance";
        RUNTIME_ERROR(error);
    }

    if( csec->GetOp() != "CO" ){
        CSmallString error;
        error << "section (" << section_name << ") is not in the correct (" << csec->GetOp()
              << ") operation, expecting (CO)";
        RUNTIME_ERROR(error);
    }

    CSmallString wsum_sname;
    wsum_sname << weight_section_name << "SUM";

    if( csec->GetMSName() != wsum_sname ){
        CSmallString error;
        error << "section (" << section_name << ") does not have the correct MS record ("
              << csec->GetMSName() << "), expected (" << wsum_sname << ")";
        RUNTIME_ERROR(error);
    }

    CPMFAccuDataPtr wsum_sec = Accu->GetSectionData(wsum_sname);
    if( ! wsum_sec ){
        CSmallString error;
        error << "unable to get weighted sum (" << wsum_sname << ") section";
        RUNTIME_ERROR(error);
    }

    double wsum = wsum_sec->GetData(ibin);

    if( wsum <= 0.0 ){
        // no data available
        return;
    }

// weighted covariance - for PMF, the covariance will be biased, 1/wsum

    if( icv >= 0 ) {
        if( csec->GetMode() == "M" ){
            cval = csec->GetData(ibin,icv) / wsum;
        } else {
            CSmallString error;
            error << "unsupported section (" << csec->GetName() << ") data mode (" << csec->GetMode()
                << "), expected (M)";
            RUNTIME_ERROR(error); 
        }
    } else {
        if( csec->GetMode() == "B" ){
            cval = csec->GetData(ibin) / wsum;
        } else {
            CSmallString error;
            error << "unsupported section (" << csec->GetName() << ") data mode (" << csec->GetMode()
                << "), expected (B)";
            RUNTIME_ERROR(error); 
        }
    }

    if( value_only ) return;

// effective number of samples

    CSmallString w2sum_sname;
    w2sum_sname << weight_section_name << "SUM2";

    CPMFAccuDataPtr w2sum_sec = Accu->GetSectionData(w2sum_sname);
    if( ! w2sum_sec ){
        CSmallString error;
        error << "unable to get weighted squared sum (" << w2sum_sname << ") section";
        RUNTIME_ERROR(error);
    }

    double w2sum = w2sum_sec->GetData(ibin);

    if( w2sum <= 0.0 ){
        // no error data available
        return;
    }

    // Kish effective sample size
    double neff = (wsum * wsum) / w2sum;

    if( neff <= 1.0 ){
        // no error data available
        return;
    }

// variance of X

    CSmallString mx_sname = csec->GetMXName();
    CSmallString m2x_sname;
    m2x_sname << "M2" << mx_sname.GetSubString(1,-1);

    CPMFAccuDataPtr m2x = Accu->GetSectionData(m2x_sname);
    if( ! m2x ){
        CSmallString error;
        error << "unable to get M2X (" << m2x_sname << ") section for covariance standard error";
        RUNTIME_ERROR(error);
    }

    if( m2x->GetOp() != "M2" ){
        CSmallString error;
        error << "section (" << m2x_sname << ") is not in the correct (" << m2x->GetOp()
              << ") operation, expecting (M2)";
        RUNTIME_ERROR(error);
    }

    if( m2x->GetMSName() != wsum_sname ){
        CSmallString error;
        error << "section (" << m2x_sname << ") does not have the correct MS record ("
              << m2x->GetMSName() << "), expected (" << wsum_sname << ")";
        RUNTIME_ERROR(error);
    }

    if( m2x->GetMXName() != mx_sname ){
        CSmallString error;
        error << "section (" << m2x_sname << ") does not have the correct MX record ("
              << m2x->GetMXName() << "), expected (" << mx_sname << ")";
        RUNTIME_ERROR(error);
    }

    double xvar = 0.0;

    bool m_mode = false;

    // m2x and m2y cannot be both in M mode, valid modes are: 1) B,B for (icv < 0) and 2) B,M or M,B for icv >= 0
    if( (m2x->GetMode() == "M") && (icv >= 0) ){  // if icv >= 0, it could be both M or B
        xvar = m2x->GetData(ibin,icv);
        m_mode = true;
    } else if( m2x->GetMode() == "B" ){
        xvar = m2x->GetData(ibin);
    } else {
        CSmallString error;
        error << "unsupported M2X data mode (" << m2x->GetMode()
              << "), expected (M) or (B)";
        RUNTIME_ERROR(error);
    }

    xvar /= wsum;

// variance of Y

    CSmallString my_sname = csec->GetMYName();
    CSmallString m2y_sname;
    m2y_sname << "M2" << my_sname.GetSubString(1,-1);

    CPMFAccuDataPtr m2y = Accu->GetSectionData(m2y_sname);
    if( ! m2y ){
        CSmallString error;
        error << "unable to get M2Y (" << m2y_sname << ") section for covariance standard error";
        RUNTIME_ERROR(error);
    }

    if( m2y->GetOp() != "M2" ){
        CSmallString error;
        error << "section (" << m2y_sname << ") is not in the correct (" << m2y->GetOp()
              << ") operation, expecting (M2)";
        RUNTIME_ERROR(error);
    }

    if( m2y->GetMSName() != wsum_sname ){
        CSmallString error;
        error << "section (" << m2y_sname << ") does not have the correct MS record ("
              << m2y->GetMSName() << "), expected (" << wsum_sname << ")";
        RUNTIME_ERROR(error);
    }

    if( m2y->GetMXName() != my_sname ){
        CSmallString error;
        error << "section (" << m2y_sname << ") does not have the correct MX record ("
              << m2y->GetMXName() << "), expected (" << my_sname << ")";
        RUNTIME_ERROR(error);
    }

    double yvar = 0.0;

    // m2x and m2y cannot be both in M mode, valid modes are: 1) B,B for (icv < 0) and 2) B,M or M,B for icv >= 0
    if( (m2y->GetMode() == "M") && (icv >= 0) && (m_mode == false)  ){ // if icv >= 0 and m_mode == false, it could be both M or B
        yvar = m2y->GetData(ibin,icv);
        m_mode = true;
    } else if( m2y->GetMode() == "B" ){
        yvar = m2y->GetData(ibin);
    } else {
        CSmallString error;
        error << "unsupported M2Y data mode (" << m2y->GetMode()
              << "), expected (M) or (B)";
        RUNTIME_ERROR(error);
    }

    if( (icv >= 0) && (m_mode == false) ){
        CSmallString error;
        error << "at least one M2 section (" << m2x->GetMode() << ", " << m2y->GetMode()
              << ") must be in the (M) mode";
        RUNTIME_ERROR(error);
    }

    yvar /= wsum;

// standard deviation and SEM of the covariance estimator

    sd = sqrt(xvar * yvar + cval * cval);

    sem = sd / sqrt(neff);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================
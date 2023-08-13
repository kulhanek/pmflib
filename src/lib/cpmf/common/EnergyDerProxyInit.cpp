// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2023 Petr Kulhanek, kulhanek@chemi.muni.cz
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

#include <EnergyDerProxyInit.hpp>
#include <ABFProxy_dG.hpp>
#include <ABFProxy_dH.hpp>
#include <ABFProxy_mTdS.hpp>
#include <CSTProxy_dG.hpp>
#include <CSTProxy_mTdS.hpp>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CEnergyDerProxyPtr CEnergyDerProxyInit::InitProxy(const CSmallString& realm,CPMFAccumulatorPtr& accu)
{
    CEnergyDerProxyPtr lproxy;

    if( realm == "dG" ){
        if( CABFProxy_dG::IsCompatible(accu) ){
            lproxy    = CABFProxy_dG_Ptr(new CABFProxy_dG);
        } else if (CCSTProxy_dG::IsCompatible(accu) ) {
            lproxy    = CCSTProxy_dG_Ptr(new CCSTProxy_dG);
        } else {
            CSmallString error;
            error << "incompatible method: " << accu->GetMethod() << " with requested realm: " <<  realm;
            RUNTIME_ERROR(error);
        }
// -----------------------------------------------
    } else if ( realm == "ICFP" ) {
        if( CABFProxy_dH::IsCompatible(accu) ){
            lproxy    = CABFProxy_dH_Ptr(new CABFProxy_dH);
        } else {
            CSmallString error;
            error << "incompatible method: " << accu->GetMethod() << " with requested realm: " <<  realm;
            RUNTIME_ERROR(error);
        }
// -----------------------------------------------
    } else if ( (realm == "-TdS") || (realm == "mTdS") ) {
        if( CABFProxy_mTdS::IsCompatible(accu) ){
            lproxy    = CABFProxy_mTdS_Ptr(new CABFProxy_mTdS);
        } else if (CCSTProxy_mTdS::IsCompatible(accu) ) {
            lproxy    = CCSTProxy_mTdS_Ptr(new CCSTProxy_mTdS);
        } else {
            CSmallString error;
            error << "incompatible method: " << accu->GetMethod() << " with requested realm: " <<  realm;
            RUNTIME_ERROR(error);
        }
// -----------------------------------------------
    } else if ( (realm == "-TdS_HP") || (realm == "mTdS_HP") ) {
        if( CABFProxy_mTdS::IsCompatible(accu) ){
            CABFProxy_mTdS_Ptr proxy = CABFProxy_mTdS_Ptr(new CABFProxy_mTdS);
            proxy->SetType(ABF_TdS_HP);
            lproxy = proxy;
        } else if (CCSTProxy_mTdS::IsCompatible(accu) ) {
            CCSTProxy_mTdS_Ptr proxy = CCSTProxy_mTdS_Ptr(new CCSTProxy_mTdS);
            proxy->SetType(CST_TdS_HP);
            lproxy = proxy;
        } else {
            CSmallString error;
            error << "incompatible method: " << accu->GetMethod() << " with requested realm: " <<  realm;
            RUNTIME_ERROR(error);
        }
// -----------------------------------------------
    } else if ( (realm == "-TdS_HR") || (realm == "mTdS_HR") ) {
        if( CABFProxy_mTdS::IsCompatible(accu) ){
            CABFProxy_mTdS_Ptr proxy = CABFProxy_mTdS_Ptr(new CABFProxy_mTdS);
            proxy->SetType(ABF_TdS_HR);
            lproxy = proxy;
        } else if (CCSTProxy_mTdS::IsCompatible(accu) ) {
            CCSTProxy_mTdS_Ptr proxy = CCSTProxy_mTdS_Ptr(new CCSTProxy_mTdS);
            proxy->SetType(CST_TdS_HR);
            lproxy = proxy;
        } else {
            CSmallString error;
            error << "incompatible method: " << accu->GetMethod() << " with requested realm: " <<  realm;
            RUNTIME_ERROR(error);
        }
// -----------------------------------------------
    } else if ( (realm == "-TdS_HK") || (realm == "mTdS_HK") ) {
        if( CABFProxy_mTdS::IsCompatible(accu) ){
            CABFProxy_mTdS_Ptr proxy = CABFProxy_mTdS_Ptr(new CABFProxy_mTdS);
            proxy->SetType(ABF_TdS_HK);
            lproxy = proxy;
        } else if (CCSTProxy_mTdS::IsCompatible(accu) ) {
            CCSTProxy_mTdS_Ptr proxy = CCSTProxy_mTdS_Ptr(new CCSTProxy_mTdS);
            proxy->SetType(CST_TdS_HK);
            lproxy = proxy;
        } else {
            CSmallString error;
            error << "incompatible method: " << accu->GetMethod() << " with requested realm: " <<  realm;
            RUNTIME_ERROR(error);
        }
// -----------------------------------------------
    } else if ( (realm == "-TdS_HV") || (realm == "mTdS_HV") ) {
        if( CABFProxy_mTdS::IsCompatible(accu) ){
            CABFProxy_mTdS_Ptr proxy = CABFProxy_mTdS_Ptr(new CABFProxy_mTdS);
            proxy->SetType(ABF_TdS_HV);
            lproxy = proxy;
        } else {
            CSmallString error;
            error << "incompatible method: " << accu->GetMethod() << " with requested realm: " <<  realm;
            RUNTIME_ERROR(error);
        }
// -----------------------------------------------
    } else if ( (realm == "-TdS_BP") || (realm == "mTdS_BP") ) {
        CABFProxy_mTdS_Ptr proxy    = CABFProxy_mTdS_Ptr(new CABFProxy_mTdS);
        proxy->SetType(ABF_TdS_BP);
        lproxy = proxy;
// -----------------------------------------------
    } else if ( (realm == "-TdS_BR") || (realm == "mTdS_BR") ) {
        CABFProxy_mTdS_Ptr proxy    = CABFProxy_mTdS_Ptr(new CABFProxy_mTdS);
        proxy->SetType(ABF_TdS_BR);
        lproxy = proxy;
// -----------------------------------------------
    } else if ( (realm == "-TdS_BK") || (realm == "mTdS_BK") ) {
        CABFProxy_mTdS_Ptr proxy    = CABFProxy_mTdS_Ptr(new CABFProxy_mTdS);
        proxy->SetType(ABF_TdS_BK);
        lproxy = proxy;
// -----------------------------------------------
    } else if ( (realm == "-TdS_BV") || (realm == "mTdS_BV") ) {
        CABFProxy_mTdS_Ptr proxy    = CABFProxy_mTdS_Ptr(new CABFProxy_mTdS);
        proxy->SetType(ABF_TdS_BV);
        lproxy = proxy;
// -----------------------------------------------
    } else {
        CSmallString error;
        error << "unsupported realm: " << realm;
        RUNTIME_ERROR(error);
    }

    return(lproxy);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================



//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

#ifndef ABPClientH
#define ABPClientH
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2020 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2007,2008 Petr Kulhanek, kulhanek@enzim.hu
//
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, Fifth Floor,
//    Boston, MA  02110-1301  USA
// =============================================================================

#include <PMFMainHeader.hpp>
#include <Client.hpp>
#include <vector>

//------------------------------------------------------------------------------

class CColVariable;
class CXMLElement;

//------------------------------------------------------------------------------

/*! \brief ABP client for fpmf library
 *
 */

class PMF_PACKAGE CABPClient : public CClient {
public:
    CABPClient(void);
    ~CABPClient(void);

// access methods -------------------------------------------------------------
    /// set number of cvs
    void SetNumberOfItems(int nitems,int ntotbins);

    /// set cv
    void SetCV(int id,const CSmallString& name,const CSmallString& type,
                  double min_value,double max_value,int nbins,double alpha,double fconv,const CSmallString& unit);

    /// get number of bins
    int GetNumOfBins(void);

// commands -------------------------------------------------------------------
    /// register client on server side
    int RegisterClient(void);

    /// unregister client on server side
    bool UnregisterClient(void);

    /// get initial data
    bool GetInitialData(int* nisamples,
                        double* idpop,
                        double* ipop);

    /// exchange data with server
    bool ExchangeData(int* nisamples,
                        double* idpop,
                        double* ipop);

// section of private data ----------------------------------------------------
private:
    int                         ClientID;       // client ID
    int                         NCVs;           // number of items
    int                         NTotBins;       // total number of bins
    std::vector<CColVariable>   CVs;            // list of CVs

    /// write data for exchange
    void WriteExchangeData(CXMLElement* p_cele,
                            int* nisamples,
                            double* idpop,
                            double* ipop);

    /// read data from exchange
    void ReadExchangeData(CXMLElement* p_rele,
                            int* nisamples,
                            double* idpop,
                            double* ipop);

    /// clear data
    void ClearExchangeData(int* nisamples,
                            double* idpop,
                            double* ipop);

    /// save cvs into XML
    void SaveCVSInfo(CXMLElement* p_tele);
};

//------------------------------------------------------------------------------

extern CABPClient      ABPClient;

//------------------------------------------------------------------------------

#endif

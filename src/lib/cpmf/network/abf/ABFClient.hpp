#ifndef ABFClientH
#define ABFClientH
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
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

/*! \brief ABF client for fpmf library
 *
 */

class PMF_PACKAGE CABFClient : public CClient {
public:
    CABFClient(void);
    ~CABFClient(void);

// access methods -------------------------------------------------------------
    /// set number of cvs
    void SetNumberOfItems(int ncvs,int ntotbins);

    /// set temperature
    void SetTemperature(double temp);

    /// set energy unit
    void SetEnergyUnit(double fconv,const CSmallString& unit);

    /// set epot enabled
    void SetEpotEnabled(bool enabled);

    /// set cv
    void SetCV(int id,const CSmallString& name,const CSmallString& type,
               double min_value,double max_value,int nbins,double fconv,const CSmallString& unit);

    /// get number of bins
    int GetNumOfBins(void);

// commands -------------------------------------------------------------------
    /// register client on server side
    int RegisterClient(void);

    /// unregister client on server side
    bool UnregisterClient(void);

    /// get initial data
    bool GetInitialData(int* nisamples,
                        double* inc_icfsum,
                        double* inc_icfsum2,
                        double* inc_epotsum,
                        double* inc_epotsum2,
                        double* inc_icfepotsum,
                        double* inc_icfepotsum2);

    /// exchange data with server
    bool ExchangeData(int* nisamples,
                        double* inc_icfsum,
                        double* inc_icfsum2,
                        double* inc_epotsum,
                        double* inc_epotsum2,
                        double* inc_icfepotsum,
                        double* inc_icfepotsum2);

// section of private data ----------------------------------------------------
private:
    int                         ClientID;       // client ID
    int                         NCVs;           // number of CVs
    int                         NTotBins;       // total number of bins
    std::vector<CColVariable>   CVs;            // list of CVs
    double                      Temperature;
    double                      EnergyFConv;
    CSmallString                EnergyUnit;
    bool                        EpotEnabled;

    /// write data for exchange
    void WriteExchangeData(CXMLElement* p_cele,
                            int* nisamples,
                            double* inc_icfsum,
                            double* inc_icfsum2,
                            double* inc_epotsum,
                            double* inc_epotsum2,
                            double* inc_icfepotsum,
                            double* inc_icfepotsum2);

    /// read data from exchange
    void ReadExchangeData(CXMLElement* p_rele,
                            int* nisamples,
                            double* inc_icfsum,
                            double* inc_icfsum2,
                            double* inc_epotsum,
                            double* inc_epotsum2,
                            double* inc_icfepotsum,
                            double* inc_icfepotsum2);

    /// clear data
    void ClearExchangeData(int* nisamples,
                            double* inc_icfsum,
                            double* inc_icfsum2,
                            double* inc_epotsum,
                            double* inc_epotsum2,
                            double* inc_icfepotsum,
                            double* inc_icfepotsum2);

    /// save cvs into XML
    void SaveCVSInfo(CXMLElement* p_tele);
};

//------------------------------------------------------------------------------

extern CABFClient      ABFClient;

//------------------------------------------------------------------------------

#endif

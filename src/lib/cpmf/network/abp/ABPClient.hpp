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
#include <PMFAccumulator.hpp>

//------------------------------------------------------------------------------

/*! \brief ABP client for fpmf library
 *
 */

class PMF_PACKAGE CABPClient : public CClient {
public:
    CABPClient(void);
    ~CABPClient(void);

// commands -------------------------------------------------------------------
    /// register client on server side
    int RegisterClient(void);

    /// unregister client on server side
    bool UnregisterClient(void);

    /// get initial data
    bool GetInitialData(double* nsamples,
                        double* dpop,
                        double* pop);

    /// exchange data with server
    bool ExchangeData(  double* inc_nsamples,
                        double* inc_dpop,
                        double* inc_pop);

// section of public data -----------------------------------------------------
public:
   CPMFAccumulatorPtr       Accu;
   int                      NumOfCVs;
   int                      NumOfBins;
   CSimpleVector<double>    Widths;

// section of private data ----------------------------------------------------
private:
    int                         ClientID;       // client ID

    void ClearExchangeData( double* inc_nsamples,
                            double* inc_dpop,
                            double* inc_pop);
};

//------------------------------------------------------------------------------

extern CABPClient      ABPClient;

//------------------------------------------------------------------------------

#endif

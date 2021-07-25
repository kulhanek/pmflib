#ifndef MTDClientH
#define MTDClientH
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
#include <PMFAccumulator.hpp>

//------------------------------------------------------------------------------

/*! \brief MTD client for fpmf library
 *
 */

class PMF_PACKAGE CMTDClient : public CClient {
public:
    CMTDClient(void);
    ~CMTDClient(void);

// commands -------------------------------------------------------------------
    /// register client on server side
    int RegisterClient(void);

    /// unregister client on server side
    bool UnregisterClient(void);

    /// get initial data
    bool GetInitialData(double* nsamples,
                        double* mtdforce,
                        double* mtdpot);

    /// exchange data with server
    bool ExchangeData(  double* inc_nsamples,
                        double* inc_mtdforce,
                        double* inc_mtdpot);

// section of public data -----------------------------------------------------
public:
   CPMFAccumulatorPtr       Accu;
   int                      NumOfCVs;
   int                      NumOfBins;
   CSimpleVector<double>    Widths;
   double                   WTemp;

// section of private data ----------------------------------------------------
private:
    int                 ClientID;

    void ClearExchangeData( double* inc_nsamples,
                            double* inc_mtdforce,
                            double* inc_mtdpot);
};

//------------------------------------------------------------------------------

extern CMTDClient      MTDClient;

//------------------------------------------------------------------------------

#endif

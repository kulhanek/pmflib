#ifndef EnergySurfaceH
#define EnergySurfaceH
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2021 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2011 Petr Kulhanek, kulhanek@chemi.muni.cz
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

#include <PMFMainHeader.hpp>
#include <SmallString.hpp>
#include <SimpleVector.hpp>
#include <PMFAccumulator.hpp>

//------------------------------------------------------------------------------

class CEnergySurface;
typedef boost::shared_ptr<CEnergySurface>    CEnergySurfacePtr;

//------------------------------------------------------------------------------

/** \brief multidimensional surface
*/

class PMF_PACKAGE CEnergySurface {
public:
// constructor and destructor -------------------------------------------------
    CEnergySurface(void);
    virtual ~CEnergySurface(void);

// setup method ---------------------------------------------------------------
    /// return number of cvs
    int GetNumOfCVs(void) const;

    /// return number of bins
    int GetNumOfBins(void) const;

    /// return coordinate definition
    const CColVariablePtr GetCV(int cv) const;

// allocate/deallocate --------------------------------------------------------
    /// allocate array
    void Allocate(CPMFAccumulatorPtr accu);

    /// allocate array
    void Allocate(CEnergySurfacePtr surf,const std::vector<bool>& enabled_cvs);

    /// deallocate array
    void Deallocate(void);

    /// set all items to zero
    void Clear(void);

// access methods -------------------------------------------------------------
    /// set energy from point of index
    void SetEnergy(int ibin,double value);

    /// get energy from point of index
    double GetEnergy(int ibin) const;

    /// get energy from point of index
    double GetEnergyRealValue(int ibin) const;

    /// set error from point of index
    void SetError(int ibin,double value);

    /// get error from point of index
    double GetError(int ibin) const;

    /// get error from point of index
    double GetErrorRealValue(int ibin) const;

    /// get error from point of index including sigma level
    double GetErrorRealValueWithSLevel(int ibin) const;

    /// get sigmaF2 from sampled area
    double GetSigmaF2(bool includeglued=false) const;

    /// get sigmaF from sampled area
    double GetSigmaF(bool includeglued=false) const;

    /// get RMS error from sampled area at slevel
    double GetRMSError(bool includeglued=false) const;

    /// get max error at slevel
    double GetMaxError(bool includeglued=false) const;

    /// get sigmaF2 from all area
    double GetSigmaF2All(void) const;


    /// get temperature in iu
    double GetTemperature(void) const;

    /// get temperature in unit
    double GetRealTemperature(void) const;

    /// get temperature unit
    const CSmallString& GetTemperatureUnit(void) const;

    /// get temperature conversion factor to given unit
    double GetTemperatureFConv(void);


    /// get energy conversion factor to given unit
    double GetEnergyFConv(void);

    /// get energy unit
    const CSmallString& GetEnergyUnit(void) const;

// access methods -------------------------------------------------------------
    /// set number of samples
    void SetNumOfSamples(int ibin,int value);

    /// get number of samples
    int GetNumOfSamples(int ibin) const;

    /// return value of global minimum
    double GetGlobalMinimumValue(void) const;

    /// add offset to the whole surface
    void ApplyOffset(double offset);

    /// set sigma-level for error analysis
    void SetSLevel(double slevel);

    /// get sigma-level for error analysis
    double GetSLevel(void) const;

    /// set unsampled region to max energy from sampled region
    void AdaptUnsampledToMaxEnergy(void);

    /// set unsampled region to provided max energy
    void AdaptUnsampledToMaxEnergy(double maxene);

// point positions -------------------------------------------------------------
    /// convert point index to point position
    void GetPoint(int ibin,CSimpleVector<double>& point) const;

    /// convert point index to point position
    void GetIPoint(int ibin,CSimpleVector<int>& point) const;

    /// return global index
    int GetGlobalIndex(const CSimpleVector<int>& position) const;

    /// reduce ipoint
    void ReduceIPoint(const std::vector<bool>& keepcvs,CSimpleVector<int>& midx,CSimpleVector<int>& ridx);

    /// convert point to bin
    int IPoint2Bin(const CSimpleVector<int>& point);

// global minimum --------------------------------------------------------------
    /// set position of global minimum - spec in real units
    void SetGlobalMin(const CSmallString& spec);

    /// set position of global minimum - in internal units
    void SetGlobalMin(const CSimpleVector<double>& pos);

    /// is global minimum set
    bool IsGlobalMinSet(void);

    /// get position of global minima - in internal units
    CSimpleVector<double> GetGlobalMinPos(void);

    /// get position of global minima - bin
    int GetGlobalMinBin(void);

    /// get global minimum energy
    double GetGlobalMinEnergy(void);

    /// find the closest bin to provided global min
    void FindGlobalMinBin(void);

    /// find global minimum
    void FindGlobalMin(void);

// operators ------------------------------------------------------------------
    /// add energy surface
    void AddFES(CEnergySurfacePtr source);

    /// divide energy surface by constant
    void DivideFES(double number);

    /// reduce FES
    CEnergySurfacePtr ReduceFES(const std::vector<bool>& keepcvs);

// section of private data ----------------------------------------------------
private:
    int                             NumOfCVs;       // number of cvs
    int                             NumOfBins;      // number of bins
    std::vector<CColVariablePtr>    CVs;            // collective variables

    double                          Temperature;
    double                          TemperatureFConv;
    CSmallString                    TemperatureUnit;

    double                          EnergyFConv;
    CSmallString                    EnergyUnit;

    CSimpleVector<double>           Energy;         // energy array
    CSimpleVector<double>           Error;          // error array
    CSimpleVector<int>              Samples;        // number of samples
    double                          SLevel;         // sigma level for error calculation

    bool                            GlobalMinSet;   // true if gpos set by SetGlobalMin()
    CSimpleVector<double>           GPos;           // global position, either detected or use
    bool                            GPosSet;        // true is gpos set by any means, either SetGlobalMin() or from FES
    int                             GPosBin;
};

//------------------------------------------------------------------------------

#endif

#ifndef EnergySurfaceH
#define EnergySurfaceH
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
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
#include <MTDHistory.hpp>
#include <ABFAccumulator.hpp>
#include <ABPAccumulator.hpp>
#include <SimpleVector.hpp>

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
    int GetNumberOfCoords(void) const;

    /// return number of points
    int GetNumberOfPoints(void) const;

    /// return coordinate definition
    const CColVariable* GetCoordinate(unsigned int cv) const;

// allocate/deallocate --------------------------------------------------------
    /// allocate array
    void Allocate(const CMTDHistory* mtd_hist);

    /// allocate array
    void Allocate(const CABFAccumulator* abf_accu);

    /// allocate array
    void Allocate(const CABPAccumulator* abp_accu);

    /// allocate array
    void Allocate(const CABFAccumulator* abf_accu,const std::vector<bool>& enabled_cvs);

    /// allocate array
    void Allocate(const CEnergySurface* p_surf,const std::vector<bool>& enabled_cvs);

    /// deallocate array
    void Deallocate(void);

    /// set all items to zero
    void Clear(void);

// access methods -------------------------------------------------------------
    /// set energy from point of index
    void SetEnergy(unsigned int index,const double& value);

    /// get energy from point of index
    const double& GetEnergy(unsigned int index) const;

    /// set error from point of index
    void SetError(unsigned int index,const double& value);

    /// get error from point of index
    const double& GetError(unsigned int index) const;

    /// get sigmaF2 from sampled area
    double GetSigmaF2(bool includeglued=false) const;

    /// get sigmaF2+ from sampled area
    double GetSigmaF2p(bool includeglued=false) const;

    /// get sigmaF2- from sampled area
    double GetSigmaF2m(bool includeglued=false) const;

// access methods -------------------------------------------------------------
    /// set number of samples
    void SetNumOfSamples(unsigned int index,const int& value);

    /// get number of samples
    const int& GetNumOfSamples(unsigned int index) const;

    /// convert point index to point position
    void GetPoint(unsigned int index,CSimpleVector<double>& point) const;

    /// convert point index to point position
    void GetIPoint(unsigned int index,CSimpleVector<int>& point) const;

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

    /// reduce FES
    bool ReduceFES(const std::vector<bool>& keepcvs,double temp,CEnergySurface* p_rsurf);

    /// reduce ipoint
    void ReduceIPoint(const std::vector<bool>& keepcvs,CSimpleVector<int>& midx,CSimpleVector<int>& ridx);

    /// convert point to bin
    int IPoint2Bin(const CSimpleVector<int>& point);

// operators ------------------------------------------------------------------
    /// add energy surface
    void operator+=(const CEnergySurface& source);

    /// divide energy surface by constant
    void operator/=(const double& number);

// section of private data ----------------------------------------------------
private:
    int                         NumOfCVs;       // number of cvs
    int                         TotNPoints;     // total energy size
    CSimpleVector<CColVariable> Sizes;          // dimensions
    CSimpleVector<double>       Energy;         // energy array
    CSimpleVector<double>       Error;          // error array
    CSimpleVector<int>          Samples;        // number of samples
    double                      SLevel;         // sigma level for error calculation

    void CalculateFES_Part(CMTDHistory& mtd_hist,
                           CSimpleVector<double>& point,unsigned int mtd_time,
                           unsigned int& loc,unsigned int cv);

    void CalculateFES_MTDParam_Part(CSimpleVector<double>& params,
                                    CSimpleVector<double>& point,unsigned int& loc,unsigned int cv);

    double CalculateValue(const CSimpleVector<double>& params,
                          const CSimpleVector<double>& point);
};

//------------------------------------------------------------------------------

#endif

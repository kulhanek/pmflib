#ifndef PMFAccuDataH
#define PMFAccuDataH
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
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

#include <PMFMainHeader.hpp>
#include <SmallString.hpp>
#include <SimpleVector.hpp>
#include <stdio.h>
#include <memory>
#include <boost/shared_ptr.hpp>

//------------------------------------------------------------------------------

class CPMFAccuData;
typedef boost::shared_ptr<CPMFAccuData>       CPMFAccuDataPtr;

//------------------------------------------------------------------------------

/** \brief PMF accumulator data segment
*/

class PMF_PACKAGE CPMFAccuData {
public:
    CPMFAccuData(int nbins, int ncvs);

// I/O operation ---------------------------------------------------------------
    /// load data section
    void Load(FILE* p_fin,const CSmallString& keyline);

    /// save data section
    void Save(FILE* p_fout);

    /// load data section
    void Load(CXMLElement* p_ele);

    /// save data section
    void Save(CXMLElement* p_ele);

// setup methods ---------------------------------------------------------------
    // reset data to zero
    void Reset(void);

// operation -------------------------------------------------------------------
    /// is compatible with right section data?
    bool CheckCompatibility(CPMFAccuDataPtr right) const;

    /// duplicate without data
    CPMFAccuDataPtr CreateTheSame(void) const;

    // combine - AD
    void CombineAD(CPMFAccuDataPtr left,CPMFAccuDataPtr right);

    // combine - WA
    void CombineWA(CPMFAccuDataPtr left,CPMFAccuDataPtr left_nsamples,CPMFAccuDataPtr right,CPMFAccuDataPtr right_nsamples);

    // combine - M2
    void CombineM2(CPMFAccuDataPtr left,CPMFAccuDataPtr left_nsamples,CPMFAccuDataPtr left_mean,
                   CPMFAccuDataPtr right,CPMFAccuDataPtr right_nsamples,CPMFAccuDataPtr right_mean);

    // combine - CO
    void CombineCO(CPMFAccuDataPtr left,CPMFAccuDataPtr left_nsamples,CPMFAccuDataPtr left_xmean,CPMFAccuDataPtr left_ymean,
                   CPMFAccuDataPtr right,CPMFAccuDataPtr right_nsamples,CPMFAccuDataPtr right_xmean,CPMFAccuDataPtr right_ymean);

// access methods --------------------------------------------------------------
    // get data section name
    const CSmallString& GetName(void) const;

    // get data length
    int GetLength(void) const;

    // get num of bins
    int GetNumOfBins(void) const;

    // get num of cvs
    int GetNumOfCVs(void) const;

    // get data section type
    const CSmallString& GetType(void) const;

    // get data section mode
    const CSmallString& GetMode(void) const;

    // get data section op
    const CSmallString& GetOp(void) const;

    // get data section name for complementary mean
    const CSmallString& GetMXName(void) const;

    // get data section name for complementary mean
    const CSmallString& GetMYName(void) const;

    /// get data
    double GetData(int ibin, int cv=0) const;

    /// set data
    void SetData(int ibin, double value);

    /// set data
    void SetData(int ibin, int cv, double value);

    /// get blob data
    void GetDataBlob(double* p_blob);

    /// set blob data
    void SetDataBlob(double* p_blob);

// section of private data -----------------------------------------------------
private:
    int                     NumOfBins;
    int                     NumOfCVs;
    CSmallString            Name;       // name of the section
    CSmallString            Op;         // data operation
                                        // AD - add
                                        // WA - weighted average
                                        // M2 - second moment
                                        // CO - co-variance

    CSmallString            Type;       // data type: R - real, I - integer
    CSmallString            Mode;       // data mode: B - per bins, C - per CVs, M - mixed per bins and cvs
    int                     Size;       // size of data
    CSimpleVector<double>   Data;       // all data are kept as real numbers

    // required for combine operation of variance and co-variance
    CSmallString            MXName;     // name of data section with mean of X
    CSmallString            MYName;     // name of data section with mean of Y

    friend class CPMFAccumulator;
};

//------------------------------------------------------------------------------

#endif

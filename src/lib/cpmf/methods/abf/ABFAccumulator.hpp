#ifndef ABFAccumulatorH
#define ABFAccumulatorH
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
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
#include <ColVariable.hpp>
#include <SmallString.hpp>
#include <stdio.h>
#include <SimpleVector.hpp>

//------------------------------------------------------------------------------

enum EABFAccuValue {
    EABF_DG_VALUE       = 0,  // dG/dksi
    EABF_DG_SIGMA       = 1,  // fluctuation of ICF
    EABF_DG_ERROR       = 2,  // error of dG/dksi
    EABF_H_VALUE        = 3,  // enthalpy
    EABF_H_SIGMA        = 4,  // fluctuation of potential energy
    EABF_H_ERROR        = 5,  // error of enthalpy
    EABF_TDS_VALUE      = 6,  // -TdS/dksi
    EABF_TDS_SIGMA      = 7,  // fluctuation of -TdS/dksi
    EABF_TDS_ERROR      = 8,  // error of -TdS/dksi
};

//------------------------------------------------------------------------------

/** \brief ABF accumulator
*/

class PMF_PACKAGE CABFAccumulator {
public:
// constructor and destructor -------------------------------------------------
    CABFAccumulator(void);
    ~CABFAccumulator(void);

// input and output methods ---------------------------------------------------
    /// load accumulator data from file with name
    void Load(const CSmallString& name);

    /// load accumulator data from file
    void Load(FILE* fin);

    /// save accumulator data from file with name
    void Save(const CSmallString& name);

    /// save accumulator data from file
    void Save(FILE* fout);

    /// save mask
    void SaveMask(FILE* fout);

    /// save mask
    void SaveMaskGnuPlot(FILE* fout);

// dimension specification ----------------------------------------------------
    /// set number of coordinates, all previous data are destroyed
    void SetNumOfCVs(int ncvs);

    /// set coordinate data
    void SetCV(int id,
                const CSmallString& name,
                const CSmallString& type,
                double min_value,double max_value,int nbins);

    /// set coordinate data
    void SetCV(int id,
                const CSmallString& name,
                const CSmallString& type,
                double min_value,double max_value,int nbins,
                double fconv, const CSmallString& unit);

    /// allocate all arrays specified via SetNumOfCVs and SetCV
    void FinalizeAllocation(void);

    /// deallocate all array
    void Deallocate(void);

// access data methods --------------------------------------------------------
    /// return number of coordinates
    int GetNumOfCVs(void) const;

    /// return coordinate definition
    const CColVariable* GetCV(int cv) const;

    /// return number of bins within limit
    int GetNumOfBins(int limit=0) const;

    /// return global index
    int GetGlobalIndex(const CSimpleVector<int>& position) const;

    /// convert point index to point position
    void GetPoint(unsigned int index,CSimpleVector<double>& point) const;

    /// convert point index to point position
    void GetPointRValues(unsigned int index,CSimpleVector<double>& point) const;

    /// convert point index to point position
    void GetIPoint(unsigned int index,CSimpleVector<int>& point) const;

    /// set temperature
    void SetTemperature(double temp);

    /// get temperature
    double GetTemperature(void);

    /// set temperature
    void SetEnergyUnit(double fconv,const CSmallString& unit);

    /// etot flag
    void SetEtotEnabled(bool enabled);

//---------------------------------------------------------------------------------

    /// return number of samples for a given bin position
    int GetNumOfSamples(const CSimpleVector<int>& position) const;

    /// return number of samples for a given bin index
    int GetNumOfSamples(int ibin) const;

    /// return mask weight for a given bin index
    double GetMaskWeight(int ibin) const;

    /// set mask weight for a given bin index
    void SetMaskWeight(int ibin,double weight);

    /// return ICF sum
    double GetMICF(int icv,int ibin) const;

    /// return ICF square sum
    double GetM2ICF(int icv,int ibin) const;

    /// return Epot sum
    double GetMEtot(int ibin) const;

    /// return Epot square sum
    double GetM2Etot(int ibin) const;

    /// return ICF*Epot sum
    double GetICFMEtot(int icv,int ibin) const;

    /// return ICF*Epot square sum
    double GetICFM2Etot(int icv,int ibin) const;

    /// set ABF force sum
    void SetMICF(int icv,int ibin,const double& value);

    /// set ABF force sum square
    void SetM2ICF(int icv,int ibin,const double& value);

    /// set ABF force sum square
    void SetNumberOfABFSamples(int ibin,const int& value);

    /// return pointer to NSamples array
    int* GetNSamplesArray(void);

    /// return pointer to ICF sum array
    double* GetMICFArray(void);

    /// return pointer to ICF square sum array
    double* GetM2ICFArray(void);

    /// return pointer to Epot sum array
    double* GetMEtotArray(void);

    /// return pointer to Epot square sum array
    double* GetM2EtotArray(void);

    /// return pointer to ICF*Epot sum array
    double* GetICFMEtotArray(void);

    /// return pointer to ICF*Epot square sum array
    double* GetICFM2EtotArray(void);

    /// set number of statistically correlated samples (it influences calculated errors of mean forces)
    void SetNCorr(double ncorr);

    /// get value - dG/dksi
    double GetValue(int icv,int ibin,EABFAccuValue realm) const;

    /// get value - H=<PotEne>
    double GetValue(int ibin,EABFAccuValue realm) const;

// mathematical operation -----------------------------------------------------
    /// set all data to zero
    void Clear(void);

// information methods --------------------------------------------------------
    /// load cvs info
    void LoadCVSInfo(CXMLElement* p_iele);

    /// check cv info
    bool CheckCVSInfo(CXMLElement* p_iele) const;

    /// check cv info between two accumulators
    bool CheckCVSInfo(const CABFAccumulator* p_accu) const;

    /// save cvs info
    void SaveCVSInfo(CXMLElement* p_tele) const;

    /// print cvs info
    void PrintCVSInfo(std::ostream& vout);

// data transfer methods ------------------------------------------------------
    /// load only abf data from XML file
    void ReadABFData(CXMLElement* p_ele);

    /// add abf data from XML file
    void AddABFData(CXMLElement* p_ele);

    /// write ABF data to XML file
    void WriteABFData(CXMLElement* p_rele);

// arithmetic operations ------------------------------------------------------
    /// add another accumulator
    void AddABFAccumulator(const CABFAccumulator* p_accu);

    /// subtract another accumulator
    void SubABFAccumulator(const CABFAccumulator* p_accu);

// section of private data ----------------------------------------------------
private:
    int                         NCVs;           // number of CVs
    double                      NCorr;          // number of correlated samples (for error evaluation)
    double                      Temperature;    // temperature
    double                      EnergyFConv;    // energy unit
    CSmallString                EnergyUnit;
    // flags
    bool                        EtotAvailable;
    bool                        EkinIncluded;

    CSimpleVector<CColVariable> CVs;            // collective variables
    int                         TotNBins;       // number of total bins

    CSimpleVector<int>          NSamples;       // number of hits into bins
    CSimpleVector<double>       Mask;

    // all values are in internal units
    CSimpleVector<double>       MICF;           // mean instantaneous collective force (ICF)
    CSimpleVector<double>       M2ICF;          // M2 moment of ICF
    CSimpleVector<double>       MICFKin;        // mean instantaneous collective force (ICF) - kinetic term
    CSimpleVector<double>       M2ICFKin;       // M2 moment of ICF

    // enthalpy/entropy
    CSimpleVector<double>       MEtot;          // mean Etot
    CSimpleVector<double>       M2Etot;         // M2 moment of Etot
    CSimpleVector<double>       MEpot;          // mean Epot
    CSimpleVector<double>       M2Epot;         // M2 moment of Epot

    CSimpleVector<double>       CDS;            // accumulated ICF*Epot
    CSimpleVector<double>       M2CDS;          // accumulated ICF*Epot force squares

// helper methods -------------------------------------------------------------
    /// return index to MICF array for particular item and bin
    int map(int item,int bin) const;

    void Load_v6(char* fline,FILE* fin);
};

//------------------------------------------------------------------------------

#endif

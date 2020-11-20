#ifndef ColVariableH
#define ColVariableH
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
#include <ostream>

//------------------------------------------------------------------------------

class CXMLElement;

//------------------------------------------------------------------------------

class PMF_PACKAGE CColVariable {
public:
    CColVariable(void);

// ----------------------------------------------------------------------------
    /// set cv name
    void SetName(const CSmallString& lname);

    /// set cv type
    void SetType(const CSmallString& ltype);

    /// set cv min vale
    void SetMinValue(double lmin);

    /// set cv max vale
    void SetMaxValue(double lmax);

    /// set max CV movement in STM
    void SetMaxMovement(double lmax);

    /// set cv data
    void SetCoord(int lid,const CSmallString& lname,const CSmallString& ltype,
                  double lmin_value,double lmax_value,unsigned int lnbins);

    /// set cv data
    void SetCoord(int lid,const CSmallString& lname,const CSmallString& ltype,
                  double lmin_value,double lmax_value,unsigned int lnbins,
                  double lfconv,const CSmallString& lunit);

    /// set cv data
    void SetCoord(int lid,const CSmallString& lname,const CSmallString& ltype,
                  double lmin_value,double lmax_value,unsigned int lnbins,
                  double alpha,
                  double lfconv,const CSmallString& lunit);

    /// set cv data
    void SetCoord(int lid,const CSmallString& lname,const CSmallString& ltype);

    /// copy from
    void CopyFrom(const CColVariable* p_coord);

// ----------------------------------------------------------------------------
    /// return number of bins
    unsigned int GetNumberOfBins(void) const;

    /// return minimal value
    double GetMinValue(void) const;

    /// return maximal value
    double GetMaxValue(void) const;

    /// return maximal CV movement in STM
    double GetMaxMovement(void) const;

    /// return coordinate range (max-min)
    double GetRange(void) const;

    /// return bin width
    double GetBinWidth(void) const;

    /// return value for particular bin
    double GetValue(unsigned int bin) const;

    /// return value in user specified unit for particular bin
    double GetRValue(unsigned int bin) const;

    /// get difference between two CVs
    double GetDifference(double left,double right) const;

    /// return coordinate type
    const CSmallString& GetType(void) const;

    /// return coordinate name
    const CSmallString& GetName(void) const;

    /// return unit name
    const CSmallString& GetUnit(void) const;

    /// check periodicity
    bool IsPeriodic(void) const;

// information mathods --------------------------------------------------------
    /// print info about cv
    void PrintInfo(std::ostream& vout);

    /// load cv info
    void LoadInfo(CXMLElement* p_ele);

    /// check cv info
    bool CheckInfo(CXMLElement* p_ele) const;

    /// check cv info
    bool CheckInfo(const CColVariable* p_coord) const;

    /// save cv info
    void SaveInfo(CXMLElement* p_ele) const;

// section of private data ----------------------------------------------------
private:
    // all values are in internal units
    int             ID;             // CV id
    CSmallString    Name;           // name of coordinate
    CSmallString    Type;           // type of coordinate
    CSmallString    Unit;           // unit
    double          FConv;          // conversion factor from i2r units
    double          MinValue;       // left boundary of coordinate
    double          MaxValue;       // right boundary of coordinate
    // MUST BE SIGNED VALUE
    int             NBins;          // number of coordinate bins
    double          BinWidth;       // (max_value-min_value)/nbins
    double          Width;          // max_value-min_value
    double          MaxMovement;    // used by STM
    double          Alpha;          // used by ABP

    friend class CMTDHistory;
    friend class CABFAccumulator;
    friend class CABPAccumulator;
    friend class CRSTAccumulator;
    friend class CBeadList;
};

//------------------------------------------------------------------------------

#endif

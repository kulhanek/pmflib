#ifndef ESPrinterH
#define ESPrinterH
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
#include <EnergySurface.hpp>

//------------------------------------------------------------------------------

enum EESPrinterFormat {
    EESPF_PLAIN,
    EESPF_GNUPLOT,
};

//------------------------------------------------------------------------------

class PMF_PACKAGE CESPrinter {
public:
// constructor and destructor -------------------------------------------------
    CESPrinter(void);
    virtual ~CESPrinter(void);

// setup methods --------------------------------------------------------------
    /// printed energy surface
    void SetPrintedES(CEnergySurfacePtr p_es);

    /// set output format
    void SetOutputFormat(EESPrinterFormat format);

    /// set XFormat
    void SetXFormat(const CSmallString& xform);

    /// set YFormat
    void SetYFormat(const CSmallString& yform);

    /// only those data with limit higher than limit will be printed
    void SetSampleLimit(int limit);

    /// include error into output
    void SetIncludeError(bool set);

    /// include bin status
    void SetIncludeBinStat(bool set);

    /// should we include glued area to energy calculation?
    void IncludeGluedAreas(bool set);

// printing methods -----------------------------------------------------------
    /// print energy surface
    void Print(const CSmallString& name);

    /// print energy surface
    void Print(FILE* fout);

// section of private data ----------------------------------------------------
private:
    CEnergySurfacePtr       EnergySurface;
    EESPrinterFormat        Format;
    CSmallString            XFormat;
    CSmallString            YFormat;
    int                     PrintLimit;
    bool                    IncludeError;
    bool                    IncludeGluedBins;
    bool                    IncludeBinStat;

    void PrintPlain(FILE* fout);
};

//------------------------------------------------------------------------------

#endif


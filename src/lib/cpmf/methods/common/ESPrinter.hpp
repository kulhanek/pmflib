#ifndef ESPrinterH
#define ESPrinterH
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
#include <SimpleVector.hpp>

//------------------------------------------------------------------------------

class CEnergySurface;

//------------------------------------------------------------------------------

enum EESPrinterFormat {
    EESPF_PLAIN,
    EESPF_GNUPLOT,
    EESPF_PMF_FES
};

//------------------------------------------------------------------------------

class PMF_PACKAGE CESPrinter {
public:
// constructor and destructor -------------------------------------------------
    CESPrinter(void);
    virtual ~CESPrinter(void);

// setup methods --------------------------------------------------------------
    /// printed energy surface
    void SetPrintedES(const CEnergySurface* p_es);

    /// set output format
    void SetOutputFormat(EESPrinterFormat format);

    /// set XFormat
    void SetXFormat(const CSmallString& xform);

    /// set YFormat
    void SetYFormat(const CSmallString& yform);

    /// only those data with limit higher than limit will be printed
    void SetSampleLimit(int limit);

    /// include errors into output
    void SetIncludeErrors(bool set);

// printing methods -----------------------------------------------------------
    /// print energy surface
    void Print(const CSmallString& name);

    /// print energy surface
    void Print(FILE* fout);

// section of private data ----------------------------------------------------
private:
    const CEnergySurface*   EnergySurface;
    EESPrinterFormat        Format;
    CSmallString            XFormat;
    CSmallString            YFormat;
    int                     PrintLimit;
    bool                    IncludeErrors;

    void PrintPlain(FILE* fout);
    void PrintPMF_FES(FILE* fout);

    void Print_Part(FILE* fout,CSimpleVector<double>& point,
                    unsigned int& loc,unsigned int cv);
};

//------------------------------------------------------------------------------

#endif


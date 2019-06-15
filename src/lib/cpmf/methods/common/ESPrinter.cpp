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

#include <errno.h>
#include <string.h>
#include <ESPrinter.hpp>
#include <EnergySurface.hpp>
#include <ErrorSystem.hpp>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CESPrinter::CESPrinter(void)
{
    EnergySurface = NULL;
    XFormat = "%f";
    YFormat = "%f";
    PrintLimit = 0;
    Format = EESPF_PLAIN;
    IncludeErrors = false;
}

//------------------------------------------------------------------------------

CESPrinter::~CESPrinter(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CESPrinter::SetPrintedES(const CEnergySurface* p_es)
{
    EnergySurface = p_es;
}

//------------------------------------------------------------------------------

void CESPrinter::SetOutputFormat(EESPrinterFormat format)
{
    Format = format;
}

//------------------------------------------------------------------------------

void CESPrinter::SetXFormat(const CSmallString& xform)
{
    XFormat = xform;
}

//------------------------------------------------------------------------------

void CESPrinter::SetYFormat(const CSmallString& yform)
{
    YFormat = yform;
}

//------------------------------------------------------------------------------

void CESPrinter::SetSampleLimit(int limit)
{
    PrintLimit = limit;
}

//------------------------------------------------------------------------------

void CESPrinter::SetIncludeErrors(bool set)
{
    IncludeErrors = set;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CESPrinter::Print(const CSmallString& name)
{
    FILE* fout = fopen(name, "w");

    if(fout == NULL) {
        CSmallString error;
        error << "unable to open file '" << name << "' (" << strerror(errno) << ")";
        RUNTIME_ERROR(error);
    }

    Print(fout);
    fclose(fout);
}

//------------------------------------------------------------------------------

void CESPrinter::Print(FILE* fout)
{
    if( fout == NULL ) {
        INVALID_ARGUMENT("stream is not open");
    }

    if( EnergySurface == NULL ) {
        RUNTIME_ERROR("no surface is associated with printer");
    }

    if(EnergySurface->GetNumberOfCoords() == 0) {
        RUNTIME_ERROR("at least one coordinate has to be present in ES");
    }

    if(EnergySurface->GetNumberOfPoints() == 0) {
        RUNTIME_ERROR("no point is present in ES");
    }

    switch(Format) {
        case EESPF_PLAIN:
            PrintPlain(fout);
            break;
        case EESPF_GNUPLOT:
            PrintPlain(fout);
            break;
        case EESPF_PMF_FES:
            PrintPMF_FES(fout);
            break;
    }
}

//------------------------------------------------------------------------------

void CESPrinter::PrintPlain(FILE* fout)
{
// allocate point
    CSimpleVector<double> point;
    point.CreateVector(EnergySurface->GetNumberOfCoords());

// print surface
    unsigned int loc = 0;
    Print_Part(fout,point,loc,0);
}

//------------------------------------------------------------------------------

void CESPrinter::PrintPMF_FES(FILE* fout)
{
// 10  format(A3,1X,I3)

// write FES header ------------------
    if(fprintf(fout,"FES V1 %3d\n",EnergySurface->GetNumberOfCoords()) <= 0) {
        CSmallString error;
        error << "unable to write header";
        RUNTIME_ERROR(error);
    }

// write coordinate specification ----------------
    for(unsigned int i=0; i < EnergySurface->GetNumberOfCoords(); i++) {
//20  format(I2,1X,A10,1X,E18.11,1X,E18.11,1X,I6)
//25  format(I2,1X,A55)
        if(fprintf(fout,"%2d %10s %18.11e %18.11E %6d\n",i+1,
                   (const char*)EnergySurface->GetCoordinate(i)->GetType(),
                   EnergySurface->GetCoordinate(i)->GetMinValue(),
                   EnergySurface->GetCoordinate(i)->GetMaxValue(),
                   EnergySurface->GetCoordinate(i)->GetNumberOfBins()) <= 0) {
            CSmallString error;
            error << "unable to write coordinate definition id: " << i+1;
            RUNTIME_ERROR(error);
        }
        if(fprintf(fout,"%2d %55s\n",i+1,
                   (const char*)EnergySurface->GetCoordinate(i)->GetName()) <= 0) {
            CSmallString error;
            error << "unable to write coordinate definition id: " << i+1;
            RUNTIME_ERROR(error);
        }
    }

// FES energies =================================================================

// 40  format(4(E19.11,1X))
    int counter = 0;
    for(unsigned int i = 0; i < EnergySurface->GetNumberOfPoints(); i++) {
        if(fprintf(fout,"%19.11E ",EnergySurface->GetEnergy(i)) <= 0) {
            CSmallString error;
            error << "unable to write energy for point " << i;
            RUNTIME_ERROR(error);
        }
        if(counter % 4 == 3) fprintf(fout,"\n");
        counter++;
    }
    if(counter % 4 != 0) fprintf(fout,"\n");
}

//------------------------------------------------------------------------------

void CESPrinter::Print_Part(FILE* fout,CSimpleVector<double>& point,
                            unsigned int& loc,unsigned int cv)
{
    if(cv >= EnergySurface->GetNumberOfCoords()) {
        // calculate value

        if(EnergySurface->GetNumOfSamples(loc) >= PrintLimit) {
            // print point position
            CSmallString xform = XFormat + " ";
            for(unsigned int i=0; i < EnergySurface->GetNumberOfCoords(); i++) {
                if(fprintf(fout,xform,point[i]) <= 0) {
                    CSmallString error;
                    error << "unable to print X data to the file (" << strerror(errno) << ")";
                    RUNTIME_ERROR(error);
                }
            }
            CSmallString yform = YFormat + " ";
            // and value
            double value = EnergySurface->GetEnergy(loc);
            if(fprintf(fout,yform,value) <= 0) {
                CSmallString error;
                error << "unable to print Y data to the file (" << strerror(errno) << ")";
                RUNTIME_ERROR(error);
            }
            if( IncludeErrors ){
                double error = EnergySurface->GetError(loc);
                if(fprintf(fout,yform,error) <= 0) {
                    CSmallString error;
                    error << "unable to print Y (error) data to the file (" << strerror(errno) << ")";
                    RUNTIME_ERROR(error);
                }
            }
            fprintf(fout,"\n");
        } else {
            if(EnergySurface->GetNumberOfCoords() == 1) fprintf(fout,"\n");
        }
        loc++;
        return;
    }

    const CColVariable* p_coord = EnergySurface->GetCoordinate(cv);

// cycle through variable
    for(unsigned int i = 0; i < p_coord->GetNumberOfBins(); i++) {
        point[cv] = p_coord->GetValue(i);
        Print_Part(fout,point,loc,cv+1);
    }

// write block delimiter - required by GNUPlot
    if(Format == EESPF_GNUPLOT) {
        if(fprintf(fout,"\n") <= 0) {
            CSmallString error;
            error << "unable to print GNUPlot delimiter to the file (" << strerror(errno) << ")";
            RUNTIME_ERROR(error);
        }
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================


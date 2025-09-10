// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2025 Petr Kulhanek, kulhanek@chemi.muni.cz
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

#include <GPREngineAUS.hpp>
#include <iomanip>

//------------------------------------------------------------------------------

using namespace std;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CGPREngineAUS::CGPREngineAUS(void)
{
    DoBalanceResiduals = false;
}

//------------------------------------------------------------------------------

CGPREngineAUS::~CGPREngineAUS(void)
{
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CGPREngineAUS::SetAccumulator(CPMFAccumulatorPtr accu)
{
}

//------------------------------------------------------------------------------

void CGPREngineAUS::SetOutputFEN(CEnergySurfacePtr p_surf)
{
    ASurface = p_surf;
}

//------------------------------------------------------------------------------

void CGPREngineAUS::SetOutputINT(CEnergySurfacePtr p_surf)
{
    USurface = p_surf;
}

//------------------------------------------------------------------------------

void CGPREngineAUS::SetOutputTDS(CEnergySurfacePtr p_surf)
{
    SSurface = p_surf;
}

//------------------------------------------------------------------------------

void CGPREngineAUS::SetOutputRES(CEnergySurfacePtr p_surf)
{
    RSurface = p_surf;
}

//------------------------------------------------------------------------------

void CGPREngineAUS::SetBalanceResiduals(bool iset)
{
    DoBalanceResiduals = iset;
}

//------------------------------------------------------------------------------

void CGPREngineAUS::PrepForMFInfo(void)
{
    NeedInv = true;
}

//------------------------------------------------------------------------------

int CGPREngineAUS::GetNumOfTasks(void)
{
    return(0);
}

//------------------------------------------------------------------------------

bool CGPREngineAUS::WriteMFInfo(const CSmallString& name,int task)
{
    return(false);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CGPREngineAUS::CalcResiduals(CVerboseStr& vout,bool balanced)
{
    double mf = 0.0;
    double m2 = 0.0;
    double n  = 0.0;
    double maxabs = 0.0;

    // https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Welford's_online_algorithm

    for(size_t ibin=0; ibin < NumOfBins; ibin++){
        if( ASurface->GetNumOfSamples(ibin) >  0 ){
            double a = ASurface->GetEnergy(ibin);
            double u = USurface->GetEnergy(ibin);
            double s = SSurface->GetEnergy(ibin);
            double res = a-(u+s);

            if( fabs(res) >  maxabs ){
                maxabs = fabs(res);
            }

            n++;
            double dx1 = (res - mf);
            mf = mf + dx1/n;
            double dx2 = (res - mf);
            m2 = m2 + dx1*dx2;

            if( RSurface ) {
                RSurface->SetNumOfSamples(ibin,ASurface->GetNumOfSamples(ibin));
                RSurface->SetEnergy(ibin,res);
                if( IncludeError ) {
                    double ea1 = ASurface->GetError(ibin);
                    double eu1 = USurface->GetError(ibin);
                    double es1 = SSurface->GetError(ibin);
                    double er = sqrt(ea1*ea1+eu1*eu1+es1*es1);
                    RSurface->SetError(ibin,er);
                }
            }
        }
    }

    double sigmaf = 0.0;
    if( n > 0 ) {
        sigmaf = sqrt(m2/n);
    }

    if( balanced ){
        vout << "      Residual AVE|B    = " << setw(10) << setprecision(5) << mf << endl;
        vout << "      Residual MaxAbs|B = " << setw(10) << setprecision(5) << maxabs << endl;
        vout << "      Residual SigmaF|B = " << setw(10) << setprecision(5) << sigmaf << endl;
    } else {
        vout << "      Residual AVE      = " << setw(10) << setprecision(5) << mf << endl;
        vout << "      Residual MaxAbs   = " << setw(10) << setprecision(5) << maxabs << endl;
        vout << "      Residual SigmaF   = " << setw(10) << setprecision(5) << sigmaf << endl;
    }
}

//------------------------------------------------------------------------------

void CGPREngineAUS::BalanceResiduals(void)
{
    double mf = 0.0;
    double n  = 0.0;

    // https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Welford's_online_algorithm

    for(size_t ibin=0; ibin < NumOfBins; ibin++){
        if( ASurface->GetNumOfSamples(ibin) >  0 ){
            double a = ASurface->GetEnergy(ibin);
            double u = USurface->GetEnergy(ibin);
            double s = SSurface->GetEnergy(ibin);
            double res = a-(u+s);
            n++;
            double dx1 = (res - mf);
            mf = mf + dx1/n;
        }
    }

// balance residuals to U and -TdS
    for(size_t ibin=0; ibin < NumOfBins; ibin++){
        if( ASurface->GetNumOfSamples(ibin) >  0 ){
            USurface->SetEnergy(ibin,USurface->GetEnergy(ibin)+mf/2.0);
            SSurface->SetEnergy(ibin,SSurface->GetEnergy(ibin)+mf/2.0);
        }
    }
}

//------------------------------------------------------------------------------

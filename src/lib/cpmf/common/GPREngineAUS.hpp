#ifndef GPREngineAUSH
#define GPREngineAUSH
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

#include <PMFMainHeader.hpp>
#include <GPREngine.hpp>
#include <EnergySurface.hpp>
#include <BaseProxy.hpp>

//------------------------------------------------------------------------------

/** \brief base for AUS GPR engines
*/

class PMF_PACKAGE CGPREngineAUS: public CGPREngine, public CBaseProxy {
public:
// constructor and destructor -------------------------------------------------
    CGPREngineAUS(void);
    virtual ~CGPREngineAUS(void);

// setup methods ---------------------------------------------------------------
    /// set accumulator
    virtual void SetAccumulator(CPMFAccumulatorPtr accu);

    /// set output energy surfaces
    virtual void SetOutputFEN(CEnergySurfacePtr p_surf);    // dA(x)
    virtual void SetOutputINT(CEnergySurfacePtr p_surf);    // dU(x)
    virtual void SetOutputTDS(CEnergySurfacePtr p_surf);    // -TdS(x)
    virtual void SetOutputRES(CEnergySurfacePtr p_surf);    // residuals: d(A)-(dU(x)-TdS(x))

    /// balance residual errors
    virtual void SetBalanceResiduals(bool iset);

    /// prepare for subsequent call WriteMFInfo
    virtual void PrepForMFInfo(void);

    /// return number of tasks
    virtual int GetNumOfTasks(void);

    /// write file with derivatives
    virtual bool WriteMFInfo(const CSmallString& name,int task);

// section of protected data ---------------------------------------------------
protected:
    CEnergySurfacePtr       ASurface;   // free energy
    CEnergySurfacePtr       USurface;   // internal energy
    CEnergySurfacePtr       SSurface;   // -TdS contribution
    CEnergySurfacePtr       RSurface;   // residuals: d(A)-(dU(x)-TdS(x))

    bool DoBalanceResiduals;

// calculate residuals
    void CalcResiduals(CVerboseStr& vout,bool balanced);
    void BalanceResiduals(void);
};

//------------------------------------------------------------------------------

typedef boost::shared_ptr<CGPREngineAUS>    CGPREngineAUSPtr;

//------------------------------------------------------------------------------

#endif

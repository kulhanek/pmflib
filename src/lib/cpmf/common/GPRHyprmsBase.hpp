#ifndef GPRHyprmsBaseH
#define GPRHyprmsBaseH
// =============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -----------------------------------------------------------------------------
//    Copyright (C) 2023 Petr Kulhanek, kulhanek@chemi.muni.cz
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
#include <vector>
#include <SimpleVector.hpp>
#include <boost/shared_ptr.hpp>

//------------------------------------------------------------------------------

/** \brief integrator of ABF accumulator employing gaussian process
*/

class PMF_PACKAGE CGPRHyprmsBase {
public:
// constructor and destructor -------------------------------------------------
    CGPRHyprmsBase(void);
    virtual ~CGPRHyprmsBase(void);

// base methods ---------------------------------------------------------------
    /// get log of Marginal Likelihood
    virtual double GetLogML(void);

    /// get derivative of logML wrt hyperparameters
    /// order sigmaf2, covar, wfac, ncorr, sigman2: only requested ders are calculated
    /// derivatives are ADDED to der
    virtual void GetLogMLDerivatives(const std::vector<bool>& flags,CSimpleVector<double>& der);

    /// get the log of pseudo-likelihood from leave-one-out cross-validation (LOO-CV)
    virtual double GetLogPL(void);

    /// get derivative of logPL wrt hyperparameters
    /// order sigmaf2, covar, wfac, ncorr, sigman2: only requested ders are calculated
    /// derivatives are ADDED to der
    virtual void GetLogPLDerivatives(const std::vector<bool>& flags,CSimpleVector<double>& der);
};

//------------------------------------------------------------------------------

typedef boost::shared_ptr<CGPRHyprmsBase>    CGPRHyprmsBasePtr;

//------------------------------------------------------------------------------

#endif

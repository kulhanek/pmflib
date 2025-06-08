#ifndef CVSplineSmoothingCubicH
#define CVSplineSmoothingCubicH
// ===============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -------------------------------------------------------------------------------
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
#include <CVSpline.hpp>
#include <SimpleVector.hpp>

// Smoothing with cubic splines
// D. S. G. Pollock
// Queen Mary and Westfield College, The University of London

// Smoothing Cubic Spline, algorithm (84)

//------------------------------------------------------------------------------

class PMF_PACKAGE CCVSplineSmoothingCubic : public CCVSpline {
public:
    CCVSplineSmoothingCubic(void);
    virtual ~CCVSplineSmoothingCubic(void);

// setup method ----------------------------------------------------------------
    /// load spline setup
    virtual bool LoadSetup(CPrmFile& prmfile,std::ostream& vout);

    /// print setup
    virtual void PrintSetup(std::ostream& vout);

// setup method ----------------------------------------------------------------
    /// clear all data
    virtual void Clear(void);

    /// allocate memory for data
    virtual void Allocate(int numofknots);

    /// register data point
    virtual void SetPoint(int knotid,double alpha,double cv);

    /// finalize spline
    virtual void BuildSpline(void);

    /// set lambda
    void SetLambda(double lam);

    /// set sigma
    void SetSigma(int knotid,double sig);

// information methods ---------------------------------------------------------
    /// get CV value for given alpha
    virtual double GetCV(double alpha);

    /// get CV first derivatives for given alpha
    virtual double GetCVFirstDer(double alpha);

// section of private data -----------------------------------------------------
private:
    int                     n;      // number of knots - 1
    CSimpleVector<double>   x;      // alphas       indexing: 0,1,...,n
    CSimpleVector<double>   y;      // CV values    indexing: 0,1,...,n
    CSimpleVector<double>   sigma;
    double                  lambda;
    double                  all_sigma;

    // spline
    CSimpleVector<double>   sa, sb, sc, sd;   //indexing: 0,1,...,n

    void Quincunx(CSimpleVector<double>& u, CSimpleVector<double> &v,
                  CSimpleVector<double>& w, CSimpleVector<double> &q);
};

//------------------------------------------------------------------------------

typedef std::shared_ptr<CCVSplineSmoothingCubic> CCVSplineSmoothingCubicPtr;

//------------------------------------------------------------------------------

#endif

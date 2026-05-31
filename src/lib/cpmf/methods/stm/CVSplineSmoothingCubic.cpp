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
// ===============================================================================

#include <CVSplineSmoothingCubic.hpp>
#include <iomanip>

//------------------------------------------------------------------------------

using namespace std;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CCVSplineSmoothingCubic::CCVSplineSmoothingCubic(void)
{
    lambda      = 0.999; // 1.0 - interpolating spline
    all_sigma   = 0.01;
}

//------------------------------------------------------------------------------

CCVSplineSmoothingCubic::~CCVSplineSmoothingCubic(void)
{
    Clear();
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CCVSplineSmoothingCubic::LoadSetup(CPrmFile& prmfile,std::ostream& vout)
{
    if( prmfile.GetDoubleByKey("lambda",lambda) == true  ) {
        vout << "Lambda (lambda)                                = " << left << setw(20) << lambda << endl;
    } else {
        vout << "Lambda (lambda)                                = " << left << setw(20) << lambda << "  (default)" << endl;
    }
    if( prmfile.GetDoubleByKey("sigma",all_sigma) == true  ) {
        vout << "All sigmas (sigma)                             = " << left << setw(20) << all_sigma << endl;
    } else {
        vout << "Lambda (lambda)                                = " << left << setw(20) << all_sigma << "  (default)" << endl;
    }

    vout << endl;
    return(true);
}

//------------------------------------------------------------------------------

void CCVSplineSmoothingCubic::PrintSetup(std::ostream& vout)
{
    vout << "Type        = smoothing cubic spline" << endl;
    vout << "Lambda      = " << lambda << endl;
    vout << "Sigma (all) = " << all_sigma << endl;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CCVSplineSmoothingCubic::LoadInfo(CXMLElement* p_ele)
{
    if( p_ele == NULL ) return(false);
    CSmallString type;

    p_ele->GetAttribute("type",type);
    if( type != "smoothing-cubic") return(false);

    lambda = 1.0; // 1.0 - interpolating spline
    all_sigma = 1.0;

    p_ele->GetAttribute("lambda",lambda);
    p_ele->GetAttribute("all_sigma",all_sigma);

    return(true);
}

//------------------------------------------------------------------------------

void CCVSplineSmoothingCubic::SaveInfo(CXMLElement* p_ele)
{
    if( p_ele == NULL ) return;
    p_ele->SetAttribute("type","smoothing-cubic");
    p_ele->SetAttribute("lambda",lambda);
    p_ele->SetAttribute("all_sigma",all_sigma);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CCVSplineSmoothingCubic::Clear(void)
{
    x.FreeVector();
    y.FreeVector();
    sa.FreeVector();
    sb.FreeVector();
    sc.FreeVector();
    sd.FreeVector();
    n = -1;
}

//------------------------------------------------------------------------------

void CCVSplineSmoothingCubic::Allocate(int numofknots)
{
    Clear();

    n = numofknots-1;
    if( n <= 0 ){
        n = 0;
        return;
    }

    x.CreateVector(n+1); // 0,1,...,n
    x.SetZero();

    y.CreateVector(n+1); // 0,1,...,n
    y.SetZero();

    sigma.CreateVector(n+1); // 0,1,...,n
    sigma.Set(all_sigma);

    sa.CreateVector(n+1); // 0,1,...,n
    sa.SetZero();

    sb.CreateVector(n+1); // 0,1,...,n
    sb.SetZero();

    sc.CreateVector(n+1); // 0,1,...,n
    sc.SetZero();

    sd.CreateVector(n+1); // 0,1,...,n
    sd.SetZero();
}

//------------------------------------------------------------------------------

void CCVSplineSmoothingCubic::SetPoint(int knotid,double alpha,double cv)
{
    if( (knotid < 0) || (knotid > n)) {
        RUNTIME_ERROR("knotid is out-of-range");
    }

    x[knotid] = alpha;
    y[knotid] = cv;
    sigma[knotid] = all_sigma;
}

//------------------------------------------------------------------------------

void CCVSplineSmoothingCubic::SetLambda(double lam)
{
    if( (lam <= 0) || (lam > 1.0) ){
        RUNTIME_ERROR("lambda out-of-range (0.0;1.0>");
    }
    lambda = lam;
}

//------------------------------------------------------------------------------

void CCVSplineSmoothingCubic::SetSigma(int knotid,double sig)
{
    if( (knotid < 0) || (knotid > n)) {
        RUNTIME_ERROR("knotid is out-of-range");
    }

    sigma[knotid] = sig;
}

//------------------------------------------------------------------------------

void CCVSplineSmoothingCubic::BuildSpline(void)
{
    sa.SetZero();
    sb.SetZero();
    sc.SetZero();
    sd.SetZero();

    if (n <= 0){
        RUNTIME_ERROR("not enough of knots");
    }

    if (n == 1) {
        // linear interpolation - between two points
        sd[0] = y[0];
        sc[0] = (y[1] - y[0]) / (x[1] - x[0]);
        return;
    }

// helpers
    CSimpleVector<double>   h, r, f, p, q, u, v, w;

    h.CreateVector(n+1); // 0,1,...,n
    h.SetZero();

    r.CreateVector(n+1);
    r.SetZero();

    f.CreateVector(n+1);
    f.SetZero();

    p.CreateVector(n+1);
    p.SetZero();

    q.CreateVector(n+1);
    q.SetZero();

    u.CreateVector(n+1);
    u.SetZero();

    v.CreateVector(n+1);
    v.SetZero();

    w.CreateVector(n+1);
    w.SetZero();

// --------------------
    double mu = 2.0*(1.0-lambda)/(3.0*lambda);

    h[0] = x[1] - x[0];
    r[0] = 3.0/h[0];
    for(int i=1; i <= n-1; i++){
        h[i] = x[i+1] - x[i];
        r[i] = 3.0/h[i];
        f[i] = -(r[i-1]+r[i]);
        p[i] = 2.0*(x[i+1]-x[i-1]);
        q[i] = 3.0*(y[i+1]-y[i])/h[i] - 3.0*(y[i]-y[i-1])/h[i-1];
    }

// --------------------
    v[0] = h[0];
    for(int i=1; i <= n-1; i++){
        u[i] = r[i-1]*r[i-1]*sigma[i-1]
             + f[i]*f[i]*sigma[i]
             + r[i]*r[i]*sigma[i+1];
        u[i] = mu*u[i] + p[i];
        v[i] = f[i]*r[i]*sigma[i] + r[i]*f[i+1]*sigma[i+1];
        v[i] = mu*v[i] + h[i];
        w[i] = mu*r[i]*r[i+1]*sigma[i+1];
    }

    // call Quincunx
    Quincunx(u,v,w,q);

// --------------------
    sd[0] = y[0] - mu*r[0]*q[1]*sigma[0];
    sd[1] = y[1] - mu*(f[1]*q[1]+r[1]*q[2])*sigma[0];
    sa[0] = q[1]/(3.0*h[0]);
    sb[0] = 0.0;
    sc[0] = (sd[1]-sd[0])/h[0] - q[1]*h[0]/3.0;
    r[0] = 0.0;
    for(int j=1; j <= n-1; j++){
        sa[j] = (q[j+1]-q[j])/(3.0*h[j]);
        sb[j] = q[j];
        sc[j] = (q[j]+q[j-1])*h[j-1]+sc[j-1];
        sd[j] = r[j-1]*q[j-1]+f[j]*q[j]+r[j]*q[j+1];
        sd[j] = y[j]-mu*sd[j]*sigma[j];
    }
}

//------------------------------------------------------------------------------

void CCVSplineSmoothingCubic::Quincunx(CSimpleVector<double>& u, CSimpleVector<double> &v,
                                       CSimpleVector<double>& w, CSimpleVector<double> &q)
{
// --------------------
    // u[-1] = 0.0
    u[0] = 0.0;

    // u[j] = u[j] - u[j-2]*w[j-2]*w[j-2] - u[j-1]*v[j-1]*v[j-1];
    // u[1] = u[1] - u[-1]*w[-1]*w[-1] - u[0]*v[0]*v[0];
    // u[1] - keep

    // v[j] = (v[j] - u[j-1]*v[j-1]*w[j-1])/u[j];
    // v[1] = (v[1] - u[0]*v[0]*w[0])/u[1];
    v[1] = v[1] / u[1];

    // w[j] = w[j]/u[j];
    // w[1] = w[1]/u[1];
    w[1] = w[1]/u[1];

    for(int j=2; j <= n-1; j++){
        u[j] = u[j] - u[j-2]*w[j-2]*w[j-2] - u[j-1]*v[j-1]*v[j-1];
        v[j] = (v[j] - u[j-1]*v[j-1]*w[j-1])/u[j];
        w[j] = w[j]/u[j];
    }

// --------------------
    // q[j] = q[j] - v[j-1]*q[j-1]-w[j-2]*q[j-2];
    // q[1] = q[1] - v[0]*q[0] - w[-1]*q[-1];
    q[1] = q[1] - v[0]*q[0];
    for(int j=2; j <= n-1; j++){
        q[j] = q[j] - v[j-1]*q[j-1]-w[j-2]*q[j-2];
    }
    for(int j=1; j <= n-1; j++){
        q[j] = q[j]/u[j];
    }

// --------------------
    //q[n+1]=0.0;
    q[n] = 0.0;
    // q[n-1] = q[n-1] - v[n-1]*q[n] - w[n-1]*q[n+1];
    // q[n-1] = q[n-1]
    // q[n-1] - keep
    for(int j=n-2; j >= 1; j--){
        q[j] = q[j] - v[j]*q[j+1] - w[j]*q[j+2];
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

double CCVSplineSmoothingCubic::GetCV(double alpha)
{
    // Handle out-of-bounds queries by clamping to the closest interval
    if( alpha <= x[0] ) {
        alpha = x[0];
    } else if( alpha >= x[n] ) {
        alpha = x[n];
    }

    // Find the right interval [x[i], x[i+1]] using linear search
    int i = 0;
    while( (i <= n - 2) && (alpha >= x[i + 1]) ) {
        i++;
    }

    double dx = alpha - x[i];
    double sp = sd[i] + sc[i] * dx + sb[i] * dx * dx + sa[i] * dx * dx * dx;

    return(sp);
}

//------------------------------------------------------------------------------

double CCVSplineSmoothingCubic::GetCVFirstDer(double alpha)
{
    // Handle out-of-bounds queries by clamping to the closest interval
    if( alpha <= x[0] ) {
        alpha = x[0];
    } else if( alpha >= x[n] ) {
        alpha = x[n];
    }

    // Find the right interval [x[i], x[i+1]] using linear search
    int i = 0;
    while( (i <= n - 2) && (alpha >= x[i + 1]) ) {
        i++;
    }

    double dx = alpha - x[i];
    double de = sc[i] + 2.0 * sb[i] * dx + 3.0 * sa[i] * dx * dx;

    return(de);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================


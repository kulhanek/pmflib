// ===============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -------------------------------------------------------------------------------
//    Copyright (C) 2011 Petr Kulhanek, kulhanek@chemi.muni.cz
//    Copyright (C) 2010 Petr Kulhanek, kulhanek@chemi.muni.cz
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

#include <Bead.hpp>
#include <ErrorSystem.hpp>
#include <XMLElement.hpp>
#include <math.h>
#include <STMPath.hpp>
#include <algorithm>


//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CBead::CBead(void)
{
    BeadList = NULL;

    NumOfCVs = 0;

    BeadID = 0;
    ClientID = -1;  // -1 - free, > 0 - registered

    Mode = BMO_UNKNOWN;
    ModeStatus = BMS_FINISHED;

    Permanent = false;
    Alpha = 0.0;
    dAdAlpha = 0.0;
    A = 0.0;
    NumOfUpdates = 0;

    beta1told   = 1.0;
    beta2told   = 1.0;
    beta1tnew   = 1.0;
    beta2tnew   = 1.0;

}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CBead::GetClientID(void)
{
    return(ClientID);
}

//------------------------------------------------------------------------------

int CBead::GetBeadID(void)
{
    return(BeadID);
}

//------------------------------------------------------------------------------

double CBead::GetPos(int cv)
{
    return(Pos[cv]);
}

//------------------------------------------------------------------------------

int CBead::GetMode(void)
{
    return(Mode);
}

//------------------------------------------------------------------------------

char CBead::GetModeString(void)
{
    switch(Mode){
        case BMO_UNKNOWN:
            return('U');
        case BMO_INITIALIZATION:
            return('I');
        case BMO_ACCUMULATION:
            return('A');
        case BMO_EQUILIBRATION:
            return('E');
        case BMO_PRODUCTION:
            return('P');
        case BMO_WAITFORRENDEZVOUS:
            return('W');
        case BMO_TERMINATE:
            return('T');
    }
    return('-');
}

//------------------------------------------------------------------------------

int CBead::GetModeStatus(void)
{
    return(ModeStatus);
}

//------------------------------------------------------------------------------

int CBead::GetModeLength(void)
{
    switch(Mode){
        case BMO_UNKNOWN:
        default:
            return(0);
        case BMO_INITIALIZATION:
            return(BeadList->InitPeriod);
        case BMO_ACCUMULATION:
            return(BeadList->AccuPeriod);
        case BMO_EQUILIBRATION:
            return(BeadList->EquiPeriod);
        case BMO_PRODUCTION:
            return(BeadList->ProdPeriod);
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CBead::InitBead(CSTMPath* p_list,int ncvs)
{
    NumOfCVs = ncvs;
    BeadList = p_list;
    OPos.CreateVector(NumOfCVs);
    OPos.SetZero();
    Pos.CreateVector(NumOfCVs);
    Pos.SetZero();
    NPos.CreateVector(NumOfCVs);
    NPos.SetZero();
    SPos.CreateVector(NumOfCVs);
    SPos.SetZero();
    FPos.CreateVector(NumOfCVs);
    FPos.SetZero();
    PPos.CreateVector(NumOfCVs);
    PPos.SetZero();
    MF.CreateVector(NumOfCVs);
    MF.SetZero();
    dCV.CreateVector(NumOfCVs);
    dCV.SetZero();
    pMF.CreateVector(NumOfCVs);
    pMF.SetZero();
    MTZ.CreateMatrix(NumOfCVs,NumOfCVs);
    MTZ.SetZero();
    P.CreateMatrix(NumOfCVs,NumOfCVs);
    P.SetZero();

    mtold.CreateVector(NumOfCVs);
    mtold.SetZero();
    mtnew.CreateVector(NumOfCVs);
    mtnew.SetZero();

    vtold.CreateVector(NumOfCVs);
    vtold.SetZero();
    vtnew.CreateVector(NumOfCVs);
    vtnew.SetZero();

    vthatold.CreateVector(NumOfCVs);
    vthatold.SetZero();
    vthatnew.CreateVector(NumOfCVs);
    vthatnew.SetZero();

    Alpha = 0;
    dAdAlpha = 0;
}

//------------------------------------------------------------------------------

void CBead::SetBeadData(int beadid,const CSimpleVector<double>& pos,bool flexible)
{
    BeadID = beadid;
    OPos = pos;
    Pos = pos;
    Permanent = !flexible;
}

//------------------------------------------------------------------------------

void CBead::SetClientID(int client_id)
{
    ClientID = client_id;
}

//------------------------------------------------------------------------------

void CBead::ReleaseBead(void)
{
    ModeStatus = BMS_PREPARED; // rollback any progress for current mode
}

//------------------------------------------------------------------------------

void CBead::MoveToNextMode(void)
{
    if( ModeStatus != BMS_FINISHED ) return; // keep current mode

    // clear spline data
    dCV.Set(0.0);
    P.SetZero();

    // clear accumulated data
    MF.Set(0.0);
    MTZ.SetZero();

    // clear derived data
    pMF.Set(0.0);
    dAdAlpha = 0.0;
    A = 0.0;

    switch(Mode){
        case BMO_UNKNOWN:
        default:
            // UN->I->(A->E->)P
            if( BeadList->InitPeriod > 0 ){
                Mode = BMO_INITIALIZATION;
            } else if( BeadList->EquiPeriod > 0 ) {
                Mode = BMO_EQUILIBRATION;
            } else if( BeadList->AccuPeriod > 0 ) {
                Mode = BMO_ACCUMULATION;
            } else if( BeadList->ProdPeriod > 0 ) {
                Mode = BMO_PRODUCTION;
            }
        break;
        case BMO_INITIALIZATION:
            if( BeadList->AccuPeriod > 0 ) {
                Mode = BMO_ACCUMULATION;
            } else {
                Mode = BMO_PRODUCTION;
            }
        break;
        case BMO_ACCUMULATION:
            if( BeadList->EquiPeriod > 0 ){
                Mode = BMO_EQUILIBRATION;
            } else {
                if( BeadList->GetSTMStatus() == ESTMS_PATH_FOUND ){
                    Mode = BMO_PRODUCTION;
                } else {
                    Mode = BMO_ACCUMULATION;
                }
            }
        break;
        case BMO_WAITFORRENDEZVOUS:
            if( BeadList->EquiPeriod > 0 ){
                Mode = BMO_EQUILIBRATION;
            } else {
                if( BeadList->GetSTMStatus() == ESTMS_PATH_FOUND ){
                    Mode = BMO_PRODUCTION;
                } else {
                    Mode = BMO_ACCUMULATION;
                }
            }
        break;
        case BMO_EQUILIBRATION:
            if( BeadList->GetSTMStatus() == ESTMS_PATH_FOUND ){
                Mode = BMO_PRODUCTION;
            } else {
                Mode = BMO_ACCUMULATION;
            }
        break;
        case BMO_PRODUCTION:
            Mode = BMO_PRODUCTION;
        break;
    }

    ModeStatus = BMS_PREPARED;
}

//------------------------------------------------------------------------------

void CBead::ResetPosUpdates(void)
{
    FPos = Pos;
    OPos = Pos;
    NPos = Pos;
    SPos = Pos;
    FPos = Pos;
}

//------------------------------------------------------------------------------

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

//------------------------------------------------------------------------------

void CBead::CalcProjector(void)
{
    if( Permanent ){
        for(int i=0; i < NumOfCVs; i++){
            pMF[i] = 0.0;
            for(int j=0; j < NumOfCVs; j++){
                P[i][j] = 0.0;
            }
        }
        return;
    }

    // calculate derivative vector length
    double slen2 = 0.0;
    for(int i=0; i < NumOfCVs; i++){
        double cvder = BeadList->CVSplines[i]->GetCVFirstDer(Alpha);
        slen2 += cvder*cvder;
    }

    if( slen2 == 0 ){
        LOGIC_ERROR("derivative segment has zero length");
    }

    // calculate projection
    for(int i=0; i < NumOfCVs; i++){
        for(int j=0; j < NumOfCVs; j++){
            if( i == j ){
                P[i][j] = 1.0;
            } else {
                P[i][j] = 0.0;
            }
            double cvder1 = BeadList->CVSplines[i]->GetCVFirstDer(Alpha);
            double cvder2 = BeadList->CVSplines[j]->GetCVFirstDer(Alpha);
            P[i][j] -= cvder1*cvder2/slen2;
        }
    }

    // projection perpendicular to the path
    for(int i=0; i < NumOfCVs; i++){
        double ps = 0;

        if( (BeadID == 1) || (BeadID == BeadList->GetNumOfBeads() ) ){
            // steepest descent movement
            for(int k=0; k < NumOfCVs; k++){
                ps += MTZ[i][k]*MF[k];
            }
        } else {
            // projection perpendicular to the path
            for(int j=0; j < NumOfCVs; j++){
                for(int k=0; k < NumOfCVs; k++){
                    ps += P[i][j]*MTZ[j][k]*MF[k];
                }
            }
        }
        pMF[i] = ps;
    }
}

// -----------------------------------------------------------------------------

void CBead::UpdatePositionGD(double step)
{
    NumOfUpdates++;

    if( Permanent ){
        for(int i=0; i < NumOfCVs; i++){
            NPos[i] = Pos[i];
        }
        return;
    }

    for(int i=0; i < NumOfCVs; i++){
        double maxmov = BeadList->CVs[i]->GetMaxMovement();
        if( (maxmov <= 0) || (fabs(pMF[i]*step) < maxmov) ){
            NPos[i] = Pos[i] - pMF[i]*step;
        } else {
            NPos[i] = Pos[i] - maxmov*sgn(pMF[i]*step);
        }
    }
}

// -----------------------------------------------------------------------------

void CBead::UpdatePositionNGD(double step,double mingnormeps)
{
    double g2 = 0.0;

    for(int i=0; i < NumOfCVs; i++){
        g2 = g2 + pMF[i]*pMF[i];
    }

    double gnorm = sqrt( g2 / (double)NumOfCVs );
    double nstep = step / (gnorm + mingnormeps);

    UpdatePositionGD(nstep);
}

// -----------------------------------------------------------------------------

void CBead::UpdatePositionNGDAuto(double step,double maxgnorm,double mingnormeps)
{
    double g2 = 0.0;

    for(int i=0; i < NumOfCVs; i++){
        g2 = g2 + pMF[i]*pMF[i];
    }

    double gnorm = sqrt( g2 / (double)NumOfCVs );
    double nstep = step / (gnorm + mingnormeps);

    if( gnorm < maxgnorm ){
        // std::cout << "GD:  " << step << " " << gnorm << " " << nstep << std::endl;
        UpdatePositionGD(step);
    } else {
        // std::cout << "NGD: " << step << " " << gnorm << " " << nstep << std::endl;
        UpdatePositionGD(nstep);
    }
}

// -----------------------------------------------------------------------------

void CBead::ResetADAM(void)
{
    beta1told   = 1.0;
    beta2told   = 1.0;
    beta1tnew   = 1.0;
    beta2tnew   = 1.0;

    mtold.SetZero();
    mtnew.SetZero();

    vtold.SetZero();
    vtnew.SetZero();

    vthatold.SetZero();
    vthatnew.SetZero();
}

// -----------------------------------------------------------------------------

void CBead::UpdatePositionADAM(double step,double beta1,double beta2,double mingnormeps)
{
    NumOfUpdates++;

    if( Permanent ){
        for(int i=0; i < NumOfCVs; i++){
            NPos[i] = Pos[i];
        }
        return;
    }

    for(int i=0; i < NumOfCVs; i++){
        mtnew[i] = beta1 * mtold[i] + (1.0 - beta1) * pMF[i];
        vtnew[i] = beta2 * vtold[i] + (1.0 - beta2) * pMF[i]*pMF[i];
    }

    // the step can be rejected later
    beta1tnew = beta1told * beta1;
    beta2tnew = beta2told * beta2;

    for(int i=0; i < NumOfCVs; i++){

        double mthat = mtnew[i]/(1.0-beta1tnew);
        double vthat = vtnew[i]/(1.0-beta2tnew);

        double dm = step*mthat/(sqrt(vthat)+mingnormeps);

        double maxmov = BeadList->CVs[i]->GetMaxMovement();

        if( (maxmov <= 0) || (fabs(dm) < maxmov) ){
            NPos[i] = Pos[i] - dm;
        } else {
            NPos[i] = Pos[i] - maxmov*sgn(dm);
        }
    }
}

// -----------------------------------------------------------------------------

void CBead::UpdatePositionAMSGrad(double step,double beta1,double beta2,double mingnormeps)
{
    NumOfUpdates++;

    if( Permanent ){
        for(int i=0; i < NumOfCVs; i++){
            NPos[i] = Pos[i];
        }
        return;
    }

    for(int i=0; i < NumOfCVs; i++){
        mtnew[i]    = beta1 * mtold[i] + (1.0 - beta1) * pMF[i];
        vtnew[i]    = beta2 * vtold[i] + (1.0 - beta2) * pMF[i]*pMF[i];
        vthatnew[i] = std::max(vthatold[i],vtnew[i]);
    //    std::cout << "pMF: " << pMF[i] << std::endl;
    }

    for(int i=0; i < NumOfCVs; i++){

     //       std::cout << "mtnew: " << mtnew[i] << " vthatnew: " << vthatnew[i] << std::endl;

        double dm = step*mtnew[i]/(sqrt(vthatnew[i])+mingnormeps);

        double maxmov = BeadList->CVs[i]->GetMaxMovement();

        if( (maxmov <= 0) || (fabs(dm) < maxmov) ){
            NPos[i] = Pos[i] - dm;
       //     std::cout << "dm: " << dm << std::endl;
        } else {
            NPos[i] = Pos[i] - maxmov*sgn(dm);
        }
    }
}

// -----------------------------------------------------------------------------

void CBead::UpdatePositionAMSGradBC(double step,double beta1,double beta2,double mingnormeps)
{
    NumOfUpdates++;

    if( Permanent ){
        for(int i=0; i < NumOfCVs; i++){
            NPos[i] = Pos[i];
        }
        return;
    }

    for(int i=0; i < NumOfCVs; i++){
        mtnew[i] = beta1 * mtold[i] + (1.0 - beta1) * pMF[i];
        vtnew[i] = beta2 * vtold[i] + (1.0 - beta2) * pMF[i]*pMF[i];
    }

    // the step can be rejected later
    beta1tnew = beta1told * beta1;
    beta2tnew = beta2told * beta2;

    for(int i=0; i < NumOfCVs; i++){

        double mthat = mtnew[i]/(1.0-beta1tnew);
        double vthat = vtnew[i]/(1.0-beta2tnew);

        vthatnew[i] = std::max(vthatold[i],vthat);

        double dm = step*mthat/(sqrt(vthatnew[i])+mingnormeps);

        double maxmov = BeadList->CVs[i]->GetMaxMovement();

        if( (maxmov <= 0) || (fabs(dm) < maxmov) ){
            NPos[i] = Pos[i] - dm;
        } else {
            NPos[i] = Pos[i] - maxmov*sgn(dm);
        }
    }
}

//------------------------------------------------------------------------------

void CBead::UpdatePositionFinalize(void)
{
    vtold = vtnew;
    mtold = mtnew;
    vthatold = vthatnew;
    beta1told = beta1tnew;
    beta2told = beta2tnew;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CBead::LoadInfo(CXMLElement* p_ele)
{
    if(p_ele == NULL) {
        INVALID_ARGUMENT("p_ele is NULL");
    }

    bool result = true;
    result &= p_ele->GetAttribute("bead_id",BeadID);
    result &= p_ele->GetAttribute("client_id",ClientID);
    result &= p_ele->GetAttribute("permanent",Permanent);
    result &= p_ele->GetAttribute("mode",Mode);
    result &= p_ele->GetAttribute("nupd",NumOfUpdates);

    if(result == false) {
        LOGIC_ERROR("unable to read some attributes");
    }

    // load positions
    CXMLBinData* p_posele = p_ele->GetFirstChildBinData("POS");
    if(p_posele == NULL) {
        LOGIC_ERROR("unable to open POS element");
    }
    Pos.Load(p_posele);

    CXMLBinData* p_pmfele = p_ele->GetFirstChildBinData("MF");
    if(p_pmfele == NULL) {
        LOGIC_ERROR("unable to open MF element");
    }
    MF.Load(p_pmfele);

    CXMLBinData* p_mtzele = p_ele->GetFirstChildBinData("MTZ");
    if(p_mtzele == NULL) {
        LOGIC_ERROR("unable to open MTZ element");
    }
    MTZ.Load(p_mtzele);

    CXMLBinData* p_rposele = p_ele->GetFirstChildBinData("OPOS");
    if(p_rposele == NULL) {
        LOGIC_ERROR("unable to open OPOS element");
    }
    OPos.Load(p_rposele);
}

//------------------------------------------------------------------------------

void CBead::SaveInfo(CXMLElement* p_ele)
{
    if(p_ele == NULL) {
        INVALID_ARGUMENT("p_ele is NULL");
    }

    p_ele->SetAttribute("bead_id",BeadID);
    p_ele->SetAttribute("client_id",ClientID);
    p_ele->SetAttribute("permanent",Permanent);
    p_ele->SetAttribute("mode",Mode);
    p_ele->SetAttribute("nupd",NumOfUpdates);

    // save position
    CXMLBinData* p_posele = p_ele->CreateChildBinData("POS");
    Pos.Save(p_posele);

    CXMLBinData* p_pmfele = p_ele->CreateChildBinData("MF");
    MF.Save(p_pmfele);

    CXMLBinData* p_mtzele = p_ele->CreateChildBinData("MTZ");
    MTZ.Save(p_mtzele);

    CXMLBinData* p_rposele = p_ele->CreateChildBinData("OPOS");
    OPos.Save(p_rposele);
}

//------------------------------------------------------------------------------

void CBead::GetProductionData(CXMLElement* p_ele)
{
    if(p_ele == NULL) {
        INVALID_ARGUMENT("p_ele is NULL");
    }

    // only if accumulation or production
    if( (GetMode() != BMO_ACCUMULATION) && (GetMode() != BMO_PRODUCTION) ) return;

    CXMLBinData* p_bposele = p_ele->GetFirstChildBinData("BPOS");
    if(p_bposele == NULL) {
        LOGIC_ERROR("unable to open BPOS element");
    }
    Pos.Load(p_bposele);

    CXMLBinData* p_pmfele = p_ele->GetFirstChildBinData("MF");
    if(p_pmfele == NULL) {
        LOGIC_ERROR("unable to open MF element");
    }
    MF.Load(p_pmfele);

    CXMLBinData* p_mtzele = p_ele->GetFirstChildBinData("MTZ");
    if(p_mtzele == NULL) {
        LOGIC_ERROR("unable to open MTZ element");
    }
    MTZ.Load(p_mtzele);

    ModeStatus = BMS_FINISHED;
}

//------------------------------------------------------------------------------

void CBead::SkipProductionData(void)
{
    ModeStatus = BMS_FINISHED;
}

//------------------------------------------------------------------------------

void CBead::WaitForRendezvous(void)
{
    Mode = BMO_WAITFORRENDEZVOUS;
    ModeStatus = BMS_FINISHED;
}

//------------------------------------------------------------------------------

void CBead::SetNextStepData(CXMLElement* p_ele)
{
    if( ModeStatus != BMS_PREPARED ){
        CSmallString error;
        error << "bead ID=" << BeadID << " is not in prepared mode, unable to set data for exchange";
        RUNTIME_ERROR(error);
    }

    if(p_ele == NULL) {
        INVALID_ARGUMENT("p_ele is NULL");
    }

    // program
    p_ele->SetAttribute("mode",GetMode());
    p_ele->SetAttribute("steps",GetModeLength());

    // and bead position
    CXMLBinData* p_bposele = p_ele->CreateChildBinData("BPOS");
    FPos.Save(p_bposele);   // use final position

    ModeStatus = BMS_RUNNING;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================


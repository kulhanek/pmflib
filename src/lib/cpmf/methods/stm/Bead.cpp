// ===============================================================================
// PMFLib - Library Supporting Potential of Mean Force Calculations
// -------------------------------------------------------------------------------
//    Copyright (C) 2025,2026 Petr Kulhanek, kulhanek@chemi.muni.cz
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

    BeadType = BTY_NORMAL;
    Alpha = 0.0;
    dAdAlpha = 0.0;
    A = 0.0;
    NumOfUpdates = 0;

    beta1told   = 1.0;
    beta2told   = 1.0;
    beta1tnew   = 1.0;
    beta2tnew   = 1.0;

    SegLength   = 0.0;
    KinkA       = 0.0;
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
    return(GetModeString(Mode));
}

//------------------------------------------------------------------------------

char CBead::GetModeString(int mode)
{
    switch(mode){
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

const CSmallString CBead::GetModeProgram(void)
{
    CSmallString mprog;

    mprog << GetModeString(Mode);

    switch(Mode){
        case BMO_UNKNOWN:
            break;
        case BMO_INITIALIZATION:
            if( BeadList->SoloInitPeriod == false ){
                if( BeadList->AccuPeriod > 0 ){
                    mprog << "+" << GetModeString(BMO_ACCUMULATION);
                } else {
                    if( BeadList->ProdPeriod > 0 ){
                        mprog << "+" << GetModeString(BMO_PRODUCTION);
                    }  
                }
            }
        case BMO_ACCUMULATION:
            break;
        case BMO_EQUILIBRATION:
            if( BeadList->SoloInitPeriod == false ){
                if( (BeadList->AccuPeriod > 0) && (BeadList->GetSTMStatus() != ESTMS_PATH_FOUND) ){
                    mprog << "+" << GetModeString(BMO_ACCUMULATION);
                } else {
                    if( BeadList->ProdPeriod > 0 ){
                        mprog << "+" << GetModeString(BMO_PRODUCTION);
                    }  
                }
            }
        case BMO_PRODUCTION:
            break;
        case BMO_WAITFORRENDEZVOUS:
            break;
        case BMO_TERMINATE:
            break;
    }

    return(mprog);
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
    BPos.CreateVector(NumOfCVs);
    BPos.SetZero();
    MF.CreateVector(NumOfCVs);
    MF.SetZero();
    cMF.CreateVector(NumOfCVs);
    cMF.SetZero();
    pMF.CreateVector(NumOfCVs);
    pMF.SetZero();
    uMF.CreateVector(NumOfCVs);
    uMF.SetZero();
    MTZ.CreateMatrix(NumOfCVs,NumOfCVs);
    MTZ.SetZero();
    dCVdAlpha.CreateVector(NumOfCVs);
    dCVdAlpha.SetZero();

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

void CBead::SetBeadData(int beadid,const CSimpleVector<double>& pos,int btype)
{
    BeadID = beadid;
    OPos = pos;
    Pos = pos;
    BeadType = btype;
}

//------------------------------------------------------------------------------

void CBead::SetClientID(int client_id)
{
    ClientID = client_id;
}

//------------------------------------------------------------------------------

void CBead::ReleaseBead(void)
{
    if( Mode == BMO_WAITFORRENDEZVOUS ) return; // this bead cannot be released

    ModeStatus = BMS_PREPARED; // rollback any progress for current mode
    
    if( Mode == BMO_ACCUMULATION ){
        if( (BeadList->SoloEquiPeriod == false) && (BeadList->EquiPeriod > 0) && (BeadList->AsynchronousMode == true) ){
            // this bead client will be effectivelly restarted from coordinates used for equi+accu
            // thus switch even further to equi mode
            Mode = BMO_EQUILIBRATION;
        }
    }
    if( Mode == BMO_PRODUCTION ){
        if( (BeadList->SoloEquiPeriod == false) && (BeadList->EquiPeriod > 0) && (BeadList->AsynchronousMode == true) ){
            // this bead client will be effectivelly restarted from coordinates used for equi+prod
            // thus switch even further to equi mode
            Mode = BMO_EQUILIBRATION;
        }
    }
}

//------------------------------------------------------------------------------

void CBead::MoveToNextMode(void)
{
    if( ModeStatus != BMS_FINISHED ) return; // keep current mode


    // clear accumulated data
    MF.Set(0.0);
    MTZ.SetZero();

    // clear derived data
    cMF.Set(0.0);
    pMF.Set(0.0);
    uMF.Set(0.0);
    dAdAlpha = 0.0;
    A = 0.0;

    switch(Mode){
        case BMO_UNKNOWN:
        default:
            // UN->I->(A->W->E->)P
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

void CBead::CalcBead(void)
{
    if( BeadType == BTY_PERMANENT ){
        cMF.SetZero();
        pMF.SetZero();
        uMF.SetZero();
        return;
    }

// corrected MF
    for(int i=0; i < NumOfCVs; i++){
        double sc = ( BeadList->CVs[i]->GetMaxValue() - BeadList->CVs[i]->GetMinValue() );
        double ps = 0;
        for(int j=0; j < NumOfCVs; j++){
            ps += sc * MTZ[i][j] * MF[j];
        }
        cMF[i] = ps;
    }

// calculate derivative vector length and dot product
    double dot = 0.0;
    double slen2 = 0.0;
    for(int i=0; i < NumOfCVs; i++){
        double cvder = BeadList->CVSplines[i]->GetCVFirstDer(Alpha);
        slen2 += cvder*cvder;
        dot += cvder * cMF[i];
    }

    if( slen2 == 0 ){
        LOGIC_ERROR("derivative segment has zero length");
    }

// MF perpendicular to path
// bead.pGrad[:] = bead.Grad[:] - bead.dCVdAlpha[:] * (np.dot(bead.dCVdAlpha, bead.Grad) / slen2)

    for(int i=0; i < NumOfCVs; i++){
        double cvder = BeadList->CVSplines[i]->GetCVFirstDer(Alpha);
        pMF[i] = cMF[i] - cvder * dot / slen2;
    }
        

// get final force
    for(int i=0; i < NumOfCVs; i++){
        switch(BeadType){
            case(BTY_FREE):
                // steepest descent movement
                uMF[i] = cMF[i];
            break;
            case(BTY_NORMAL):
                // movement perpendicular to path
                uMF[i] = pMF[i];
            break;
            case(BTY_PERMANENT):
            default:
                // steepest descent movement
                uMF[i] = 0.0;
            break;   
        }
    }
}

// -----------------------------------------------------------------------------

void CBead::UpdatePositionGD(double step)
{
    NumOfUpdates++;

    if( BeadType == BTY_PERMANENT ){
        for(int i=0; i < NumOfCVs; i++){
            NPos[i] = Pos[i];
        }
        return;
    }

    for(int i=0; i < NumOfCVs; i++){
        double maxmov = BeadList->CVs[i]->GetMaxMovement();
        if( (maxmov <= 0) || (fabs(uMF[i]*step) < maxmov) ){
            NPos[i] = Pos[i] - uMF[i]*step;
        } else {
            NPos[i] = Pos[i] - maxmov*sgn(uMF[i]*step);
        }
    }
}

// -----------------------------------------------------------------------------

void CBead::UpdatePositionNGD(double step,double mingnormeps)
{
    double g2 = 0.0;

    for(int i=0; i < NumOfCVs; i++){
        g2 = g2 + uMF[i]*uMF[i];
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
        g2 = g2 + uMF[i]*uMF[i];
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

    if( BeadType == BTY_PERMANENT ){
        for(int i=0; i < NumOfCVs; i++){
            NPos[i] = Pos[i];
        }
        return;
    }

    for(int i=0; i < NumOfCVs; i++){
        mtnew[i] = beta1 * mtold[i] + (1.0 - beta1) * uMF[i];
        vtnew[i] = beta2 * vtold[i] + (1.0 - beta2) * uMF[i]*uMF[i];
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

void CBead::UpdatePositionADABelief(double step,double beta1,double beta2,double mingnormeps)
{
    NumOfUpdates++;

    if( BeadType == BTY_PERMANENT ){
        for(int i=0; i < NumOfCVs; i++){
            NPos[i] = Pos[i];
        }
        return;
    }

    for(int i=0; i < NumOfCVs; i++){
        mtnew[i] = beta1 * mtold[i] + (1.0 - beta1) * uMF[i];
        vtnew[i] = beta2 * vtold[i] + (1.0 - beta2) * ((uMF[i] - mtnew[i])*(uMF[i] - mtnew[i]) + mingnormeps);
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

    if( BeadType == BTY_PERMANENT ){
        for(int i=0; i < NumOfCVs; i++){
            NPos[i] = Pos[i];
        }
        return;
    }

    for(int i=0; i < NumOfCVs; i++){
        mtnew[i]    = beta1 * mtold[i] + (1.0 - beta1) * uMF[i];
        vtnew[i]    = beta2 * vtold[i] + (1.0 - beta2) * uMF[i]*uMF[i];
        vthatnew[i] = std::max(vthatold[i],vtnew[i]);
    //    std::cout << "uMF: " << uMF[i] << std::endl;
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

    if( BeadType == BTY_PERMANENT ){
        for(int i=0; i < NumOfCVs; i++){
            NPos[i] = Pos[i];
        }
        return;
    }

    for(int i=0; i < NumOfCVs; i++){
        mtnew[i] = beta1 * mtold[i] + (1.0 - beta1) * uMF[i];
        vtnew[i] = beta2 * vtold[i] + (1.0 - beta2) * uMF[i]*uMF[i];
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
    result &= p_ele->GetAttribute("btype",BeadType);
    result &= p_ele->GetAttribute("mode",Mode);
    result &= p_ele->GetAttribute("status",ModeStatus);
    result &= p_ele->GetAttribute("nupd",NumOfUpdates);

    if(result == false) {
        LOGIC_ERROR("unable to read some attributes");
    }

    // load bead

    // here we use scaled positions
    CXMLBinData* p_posele = p_ele->GetFirstChildBinData("S-POS");
    if(p_posele == NULL) {
        LOGIC_ERROR("unable to open POS element");
    }
    Pos.Load(p_posele);

    // here we use scaled positions
    CXMLBinData* p_rposele = p_ele->GetFirstChildBinData("S-OPOS");
    if(p_rposele == NULL) {
        LOGIC_ERROR("unable to open OPOS element");
    }
    OPos.Load(p_rposele);


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
}

//------------------------------------------------------------------------------

void CBead::SaveInfo(CXMLElement* p_ele)
{
    if(p_ele == NULL) {
        INVALID_ARGUMENT("p_ele is NULL");
    }

    p_ele->SetAttribute("bead_id",BeadID);
    p_ele->SetAttribute("client_id",ClientID);
    p_ele->SetAttribute("btype",BeadType);
    p_ele->SetAttribute("mode",Mode);
    p_ele->SetAttribute("status",ModeStatus);
    p_ele->SetAttribute("nupd",NumOfUpdates);

    // save bead

    // here we use scaled positions
    CXMLBinData* p_posele = p_ele->CreateChildBinData("S-POS");
    Pos.Save(p_posele);

    // here we use scaled positions
    CXMLBinData* p_rposele = p_ele->CreateChildBinData("S-OPOS");
    OPos.Save(p_rposele);

    CXMLBinData* p_pmfele = p_ele->CreateChildBinData("MF");
    MF.Save(p_pmfele);

    CXMLBinData* p_mtzele = p_ele->CreateChildBinData("MTZ");
    MTZ.Save(p_mtzele);
}

//------------------------------------------------------------------------------

void CBead::GetProductionData(CXMLElement* p_ele)
{
    if(p_ele == NULL) {
        INVALID_ARGUMENT("p_ele is NULL");
    }

    // only if accumulation or production
    if( (GetMode() != BMO_ACCUMULATION) && (GetMode() != BMO_PRODUCTION) ) return;

    // load unscaled BPos value
    CXMLBinData* p_bposele = p_ele->GetFirstChildBinData("BPOS");
    if(p_bposele == NULL) {
        LOGIC_ERROR("unable to open BPOS element");
    }
    BPos.Load(p_bposele);

    // convert to scaled
    for(int i=0; i < NumOfCVs; i++){
        Pos[i] = BeadList->CVs[i]->GetScaledValue(BPos[i]);
    }

    // load unscaled data
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

void CBead::SetWaitForRendezvous(void)
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

    // convert to unscaled
    for(int i=0; i < NumOfCVs; i++){
        BPos[i] = BeadList->CVs[i]->GetUnscaledValue(FPos[i]); // use final position
    }

    // and bead position
    CXMLBinData* p_bposele = p_ele->CreateChildBinData("BPOS");
    BPos.Save(p_bposele);   // use final position

    ModeStatus = BMS_RUNNING;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================


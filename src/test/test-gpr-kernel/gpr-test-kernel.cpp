#include <stdio.h>
#include <SimpleVector.hpp>
#include <GPRKernel.hpp>
#include <boost/format.hpp>

//------------------------------------------------------------------------------

using namespace std;
using namespace boost;

//------------------------------------------------------------------------------

void SetRandomPos(CSimpleVector<double>& vec);
void SetRandomPositive(CSimpleVector<double>& vec);
void TestKDer(CSimpleVector<double>& akder,CSimpleVector<double>& nkder);
void TestKBlock(CFortranMatrix& akblock,CFortranMatrix& nkblock);

void TestSE(void);
void TestRQ(void);

//------------------------------------------------------------------------------

CSimpleVector<double> ipos;
CSimpleVector<double> jpos;
CSimpleVector<double> wfac;
CSimpleVector<double> nkder;
CSimpleVector<double> akder;
CFortranMatrix nkblock;
CFortranMatrix akblock;

CPMFAccumulatorPtr accu;

int nfails = 0;
double trh = 1e-3;

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int main(void)
{
    accu = CPMFAccumulatorPtr(new CPMFAccumulator);
    accu->SetNumOfCVs(2);
    accu->SetCV(0,"d1","DIS",-1.0,1.0,100);
    accu->SetCV(1,"t1","DIH",10.0,20.0,20);

    ipos.CreateVector(accu->GetNumOfCVs());
    jpos.CreateVector(accu->GetNumOfCVs());
    wfac.CreateVector(accu->GetNumOfCVs());
    nkder.CreateVector(accu->GetNumOfCVs());
    akder.CreateVector(accu->GetNumOfCVs());

    nkblock.CreateMatrix(accu->GetNumOfCVs(),accu->GetNumOfCVs());
    akblock.CreateMatrix(accu->GetNumOfCVs(),accu->GetNumOfCVs());

    for(int k=0; k < 100; k++){

    SetRandomPos(ipos);
    cout << "# I positions:" << endl;
    for(size_t i=0; i < ipos.GetLength(); i++){
        cout << format("%15.7f ")%ipos[i];
    }
    cout << endl;

    SetRandomPos(jpos);
    cout << "# J positions:" << endl;
    for(size_t i=0; i < jpos.GetLength(); i++){
        cout << format("%15.7f ")%jpos[i];
    }
    cout << endl;

    SetRandomPositive(wfac);
    cout << "# WFac:" << endl;
    for(size_t i=0; i < wfac.GetLength(); i++){
        cout << format("%15.7f ")%wfac[i];
    }
    cout << endl;

    TestSE();
    TestRQ();

    }

    cout << endl;
    cout << "# number of fails: " <<  nfails <<  endl;
}

//------------------------------------------------------------------------------

void TestSE(void)
{
    cout << endl;
    cout << "# ARDSE" << endl;

    CGPRKernel kernel;
    kernel.SetKernel("ardse");
    kernel.SetAccumulator(accu);
    kernel.SetWFac(wfac);
    kernel.SetupKernel();

    cout << "<<< I" << endl;
    kernel.GetKernelDerIAna(ipos,jpos,akder);
    kernel.GetKernelDerINum(ipos,jpos,nkder);

    TestKDer(akder,nkder);

    cout << "<<< J" << endl;
    kernel.GetKernelDerJAna(ipos,jpos,akder);
    kernel.GetKernelDerJNum(ipos,jpos,nkder);

    TestKDer(akder,nkder);

    cout << "<<< IJ" << endl;
    kernel.GetKernelDerIJAna(ipos,jpos,akblock);
    kernel.GetKernelDerIJNum(ipos,jpos,nkblock);

    TestKBlock(akblock,nkblock);

    cout << "<<< wfac" << endl;
    for(int i=0; i < accu->GetNumOfCVs(); i++){
        double ana = kernel.GetKernelValueWFacDerAna(ipos,jpos,i);
        double num = kernel.GetKernelValueWFacDerAna(ipos,jpos,i);
        double diff = fabs(ana-num);
        cout << format("%15.7f %15.7f %15.7f")%ana%num%diff;
        if( diff > trh ){
            cout << " FAIL";
            nfails++;
        } else {
            cout << " OK";
        }
        cout << endl;
    }

    cout << "<<< wfac I" << endl;
    for(int i=0; i < accu->GetNumOfCVs(); i++){
        kernel.GetKernelDerIWFacDerAna(ipos,jpos,i,akder);
        kernel.GetKernelDerIWFacDerNum(ipos,jpos,i,nkder);

        TestKDer(akder,nkder);
    }

    cout << "<<< wfac J" << endl;
    for(int i=0; i < accu->GetNumOfCVs(); i++){
        kernel.GetKernelDerJWFacDerAna(ipos,jpos,i,akder);
        kernel.GetKernelDerJWFacDerNum(ipos,jpos,i,nkder);

        TestKDer(akder,nkder);
    }

    cout << "<<< wfac IJ" << endl;
    for(int i=0; i < accu->GetNumOfCVs(); i++){
        kernel.GetKernelDerIJWFacDerAna(ipos,jpos,i,akblock);
        kernel.GetKernelDerIJWFacDerNum(ipos,jpos,i,nkblock);

        TestKBlock(akblock,nkblock);
    }
}

//------------------------------------------------------------------------------

void TestRQ(void)
{
    cout << endl;
    cout << "# ARDRQ" << endl;

    CGPRKernel kernel;
    kernel.SetKernel("ardrq");
    kernel.SetAccumulator(accu);
    kernel.SetWFac(wfac);
    kernel.SetupKernel();

    cout << "<<< I" << endl;
    kernel.GetKernelDerIAna(ipos,jpos,akder);
    kernel.GetKernelDerINum(ipos,jpos,nkder);

    TestKDer(akder,nkder);

    cout << "<<< J" << endl;
    kernel.GetKernelDerJAna(ipos,jpos,akder);
    kernel.GetKernelDerJNum(ipos,jpos,nkder);

    TestKDer(akder,nkder);

    cout << "<<< IJ" << endl;
    kernel.GetKernelDerIJAna(ipos,jpos,akblock);
    kernel.GetKernelDerIJNum(ipos,jpos,nkblock);

    TestKBlock(akblock,nkblock);

    cout << "<<< wfac" << endl;
    for(int i=0; i < accu->GetNumOfCVs(); i++){
        double ana = kernel.GetKernelValueWFacDerAna(ipos,jpos,i);
        double num = kernel.GetKernelValueWFacDerAna(ipos,jpos,i);
        double diff = fabs(ana-num);
        cout << format("%15.7f %15.7f %15.7f")%ana%num%diff;
        if( diff > trh ){
            cout << " FAIL";
            nfails++;
        } else {
            cout << " OK";
        }
        cout << endl;
    }

    cout << "<<< wfac I" << endl;
    for(int i=0; i < accu->GetNumOfCVs(); i++){
        kernel.GetKernelDerIWFacDerAna(ipos,jpos,i,akder);
        kernel.GetKernelDerIWFacDerNum(ipos,jpos,i,nkder);

        TestKDer(akder,nkder);
    }

    cout << "<<< wfac J" << endl;
    for(int i=0; i < accu->GetNumOfCVs(); i++){
        kernel.GetKernelDerJWFacDerAna(ipos,jpos,i,akder);
        kernel.GetKernelDerJWFacDerNum(ipos,jpos,i,nkder);

        TestKDer(akder,nkder);
    }

    cout << "<<< wfac IJ" << endl;
    for(int i=0; i < accu->GetNumOfCVs(); i++){
        kernel.GetKernelDerIJWFacDerAna(ipos,jpos,i,akblock);
        kernel.GetKernelDerIJWFacDerNum(ipos,jpos,i,nkblock);

        TestKBlock(akblock,nkblock);
    }
}

//------------------------------------------------------------------------------

void TestKDer(CSimpleVector<double>& akder,CSimpleVector<double>& nkder)
{
    for(size_t i=0; i < akder.GetLength(); i++){
        double diff = fabs(akder[i]-nkder[i]);
        cout << format("%15.7f %15.7f %15.7f")%akder[i]%nkder[i]%diff;
        if( diff > trh ){
            cout << " FAIL";
            nfails++;
        } else {
            cout << " OK";
        }
        cout << endl;
    }
}

//------------------------------------------------------------------------------

void TestKBlock(CFortranMatrix& akblock,CFortranMatrix& nkblock)
{
    for(size_t i=0; i < akblock.GetNumberOfRows(); i++){
        for(size_t j=0; j < akblock.GetNumberOfColumns(); j++){
            double diff = fabs(akblock[i][j]-nkblock[i][j]);
            cout << format("%15.7f %15.7f %15.7f")%akblock[i][j]%nkblock[i][j]%diff;
            if( diff > trh ){
                cout << " FAIL";
                nfails++;
            } else {
                cout << " OK";
            }
            cout << endl;
        }
    }
}

//------------------------------------------------------------------------------

void SetRandomPos(CSimpleVector<double>& vec)
{
   for(int i=0; i < accu->GetNumOfCVs(); i++){
        double r = (double)rand() / (double)RAND_MAX;
        CColVariablePtr cv = accu->GetCV(i);
        double v = r*cv->GetRange() + cv->GetMinValue();
        vec[i] = v;
   }
}

//------------------------------------------------------------------------------

void SetRandomPositive(CSimpleVector<double>& vec)
{
   for(size_t i=0; i < vec.GetLength(); i++){
        double f = (double)rand() / (double)RAND_MAX;
        double val = 20.0 * f + 0.1;
        vec[i] = val;
   }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

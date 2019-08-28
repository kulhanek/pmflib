#include <stdio.h>
#include <vector>
#include <iostream>

using namespace std;

// -----------------------------------------------------------------

int NCoords;
int NBins[10];
int TotNBins;

// -----------------------------------------------------------------

void GetIPoint(unsigned int index,std::vector<int>& point)
{
    for(int k=NCoords-1; k >= 0; k--) {
        int ibin = index % NBins[k];
        point[k] = ibin;
        index = index / NBins[k];
    }
}

// -----------------------------------------------------------------

int GetFBinIndex(const std::vector<int>& position)
{
    int glbindex = 0;
    for(int i=0; i < NCoords; i++) {
        int nbins = NBins[i];
        int pos = position[i];
        glbindex = glbindex*nbins + pos;
    }

    return(glbindex);
}

// -----------------------------------------------------------------

int main(void)
{
    NCoords = 2;
    NBins[0] = 140;
    NBins[1] = 33;
    TotNBins = NBins[0]*NBins[1];

    std::vector<int>    ipoint;
    ipoint.resize(NCoords);

    for(int glbidx = 0; glbidx < TotNBins; glbidx++){
        GetIPoint(glbidx,ipoint);
        int rglbindx = GetFBinIndex(ipoint);
        cout << "glbidx = " << glbidx << " cv1 = " << ipoint[0] << " cv2 = " << ipoint[1] << endl;
        if( glbidx != rglbindx ){
            cout << "index violation " << glbidx << " " << rglbindx << endl;
            return(1);
        }
    }

    return(0);
}


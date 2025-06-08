
#include <iostream>
#include <CVSplineInterpolatingCubic.hpp>
#include <CVSplineSmoothingCubic.hpp>
#include <vector>

int main() {
    std::vector<double> x = {0, 1, 2, 3, 4, 5};
    std::vector<double> y = {0, 1.2, 0.9, 2.5, 2.0, 3.1};
//    double lambda = 1.0;  // adjust to control smoothness


    for(size_t i=0; i < x.size(); i++ ){
        printf("%10.2f %10.5f\n",x[i],y[i]);
    }

    printf("\n");

// natural cubic spline
    CCVSplineInterpolatingCubic spline_nc;
    spline_nc.Allocate(x.size());
    for(size_t i=0; i< x.size(); i++){
        spline_nc.SetPoint(i,x[i],y[i]);
    }
    spline_nc.BuildSpline();

    for (double xi = 0; xi <= 5.0; xi += 0.1) {
        double yi = spline_nc.GetCV(xi);
        printf("%10.2f %10.5f\n",xi,yi);
    }

    printf("\n");

// smoothing cubic spline
    CCVSplineSmoothingCubic spline_sc;
    spline_sc.Allocate(x.size());
    for(size_t i=0; i< x.size(); i++){
        spline_sc.SetPoint(i,x[i],y[i]);
    }
    spline_sc.BuildSpline();

    for (double xi = 0; xi <= 5.0; xi += 0.1) {
        double yi = spline_sc.GetCV(xi);
        printf("%10.2f %10.5f\n",xi,yi);
    }

    return 0;
}

#ifndef BIN_H
#define BIN_H

#include "Wrapper.hpp"
#include <vector>

using namespace std;

#define BINSIZE 16
#define DENSITYUPPER 0.97
#define BETA 0.0015
#define BINNUM 10

class Bin
{
public:
    Bin(){};
    Bin(double x, double y, unsigned int size, double width, double height);
    double x();
    double y();

    void computeO(wrapper::Placement &placement);

    vector <long double> _Ox, _Oy;
    long double density;
    vector<size_t> moduleIdx;
private:
    // center position
    double _x;
    double _y;
    long double binWidth, binHeight;
};

class Bins
{
public:
    Bins(){};
    Bins(unsigned int numModules, vector<double> boundry);

    void computeBinO(wrapper::Placement &placement);
    void computeC(wrapper::Placement &placement);
    void computeBinDensity(wrapper::Placement &placement);
    void computeGrad(wrapper::Placement &placement);

    unsigned int binNum();
    long double binWidth();
    long double binHeight();

    double maxDensity();

    void showDensity();
    void showO(wrapper::Placement &placement);
    void showC(wrapper::Placement &placement);

    void updateBins(wrapper::Placement &placement);

    unsigned int iterNum = 0, iter = 0;
    vector<double> C;
    vector<long double> gradX, gradY;
    double beta = BETA;
    int cnt = 0;
private:
    vector<vector<Bin>> b;
    long double chipWidth, chipHeight;
    long double _binWidth, _binHeight;

    pair<pair<size_t, size_t>,pair<size_t, size_t>> findModuleBin(wrapper::Placement &placement, size_t i);
};

#endif // BIN_H
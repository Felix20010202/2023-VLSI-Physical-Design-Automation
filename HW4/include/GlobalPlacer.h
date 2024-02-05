#ifndef GLOBALPLACER_H
#define GLOBALPLACER_H

#include <vector>
#include "Wrapper.hpp"
#include "LSE.h"
#include "Bin.h"

#define LEARNRATE 0.01

class GlobalPlacer
{
public:
    GlobalPlacer(wrapper::Placement &placement);

    void randomPlace(); // An example of random placement implemented by TA
    void place();

    // LSE
    void getGrad();
    void updatePos();

    void initPlace();
    void LSEPlace();
    void LSEdifferentiate(const vector<long double> &x, long double &g);
    void LSEevaluateFG(const vector<long double> &x, long double &f, vector<long double> &g);
    void LSEevaluateF(const vector<long double> &x, long double &f);
    unsigned LSEdimension();
    void setDimension();

    void testBest(double density, double hpwl);
    void restoreBest();
    


    // show data

    double alpha = 0;
    unsigned int iterNum = 10000;
    unsigned int iter = 0;

    vector<long double> historyX;
    vector<long double> historyY;
    vector<long double> GradientX;
    vector<long double> GradientY;
private:
    wrapper::Placement &_placement;
    Bins _bin;
    bool best = false;
    double bestDensity = 100;
    double bestHpwl = 0;
    vector<double> bestX;
    vector<double> bestY;
};

#endif // GLOBALPLACER_H

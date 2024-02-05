#ifndef TARGETFUNCTION_H
#define TARGETFUNCTION_H

#include "NumericalOptimizerInterface.h"
#include <cmath>

class TargetFunction : public NumericalOptimizerInterface
{
public:
    TargetFunction();

    void evaluateFG(const vector<double> &x, double &f, vector<double> &g);
    void evaluateF(const vector<double> &x, double &f);
    unsigned int dimension();
    
    void setDimension(unsigned int dimension);

    unsigned int _dimension = 0;
};

#endif // TARGETFUNCTION_H

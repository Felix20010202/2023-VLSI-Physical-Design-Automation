#ifndef LSE_H
#define LSE_H

#include "NumericalOptimizerInterface.h"
#include <cmath>

# define LSECOEF 0.0001

class LSE : public NumericalOptimizerInterface
{
public:

    int _netNum = 0;

    // void differentiate(const vector<long double> &x, const vector<long double> &y, long double &g);
    static void evaluateFG(const vector<long double> &x, long double &f, vector<long double> &g);
    static void evaluateF(const vector<long double> &x, long double &f);
    void setDimension(unsigned d);

private:
    unsigned _dimension = 0;
};

#endif // LSE_H

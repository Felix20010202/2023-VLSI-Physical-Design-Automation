#include "LSE.h"
#include <iostream>

using namespace std;

// LSE::LSE(wrapper::Placement &placement)
//     : _placement(placement)
// {
// }

void LSE::evaluateFG(const vector<long double> &x, long double &f, vector<long double> &g)
{
    size_t dim = x.size();
    long double y1 = 0, y2 = 0;
    LSE::evaluateF(x, f);

    for(size_t i = 0; i < dim; i++){
        y1 += exp(x[i]/LSECOEF);
        y2 += exp(-1*x[i]/LSECOEF);
    }

    g.resize(dim);
    for(size_t i = 0; i < dim; i++){
        g[i] = (exp(x[i]/LSECOEF)/LSECOEF)/y1 - (exp(-1*x[i]/LSECOEF)/LSECOEF)/y2;
        g[i] *= LSECOEF;
    }
}

void LSE::evaluateF(const vector<long double> &x, long double &f)
{
    size_t dim = x.size();

    f = 0;
    long double temp = 0;
    for(size_t i = 0; i < dim; i++){
        temp += exp(x[i]/LSECOEF);
    }
    f += log(temp);
    temp = 0;
    for(size_t i = 0; i < dim; i++){
        temp += exp(-1*x[i]/LSECOEF);
    }
    f += log(temp);
    f = LSECOEF*f;
}

void LSE::setDimension(unsigned d){
    _dimension = d;
}

// vector<vector<vector>> LSE;
// net >> xy >> 1~n
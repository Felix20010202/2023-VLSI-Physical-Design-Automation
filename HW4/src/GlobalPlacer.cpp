#include "GlobalPlacer.h"
#include "ExampleFunction.h"
#include "NumericalOptimizer.h"
#include "LSE.h"
#include "Bin.h"
#include <cstdlib>
#include <iostream>
#include <unistd.h>
#include<climits>
#include<signal.h>
#include<unistd.h>

using namespace std;

bool timeOut = false;
void sigalrm_handler(int)
{
    timeOut = true;
}

GlobalPlacer::GlobalPlacer(wrapper::Placement &placement)
    : _placement(placement), _bin()
{
    historyX = vector<long double>(placement.numModules(), 0);
    historyY = vector<long double>(placement.numModules(), 0);

    // set process time
    signal(SIGALRM, &sigalrm_handler);
    alarm(590);
}

void GlobalPlacer::randomPlace()
{
    srand(7777);
    double coreWidth = _placement.boundryRight() - _placement.boundryLeft();
    double coreHeight = _placement.boundryTop() - _placement.boundryBottom();
    for (size_t i = 0; i < _placement.numModules(); ++i)
    {
        if (_placement.module(i).isFixed())
            continue;

        double width = _placement.module(i).width();
        double height = _placement.module(i).height();
        double x = rand() % static_cast<int>(coreWidth - width) + _placement.boundryLeft();
        double y = rand() % static_cast<int>(coreHeight - height) + _placement.boundryBottom();
        _placement.module(i).setPosition(x, y);
    }
}

void GlobalPlacer::getGrad()
{
    double width = _placement.boundryRight() - _placement.boundryLeft();
    double heigth = _placement.boundryTop() - _placement.boundryBottom();
    unsigned int moduleNum = _placement.numModules();
    GradientX = vector<long double>(moduleNum, 0);
    GradientY = vector<long double>(moduleNum, 0);
    
    for (size_t i = 0; i < _placement.numNets(); ++i)
    {
        vector<long double> x, y;
        long double maxX = INT_MIN, minX = INT_MAX, maxY = INT_MIN, minY = INT_MAX;
        long double fx, fy;
        vector<long double> gx, gy;
        for(size_t j = 0; j < _placement.net(i).numPins(); ++j){
            double _x = _placement.module(_placement.net(i).pin(j).moduleId()).x();
            double _y = _placement.module(_placement.net(i).pin(j).moduleId()).y();
            if(_x > maxX) maxX = _x;
            if(_x < minX) minX = _x;
            if(_y > maxY) maxY = _y;
            if(_y < minY) minY = _y;
            x.push_back((_x-_placement.boundryLeft())/width);
            y.push_back((_y-_placement.boundryBottom())/heigth);
        
        }
        // LSE::evaluateF(x, fx);
        // LSE::evaluateF(y, fy);
        LSE::evaluateFG(x, fx, gx);
        LSE::evaluateFG(y, fy, gy);
        
        for(size_t j = 0; j < _placement.net(i).numPins(); ++j){
            GradientX[_placement.net(i).pin(j).moduleId()] += gx[j];
            GradientY[_placement.net(i).pin(j).moduleId()] += gy[j];
        }
        
    }

    _bin.computeBinO(_placement);
    _bin.computeC(_placement);
    _bin.computeBinDensity(_placement);
    _bin.computeGrad(_placement);

    double maxGX = 0, maxGY = 0;
    for(size_t i = 0; i < moduleNum; ++i){
        if(abs(GradientX[i]) > maxGX) maxGX = abs(GradientX[i]);
        if(abs(GradientY[i]) > maxGY) maxGY = abs(GradientY[i]);
    }

    for(size_t i = 0; i < moduleNum; ++i)
    {
        GradientX[i] /= maxGX;
        GradientY[i] /= maxGY;
        // GradientX[i] /= _placement.module(i).numPins();
        GradientX[i] *= width;
        // GradientY[i] /= _placement.module(i).numPins();
        GradientY[i] *= heigth;
    }
}

void GlobalPlacer::updatePos(){
    alpha = min(0.3 - 0.00012*iter, 0.25);
    // if(alpha < 0.15) alpha = 0.15;
    // alpha = 0.5;
    // alpha = 0.00 + (0.3/(1+exp(0.03*((double)iterNum/2-iter)))) + (0.3/(1+exp(-0.03*((double)iterNum/4-iter))));
    // alpha = 0.1+(0.2/(1+exp(0.03*((double)iterNum/2-iter))));
    // cout<<"alpha: "<<alpha<<endl;
    double chipWidth = _placement.boundryRight() - _placement.boundryLeft();
    double chipHeigth = _placement.boundryTop() - _placement.boundryBottom();
    for (size_t i = 0; i < _placement.numModules(); ++i){
        if(_placement.module(i).isFixed() == true){
            continue;
        }

        double GX = (alpha*GradientX[i] + _bin.gradX[i]*chipWidth)*0.8 + historyX[i]*0.2;
        double GY = (alpha*GradientY[i] + _bin.gradY[i]*chipHeigth)*0.8 + historyY[i]*0.2;
        historyX[i] = GX;
        historyY[i] = GY;
        double newX = _placement.module(i).x() - GX;
        double newY = _placement.module(i).y() - GY;

        double mw = _placement.module(i).area()/_placement.module(i).height();
        double mh = _placement.module(i).height();

        if(newX + mw > _placement.boundryRight()){
            // newX = _placement.boundryRight() - (_bin.binWidth()/2) - mw;
            // newX = (_placement.boundryRight()+_placement.boundryLeft())/2;
            // newX = _placement.module(i).x() - 0.3*GX;
            // if(newX + mw > _placement.boundryRight()){
            //     newX = _placement.module(i).x();
            // }
            // newX = _placement.module(i).x();
            newX = ((_placement.boundryRight()+_placement.boundryLeft())/2 + _placement.module(i).x())/2;
        }
        else if(newX < _placement.boundryLeft()){
            // newX = _placement.boundryLeft() + (_bin.binWidth()/2);
            // newX = (_placement.boundryRight()+_placement.boundryLeft())/2;
            // newX = _placement.module(i).x() - 0.3*GX;
            // if(newX < _placement.boundryLeft()){
            //     newX = _placement.module(i).x();
            // }
            // newX = _placement.module(i).x();
            newX = ((_placement.boundryRight()+_placement.boundryLeft())/2 + _placement.module(i).x())/2;
        }

        if(newY + mh > _placement.boundryTop()){
            // newY = _placement.boundryTop() - (_bin.binHeight()/2) - mh;
            // newY = (_placement.boundryBottom()+_placement.boundryTop())/2;
            // newY = _placement.module(i).y() - 0.3*GY;
            // if(newY + mh > _placement.boundryTop()){
            //     newY = _placement.module(i).y();
            // }
            // newY = _placement.module(i).y();
            newY = ((_placement.boundryBottom()+_placement.boundryTop())/2 + _placement.module(i).y())/2;
        }
        else if(newY < _placement.boundryBottom()){
            // newY = _placement.boundryBottom() + (_bin.binHeight()/2);
            // newY = (_placement.boundryBottom()+_placement.boundryTop())/2;
            // newY = _placement.module(i).y() - 0.3*GY;
            // if(newY < _placement.boundryBottom()){
            //     newY = _placement.module(i).y();
            // }
            // newY = _placement.module(i).y();
            newY = ((_placement.boundryBottom()+_placement.boundryTop())/2 + _placement.module(i).y())/2;
        }

        if(_placement.module(i).isFixed() == false){
            _placement.module(i).setPosition(newX, newY);
        }
    }
}

void GlobalPlacer::initPlace()
{
    // double centerX = (_placement.boundryLeft() + _placement.boundryRight())/2;
    // double centerY = (_placement.boundryBottom() + _placement.boundryTop())/2;
    // double hbw = _bin.binWidth()/2, hbh = _bin.binHeight()/2;

    // for(size_t i = 0; i < _placement.numModules(); ++i)
    // {
    //     if (_placement.module(i).isFixed())
    //         continue;

    //     double x = hbw*((double) rand() / (RAND_MAX));
    //     double y = hbh*((double) rand() / (RAND_MAX));
    //     if(rand() % 2 == 1) x *= -1;
    //     if(rand() % 2 == 1) y *= -1;
    //     double biasX = (_placement.module(i).area()/_placement.module(i).height())/2;
    //     double biasY = _placement.module(i).height()/2;
    //     _placement.module(i).setPosition(centerX  - biasX, centerY  - biasY);
    // }
    
    // exit(1);
    vector<double> boundry;
    boundry.push_back(_placement.boundryLeft());
    boundry.push_back(_placement.boundryRight());
    boundry.push_back(_placement.boundryBottom());
    boundry.push_back(_placement.boundryTop());
    _bin = Bins(_placement.numModules(), boundry);
    
}

void GlobalPlacer::LSEPlace(){
    initPlace();

    _bin.iterNum = iterNum;
    for(iter = 0; iter < iterNum; iter++){
        if(timeOut == true){
            break;
        }
        getGrad();
        updatePos();
        // _bin.showDensity();
        _bin.iter++;
        double nowHpwl = _placement.computeHpwl();
        // cout<<iter<<" HPWL: "<<nowHpwl<<endl;
        double max = _bin.maxDensity();
        // cout<<max<<endl;
        if(iter > 30){
            testBest(max, nowHpwl);
        }
    }

    restoreBest();
    cout<<_bin.maxDensity()<<" "<<_placement.computeHpwl()<<endl;
}

void GlobalPlacer::testBest(double density, double hpwl){
    if(best == false){
        best = true;
        bestDensity = density;
        bestHpwl = hpwl;
    }
    else{
        if(density > bestDensity){
            if(density >= 1.3){
                return;
            }
            // return;
        }
    }
    if(density - 0.01 > bestDensity || (hpwl > bestHpwl && bestDensity < 1.3)){
        return;
    }
    cout<<iter<<" Update best, max density: "<<density<<", HPWL: "<<hpwl<<endl;
    bestDensity = density;
    bestHpwl = hpwl;
    bestX.resize(_placement.numModules());
    bestY.resize(_placement.numModules());
    for(size_t i = 0; i < _placement.numModules(); i++){
        bestX[i] = _placement.module(i).centerX();
        bestY[i] = _placement.module(i).centerY();
    }
}

void GlobalPlacer::restoreBest(){
    if(best == false){
        return;
    }
    for(size_t i = 0; i < _placement.numModules(); i++){
        if(_placement.module(i).isFixed() == true){
            continue;
        }
        _placement.module(i).setCenterPosition(bestX[i], bestY[i]);
    }
}

void GlobalPlacer::place()
{
    randomPlace(); // as init place for LSE place
    LSEPlace();

    /* @@@ TODO
     * 1. Understand above example and modify ExampleFunction.cpp to implement the analytical placement
     * 2. You can choose LSE or WA as the wirelength model, the former is easier to calculate the gradient
     * 3. For the bin density model, you could refer to the lecture notes
     * 4. You should first calculate the form of wirelength model and bin density model and the forms of their gradients ON YOUR OWN
     * 5. Replace the value of f in evaluateF() by the form like "f = alpha*WL() + beta*BinDensity()"
     * 6. Replace the form of g[] in evaluateG() by the form like "g = grad(WL()) + grad(BinDensity())"
     * 7. Set the initial vector x in place(), set step size, set #iteration, and call the solver like above example
     * */
}
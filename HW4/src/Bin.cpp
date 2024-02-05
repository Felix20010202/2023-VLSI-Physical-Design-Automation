#include "Bin.h"
#include <iostream>
#include <cmath>

using namespace std;

Bins::Bins(unsigned int numModules, vector<double> boundry)
    {
    // unsigned int binSegmentNum = sqrt(numModules) / BINSIZE;
    unsigned int binSegmentNum = BINNUM;
    chipWidth = boundry[1] - boundry[0];
    chipHeight = boundry[3] - boundry[2];

    _binWidth = chipWidth/binSegmentNum;
    _binHeight = chipHeight/binSegmentNum;
    
    double biasX = _binWidth/2, biasY = _binHeight/2;

    b.resize(binSegmentNum);
    double x = boundry[0], y = boundry[2];
    for(size_t i = 0; i < binSegmentNum; ++i){
        b[i].resize(binSegmentNum);
        for(size_t j = 0; j < binSegmentNum; ++j){
            b[i][j] = Bin(x+biasX, y+biasY, numModules, _binWidth, _binHeight);
            x += _binWidth;
        }
        x = boundry[0];
        y += _binHeight;
    }
}

void Bins::computeBinO(wrapper::Placement &placement){
    for(size_t i = 0; i < b.size(); i++){
        for(size_t j = 0; j < b.size(); j++){
            b[i][j].computeO(placement);
        }
    }
}

void Bins::computeC(wrapper::Placement &placement){
    unsigned int size = placement.numModules();
    C.resize(size);
    for(size_t i = 0; i < size; ++i){
        long double temp = 0;
        for(size_t j = 0; j < b.size(); ++j){
            for(size_t k = 0; k < b.size(); ++k){
                temp += (b[j][k]._Ox[i])*(b[j][k]._Oy[i]);
            }
        }
        C[i] = placement.module(i).area()/temp;
        // cout<<temp<<endl;
    }
}

void Bins::computeBinDensity(wrapper::Placement &placement){
    for(size_t i = 0; i < b.size(); ++i){
        for(size_t j = 0; j < b.size(); ++j){
            long double temp = 0;
            for(size_t k = 0; k < placement.numModules(); ++k){
                temp += C[k] * b[i][j]._Ox[k] * b[i][j]._Oy[k];
            }
            b[i][j].density = temp;
            // cout<<b[i][j].density/(binWidth*binHeight)<<" ";
            // cout<<temp<<"   ";
        }
        // cout<<endl;
    }
}

void Bins::computeGrad(wrapper::Placement &placement){

    updateBins(placement);

    unsigned int size = placement.numModules();
    gradX = vector<long double>(size, 0);
    gradY = vector<long double>(size, 0);

    long double maxX = 0;
    long double maxY = 0;

    vector<vector<long double>> C2Temp;
    C2Temp.resize(b.size());
    for(size_t i = 0; i < b.size(); i++){
        C2Temp[i].resize(b.size());
        for(size_t j = 0; j < b.size(); j++){
            C2Temp[i][j] = 0;
            // for(size_t k = 0; k < size; k++){
            //     C2Temp[i][j] += C[k] * b[i][j]._Ox[k] * b[i][j]._Oy[k];
            // }
            // long double temp = 0;
            for(size_t k = 0; k < b[i][j].moduleIdx.size(); k++){
                size_t idx = b[i][j].moduleIdx[k];
                C2Temp[i][j] += C[idx] * b[i][j]._Ox[idx] * b[i][j]._Oy[idx];
            }
        }    
    }


    // speedup
    for(size_t i = 0; i < b.size(); i++){
        for(size_t j = 0; j < b[i].size(); j++){
            for(size_t k = 0; k < b[i][j].moduleIdx.size(); k++){
                long double oneBinGrad = 0;
                size_t idx = b[i][j].moduleIdx[k];
                double xi = placement.module(idx).centerX();
                double xb = b[i][j].x();
                double wb = _binWidth;
                double wi = placement.module(idx).area()/placement.module(idx).height();
                double dx = abs(xi - xb);

                double C1 = C[idx] * b[i][j]._Oy[idx];

                long double C2 = C2Temp[i][j] - C[idx] * b[i][j]._Ox[idx] * b[i][j]._Oy[idx];
                C2 -= DENSITYUPPER*_binWidth*_binHeight;

                if(0 <= dx && dx <= wb/2 + wi/2){
                    double a = 4/((wb+wi)*(2*wb+wi));

                    // oneBinGrad = -4*a*(C1*(xi-xb)*(a*pow((xi-xb),2)-1)+C2);
                    oneBinGrad = -4*a*C1*(xi-xb)*(C1*(1-a*pow((xi-xb),2)) +C2);
                    if(xi < xb){
                        oneBinGrad *= -1;
                    }
                }
                else if(wb/2 + wi/2 <= dx && dx <= wb+wi/2){
                    double b = 4/(wb*(2*wb+wi));
                    double C3 = xb+wb+wi/2;

                    oneBinGrad = 4*b*C1*(xi-C3)*(C1*( b*pow((xi-C3),2) ) + C2);
                    if(xi < xb){
                        oneBinGrad *= -1;
                    }
                }
                else if(wb+wi/2 <= dx){
                    oneBinGrad = 0;
                }
                else{
                    cerr<<"error: "<<endl;
                    exit(1);
                }
                gradX[idx] += oneBinGrad;
            }
        }
    }

    for(size_t i = 0; i < placement.numModules(); i++){
        if(abs(gradX[i]) > maxX) maxX = abs(gradX[i]);
    }

    for(size_t i = 0; i < b.size(); i++){
        for(size_t j = 0; j < b[i].size(); j++){
            for(size_t k = 0; k < b[i][j].moduleIdx.size(); k++){
                long double oneBinGrad = 0;
                size_t idx = b[i][j].moduleIdx[k];
                double yi = placement.module(idx).centerY();
                double yb = b[i][j].y();
                double hb = _binHeight;
                double hi = placement.module(idx).height();
                double dy = abs(yi - yb);

                double C1 = C[idx] * b[i][j]._Ox[idx];

                long double C2 = C2Temp[i][j] - C[idx] * b[i][j]._Ox[idx] * b[i][j]._Oy[idx];
                C2 -= DENSITYUPPER*_binWidth*_binHeight;


                if(0 <= dy && dy <= hb/2 + hi/2){
                    double a = 4/((hb+hi)*(2*hb+hi));

                    // oneBinGrad = -4*a*(C1*(xi-xb)*(a*pow((xi-xb),2)-1)+C2);
                    oneBinGrad = -4*a*C1*(yi-yb)*(C1*(1-a*pow((yi-yb),2)) +C2);
                    if(yi < yb){
                        oneBinGrad *= -1;
                    }
                }
                else if(hb/2 + hi/2 <= dy && dy <= hb+hi/2){
                    double b = 4/(hb*(2*hb+hi));
                    double C3 = yb+hb+hi/2;

                    oneBinGrad = 4*b*C1*(yi-C3)*(C1*( b*pow((yi-C3),2) ) + C2);
                    if(yi < yb){
                        oneBinGrad *= -1;
                    }
                }
                else if(hb+hi/2 <= dy){
                    // cerr<<"error: "<<yi<<" "<<yb-hb/2<<" "<<yb+hb/2<<endl;
                    // exit(1);
                }

                gradY[idx] += oneBinGrad;
            }
        }
    }

    for(size_t i = 0; i < placement.numModules(); i++){
        if(abs(gradY[i]) > maxY) maxY = abs(gradY[i]);
    }

    // -------------

    // for(size_t i = 0; i < size; i++){
    //     long double temp = 0;
    //     for(size_t j = 0; j < b.size(); j++){
    //         for(size_t k = 0; k < b.size(); k++){
    //             long double oneBinGrad = 0;
    //             double xi = placement.module(i).centerX();
    //             double xb = b[j][k].x();
    //             double wb = _binWidth;
    //             double wi = placement.module(i).area()/placement.module(i).height();
    //             double dx = abs(xi - xb);

    //             double C1 = C[i] * b[j][k]._Oy[i];
                
    //             long double C2 = C2Temp[j][k] - C[i] * b[j][k]._Ox[i] * b[j][k]._Oy[i];
    //             C2 -= DENSITYUPPER*_binWidth*_binHeight;
                
    //             if(0 <= dx && dx <= wb/2 + wi/2){
    //                 double a = 4/((wb+wi)*(2*wb+wi));

    //                 // oneBinGrad = -4*a*(C1*(xi-xb)*(a*pow((xi-xb),2)-1)+C2);
    //                 oneBinGrad = -4*a*C1*(xi-xb)*(C1*(1-a*pow((xi-xb),2)) +C2);
    //                 if(xi < xb){
    //                     oneBinGrad *= -1;
    //                 }
    //             }
    //             else if(wb/2 + wi/2 <= dx && dx <= wb+wi/2){
    //                 double b = 4/(wb*(2*wb+wi));
    //                 double C3 = xb+wb+wi/2;

    //                 oneBinGrad = 4*b*C1*(xi-C3)*(C1*( b*pow((xi-C3),2) ) + C2);
    //                 if(xi < xb){
    //                     oneBinGrad *= -1;
    //                 }
    //             }
    //             else if(wb+wi/2 <= dx){
    //                 oneBinGrad = 0;
    //             }

    //             temp += oneBinGrad;
    //         }
    //     }

    //     if(abs(temp) > maxX) maxX = abs(temp);
    //     gradX[i] = temp;
    // }

    // for(size_t i = 0; i < size; i++){
    //     long double temp = 0;
    //     for(size_t j = 0; j < b.size(); j++){
    //         for(size_t k = 0; k < b.size(); k++){
    //             long double oneBinGrad = 0;
    //             double yi = placement.module(i).centerY();
    //             double yb = b[j][k].y();
    //             double hb = _binHeight;
    //             double hi = placement.module(i).height();
    //             double dy = abs(yi - yb);

    //             double C1 = C[i] * b[j][k]._Ox[i];
                
    //             long double C2 = C2Temp[j][k] - C[i] * b[j][k]._Ox[i] * b[j][k]._Oy[i];
    //             C2 -= DENSITYUPPER*_binWidth*_binHeight;
                
    //             if(0 <= dy && dy <= hb/2 + hi/2){
    //                 double a = 4/((hb+hi)*(2*hb+hi));

    //                 // oneBinGrad = -4*a*(C1*(xi-xb)*(a*pow((xi-xb),2)-1)+C2);
    //                 oneBinGrad = -4*a*C1*(yi-yb)*(C1*(1-a*pow((yi-yb),2)) +C2);
    //                 if(yi < yb){
    //                     oneBinGrad *= -1;
    //                 }
    //             }
    //             else if(hb/2 + hi/2 <= dy && dy <= hb+hi/2){
    //                 double b = 4/(hb*(2*hb+hi));
    //                 double C3 = yb+hb+hi/2;

    //                 oneBinGrad = 4*b*C1*(yi-C3)*(C1*( b*pow((yi-C3),2) ) + C2);
    //                 if(yi < yb){
    //                     oneBinGrad *= -1;
    //                 }
    //             }
    //             else if(hb+hi/2 <= dy){
    //                 oneBinGrad = 0;
    //             }

    //             temp += oneBinGrad;
    //         }
    //     }

    //     if(abs(temp) > maxY) maxY = abs(temp);
    //     gradY[i] = temp;

    // }


    // cout<<"Height: "<<(placement.boundryTop()-placement.boundryBottom())<<endl;

    for(size_t i = 0; i < size; i++){
        gradX[i] /= maxX;
        gradX[i] *= beta;
        gradY[i] /= maxY;
        gradY[i] *= beta;
    }

    // beta = 0.35+(0.4/(1+exp(0.001*(((double)10000/2-min(iter,(unsigned int)10000))))));
    if(iter < 20){
        beta = 0.001 + iter*0.00005;
    }
    else{
        beta = min(0.01 + (iter-20)*0.003, 0.51);
    }
    // cout<<"beta: "<<beta<<endl;
}

unsigned int Bins::binNum(){
    // return pow(b.size(), 2);
    return BINNUM*BINNUM;
}

void Bins::showDensity(){
    long double max = 0;
    cout<<endl;
    for(size_t i = 0; i < b.size(); i++){
        for(size_t j = 0; j < b.size(); j++){
            if(b[i][j].density/(_binWidth*_binHeight) > max){
                max = b[i][j].density/(_binWidth*_binHeight);
            }
            cout<<b[i][j].density/(_binWidth*_binHeight)<<" ";
        }
        cout<<endl;
    }
    cout<<"Max density value: "<<max<<endl;
}

double Bins::maxDensity(){
    double max = 0;
    for(size_t i = 0; i < b.size(); i++){
        for(size_t j = 0; j < b.size(); j++){
            if(b[i][j].density/(_binWidth*_binHeight) > max){
                max = b[i][j].density/(_binWidth*_binHeight);
            }
        }
    }
    return max;
}

Bin::Bin(double x, double y, unsigned int size, double width, double height){
    _x = x;
    _y = y;
    _Ox.resize(size);
    _Oy.resize(size);
    binWidth = width;
    binHeight = height;
}

double Bin::x(){
    return _x;
}

double Bin::y(){
    return _y;
}

void Bin::computeO(wrapper::Placement &placement){
    unsigned int size = placement.numModules();
    
    _Ox.resize(size);
    _Oy.resize(size);

    // Ox
    for(size_t i = 0; i < size; ++i)
    {
        double dx = abs(placement.module(i).centerX() - _x);
        double wi = placement.module(i).area()/placement.module(i).height();
        if(0 <= dx && dx <= binWidth/2 + wi/2){
            long double a = 4/((binWidth + wi)*(2*binWidth + wi));
            _Ox[i] = 1 - a*dx*dx;
        }
        else if((binWidth/2 + wi/2) <= dx && dx <= binWidth + wi/2){
            long double b = 4/(binWidth*(2*binWidth + wi));
            _Ox[i] = b*pow( (dx-binWidth-(wi/2)), 2);
        }
        else if(binWidth + wi/2 <= dx){
            _Ox[i] = 0;
        }
    }

    // Oy
    for(size_t i = 0; i < size; ++i)
    {
        double dy = abs(placement.module(i).centerY() - _y);
        double hi = placement.module(i).height();
        if(0 <= dy && dy <= binHeight/2 + hi/2){
            long double a = 4/((binHeight + hi)*(2*binHeight + hi));
            _Oy[i] = 1 - a*dy*dy;
        }
        else if((binHeight/2 + hi/2) <= dy && dy <= binHeight + hi/2){
            long double b = 4/(binHeight*(2*binHeight + hi));
            _Oy[i] = b*pow( (dy-binHeight-(hi/2)) , 2);
        }
        else if(binHeight + hi/2 <= dy){
            _Oy[i] = 0;
        }
    }
}

long double Bins::binWidth(){
    return _binWidth;
}

long double Bins::binHeight(){
    return _binHeight;
}

void Bins::showO(wrapper::Placement &placement){
    unsigned int size = placement.numModules();
    for(size_t i = 0; i < b.size(); i++){
        for(size_t j = 0; j < b.size(); j++){
            for(size_t k = 0; k < size; k++){
                if(b[i][j]._Ox[k] > 1){
                    cout<<b[i][j]._Ox[k]<<endl;
                }
                if(b[i][j]._Oy[k] > 1){
                    cout<<b[i][j]._Oy[k]<<endl;
                }
            }
        }
    }
}

void Bins::showC(wrapper::Placement &placement){
    for(size_t i = 0; i < C.size(); i++){
        long double temp = 0;
        for(size_t j = 0; j < b.size(); j++){
            for(size_t k = 0; k < b.size(); k++){
                temp += b[j][k]._Ox[i] * b[j][k]._Oy[i];
            }
        }
        cout<<temp*C[i]<<" "<<placement.module(i).area()<<endl;
        // cout<<C[i]<<endl;
    }
}

void Bins::updateBins(wrapper::Placement &p){
    for(size_t i = 0; i < b.size(); i++){
        for(size_t j = 0; j < b[i].size(); j++){
            b[i][j].moduleIdx.clear();
        }
    }

    for(size_t i = 0; i < p.numModules(); i++){
        pair<pair<size_t, size_t>,pair<size_t, size_t>> box = findModuleBin(p, i);
        // cout<<"box: "<<box.first.second<<" "<<box.second.second<<endl;
        for(size_t j = box.first.first; j <= box.second.first; j++){
            for(size_t k = box.first.second; k <= box.second.second; k++){
                b[k][j].moduleIdx.push_back(i);
                // cout<<j<<" "<<k<<" : "<<i<<endl;
                // double xi = p.module(i).centerX();
                // double xb = b[j][k].x();
                // double wb = _binWidth;
                // double wi = p.module(i).area()/p.module(i).height();
                // double dx = abs(xi - xb);
                // if(wb+wi/2 <= dx){
                //     // cerr<<"error"<<endl;
                // }
            }
        }
    }

    //test 
    // for(size_t i = 0; i < p.numModules(); i++){
    //     pair<pair<size_t, size_t>,pair<size_t, size_t>> box = findModuleBin(p, i);
    //     for(size_t j = 0; j < b.size(); j++){
    //         for(size_t k = 0; k < b.size(); k++){
    //             // if(b[j][k].x() - _binWidth/2 << module)
    //         }
    //     }
    // }
}

pair<pair<size_t, size_t>,pair<size_t, size_t>> Bins::findModuleBin(wrapper::Placement &placement, size_t idx){

    double xbegin = placement.module(idx).centerX() - (placement.module(idx).area()/placement.module(idx).height())/2;
    double ybegin = placement.module(idx).centerY() - (placement.module(idx).height())/2;
    double xend = placement.module(idx).centerX() + (placement.module(idx).area()/placement.module(idx).height())/2;
    double yend = placement.module(idx).centerY() + (placement.module(idx).height())/2;

    size_t beginX = (size_t)((xbegin-placement.boundryLeft())/binWidth());
    if(beginX > 0) beginX--;
    size_t beginY = (size_t)((ybegin-placement.boundryBottom())/binHeight());
    if(beginY > 0) beginY--;
    pair<size_t, size_t> begin(beginX,beginY);

    size_t endX = min((size_t)((xend-placement.boundryLeft())/binWidth())+2, b[0].size()-1);
    size_t endY = min((size_t)((yend-placement.boundryBottom())/binHeight())+2, b[0].size()-1);
    pair<size_t, size_t> end(endX,endY);

    // cout<<"being: "<<beginX<<" "<<beginY<<endl;
    // cout<<"end: "<<endX<<" "<<endY<<endl;

    // if((b[beginX][beginY].x() - binWidth()/2) <= xbegin || (b[beginX][beginY].x() + binWidth()/2) >= xbegin && (b[beginX][beginY].y() - binHeight()/2) <= ybegin || (b[beginX][beginY].y() + binHeight()/2) >= ybegin){

    // }
    // else{
    //     cerr<<"err"<<endl;
    //     exit(1);
    // }

    // if((b[endX][endY].x() - binWidth()/2) <= xend || (b[endX][endY].x() + binWidth()/2) >= xend && (b[endX][endY].y() - binHeight()/2) <= yend || (b[endX][endY].y() + binHeight()/2) >= yend){

    // }
    // else{
    //     cerr<<"err"<<endl;
    //     exit(1);
    // }

    // if(!(b[beginY][beginX].x() - _binWidth/2 <= xbegin && xbegin <= b[beginY][beginX].x() + _binWidth/2)){
    //     cout<<b[beginY][beginX].x() - _binWidth/2<<" "<<b[beginY][beginX].x() + _binWidth/2<<" "<<xbegin<<endl;
    //     cout<<xbegin<<endl;
    //     cout<<xbegin<<" - "<<placement.boundryLeft()<<" = "<<(xbegin-placement.boundryLeft())<<endl;
    //     cout<<binWidth()<<endl;
    //     cout<<(size_t)((xbegin-placement.boundryLeft())/binWidth())<<" "<<beginX<<endl;
    //     exit(1);
    // }

    // if(!(b[beginY][beginX].y() - _binHeight/2 <= ybegin && ybegin <= b[beginY][beginX].y() + _binHeight/2)){
    //     cout<<b[beginY][beginX].y() - _binHeight/2<<" "<<b[beginY][beginX].y() + _binHeight/2<<" "<<ybegin<<endl;
    //     exit(1);
    // }


    // if(!(b[endY][endX].x() - _binWidth/2 <= xend && xend <= b[endY][endX].x() + _binWidth/2)){
    //     cout<<"xend"<<endl;
    //     cout<<b[endY][endX].x() - _binWidth/2<<" "<<b[endY][endX].x() + _binWidth/2<<" "<<xend<<endl;
    //     cout<<placement.module(idx).centerX()<<" "<<(placement.module(idx).width())/2<<endl;
    //     exit(1);
    // }


    // if(!(b[endY][endX].y() - _binHeight/2 <= yend && yend <= b[endY][endX].y() + _binHeight/2)){
    //     cout<<b[endY][endX].y() - _binHeight/2<<" "<<b[endY][endX].y() + _binHeight/2<<" "<<yend<<endl;
    //     cout<<placement.boundryTop()<<endl;
    //     cout<<placement.module(idx).centerY()<<" "<<(placement.module(idx).height())/2<<endl;
    //     exit(1);
    // }


    pair<pair<size_t, size_t>,pair<size_t, size_t>> res = pair<pair<size_t, size_t>,pair<size_t, size_t>>(begin, end);

    return res;
}

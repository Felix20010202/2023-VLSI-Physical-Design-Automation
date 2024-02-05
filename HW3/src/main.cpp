#include"SA.h"
#include<iostream>
#include<signal.h>
#include<unistd.h>
#include<cmath>

using namespace std;

bool timeOut = false;
void sigalrm_handler(int sig)
{
    timeOut = true;
}

int main(int argc, char *argv[]){
    if(argc != 3){
        cerr<<"Usage: hw3 <input> <output>";
        exit(1);
    }

    // set process time
    signal(SIGALRM, &sigalrm_handler);
    alarm(595);

    srand(SEED);

    SA sa;
    sa.fileInput(argv[1]);
    sa.softOnlyInitPlace();
    long long T = sa.getAvgCost();

    sa.process(false, T, 0.001*T, 0.5, sa.TE,timeOut);
    sa.process(true, 0.01*T, 0.0001*T, 0.9, sa.TE,timeOut);

    if(timeOut == false){
        sa.copyTree(sa.bestR);
        T = sa.getAvgCost();
        sa.TE *= 10;
        sa.process(false, T, 0.001*T, 0.7, sa.TE,timeOut);
        sa.process(true, 0.01*T, 0.0001*T, 0.9, sa.TE,timeOut);
        while(timeOut == false){
            sa.copyTree(sa.bestR);
            T = sa.getAvgCost();
            sa.TE *= 1.1;
            sa.process(true, 0.015*T, 0.0001*T, 0.9, sa.TE,timeOut);
        }
    }
    sa.outBest(argv[2]);
    
    return 0;
}
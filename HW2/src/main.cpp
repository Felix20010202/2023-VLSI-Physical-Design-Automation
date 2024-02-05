#include"FM.h"

#include<iostream>
#include<signal.h>
#include <unistd.h>
using namespace std;

// set process time
bool timeOut = false;
void sigalrm_handler(int sig)
{
    timeOut = true;
}

int main(int argc, char *argv[]){
    if(argc != 3){
        cerr<<"Usage: hw2 [input] [output]"<<endl;
        return 1;
    }
    
    // set process time
    signal(SIGALRM, &sigalrm_handler);
    alarm(250);

    FM fm;
    fm.fileInput(argv[1]);
    int mode = 0;

    while(mode < MODE_COUNT){
        if(fm.cellSimplePlace(mode) == -1){
            cout<<mode<<" fail!"<<endl;
            fm.placeReset();
            mode++;
        }
        else{
            break;
        }
    }
    fm.mkBucketList();
    fm.process(timeOut);
    fm.output(argv[2]);
    return 0;
}
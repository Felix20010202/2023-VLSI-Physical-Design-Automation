#include"SA.h"
#include"DList.h"
#include<fstream>
#include<iostream>
#include<string>
#include<cmath>
#include <algorithm>
#include<stack>
#include<unistd.h>

using namespace std;

void SA::fileInput(char * file){
    std::ios::sync_with_stdio(false);
    ifstream fin(file);
    if(!fin){
        cerr<<"error: Can't open file"<<endl;
        exit(1);
    }

    string skip, s, s2;
    fin>>skip>>chipSize[0]>>chipSize[1];
    fin>>skip>>S;
    modules.resize(S);
    for(int i = 0; i < S; i++){
        fin>>skip>>s>>modules[i].area;
        modules[i].ptr->num = i;
        mdMap[s] = i;
        modules[i].name = s;
    }
    fin>>skip>>F;
    modules.resize(S+F);
    for(int i = S; i < S+F; i++){
        fin>>skip>>s>>modules[i].pos[0]>>modules[i].pos[1]>>modules[i].WH[0]>>modules[i].WH[1];
        modules[i].isFix = false;
        modules[i].ptr->num = i;
        mdMap[s] = i;
        modules[i].name = s;
    }
    fin>>skip>>N;
    Net nTemp;
    for(int i = 0; i < N; i++){
        fin>>skip>>s>>s2>>nTemp.weight;
        nTemp.md.first = mdMap[s];
        nTemp.md.second = mdMap[s2];
        modules[mdMap[s]].nets.push_back(i);
        modules[mdMap[s2]].nets.push_back(i);
        nets.push_back(nTemp);
    }
    fin.close();
}

void SA::softOnlyInitPlace(){
    for(int i = 0; i < (S+1)/2; i++){
        if(i*2+1 < S){
            modules[i].ptr->left = modules[i*2+1].ptr;
            modules[i*2+1].ptr->parent = modules[i].ptr;
        }
        if(i*2+2 < S){
            modules[i].ptr->right = modules[i*2+2].ptr;
            modules[i*2+2].ptr->inParentR = true;
            modules[i*2+2].ptr->parent = modules[i].ptr;
        }
    }
    rootNum = 0;

    for(int i = (S+1)/2+1; i < S; i++){
        modules[i].ptr->left = nullptr;
        modules[i].ptr->right = nullptr;
    }
    // init soft width height
    for(int i = 0; i < S; i++){
        int size = std::ceil(std::sqrt(modules[i].area));
        modules[i].WH[0] = modules[i].WH[1] = size;
    }

    // get position
    initList(head, true);
    initList(vhead, false);
    SA::getPos(rootNum);

    wtSolution(best);
    now = best;
}

void SA::outResult(char *file){
    ofstream fout;
    fout.open(file);
    if(fout.is_open() == false){
        cerr<<"error: Open .out file fail"<<endl;
        exit(1);
    };
    fout<<"Wirelength "<<getWireLen()<<endl;
    fout<<"NumSoftModules "<<S<<endl;
    for(int i = 0; i < S; i++){
        fout<<modules[i].name<<" "<<modules[i].pos[0]<<" "
        <<modules[i].pos[1]<<" "<<modules[i].WH[0]<<" "<<modules[i].WH[1]<<endl;        
    }

    fout.close();
}

void SA::outBest(char *file){
    ofstream fout;
    fout.open(file);
    if(fout.is_open() == false){
        cerr<<"error: Open .out file fail"<<endl;
        exit(1);
    };
    fout<<"Wirelength "<<best.wireLen<<endl;
    fout<<"NumSoftModules "<<S<<endl;
    for(int i = 0; i < S; i++){
        fout<<modules[i].name<<" "<<best.X[i]<<" "<<best.Y[i]<<" "<<best.W[i]<<" "<<best.H[i]<<endl;        
    }

    fout.close();
}

void Btree::inorder(Btree *root){
    if(root != nullptr){
        cout<<root->num<<" ";
        inorder(root->left);
        inorder(root->right);
    }
}

long long SA::getWireLen(){
    int res = 0;
    for(int i = 0; i < N; i++){
        int m1 = nets[i].md.first;
        int m2 = nets[i].md.second;
        int x = std::abs((modules[m1].pos[0] + modules[m1].WH[0]/2) - (modules[m2].pos[0] + modules[m2].WH[0]/2));
        int y = std::abs((modules[m1].pos[1] + modules[m1].WH[1]/2) - (modules[m2].pos[1] + modules[m2].WH[1]/2));
        res += (x+y)*nets[i].weight;
    }
    return res;
}

bool SA::legal(int idx){
    bool res = true;
    for(int i = S; i < S+F; i++){
        if(isOverlap(i, idx) == true){
            modules[idx].pos[0] = modules[i].pos[0] + modules[i].WH[0];
            i = S;
            res = false;
        }
    }

    return res;
}

void SA::getPos(int idx){
    if(modules[idx].ptr->parent == nullptr){ // root
        modules[idx].pos[0] = 0;
        modules[idx].pos[1] = 0;
    }
    else{
        int p = modules[idx].ptr->parent->num;
        if(modules[idx].ptr->inParentR){
            modules[idx].pos[0] = modules[p].pos[0];
        }
        else{
            modules[idx].pos[0] = modules[p].pos[0] + modules[p].WH[0];
        }
    }

    copyList(head, headCopy);
    modules[idx].pos[1] = update(modules[idx].pos[0], modules[idx].WH[0], modules[idx].WH[1]);

    while(legal(idx) == false){
        copyList(headCopy, head);
        modules[idx].pos[1] = update(modules[idx].pos[0], modules[idx].WH[0], modules[idx].WH[1]);
    }

    if(modules[idx].pos[0] + modules[idx].WH[0] > maxX){
        maxX = modules[idx].pos[0] + modules[idx].WH[0];
    }
    if(modules[idx].pos[1] + modules[idx].WH[1] > maxY){
        maxY = modules[idx].pos[1] + modules[idx].WH[1];
    }

    if(modules[idx].ptr->left != nullptr){
        SA::getPos(modules[idx].ptr->left->num);
    }
    if(modules[idx].ptr->right != nullptr){
        SA::getPos(modules[idx].ptr->right->num);
    }

    modules[idx].isPlace = true;
}

// default update head list
int SA::update(int pos, int len, int v){
    if(head == nullptr){
        cerr<<"error: SA::update"<<endl;
        exit(1);
    }
    doublyList *ptr1 = head;

    // node in one segment
    while(ptr1 != nullptr){
        if(ptr1->start <= pos && ptr1->end >= pos+len){
            int high = ptr1->value;
            if(ptr1->start == pos && ptr1->end == pos+len){
                ptr1->value = ptr1->value + v;
            }
            else if(ptr1->start == pos && ptr1->end > pos+len){
                doublyList *newptr = new doublyList(pos+len, ptr1->end, ptr1->value);
                newptr->prev = ptr1;
                newptr->next = ptr1->next;
                if(ptr1->next != nullptr){
                    ptr1->next->prev = newptr;
                }
                ptr1->next = newptr;
                ptr1->value = ptr1->value + v;
                ptr1->end = pos+len;
            }
            else if(ptr1->start < pos && ptr1->end == pos+len){
                doublyList *newptr = new doublyList(pos, ptr1->end, ptr1->value + v);
                newptr->prev = ptr1;
                newptr->next = ptr1->next;
                if(ptr1->next != nullptr){
                    ptr1->next->prev = newptr;
                }
                ptr1->next = newptr;
                ptr1->end = pos;
            }
            else if(ptr1->start < pos && ptr1->end > pos+len){
                doublyList *newptr1 = new doublyList(pos, pos+len, ptr1->value+v);
                doublyList *newptr2 = new doublyList(pos+len, ptr1->end, ptr1->value);
                newptr2->next = ptr1->next;
                if(ptr1->next != nullptr){
                    ptr1->next->prev = newptr2;
                }
                newptr2->prev = newptr1;
                newptr1->next = newptr2;
                newptr1->prev = ptr1;
                ptr1->next = newptr1;
                ptr1->end = pos;
            }
            else{
                cerr<<"error: SA::update"<<endl;
                exit(1);
            }
            return high;
        }
        ptr1 = ptr1->next;
    }
    
    // node crocess many segments
    ptr1 = head;
    while(ptr1 != nullptr){
        if(ptr1->start <= pos && ptr1->end > pos){
            doublyList *ptr2 = ptr1->next;
            int high = ptr1->value;
            while(ptr2 != nullptr){
                if(ptr2->value > high){
                    high = ptr2->value;
                }
                if(ptr2->end >= pos+len){
                    doublyList *ptr3 = ptr1->next;
                    while(ptr3 != ptr2){
                        ptr3 = ptr3->next;
                        delete ptr3->prev;
                    }
                    doublyList *newptr1 = new doublyList(pos, pos+len, high+v);
                    newptr1->next = ptr2;
                    newptr1->prev = ptr1;
                    ptr1->next = newptr1;
                    ptr2->prev = newptr1;
                    ptr1->end = pos;
                    if(ptr1->start == ptr1->end){
                        newptr1->prev = ptr1->prev;
                        if(ptr1->prev != nullptr){
                            ptr1->prev->next = newptr1;
                        }
                        if(head == ptr1){
                            head = ptr1->next;
                        }
                        delete ptr1;
                    }
                    ptr2->start = pos+len;
                    if(ptr2->start == ptr2->end){
                        newptr1->next = ptr2->next;
                        if(ptr2->next != nullptr){
                            ptr2->next->prev = newptr1;
                        }
                        delete ptr2;
                    }
                    return high;
                }
                ptr2 = ptr2->next;            
            }
        }
        ptr1 = ptr1->next;
    }
    
    cerr<<"error: SA::update(int, int, int), end"<<endl;
    exit(1);
    return -1;
}

int SA::update(doublyList *head, doublyList **headPtr, int pos, int len, int v){
    if(head == nullptr){
        cerr<<"error: SA::update"<<endl;
        exit(1);
    }
    doublyList *ptr1 = head;
    // node in one segment
    while(ptr1 != nullptr){
        if(ptr1->start <= pos && ptr1->end >= pos+len){
            int high = ptr1->value;
            if(ptr1->start == pos && ptr1->end == pos+len){
                ptr1->value = ptr1->value + v;
            }
            else if(ptr1->start == pos && ptr1->end > pos+len){
                doublyList *newptr = new doublyList(pos+len, ptr1->end, ptr1->value);
                newptr->prev = ptr1;
                newptr->next = ptr1->next;
                if(ptr1->next != nullptr){
                    ptr1->next->prev = newptr;
                }
                ptr1->next = newptr;
                ptr1->value = ptr1->value + v;
                ptr1->end = pos+len;
            }
            else if(ptr1->start < pos && ptr1->end == pos+len){
                doublyList *newptr = new doublyList(pos, ptr1->end, ptr1->value + v);
                newptr->prev = ptr1;
                newptr->next = ptr1->next;
                if(ptr1->next != nullptr){
                    ptr1->next->prev = newptr;
                }
                ptr1->next = newptr;
                ptr1->end = pos;
            }
            else if(ptr1->start < pos && ptr1->end > pos+len){
                doublyList *newptr1 = new doublyList(pos, pos+len, ptr1->value+v);
                doublyList *newptr2 = new doublyList(pos+len, ptr1->end, ptr1->value);
                newptr2->next = ptr1->next;
                if(ptr1->next != nullptr){
                    ptr1->next->prev = newptr2;
                }
                newptr2->prev = newptr1;
                newptr1->next = newptr2;
                newptr1->prev = ptr1;
                ptr1->next = newptr1;
                ptr1->end = pos;
            }
            else{
                cerr<<"error: SA::update"<<endl;
                exit(1);
            }
            return high;
        }
        ptr1 = ptr1->next;
    }
    
    // node crocess many segments
    ptr1 = head;
    while(ptr1 != nullptr){
        if(ptr1->start <= pos && ptr1->end > pos){
            doublyList *ptr2 = ptr1->next;
            int high = ptr1->value;
            while(ptr2 != nullptr){
                if(ptr2->value > high){
                    high = ptr2->value;
                }
                if(ptr2->end >= pos+len){
                    doublyList *ptr3 = ptr1->next;
                    while(ptr3 != ptr2){
                        ptr3 = ptr3->next;
                        delete ptr3->prev;
                    }
                    doublyList *newptr1 = new doublyList(pos, pos+len, high+v);
                    newptr1->next = ptr2;
                    newptr1->prev = ptr1;
                    ptr1->next = newptr1;
                    ptr2->prev = newptr1;
                    ptr1->end = pos;
                    if(ptr1->start == ptr1->end){
                        newptr1->prev = ptr1->prev;
                        if(ptr1->prev != nullptr){
                            ptr1->prev->next = newptr1;
                        }
                        if(head == ptr1){
                            *headPtr = (*headPtr)->next;
                        }
                        delete ptr1;
                    }
                    ptr2->start = pos+len;
                    if(ptr2->start == ptr2->end){
                        newptr1->next = ptr2->next;
                        if(ptr2->next != nullptr){
                            ptr2->next->prev = newptr1;
                        }
                        delete ptr2;
                    }
                    return high;
                }
                ptr2 = ptr2->next;            
            }
        }
        ptr1 = ptr1->next;
    }
    
    cerr<<"error: SA::update, end"<<endl;
    exit(1);
    return -1;
}

void SA::update(doublyList *head,  doublyList **headPtr, int pos, int len, int v, int vStart){
    if(head == nullptr){
        cerr<<"error: SA::update"<<endl;
        exit(1);
    }
    doublyList *ptr1 = head;
    // node in one segment
    while(ptr1 != nullptr){
        if(ptr1->start <= pos && ptr1->end >= pos+len){
            if(ptr1->start == pos && ptr1->end == pos+len){
                ptr1->value = vStart + v;
            }
            else if(ptr1->start == pos && ptr1->end > pos+len){
                doublyList *newptr = new doublyList(pos+len, ptr1->end, ptr1->value);
                newptr->prev = ptr1;
                newptr->next = ptr1->next;
                if(ptr1->next != nullptr){
                    ptr1->next->prev = newptr;
                }
                ptr1->next = newptr;
                ptr1->value = vStart + v;
                ptr1->end = pos+len;
            }
            else if(ptr1->start < pos && ptr1->end == pos+len){
                doublyList *newptr = new doublyList(pos, ptr1->end, vStart + v);
                newptr->prev = ptr1;
                newptr->next = ptr1->next;
                if(ptr1->next != nullptr){
                    ptr1->next->prev = newptr;
                }
                ptr1->next = newptr;
                ptr1->end = pos;
            }
            else if(ptr1->start < pos && ptr1->end > pos+len){
                doublyList *newptr1 = new doublyList(pos, pos+len, vStart+v);
                doublyList *newptr2 = new doublyList(pos+len, ptr1->end, ptr1->value);
                newptr2->next = ptr1->next;
                if(ptr1->next != nullptr){
                    ptr1->next->prev = newptr2;
                }
                newptr2->prev = newptr1;
                newptr1->next = newptr2;
                newptr1->prev = ptr1;
                ptr1->next = newptr1;
                ptr1->end = pos;
            }
            else{
                cerr<<"error: SA::update"<<endl;
                exit(1);
            }
            return;
        }
        ptr1 = ptr1->next;
    }

    // node crocess many segments
    ptr1 = head;
    while(ptr1 != nullptr){
        if(ptr1->start <= pos && ptr1->end > pos){
            doublyList *ptr2 = ptr1->next;
            while(ptr2 != nullptr){
                if(ptr2->end >= pos+len){
                    doublyList *ptr3 = ptr1->next;
                    while(ptr3 != ptr2){
                        ptr3 = ptr3->next;
                        delete ptr3->prev;
                    }
                    doublyList *newptr1 = new doublyList(pos, pos+len, vStart+v);
                    newptr1->next = ptr2;
                    newptr1->prev = ptr1;
                    ptr1->next = newptr1;
                    ptr2->prev = newptr1;
                    ptr1->end = pos;
                    if(ptr1->start == ptr1->end){
                        newptr1->prev = ptr1->prev;
                        if(ptr1->prev != nullptr){
                            ptr1->prev->next = newptr1;
                        }
                        if(head == ptr1){
                            *headPtr = (*headPtr)->next;
                        }
                        delete ptr1;
                    }
                    ptr2->start = pos+len;
                    if(ptr2->start == ptr2->end){
                        newptr1->next = ptr2->next;
                        if(ptr2->next != nullptr){
                            ptr2->next->prev = newptr1;
                        }
                        delete ptr2;
                    }
                    return;
                }
                ptr2 = ptr2->next;            
            }
        }
        ptr1 = ptr1->next;
    }
    
    cerr<<"error: SA::update, end"<<endl;
    exit(1);
}

long long SA::getAvgCost(){
    T = INT_MAX;
    endingT = 0;
    int sampleNum = 1000;
    vector<long long> c(sampleNum, 0);
    int cnt = 0;
    while (T > endingT){
        while(!thermalEquilibrium()){
            initRound();
            perturb();
            getPos(rootNum);
            wtSolution(next);
            long long costNext = cost(next), costNow = cost(now);
            if (costNext < costNow){
                now = next;
                if (costNow < cost(best)){
                    best = now;
                }   
            }
            else if (accept(T, costNext-costNow))
                now = next;
            if(costNext-costNow > 0){
                c[cnt] = costNext;
                cnt++;
                if(cnt == sampleNum){
                    long long avg = 0;
                    for(int i = 0; i < sampleNum; i++){
                        avg += c[i];
                    }
                    return avg/sampleNum;
                }
            }
        }
        T *= decreaseRate;
    }
    return -1;
}

void SA::process(bool loadBest, long long t, long long et, double d, long long argTE, bool &timeOut){
    if(loadBest){
        copyTree(bestR);
    }
    TE = argTE;
    T = t;
    endingT = et;
    decreaseRate = d;
    while (T > endingT && timeOut == false){
        while(!thermalEquilibrium() && timeOut == false){
            initRound();
            copyTree(nowR);
            perturb();
            getPos(rootNum);
            moveHorizontally(rootNum);
            // isPlaceLegal();
            wtSolution(next);
            if(cost(next) < cost(now)){
                now = next;
                if (cost(now) < cost(best)){
                    // cout<<"update: best"<<endl;
                    best = now;
                    copyTree(bestR);
                }   
            }
            else if(accept(T, (cost(next)-cost(now)))){
                now = next;
            }
            else{
                revertTree(nowR);
            }
        }
        T *= decreaseRate;
        // cout<<"T = "<<T<<endl;
    }
    now = best;
}

void SA::perturb(){
    double r = ((double) rand() / (RAND_MAX));
    // change height width 
    if(r < 0.6){
        int idx = rand() % S;
        // rate between 0.5 and 2
        double rate = (((double) rand() / (RAND_MAX))*1.5) + 0.5;
        long long newW = std::ceil(std::sqrt((double)modules[idx].area/rate));
        long long newH = std::ceil((double)modules[idx].area/newW);
        modules[idx].WH[0] = newW;
        modules[idx].WH[1] = newH;
    }
    // Move a module to another place.
    else if(r < 0.8){
        // delete A, move to be B's child
        int A = rand() % S, B = rand() % S;
        while(A == B){
            B = rand() % S;
        }
        deleteNode(A);
        if(A == rootNum){
            Btree *ptr = modules[B].ptr;
            while(ptr != nullptr){
                rootNum = ptr->num;
                ptr = ptr->parent;
            }
        }
        insertNode(B, A);
    }
    // Swap two modules
    else{
        int A = rand() % S, B = rand() % S;
        Btree *temp = modules[A].ptr;
        modules[A].ptr = modules[B].ptr;
        modules[B].ptr = temp;
        modules[A].ptr->num = A;
        modules[B].ptr->num = B;
        if(A == rootNum){
            rootNum = B;
        }
        else if(B == rootNum){
            rootNum = A;
        }
    }
}

bool SA::accept(int T, int cost){
    double r = ((double) rand() / (RAND_MAX));
    // cout<<"p = "<<pow(exp(1), -1*((double)cost/T))<<endl;
    return pow(exp(1), -1*((double)cost/T)) > r;
}

long long SA::cost(Solution &s){
    // return 0.8*s.wireLen;
    // if(s.maxX > chipSize[0] || s.maxY > chipSize[1]){
    //     // return 0.8*s.wireLen + 0.2*std::abs(s.maxX*s.maxY - chipSize[0]*chipSize[1]);
    //     return 100*s.wireLen;
    // }
    long long X = 0,Y = 0;
    if(s.maxX > chipSize[0]) X = s.maxX;
    if(s.maxY > chipSize[1]) Y = s.maxY;
    // if(X != 0 || Y != 0){
    //     return 10*s.wireLen;
    // }
    long long dieSpace = getDieSpace();
    // cout<<"die space: "<<dieSpace<<endl;
    return 0.4*s.wireLen + 0.3*dieSpace + 0.3*((X+Y)*(X+Y));
}

bool SA::thermalEquilibrium(){
    if(nowTE < TE){
        nowTE++;
        return false;
    }
    nowTE = 0;
    return true;
}

void SA::initList(doublyList *head, bool h){
    doublyList *ptr = head->next, *qtr = head->next;
    while(ptr != nullptr){
        qtr = ptr;
        ptr = ptr->next;
        delete qtr;
    }
    head->end = INT_MAX;
    head->start = head->value = 0;
    head->next = head->prev = nullptr;
    if(h == true){
        for(int i = S; i < S+F; i++){
           if(modules[i].pos[1] == 0){
                update(modules[i].pos[0], modules[i].WH[0], modules[i].WH[1]);
            }
        }
    }
    else{
        for(int i = S; i < S+F; i++){
           if(modules[i].pos[0] == 0){
                modules[i].isUpdate = true;
                update(vhead, &vhead, modules[i].pos[1], modules[i].WH[1], modules[i].WH[0]);
            }
        }
    }
}

void SA::copyList(doublyList *head, doublyList *headCopy){
    doublyList *ptr1 = head;
    if(headCopy == nullptr){
        doublyList *ptr2 = new doublyList(ptr1);
        while(ptr1 != nullptr){
            doublyList *newptr = new doublyList(ptr1);
            ptr2->next = newptr;
            newptr->prev = ptr2;
            ptr2 = ptr2->next;
            ptr1 = ptr1->next;
        }
        return;
    }
    doublyList *ptr2 = headCopy, *qtr = nullptr;
    while(ptr1 != nullptr){
        if(ptr2 != nullptr){
            ptr2->start = ptr1->start;
            ptr2->end = ptr1->end;
            ptr2->value = ptr1->value;
            qtr = ptr2;
            ptr2 = ptr2->next;
        }
        else{
            qtr->next = new doublyList(ptr1);
            qtr->next->prev = qtr;
            qtr = qtr->next;
        }
        ptr1 = ptr1->next;
    }
    if(ptr2 != nullptr){
        ptr2->prev->next = nullptr;
        while(ptr2 != nullptr){
            qtr = ptr2;
            ptr2 = ptr2->next;
            delete qtr;
        }
    }
}

bool SA::isOverlap(int A, int B){
    if(modules[A].pos[0] >= modules[B].pos[0] + modules[B].WH[0] || modules[B].pos[0] >=  modules[A].pos[0] + modules[A].WH[0]){
        return false;
    }

    if(modules[A].pos[1] >= modules[B].pos[1] + modules[B].WH[1] || modules[B].pos[1] >=  modules[A].pos[1] + modules[A].WH[1]){
        return false;
    }

    return true;
}

void SA::initRound(){
    maxX = maxY = 0;
    for(int i = S; i < S+F; i++){
        modules[i].isUpdate = false;
    }
    for(int i = 0; i < S; i++){
        modules[i].pos[0] = modules[i].pos[1] = -1;
        modules[i].isPlace = false;
    }
    initList(head, true);
    initList(vhead, false);
}

void SA::wtSolution(Solution &s){
    s.X.clear();
    s.Y.clear();
    s.W.clear();
    s.H.clear();
    s.wireLen = getWireLen();
    s.maxX = s.maxY = 0;
    for(int i = 0; i < S; i++){
        if(modules[i].pos[0] + modules[i].WH[0] > s.maxX){
            s.maxX = modules[i].pos[0] + modules[i].WH[0];
        }
        if(modules[i].pos[1] + modules[i].WH[1] > s.maxY){
            s.maxY = modules[i].pos[1] + modules[i].WH[1];
        }
        s.X.push_back(modules[i].pos[0]);
        s.Y.push_back(modules[i].pos[1]);
        s.W.push_back(modules[i].WH[0]);
        s.H.push_back(modules[i].WH[1]);
    }
}

void SA::showTree(){
    for(int i = 0; i < S; i++){
        cout<<i<<" left: ";
        if(modules[i].ptr->left != nullptr){
            cout<<modules[i].ptr->left->num;
        }
        else{
            cout<<"-1";
        }
        cout<<" right: ";
        if(modules[i].ptr->right != nullptr){
            cout<<modules[i].ptr->right->num;
        }
        else{
            cout<<"-1";
        }
        cout<<" parent: ";
        if(modules[i].ptr->parent != nullptr){
            cout<<modules[i].ptr->parent->num;
        }
        else{
            cout<<"-1";
        }
        cout<<endl;
    }
    isTreeLegal();
}

void SA::showList(doublyList *head){
    doublyList *ptr = head;
    while(ptr != nullptr){
        cout<<"start: "<<ptr->start<<", end: "<<ptr->end<<", val: "<<ptr->value<<" | ";
        ptr = ptr->next;
    }
    cout<<endl;
}

void SA::isTreeLegal(){
    
    for(int i = 0; i < S; i++){
        if(modules[i].ptr->right != nullptr && modules[i].ptr->right->parent != modules[i].ptr){
            cerr<<"error: SA::isTreeLegal "<<modules[i].ptr->right->parent->num<<" "<<modules[i].ptr->num<<endl;
            exit(1);
        }
        if(modules[i].ptr->left != nullptr && modules[i].ptr->left->parent != modules[i].ptr){
            cerr<<"error: SA::isTreeLegal "<<modules[i].ptr->left->parent->num<<" "<<modules[i].ptr->num<<endl;
            exit(1);
        }
    }
    for(int i = 0; i < S; i++){
        if(modules[i].ptr->parent != nullptr && modules[i].ptr->inParentR == true && modules[i].ptr->parent->right != modules[i].ptr){
            cerr<<"error: SA::isTreeLegal, ParentR"<<modules[i].ptr->num<<endl;
            exit(1);
        }
        if(modules[i].ptr->parent != nullptr && modules[i].ptr->inParentR == false && modules[i].ptr->parent->left != modules[i].ptr){
            cerr<<"error: SA::isTreeLegal, ParentR"<<modules[i].ptr->num<<endl;
            exit(1);
        }
    }
}

void SA::isPlaceLegal(){
    bool e = false;
    for(int i = 0; i < S; i++){
        for(int j = i+1; j < S+F; j++){
            if(modules[i].pos[0] >= modules[j].pos[0] + modules[j].WH[0] || modules[j].pos[0] >=  modules[i].pos[0] + modules[i].WH[0]){
                break;
            }

            if(modules[i].pos[1] >= modules[j].pos[1] + modules[j].WH[1] || modules[j].pos[1] >=  modules[i].pos[1] + modules[i].WH[1]){
                break;
            }
            cerr<<"overlap: "<<modules[i].name<<" "<<modules[j].name<<endl;
            cerr<<modules[i].name<<" xy: "<<modules[i].pos[0]<<" "<<modules[i].pos[1]<<" wh: "<<modules[i].WH[0]<<" "<<modules[i].WH[1]<<endl;
            cerr<<modules[j].name<<" xy: "<<modules[j].pos[0]<<" "<<modules[j].pos[1]<<" wh: "<<modules[j].WH[0]<<" "<<modules[j].WH[1]<<endl<<endl;
            e = true;
        }
    }
    if(e == true){
        exit(1);
    }
}

void SA::isListLegal(doublyList *head){
    doublyList *ptr = head;
    while(ptr != nullptr){
        if(ptr->start == ptr->end){
            cerr<<"error: SA::isListLegal, start == end"<<endl;
            exit(1);
        }
        if(ptr->next != nullptr && ptr->next->start < ptr->end){
            cerr<<"error: SA::isListLegal, interval"<<endl;
            exit(1);
        }
        if(ptr->next != nullptr && ptr->next->prev != ptr){
            cerr<<"error: SA::isListLegal, back link"<<endl;
            exit(1);
        }
        ptr = ptr->next;
    }
}

void SA::deleteNode(int idx){
    Btree *root = modules[idx].ptr;
    if(root == nullptr){
        cout<<"error: Btree::deletion"<<endl;
        exit(1);
    }
    // case 1
    if(root->left == nullptr && root->right == nullptr){
        if(root->parent != nullptr){
            if(root->inParentR == true){
                root->parent->right = nullptr;
            }
            else{
                root->parent->left = nullptr;
            }
        }
        root->left = root->right = root->parent = nullptr;
    }
    // case 2
    else if(root->left != nullptr && root->right == nullptr){
        if(root->inParentR == true){
            if(root->parent != nullptr){
                root->parent->right = root->left;
            }
            root->left->parent = root->parent;
            root->left->inParentR = true;
        }
        else{
            if(root->parent != nullptr){
                root->parent->left = root->left;
            }
            root->left->parent = root->parent;
            root->left->inParentR = false;
        }
        root->left = root->parent = nullptr;
    }
    else if(root->left == nullptr && root->right != nullptr){
        if(root->inParentR == true){
            if(root->parent != nullptr){
                root->parent->right = root->right;
            }
            root->right->parent = root->parent;
            root->right->inParentR = true;
        }
        else{
            if(root->parent != nullptr){
                root->parent->left = root->right;
            }
            root->right->parent = root->parent;
            root->right->inParentR = false;
        }
        root->right = root->parent = nullptr;
    }
    // case 3, right child replace
    else{
        Btree *ptr = nullptr, *qtr = root->right;
        while(qtr != nullptr){
            ptr = qtr;
            qtr = qtr->left;
        }
        
        int swapNode = ptr->num;
        modules[idx].ptr = modules[swapNode].ptr;
        modules[swapNode].ptr = root;
        modules[idx].ptr->num = idx;
        modules[swapNode].ptr->num = swapNode;

        if(idx == rootNum){
            rootNum = swapNode;
        }

        deleteNode(idx);
    }
}

void SA::insertNode(int B, int A){
    Btree *p = modules[B].ptr;
    Btree *iNode = modules[A].ptr;
    int r1 = rand()%2, r2 = rand()%2;
    if(r1 == 1){
        if(p->right != nullptr){
            if(r2 == 1){
                iNode->right = p->right;
                p->right->parent = iNode;
                p->right->inParentR = true;
            }
            else{
                iNode->left = p->right;
                p->right->parent = iNode;
                p->right->inParentR = false;
            }
        }
        p->right = iNode;
        iNode->inParentR = true;
    }
    else{
        if(p->left != nullptr){
            if(r2 == 1){
                iNode->right = p->left;
                p->left->parent = iNode;
                p->left->inParentR = true;
            }
            else{
                iNode->left = p->left;
                p->left->parent = iNode;
                p->left->inParentR = false;
            }
        }
        p->left = iNode;
        iNode->inParentR = false;
    }
    iNode->parent = p;
}

void SA::moveHorizontally(int idx){
    long long high = 0;
    for(int i = 0; i < S+F; i++){
        if(i == idx){
            continue;
        }
        if((modules[idx].pos[1] < modules[i].pos[1] + modules[i].WH[1] && modules[i].pos[1] < modules[idx].pos[1] + modules[idx].WH[1]) && modules[idx].pos[0] >= modules[i].pos[0]){
            if(modules[i].pos[0] + modules[i].WH[0] > high){
                high = modules[i].pos[0] + modules[i].WH[0];
            }
        }
    }

    modules[idx].pos[0] = high;
    if(modules[idx].ptr->left != nullptr){
        moveHorizontally(modules[idx].ptr->left->num);
    }
    if(modules[idx].ptr->right != nullptr){
        moveHorizontally(modules[idx].ptr->right->num);
    }
}

void SA::moveLeft(int idx){
    long long high = 0;
    for(int i = 0; i < S+F; i++){
        if(modules[i].isPlace == false){
            continue;
        }
        if(i == idx){
            continue;
        }
        if((modules[idx].pos[1] < modules[i].pos[1] + modules[i].WH[1] && modules[i].pos[1] < modules[idx].pos[1] + modules[idx].WH[1]) && modules[idx].pos[0] >= modules[i].pos[0]){
            if(modules[i].pos[0] + modules[i].WH[0] > high){
                high = modules[i].pos[0] + modules[i].WH[0];
            }
        }
    }

    modules[idx].pos[0] = high;
}

void SA::moveUp(){
    int idx = 0; 
    vector<pair<int, long long>> order;
    for(int i = 0; i < S; i++){
        order.push_back(pair<int, long long>(i, modules[i].pos[1]));
    }
    std::sort(order.begin(), order.end(), [](pair<int, long long> A, pair<int, long long> B){
        return A.second > B.second;
    });

    long long low = chipSize[1];
    for(int j = 0; j < S; j++){
        idx = order[j].first;
        low = chipSize[1] + modules[idx].pos[1];

        for(int i = 0; i < S+F; i++){
            if(i == idx){
                continue;
            }
            if((modules[idx].pos[0] < modules[i].pos[0] + modules[i].WH[0] && modules[i].pos[0] < modules[idx].pos[0] + modules[idx].WH[0]) && modules[idx].pos[1] <= modules[i].pos[1]){
                if(modules[i].pos[1] < low){
                    low = modules[i].pos[1];
                }
            }
        }

        modules[idx].pos[1] = low - modules[idx].WH[1];
    }
}

long long SA::getDieSpace(){
    long long res = chipSize[0]*chipSize[1];
    for(int i = S; i < S+F; i++){
        res -= modules[i].WH[0]*modules[i].WH[1];
    }
    for(int i = 0; i < S; i++){
        if(modules[i].pos[0] + modules[i].WH[0] <= chipSize[0] && modules[i].pos[1] + modules[i].WH[1] <= chipSize[1]){
            res -= modules[i].WH[0]*modules[i].WH[1];
        }
        else{
            long long W = 0, H = 0;
            if(chipSize[0] - modules[i].pos[0] > 0){
                W = chipSize[0] - modules[i].pos[0];
            }
            if(chipSize[1] - modules[i].pos[1] > 0){
                H = chipSize[1] - modules[i].pos[1];
            }
            res -= W*H;
        }
    }

    return res;
}


// discard
// void SA::fillUpdirectly(){
//     stack<int> m;
//     // <index, area>
//     vector<pair<int, long long>> order;
//     for(int i = 0; i < S; i++){
//         modules[i].isPlace = false;
//         order.push_back(pair<int, long long>(i, modules[i].area));
//     }
//     sort(order.begin(), order.end(), [](pair<int, long long> A, pair<int, long long> B){
//         return A.second > B.second;
//     });
//     initList(head, true);
//     for(int i = 0; i < S; i++){
//         long long start = 0, end = 0, high = 0;
//         doublyList *ptr = head;
//         long long low = chipSize[1];
//         while(ptr != nullptr){
//             if(ptr->start < chipSize[0] && ptr->end > chipSize[0] && ptr->value < low){
//                 start = ptr->start;
//                 end = chipSize[0];
//                 low = ptr->value;
//             }
//             else if(ptr->value < low){
//                 start = ptr->start;
//                 end = ptr->end;
//                 low = ptr->value;
//             }
//             ptr = ptr->next;
//         }
//         for(int j = 0; j < S; j++){
//             if(modules[order[j].first].isPlace){
//                 continue;
//             }
//             // long long newW = std::ceil(std::sqrt((double)modules[idx].area/rate));
//             // long long newH = low; //std::ceil((double)modules[idx].area/newW);
//             // if(modules[order[j].first].)
//         }
//     }
// }

bool SA::isResultInBound(){
    if(maxX > chipSize[0] || maxY > chipSize[1]){
        return false;
    }
    return true;
}

bool SA::isBestInBound(){
    if(best.maxX > chipSize[0] || best.maxY > chipSize[1]){
        return false;
    }
    return true;
}

void SA::copyTree(Record &r){
    r.copyT.resize(S);
    r.copyWH.resize(S);
    r.copyPtr.resize(S);
    r.copyRootNum = rootNum;
    
    for(int i = 0; i < S; i++){
        r.copyWH[i].push_back(-1);
        r.copyWH[i].push_back(-1);
        r.copyPtr[i] = modules[i].ptr;
        r.copyT[i].inParentR = modules[i].ptr->inParentR;
        r.copyT[i].left = modules[i].ptr->left;
        r.copyT[i].right = modules[i].ptr->right;
        r.copyT[i].parent = modules[i].ptr->parent;
        r.copyT[i].num = i;
        r.copyWH[i][0] = modules[i].WH[0];
        r.copyWH[i][1] = modules[i].WH[1];
    }
}

void SA::revertTree(Record &r){
    rootNum = r.copyRootNum;
    for(int i = 0; i < S; i++){
        modules[i].ptr = r.copyPtr[i];
        modules[i].WH[0] = r.copyWH[i][0];
        modules[i].WH[1] = r.copyWH[i][1];
        modules[i].ptr->inParentR = r.copyT[i].inParentR;
        modules[i].ptr->left = r.copyT[i].left;
        modules[i].ptr->right = r.copyT[i].right;
        modules[i].ptr->parent = r.copyT[i].parent;
        modules[i].ptr->num = i;
    }
}
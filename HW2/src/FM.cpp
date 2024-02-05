#include"FM.h"
#include <fstream>
#include<string>
#include<iostream>
#include <cstdlib>
#include <climits>
#include<algorithm>
using namespace std;

int globalMaxGain = 0;
char buffer[1024*1024];

void bucketList::add(int idx, doubleList &node){
    if(list.size() < idx){
        cerr<<"Memory size: "<<list.size()<<", index: "<<idx<<endl;
        cerr<<"error: bucketList out of range"<<endl;
        exit(1);
    }
    if(idx > maxPtr){
        maxPtr = idx;
    }
    doubleList *ptr = list[idx].right;
    list[idx].right = &node;
    node.right = ptr;
    node.left = &(list[idx]);
    if(ptr != nullptr){
        ptr->left = &node;
    }
    cellCnt++;
}

// output: cell number
int bucketList::pop(int idx){
    if(idx < 0){
        cerr<<"error: bucketList::pop, out of bound."<<endl;
        exit(1);
    }
    if(maxPtr == -1){
        return -1;
    }
    doubleList *ptr = list[idx].right;
    if(ptr == nullptr){
        cerr<<"error: bucketList::pop, no node in index "<<idx<<endl;
        return -1;
        exit(1);
    }
    list[idx].right = ptr->right;
    if(ptr->right != nullptr)
        ptr->right->left = &(list[idx]);
    int res = ptr->cellNum;
    cellCnt--;
    ptr->left = nullptr;
    ptr->right = nullptr;
    //update max gain ptr
    if(list[idx].right == nullptr){
        if(cellCnt == 0){
            maxPtr = -1;
        }
        for(int i = globalMaxGain*2; i >= 0; i--){
            if(list[i].right != nullptr){
                maxPtr = i;
                break;
            }
        }
    }
    return res;
}

void Cell::updateList(doubleList &head){
    if(lock == true || ptr == nullptr){
        cerr<<"error: Cell::deleteList"<<endl;
        exit(1);
    }
    else if(ptr->left == nullptr){
        cerr<<"error: Cell::deleteList left child is null"<<endl;
        exit(1);
    }
    // It will be tail
    if(ptr->right != nullptr){
        ptr->right->left = ptr->left;
    }
    ptr->left->right = ptr->right;
    // insert
    doubleList *temp = head.right;
    ptr->right = temp;
    ptr->left = &head;
    head.right = ptr;
    if(temp != nullptr){
        temp->left = ptr;
    }
}

void FM::fileInput(char *path){
    std::ios::sync_with_stdio(false);
    ifstream fin(path);
    // std::streambuf *buf = fin.rdbuf();
    if(!fin){
        cerr<<"error: Can't open file"<<endl;
        exit(1);
    }
 
    string skip, s, s2;
    int num;
    fin>>skip>>tech.techNum;
    // cout<<"reading Tech ..."<<endl;
    for(int i = 0; i < tech.techNum; i++){
        fin>>skip>>s>>num;
        tech.techMap[s] = i;
        tech.techSize.push_back(num);

        unordered_map<string, pair<int, int>> Cell;
        for(int j = 0; j < num; j++){
            pair<int,int> H_W;
            fin>>skip>>s>>H_W.first>>H_W.second;
            Cell[s] = H_W;
        }
        tech.libCell.push_back(Cell);
    }

    fin>>skip>>dieH_W.first>>dieH_W.second;
    fin>>skip>>s>>dieUtil[0];
    dieTech.first = tech.techMap[s];
    fin>>skip>>s>>dieUtil[1];
    dieTech.second = tech.techMap[s];
    dieArea = (dieH_W.first * dieH_W.second);
    dieMax[0] = dieArea*((float)dieUtil[0]/100);
    dieMax[1] = dieArea*((float)dieUtil[1]/100);

    // cout<<"reading Cell..."<<endl;
    fin>>skip>>C;
    for(int i = 0; i < C; i++){
        fin>>skip>>s>>s2;
        cellMap[s] = i;
        cells.push_back(Cell(s, s2, i));
    }
    // cout<<"reading Net..."<<endl;
    fin>>skip>>N;
    for(int i = 0; i < N; i++){
        fin>>skip>>s>>num;
        netMap[s] = i;
        nets.push_back(Net(num));
        if(num > pinMax) pinMax = num;
        for(int j = 0; j < num; j++){
            fin>>skip>>s;
            nets[i].cell.push_back(cellMap[s]);
            cells[cellMap[s]].net.push_back(i);
        }
    }
    fin.close();
}

int FM::initCellGain(){
    // whether a net is critical
    for(int i = 0; i < N; i++){
        if(nets[i].Distribution[0] == 1 || nets[i].Distribution[0] == 0 ||
           nets[i].Distribution[1] == 1 || nets[i].Distribution[1] == 0){
            nets[i].critical = true;
        }
    }

    for(int i = 0; i < C; i++){
        cells[i].gain = 0;
    }

    // computing cell gain
    for(int i = 0; i < C; i++){
        vector<int> net(cells[i].net);

        int idx = 0;
        if(cells[i].inA == false) idx = 1;

        for(const int j : net){
            if(nets[j].Distribution[idx] == 1){
                cells[i].gain++;
            }
            if(nets[j].Distribution[!idx] == 0){
                cells[i].gain--;
            }
        }
    }
    return 0;
}

void FM::addCell2List(){
    for(int i = 0; i <= globalMaxGain*2; i++){
        listA->list[i].right = nullptr;
        listB->list[i].right = nullptr;
    }
    for(int i = 0; i < C; i++){
        int idx = cells[i].gain + maxGain;
        if(idx < 0 || idx > 2*maxGain){
            showBucketList(0);
            cerr<<"error: doubleList index out of range "<<idx<<endl;
            exit(1);
        }
        if(cells[i].inA){
            listA->add(idx, *(cells[i].ptr));
        }
        else{
            listB->add(idx, *(cells[i].ptr));
        }
    }
}

void FM::mkBucketList(){
    for(int i = 0; i < C; i++){
        if(cells[i].net.size() > maxGain){
            maxGain = cells[i].net.size();
        }
    }

    globalMaxGain = maxGain; 
    listA = new bucketList(2*maxGain+1);
    listB = new bucketList(2*maxGain+1);
}

void FM::lock(int base){
    cells[base].lock = true;
    for(int i = 0; i < cells[base].net.size(); i++){
        nets[cells[base].net[i]].lockPin++;
    }
    lockCellNum++;
}

void FM::initRound(){
    lockCellNum = 0;
    for(int i = 0; i < C; i++){
        cells[i].lock = false;
    }
    for(int i = 0; i < N; i++){
        nets[i].lockPin = 0;
    }
    initCellGain();
    addCell2List();
}

bool FM::test(int base){
    pair<int, int> areaA = tech.libCell[dieTech.first][cells[base].libCell];
    pair<int, int> areaB = tech.libCell[dieTech.second][cells[base].libCell];
    long long putASize = ((long long)areaA.first*areaA.second);
    long long putBSize = ((long long)areaB.first*areaB.second);

    if(cells[base].inA && dieUsed[1] + putBSize <= dieMax[1]){
        return true;
    }
    else if(cells[base].inA == false && dieUsed[0] + putASize <= dieMax[0]){
        return true;
    }
    return false;
}

void FM::isNetCritical(int idx){
    if(nets[idx].lockPin > 4){
        nets[idx].critical = false;
        return;
    }
    if(nets[idx].Distribution[0] == 1 || nets[idx].Distribution[0] == 0 ||
       nets[idx].Distribution[1] == 1 || nets[idx].Distribution[1] == 1){
        nets[idx].critical = true;
        return;
    }
    nets[idx].critical = false;
}

bool FM::testAndSwap(int base){
    pair<int, int> areaA = tech.libCell[dieTech.first][cells[base].libCell];
    pair<int, int> areaB = tech.libCell[dieTech.second][cells[base].libCell];
    long long putASize = ((long long)areaA.first*areaA.second);
    long long putBSize = ((long long)areaB.first*areaB.second);

    // A swap to B
    if(cells[base].inA && dieUsed[1] + putBSize <= dieMax[1]){
        dieUsed[1] += putBSize;
        dieUsed[0] -= putASize;
        cellNum[0]--;
        cellNum[1]++;
        return true;
    }
    // B swap to A
    else if(cells[base].inA == false && dieUsed[0] + putASize <= dieMax[0]){
        dieUsed[0] += putASize;
        dieUsed[1] -= putBSize;
        cellNum[1]--;
        cellNum[0]++;
        return true;
    }
    // else no swap
    return false;
}

void FM::testNet(int base){
    cout<<"base: "<<base<<" ";
    for(const int &n : cells[base].net){
        if(nets[n].Distribution[0] > 0 && nets[n].Distribution[1] > 0)
            cout<<n<<" ";
    }
    cout<<endl;
}

void FM::showBucketList(int mode){
    ofstream fout;
    if(mode == 0){
        fout.open("./BucketList.txt");
    }
    else{
        fout.open("./BucketList.txt", std::ios_base::app);
    }
    if(fout.is_open() == false){
        cerr<<"error: FM::showBucketList, open file fail"<<endl;
        exit(1);
    };
    int cnt = 0;
    fout<<"List A: "<<endl;
    for(int i = 2*maxGain; i >= 0; i--){
        fout<<i<<": ";
        doubleList *ptr = listA->list[i].right;
        while(ptr != nullptr){
            fout<<ptr->cellNum<<" ";
            ptr = ptr->right;
            cnt++;
        }
        fout<<endl;
    }
    fout<<"List B: "<<endl;
    for(int i = 2*maxGain; i >= 0; i--){
        fout<<i<<": ";
        doubleList *ptr = listB->list[i].right;
        while(ptr != nullptr){
            fout<<ptr->cellNum<<" ";
            ptr = ptr->right;
            cnt++;
        }
        fout<<endl;
    }
    fout<<"Total cells "<<C<<", in list "<<cnt<<endl<<endl;
    fout.close();
}

int FM::chooseCell(){
    int gainA = listA->maxPtr, gainB = listB->maxPtr;
    if(gainA >= gainB && gainA != -1){
        return listA->pop(gainA);
    }
    else if(gainB != -1){
        return listB->pop(gainB);
    }
    return -1;
}

void FM::testBucketList(){
    int cnt = 0;
    for(int i = 2*maxGain; i >= 0; i--){
        doubleList *ptr = listA->list[i].right;
        while(ptr != nullptr){
            ptr = ptr->right;
            cnt++;
        }
    }
    for(int i = 2*maxGain; i >= 0; i--){
        doubleList *ptr = listB->list[i].right;
        while(ptr != nullptr){
            ptr = ptr->right;
            cnt++;
        }
    }
    if(cnt != C-lockCellNum){
        cerr<<"C: "<<C<<", lock "<<lockCellNum<<endl;
        cerr<<"number of nodes in List: "<<cnt<<endl;
        exit(1);
    }
}

void FM::testGain(){
    for(int i = 0; i < C; i++){
        cells[i].testGain = 0;
    }
    for(int i = 0; i < C; i++){
        int idx = 0;
        if(cells[i].inA == false) idx = 1;
        vector<int> net(cells[i].net);

        for(const int j : net){
            if(nets[j].Distribution[idx] == 1){
                cells[i].testGain++;
            }
            if(nets[j].Distribution[!idx] == 0){
                cells[i].testGain--;
            }
        }
    }
    for(int i = 0; i < C; i++){
        if(cells[i].testGain != cells[i].gain && cells[i].lock == false){
            cerr<<"error: FM::testGain cell "<<i<<", need "<<cells[i].testGain<<", val "<<cells[i].gain<<endl;
            exit(1);
        }
    }
}


void FM::updateGain(int base){
    unordered_set<int> modifyCell;
    int F = 0, T = 1;
    if(cells[base].inA == true){
        F = 0; T = 1;
    }
    else{
        F = 1; T = 0;
    }
    cells[base].inA = !cells[base].inA;

    for(int i = 0; i < cells[base].net.size(); i++){
        int netIdx = cells[base].net[i];
        if(nets[netIdx].Distribution[T] == 0){
            for (const int &s : nets[netIdx].set[F]) {
                if(cells[s].lock == false){
                    cells[s].gain++;
                    modifyCell.insert(s);
                }
            }
        }
        else if(nets[netIdx].Distribution[T] == 1){
            for (const int &s : nets[netIdx].set[T]) {
                if(cells[s].lock == false){
                    cells[s].gain--;
                    modifyCell.insert(s);
                }
            }
        }
        nets[netIdx].set[F].erase(base);
        nets[netIdx].set[T].insert(base);
        nets[netIdx].Distribution[F]--;
        nets[netIdx].Distribution[T]++;


        int A = 0, B = 0; 
        for(int j = 0;j < nets[netIdx].cell.size(); j++){
            if(cells[nets[netIdx].cell[j]].inA){
                A++;
            }
            else B++;
        }
        if(nets[netIdx].Distribution[0] != A || nets[netIdx].Distribution[1] != B){
            cerr<< nets[netIdx].Distribution[0]<<" "<<A<<" "<<nets[netIdx].Distribution[1]<<" "<<B<<endl;
            nets[netIdx].Distribution[0] = A;
            nets[netIdx].Distribution[1] = B;
            exit(1);
        }

        if(nets[netIdx].Distribution[F] == 0){
            for (const int &s : nets[netIdx].set[T]) {
                if(cells[s].lock == false){
                    cells[s].gain--;
                    modifyCell.insert(s);
                }
            }
        }
        else if(nets[netIdx].Distribution[F] == 1){
            for (const int &s : nets[netIdx].set[F]) {
                if(cells[s].lock == false){
                    cells[s].gain++;
                    modifyCell.insert(s);
                }
            }
        }
    }

    // testGain();
    for(const int &s : modifyCell){
        int idx = cells[s].gain + maxGain;
        // cout<<"update "<<s<<", new gain "<<cells[s].gain<<endl;
        if(idx < 0){
            cerr<<"error: FM::updateGain, out of range "<<idx<<endl;
            exit(1);
        }
        if(cells[s].lock == false){
            if(cells[s].inA){
                cells[s].updateList(listA->list[idx]);
                if(idx > listA->maxPtr){
                    listA->maxPtr = idx; 
                }
                if(listA->list[listA->maxPtr].right == nullptr){
                    for(int i = globalMaxGain*2; i >= 0; i--){
                        if(listA->list[i].right != nullptr){
                            listA->maxPtr = i;
                            break;
                        }
                    }
                }
            }
            else{
                cells[s].updateList(listB->list[idx]);
                if(idx > listB->maxPtr){
                    listB->maxPtr = idx;
                }
                if(listB->list[listB->maxPtr].right == nullptr){
                    for(int i = globalMaxGain*2; i >= 0; i--){
                        if(listB->list[i].right != nullptr){
                            listB->maxPtr = i;
                            break;
                        }
                    }
                }
            }
        }
    }
}

void FM::recovery(int base){
    int F, T;
    if(cells[base].inA){
        F = 0; T = 1;
    }
    else{
        F = 1; T = 0;
    }
    cells[base].inA = !cells[base].inA;
    for(int i = 0; i < cells[base].net.size(); i++){
        int netIdx = cells[base].net[i];
        
        nets[netIdx].set[F].erase(base);
        nets[netIdx].set[T].insert(base);
        nets[netIdx].Distribution[F]--;
        nets[netIdx].Distribution[T]++;
    }
}

void FM::process(bool &timeOut){
    int round = 0;
    int minCutSize = cutSize;
    // bool retry = false;
    while(1){
        // cout<<"Round "<<++round<<endl;
        initRound();
        vector<int> gain;
        vector<int> swapIdx;
        for(int i = 0; i < C; i++){
            int baseCell = chooseCell();
            // cout<<"baseCell: "<<baseCell<<endl;
            if(baseCell == -1){
                showBucketList(1);
                cerr<<listA->maxPtr<<" "<<listB->maxPtr<<endl;
                cerr<<"error: Bucket list has problem."<<endl;
                cerr<<"error: Find "<<i<<" cell"<<endl;
                cerr<<"error: Can't find a cell, bucket list don't have enough cells."<<endl;
                exit(1);
            }
            lock(baseCell);
            if(testAndSwap(baseCell) == true){
                gain.push_back(cells[baseCell].gain);
                swapIdx.push_back(baseCell);
                updateGain(baseCell);
            }
        }
        long long accMax = 0, acc = 0;
        int maxIdx = -1;
        for(int i = 0; i < gain.size(); i++){
            acc += gain[i];
            if(acc > accMax){
                accMax = acc;
                maxIdx = i;
            }
        }
        cutSize -= accMax;
        // cout<<"Swap "<<maxIdx+1<<" nodes, total gain "<<accMax<<endl<<"New cut size "<<cutSize<<endl<<endl;
        // Restore previous swap
        for(int i = gain.size()-1; i > maxIdx; i--){
            if(testAndSwap(swapIdx[i]) == false){
                cerr<<"error: FM::process restore"<<endl;
                exit(1);
            }
            recovery(swapIdx[i]);
        }
        // if(retry == false && accMax <= 0){
        //     retry = 
        // }
        // if(cutSize < minCutSize){
        //     minCutSize = cutSize;
        //     cout<<"New cut size "<<minCutSize<<endl;
        // }

        if(timeOut == true || accMax <= 0){
            break;
        }
    }
}

void FM::placeReset(){
    dieUsed[0] = 0;
    dieUsed[1] = 0;
    cellNum[0] = 0;
    cellNum[1] = 0;
    for(int i = 0; i < N; i++){
        nets[i].Distribution[0] = 0;
        nets[i].Distribution[1] = 0;
    }
}

int FM::getCutSize(){
    int cutNum = 0;
    for(int i = 0; i < N; i++){
        if(nets[i].Distribution[0] < 0 || nets[i].Distribution[1] < 0){
            cerr<<"error: FM::getCutSize"<<endl;
            exit(1);
        }
        if(nets[i].Distribution[0] > 0 && nets[i].Distribution[1] > 0){
            cutNum++;
        }
    }
    cutSize = cutNum;
    return cutNum;
}

void FM::output(char *path){
    ofstream fout;
    fout.open(path);
    if(fout.is_open() == false){
        cerr<<"error: Open .out file fail"<<endl;
        exit(1);
    };
    fout<<"CutSize "<<getCutSize()<<endl;
    fout<<"DieA "<<cellNum[0]<<endl;
    for(int i = 0; i < C; i++){
        if(cells[i].inA == true){
            fout<<cells[i].name<<endl;
        }
    }
    fout<<"DieB "<<cellNum[1]<<endl;
    for(int i = 0; i < C; i++){
        if(cells[i].inA == false){
            fout<<cells[i].name<<endl;
        }
    }
    fout.close();
}

int FM::cellSimplePlace(int mode){
    if(C == 0 || N == 0){
        cerr<<"error: No cell input"<<endl;
        exit(1);
    }

    for(int i = 0; i < C ; i++){
        pair<int, int> areaA = tech.libCell[dieTech.first][cells[i].libCell];
        pair<int, int> areaB = tech.libCell[dieTech.second][cells[i].libCell];
        cells[i].putASize = ((long long)areaA.first*areaA.second);
        cells[i].putBSize = ((long long)areaB.first*areaB.second);
        long long temp = std::max(cells[i].putASize, cells[i].putBSize);
        if(temp > maxCellSize && dieArea*0.01 > temp){
            maxCellSize = temp;
        }
    }

    
    // select micro O(n);
    if(mode == 0){
        unordered_set<int> microSet;
        vector<pair<int, long long>> micro;
        for(int i = 0; i < C; i++){
            if( std::min(cells[i].putASize, cells[i].putBSize) > dieArea*0.01){
                microSet.insert(i);
                micro.push_back(pair<int, long long>(i, std::max(cells[i].putASize, cells[i].putBSize)));
            }
        }
        
        // if(micro.size() == 0){
        //     return -1;
        // }

        sort(micro.begin(), micro.end(), 
        [](pair<int, long long> a, pair<int, long long> b){
            return a.second > b.second; // 降序排列
        });

        for(const auto &m : micro){
            int i = m.first;
            if(cells[i].putASize <= cells[i].putBSize && dieUsed[0] + cells[i].putASize < dieMax[0]*0.7){
                for(int j = 0; j < cells[i].net.size(); j++){
                    nets[cells[i].net[j]].Distribution[0]++;
                }
                dieUsed[0] += cells[i].putASize;
                cells[i].inA = true;
                cellNum[0]++;
            }
            else if(dieUsed[1] + cells[i].putBSize < dieMax[1]){ //put micro to B
                for(int j = 0; j < cells[i].net.size(); j++){
                    nets[cells[i].net[j]].Distribution[1]++;
                }
                dieUsed[1] += cells[i].putBSize;
                cells[i].inA = false;
                cellNum[1]++;
            }
            else{
                for(int j = 0; j < cells[i].net.size(); j++){
                    nets[cells[i].net[j]].Distribution[0]++;
                }
                dieUsed[0] += cells[i].putASize;
                cells[i].inA = true;
                cellNum[0]++;
            }
        }
        for(int i = 0; i < C; i++){
            if(microSet.count(i) == 0){
                if(cells[i].putASize <= cells[i].putBSize && dieUsed[0] + cells[i].putASize <= dieMax[0] - maxCellSize){
                    for(int j = 0; j < cells[i].net.size(); j++){
                        nets[cells[i].net[j]].Distribution[0]++;
                    }
                    dieUsed[0] += cells[i].putASize;
                    cells[i].inA = true;
                    cellNum[0]++;
                }
                // put B
                else if(dieUsed[1] + cells[i].putBSize <= dieMax[1]){
                    for(int j = 0; j < cells[i].net.size(); j++){
                        nets[cells[i].net[j]].Distribution[1]++;
                    }
                    dieUsed[1] += cells[i].putBSize;
                    cells[i].inA = false;
                    cellNum[1]++;
                }
                else{
                    cerr<<"Usage limits: 0."<<dieUtil[0]<<", 0."<<dieUtil[1]<<endl;
                    cerr<<"A Usage rate : "<<((float)dieUsed[0]/dieArea)<<endl;
                    cerr<<"B Usage rate : "<<((float)dieUsed[1]/dieArea)<<endl;
                    cerr<<"Total cell count : "<<C<<", Number of non-placement cells is "<<(C-i)<<endl;
                    cerr<<"error: cellSimplePlace fail!"<<endl;
                    return -1;
                }
            }
        }
    }
    else if(mode == 1){
        for(int i = 0; i < C; i++){
            // put A
            if(cells[i].putASize <= cells[i].putBSize && dieUsed[0] + cells[i].putASize <= dieMax[0] - maxCellSize){
                for(int j = 0; j < cells[i].net.size(); j++){
                    nets[cells[i].net[j]].Distribution[0]++;
                }
                dieUsed[0] += cells[i].putASize;
                cells[i].inA = true;
                cellNum[0]++;
            }
            // put B
            else if(dieUsed[1] + cells[i].putBSize <= dieMax[1]){
                for(int j = 0; j < cells[i].net.size(); j++){
                    nets[cells[i].net[j]].Distribution[1]++;
                }
                dieUsed[1] += cells[i].putBSize;
                cells[i].inA = false;
                cellNum[1]++;
            }
            else{
                cerr<<"Usage limits: 0."<<dieUtil[0]<<", 0."<<dieUtil[1]<<endl;
                cerr<<"A Usage rate : "<<((float)dieUsed[0]/dieArea)<<endl;
                cerr<<"B Usage rate : "<<((float)dieUsed[1]/dieArea)<<endl;
                cerr<<"Total cell count : "<<C<<", Number of non-placement cells is "<<(C-i)<<endl;
                cerr<<"error: cellSimplePlace fail!"<<endl;
                return -1;
            }
        }
    }
    else if(mode == 2){
        unordered_set<int> microSet;
        vector<pair<int, long long>> micro;
        for(int i = 0; i < C; i++){
            if( std::min(cells[i].putASize, cells[i].putBSize) > dieArea*0.01){
                microSet.insert(i);
                micro.push_back(pair<int, long long>(i, std::max(cells[i].putASize, cells[i].putBSize)));
            }
        }
        
        // if(micro.size() == 0){
        //     return -1;
        // }

        sort(micro.begin(), micro.end(), 
        [](pair<int, long long> a, pair<int, long long> b){
            return a.second > b.second; // 降序排列
        });

        for(const auto &m : micro){
            int i = m.first;
            if(cells[i].putASize <= cells[i].putBSize && dieUsed[0] + cells[i].putASize < dieMax[0]*0.7){
                for(int j = 0; j < cells[i].net.size(); j++){
                    nets[cells[i].net[j]].Distribution[0]++;
                }
                dieUsed[0] += cells[i].putASize;
                cells[i].inA = true;
                cellNum[0]++;
            }
            else if(dieUsed[1] + cells[i].putBSize < dieMax[1]){ //put micro to B
                for(int j = 0; j < cells[i].net.size(); j++){
                    nets[cells[i].net[j]].Distribution[1]++;
                }
                dieUsed[1] += cells[i].putBSize;
                cells[i].inA = false;
                cellNum[1]++;
            }
            else{
                for(int j = 0; j < cells[i].net.size(); j++){
                    nets[cells[i].net[j]].Distribution[0]++;
                }
                dieUsed[0] += cells[i].putASize;
                cells[i].inA = true;
                cellNum[0]++;
            }
        }
        for(int i = 0; i < C; i++){
            if(microSet.count(i) == 0){
                if(cells[i].putASize < cells[i].putBSize && dieUsed[0] + cells[i].putASize <= dieMax[0] - maxCellSize){
                    for(int j = 0; j < cells[i].net.size(); j++){
                        nets[cells[i].net[j]].Distribution[0]++;
                    }
                    dieUsed[0] += cells[i].putASize;
                    cells[i].inA = true;
                    cellNum[0]++;
                }
                // put B
                else if(dieUsed[1] + cells[i].putBSize <= dieMax[1]){
                    for(int j = 0; j < cells[i].net.size(); j++){
                        nets[cells[i].net[j]].Distribution[1]++;
                    }
                    dieUsed[1] += cells[i].putBSize;
                    cells[i].inA = false;
                    cellNum[1]++;
                }
                else if(dieUsed[0] + cells[i].putASize <= dieMax[0]){
                    for(int j = 0; j < cells[i].net.size(); j++){
                        nets[cells[i].net[j]].Distribution[0]++;
                    }
                    dieUsed[0] += cells[i].putASize;
                    cells[i].inA = true;
                    cellNum[0]++;
                }
                else{
                    cerr<<"Usage limits: 0."<<dieUtil[0]<<", 0."<<dieUtil[1]<<endl;
                    cerr<<"A Usage rate : "<<((float)dieUsed[0]/dieArea)<<endl;
                    cerr<<"B Usage rate : "<<((float)dieUsed[1]/dieArea)<<endl;
                    cerr<<"Total cell count : "<<C<<", Number of non-placement cells is "<<(C-i)<<endl;
                    cerr<<"error: cellSimplePlace fail!"<<endl;
                    return -1;
                }
            }
        }
    }
    else if(mode == 3){
        vector<pair<int, long long>> sortSize;
        for(int i = 0; i < C; i++){
            sortSize.push_back(pair<int, long long>(i, std::max(cells[i].putASize , cells[i].putASize)));
        }
        sort(sortSize.begin(), sortSize.end(), 
            [](pair<int, long long> a, pair<int, long long> b){
                return a.second > b.second; // 降序排列
            });
        for(int idx = 0; idx < C; idx++){
            int i = sortSize[idx].first;
            // put A
            if(cells[i].putASize <= cells[i].putBSize && dieUsed[0] + cells[i].putASize <= dieMax[0]*0.95){
                for(int j = 0; j < cells[i].net.size(); j++){
                    nets[cells[i].net[j]].Distribution[0]++;
                }
                dieUsed[0] += cells[i].putASize;
                cells[i].inA = true;
                cellNum[0]++;
            }
            // put B
            else if(dieUsed[1] + cells[i].putBSize <= dieMax[1]*0.95 - maxCellSize){
                for(int j = 0; j < cells[i].net.size(); j++){
                    nets[cells[i].net[j]].Distribution[1]++;
                }
                dieUsed[1] += cells[i].putBSize;
                cells[i].inA = false;
                cellNum[1]++;
            }
            else if(dieUsed[0] + cells[i].putASize <= dieMax[0] - maxCellSize){
                for(int j = 0; j < cells[i].net.size(); j++){
                    nets[cells[i].net[j]].Distribution[0]++;
                }
                dieUsed[0] += cells[i].putASize;
                cells[i].inA = true;
                cellNum[0]++;
            }
            else if(dieUsed[1] + cells[i].putBSize <= dieMax[1] - maxCellSize){
                for(int j = 0; j < cells[i].net.size(); j++){
                    nets[cells[i].net[j]].Distribution[1]++;
                }
                dieUsed[1] += cells[i].putBSize;
                cells[i].inA = false;
                cellNum[1]++;
            }
            else{
                cerr<<"Usage limits: 0."<<dieUtil[0]<<", 0."<<dieUtil[1]<<endl;
                cerr<<"A Usage rate : "<<((float)dieUsed[0]/dieArea)<<endl;
                cerr<<"B Usage rate : "<<((float)dieUsed[1]/dieArea)<<endl;
                cerr<<"Total cell count : "<<C<<", Number of non-placement cells is "<<(C-i)<<endl;
                cerr<<"error: cellSimplePlace fail!"<<endl;
                return -1;
            }
        }
    }
    else if(mode == 4){
        unordered_set<int> microSet;
        vector<pair<int, long long>> micro;
        for(int i = 0; i < C; i++){
            if( std::min(cells[i].putASize, cells[i].putBSize) > dieArea*0.01){
                microSet.insert(i);
                micro.push_back(pair<int, long long>(i, std::max(cells[i].putASize, cells[i].putBSize)));
            }
        }
        
        // if(micro.size() == 0){
        //     return -1;
        // }

        sort(micro.begin(), micro.end(), 
        [](pair<int, long long> a, pair<int, long long> b){
            return a.second < b.second; // 降序排列
        });

        for(const auto &m : micro){
            int i = m.first;
            if(cells[i].putASize <= cells[i].putBSize && dieUsed[0] + cells[i].putASize < dieMax[0]*0.7){
                for(int j = 0; j < cells[i].net.size(); j++){
                    nets[cells[i].net[j]].Distribution[0]++;
                }
                dieUsed[0] += cells[i].putASize;
                cells[i].inA = true;
                cellNum[0]++;
            }
            else if(dieUsed[1] + cells[i].putBSize < dieMax[1]){ //put micro to B
                for(int j = 0; j < cells[i].net.size(); j++){
                    nets[cells[i].net[j]].Distribution[1]++;
                }
                dieUsed[1] += cells[i].putBSize;
                cells[i].inA = false;
                cellNum[1]++;
            }
            else{
                for(int j = 0; j < cells[i].net.size(); j++){
                    nets[cells[i].net[j]].Distribution[0]++;
                }
                dieUsed[0] += cells[i].putASize;
                cells[i].inA = true;
                cellNum[0]++;
            }
        }
        for(int i = 0; i < C; i++){
            if(microSet.count(i) == 0){
                if(cells[i].putASize <= cells[i].putBSize && dieUsed[0] + cells[i].putASize <= dieMax[0] - maxCellSize){
                    for(int j = 0; j < cells[i].net.size(); j++){
                        nets[cells[i].net[j]].Distribution[0]++;
                    }
                    dieUsed[0] += cells[i].putASize;
                    cells[i].inA = true;
                    cellNum[0]++;
                }
                // put B
                else if(dieUsed[1] + cells[i].putBSize <= dieMax[1]){
                    for(int j = 0; j < cells[i].net.size(); j++){
                        nets[cells[i].net[j]].Distribution[1]++;
                    }
                    dieUsed[1] += cells[i].putBSize;
                    cells[i].inA = false;
                    cellNum[1]++;
                }
                else{
                    cerr<<"Usage limits: 0."<<dieUtil[0]<<", 0."<<dieUtil[1]<<endl;
                    cerr<<"A Usage rate : "<<((float)dieUsed[0]/dieArea)<<endl;
                    cerr<<"B Usage rate : "<<((float)dieUsed[1]/dieArea)<<endl;
                    cerr<<"Total cell count : "<<C<<", Number of non-placement cells is "<<(C-i)<<endl;
                    cerr<<"error: cellSimplePlace fail!"<<endl;
                    return -1;
                }
            }
        }
    }
    else if(mode == 5){
        vector<pair<int, long long>> sortSize;
        for(int i = 0; i < C; i++){
            sortSize.push_back(pair<int, long long>(i, std::max(cells[i].putASize , cells[i].putASize)));
        }
        sort(sortSize.begin(), sortSize.end(), 
            [](pair<int, long long> a, pair<int, long long> b){
                return a.second < b.second; // 降序排列
            });
        for(int idx = 0; idx < C; idx++){
            int i = sortSize[idx].first;
            // put A
            if(cells[i].putASize <= cells[i].putBSize && dieUsed[0] + cells[i].putASize <= dieMax[0]*0.8){
                for(int j = 0; j < cells[i].net.size(); j++){
                    nets[cells[i].net[j]].Distribution[0]++;
                }
                dieUsed[0] += cells[i].putASize;
                cells[i].inA = true;
                cellNum[0]++;
            }
            // put B
            else if(dieUsed[1] + cells[i].putBSize <= dieMax[1]*0.8 - maxCellSize){
                for(int j = 0; j < cells[i].net.size(); j++){
                    nets[cells[i].net[j]].Distribution[1]++;
                }
                dieUsed[1] += cells[i].putBSize;
                cells[i].inA = false;
                cellNum[1]++;
            }
            else if(dieUsed[0] + cells[i].putASize <= dieMax[0] - maxCellSize){
                for(int j = 0; j < cells[i].net.size(); j++){
                    nets[cells[i].net[j]].Distribution[0]++;
                }
                dieUsed[0] += cells[i].putASize;
                cells[i].inA = true;
                cellNum[0]++;
            }
            else if(dieUsed[1] + cells[i].putBSize <= dieMax[1] - maxCellSize){
                for(int j = 0; j < cells[i].net.size(); j++){
                    nets[cells[i].net[j]].Distribution[1]++;
                }
                dieUsed[1] += cells[i].putBSize;
                cells[i].inA = false;
                cellNum[1]++;
            }
            else{
                cerr<<"Usage limits: 0."<<dieUtil[0]<<", 0."<<dieUtil[1]<<endl;
                cerr<<"A Usage rate : "<<((float)dieUsed[0]/dieArea)<<endl;
                cerr<<"B Usage rate : "<<((float)dieUsed[1]/dieArea)<<endl;
                cerr<<"Total cell count : "<<C<<", Number of non-placement cells is "<<(C-i)<<endl;
                cerr<<"error: cellSimplePlace fail!"<<endl;
                return -1;
            }
        }
    }

    for(int i = 0; i < N; i++){
        for(int j = 0; j < nets[i].cell.size(); j++){
            int idx = nets[i].cell[j];
            if(cells[idx].inA){
                nets[i].set[0].insert(idx);
            }
            else{
                nets[i].set[1].insert(idx);
            }
        }
    }
    cout<<"Init cut size "<<getCutSize()<<endl;
    return 0;
}


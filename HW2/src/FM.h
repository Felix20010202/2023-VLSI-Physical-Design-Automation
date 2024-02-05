#ifndef FM_H
#define FM_H

#define MODE_COUNT 6

#include<vector>
#include<unordered_map>
#include<unordered_set>
using namespace std;

class doubleList{
public:
    int cellNum;
    doubleList *right = nullptr, *left = nullptr;
    
    doubleList(int num){ cellNum = num; }
    doubleList(){ cellNum = -1; }
};

class bucketList{
public:
    vector<doubleList> list;
    int maxPtr = -1;
    int cellCnt = 0;
    bucketList(int sz){
        list.resize(sz,doubleList());
    };

    void add(int idx, doubleList &node);
    int pop(int idx);
};

class Tech{
public:
    int techNum;
    unordered_map<string, int> techMap;
    vector<int> techSize;
    vector<unordered_map<string, pair<int,int>>> libCell;
};

class Cell{
public:
    bool lock = false;
    string name;
    string libCell;
    vector<int> net;
    bool inA = false;
    int gain = 0;
    int testGain = 0;
    doubleList *ptr = nullptr;
    int cellIdx = 0;
    long long putASize = 0;
    long long putBSize = 0;

    Cell(string &n,string &ln, int ci){
        lock = false;
        name = n;
        libCell = ln;
        cellIdx = ci;
        ptr = new doubleList(cellIdx);
    }
    void updateList(doubleList &head);
};

class Net{
public:
    int pinNum = 0;
    int Distribution[2] = {0,0};
    bool critical = false;
    int lockPin = 0;
    int gain = 0;
    unordered_set<int> set[2];
    vector<int> cell;
    Net(int num){
        pinNum = num;
    }
};


class FM{
public:
    Tech tech;
    int C = 0; // total of cells.
    int N = 0; // total of nets.
    bucketList *listA, *listB;
    pair<long long, long long> dieH_W;
    long long dieArea = 0;
    pair<int, int> dieTech;
    int dieUtil[2];
    long long dieMax[2] = {0,0};
    long long dieUsed[2] = {0,0};
    unordered_map<string, int> cellMap;
    unordered_map<string, int> netMap;
    vector<Cell> cells;
    vector<Net> nets;
    int lockCellNum = 0;
    int pinMax = 0;
    int cellNum[2] = {0, 0};
    int cutSize = -1;
    int maxGain = 1; // cell which connect max number of net
    long long maxCellSize = 0; // except micro

    void fileInput(char *path);
    int cellSimplePlace(int mode);
    int initCellGain();
    void mkBucketList();
    void addCell2List();
    int getCutSize();
    void output(char *path);
    void process(bool &timeOut);
    void updateGain(int base);
    int chooseCell();
    void lock(int base);
    bool test(int base);
    bool testAndSwap(int base);
    void initRound();
    void isNetCritical(int idx);
    void recovery(int base);
    void placeReset();

    // debug
    void showBucketList(int mode);
    void testBucketList();
    void testGain();
    void testNet(int base);

    FM(){};
};

# endif
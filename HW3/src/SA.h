#ifndef SA_H
#define SA_H

#include"DList.h"
#include<vector>
#include<unordered_map>
#include<climits>
#include <cstdlib>
#include <ctime>
using namespace std;

#define SEED 7777
#define TEMPERATURE 5000000000
#define ENDINGT 1
#define DRATE 0.75
#define TENUM 3000;

class Btree{
public:
    int num = -1;
    Btree *left = nullptr, *right = nullptr;
    Btree *parent = nullptr;
    bool inParentR = false;
    Btree(){};
    
    // debug
    static void inorder(Btree *root);
};

class Module{
public:
    string name;
    bool isFix = false;
    long long pos[2] = {-1, -1}; // [0] = x, [1] = y
    long long WH[2] = {-1, -1}; 
    long long area = -1;
    vector<int> nets;
    Btree *ptr;
    bool isUpdate = false;
    bool isPlace = true;

    Module(){
        ptr = new Btree;
    };
};

class Net{
public:
    pair<int, int> md;
    int weight = -1;
    Net(){};
};

class Solution{
public:
    long long wireLen = 0;
    vector<int> X; 
    vector<int> Y;
    vector<int> W;
    vector<int> H;
    long long maxX = 0;
    long long maxY = 0;
    Solution(){};
};

class Record{
public:
    vector<Btree> copyT;
    vector<Btree *> copyPtr;
    vector<vector<long long>> copyWH;
    int copyRootNum;
};

class SA{
public:
    // Input
    int S = 0;
    int F = 0;
    int N = 0;
    int rootNum = 0;
    long long chipSize[2] = {0, 0}; //[0] width, [1] height
    unordered_map<string, int> mdMap;
    vector<Module> modules;
    vector<Net> nets;
    doublyList *head = new doublyList(0, INT_MAX);
    doublyList *headCopy = new doublyList(0, INT_MAX);
    doublyList *vhead = new doublyList(0, INT_MAX);
    doublyList *vheadCopy = new doublyList(0, INT_MAX);
    
    // revert variable
    Record nowR, bestR;
    
    // SA parameter
    long long T = TEMPERATURE;
    long long endingT = TEMPERATURE*0.4;
    double decreaseRate = DRATE;
    long long TE = TENUM;
    long long nowTE = 0;
    long long maxX = 0, maxY = 0;

    Solution best, now, next;

    void fileInput(char *file);
    void outResult(char *file);
    void outBest(char *file);
    void softOnlyInitPlace();
    void process(bool loadBest, long long t, long long et, double d, long long argTE,bool &timeOut);
    bool thermalEquilibrium();
    void copyTree(Record &r);
    void revertTree(Record &r);
    void perturb();
    bool accept(int T, int cost);
    long long cost(Solution &s);
    void getPos(int idx);
    long long getWireLen();
    int update(int pos, int len, int v);
    int update(doublyList *head,  doublyList **headPtr, int pos, int len, int v);
    void update(doublyList *head,  doublyList **headPtr, int pos, int len, int v, int vStart);
    void initList(doublyList *head, bool h);
    void copyList(doublyList *head, doublyList *headCopy);
    bool legal(int idx);
    bool isOverlap(int A, int B);
    void initRound();
    void wtSolution(Solution &s);
    void deleteNode(int idx);
    void insertNode(int B, int A);
    void moveHorizontally(int idx);
    void moveLeft(int idx);
    void moveUp();
    long long getDieSpace();
    void fillUpdirectly();
    bool isResultInBound();
    bool isBestInBound();
    long long getAvgCost();

    // debug
    void showTree();
    void showList(doublyList *head);

    // test
    void isTreeLegal();
    void isPlaceLegal();
    void isListLegal(doublyList *head);

    SA(){};
};

# endif
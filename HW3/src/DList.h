#ifndef DLIST_H
#define DLIST_H

class doublyList{
public:
    int start, end, value;
    doublyList *prev = nullptr, *next = nullptr;
    void insert(int start, int end);
    static void show(doublyList *head);

    doublyList(int s, int e, int v){
        start = s;
        end = e;
        value = v;
    };
    doublyList(int s, int e){
        start = s;
        end = e;
    };
    doublyList(doublyList *l){
        start = l->start;
        end = l->end;
        value = l->value;
    }
};

# endif
#include"DList.h"
#include<iostream>

using namespace std;

void doublyList::show(doublyList *head){
    doublyList *ptr = head;
    while(ptr != nullptr){
        cout<<ptr->value<<" "<<ptr->start<<" "<<ptr->end<<" : ";
        ptr = ptr->next;
    }
    cout<<endl;
    return;
}
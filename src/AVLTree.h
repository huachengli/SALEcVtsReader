//
// Created by huacheng on 3/24/22.
//

#ifndef SALECVTSREADER_AVLTREE_H
#define SALECVTSREADER_AVLTREE_H

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

struct _Node;

struct _Node
{
    struct _Node * lc;
    struct _Node * rc;
    uintptr_t attach;
    uint32_t balance;
};

typedef struct _Node AVLNode;


int BalanceInsert(AVLNode ** root_ptr,AVLNode * item, int (*cmp)(AVLNode *,AVLNode *));
int NodeBalance(AVLNode * node);

void LRotate(AVLNode **);
void RRotate(AVLNode **);
void LRRotate(AVLNode **);
void RLRotate(AVLNode **);
void BalanceLeftSub(AVLNode ** node_ptr);
void BalanceRightSub(AVLNode ** node_ptr);
void AVLDelete(AVLNode ** root_ptr, void (*back)(AVLNode*));

#define LH 0
#define EH 1
#define RH 2

#define TALLER 1
#define BALANCE 0

int CmpInt(AVLNode * a,AVLNode *b);
void BFSTraverse(AVLNode ** root_ptr, int maxnode);
AVLNode * AVLQuery(AVLNode * root,AVLNode * item,int (*cmp)(AVLNode *,AVLNode *));
#endif //SALECVTSREADER_AVLTREE_H

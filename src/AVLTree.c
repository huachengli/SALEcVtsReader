//
// Created by huacheng on 3/24/22.
//

#include "AVLTree.h"

int NodeBalance(AVLNode * node)
{
    return (int ) (node->balance);
}

void RRotate(AVLNode ** node_ptr)
{
   /*
    *                A                       B
    *              /   \                   /  \
    *             B    C      ==>         D    A
    *           /  \                    /  \  / \
    *          D    E                  F   G E  C
    *        /  \
    *       F    G
    */

   AVLNode * node = *node_ptr;
   AVLNode * nA = node;
   AVLNode * nB = node->lc;
   AVLNode * nE = nB->rc;

   nA->lc = nE;
   nB->rc = nA;
   *node_ptr = nB;
}

void LRotate(AVLNode ** node_ptr)
{
    /*
     *                A                       C
     *              /   \                   /  \
     *             B    C      ==>         A    E
     *                 /  \              /  \  / \
     *                D    E            B   D F  G
     *                    /  \
     *                   F    G
     */

    AVLNode * node = *node_ptr;
    AVLNode * nA = node;
    AVLNode * nC = node->rc;
    AVLNode * nD = nC->lc;

    nA->rc = nD;
    nC->lc = nA;
    *node_ptr = nC;
}

void LRRotate(AVLNode ** node_ptr)
{
    /*
     *                A                       A               E
     *              /   \                   /  \            /  \
     *             B    C      ==>         E    C   ==>    B    A
     *           /  \                    /  \             / \  / \
     *          D    E                  B   G            D  F G  C
     *              /  \               / \
     *             F    G             D   F
     */

    // left rotate on node(B)
    LRotate(&(*node_ptr)->lc);
    // right rotate on node(E)
    RRotate(node_ptr);
}

void RLRotate(AVLNode ** node_ptr)
{
    RRotate(&(*node_ptr)->rc);
    LRotate(node_ptr);
}

void BalanceLeftSub(AVLNode ** node_ptr)
{
    /*
     *              A
     *            /  \
     *           B    C
     *          /  \
     *        ?(D) ?(E)
     */
    AVLNode * nA = *node_ptr;
    AVLNode * nB = (*node_ptr)->lc;
    AVLNode * nE =  nB->rc;
    switch(nB->balance)
    {
        case LH:
            nA->balance = nB->balance = EH;
            RRotate(node_ptr);
            break;
        case RH:
            switch (nE->balance) {
                case LH:
                    nA->balance = RH;
                    nB->balance = EH;
                    break;
                case EH:
                    nA->balance = EH;
                    nB->balance = EH;
                    break;
                case RH:
                    nA->balance = EH;
                    nB->balance = LH;
                    break;
            }
            nE->balance = EH;
            LRRotate(node_ptr);
            break;
        case EH:
            break;
    }
}

void BalanceRightSub(AVLNode ** node_ptr)
{
    /*
     *             A
     *           /   \
     *          B     C
     *               /  \
     *              ?D  ?E
     */
    AVLNode * nA = *node_ptr;
    AVLNode * nC = (*node_ptr)->rc;
    AVLNode * nD = nC->lc;
    switch (nC->balance)
    {
        case LH:
            switch (nD->balance)
            {
                case LH:
                    nA->balance = EH;
                    nC->balance = RH;
                    break;
                case EH:
                    nA->balance = EH;
                    nC->balance = EH;
                    break;
                case RH:
                    nA->balance = LH;
                    nC->balance = EH;
                    break;
            }
            nD->balance = EH;
            RLRotate(node_ptr);
            break;
        case RH:
            nA->balance = nC->balance = EH;
            LRotate(node_ptr);
            break;
        case EH:
            break;
    }
}



int BalanceInsert(AVLNode ** root_ptr,AVLNode * item, int (*cmp)(AVLNode *,AVLNode *))
{
    if(NULL == (*root_ptr))
    {
        *root_ptr = item;
        item->lc = item->rc = NULL;
        item->balance = EH;
        return TALLER;
    }

    int res = cmp(*root_ptr,item);
    if(res < 0)
    {
        // insert to left sub tree
        int status = BalanceInsert(&(*root_ptr)->lc,item,cmp);
        int rtn_status = BALANCE;
        if(TALLER == status)
        {
            switch ((*root_ptr)->balance)
            {
                case LH:
                    BalanceLeftSub(root_ptr);
                    rtn_status = BALANCE;
                    break;
                case EH:
                    (*root_ptr)->balance = LH;
                    rtn_status = TALLER;
                    break;
                case RH:
                    (*root_ptr)->balance = EH;
                    rtn_status = BALANCE;
                    break;
            }
        }

        return rtn_status;
    } else if(res >0)
    {
        // insert to right sub tree
        int status = BalanceInsert(&(*root_ptr)->rc,item,cmp);
        int rtn_status = BALANCE;
        if(TALLER == status)
        {
            switch ((*root_ptr)->balance)
            {
                case LH:
                    (*root_ptr)->balance = EH;
                    rtn_status = BALANCE;
                    break;
                case EH:
                    (*root_ptr)->balance = RH;
                    rtn_status = TALLER;
                    break;
                case RH:
                    BalanceRightSub(root_ptr);
                    rtn_status = BALANCE;
                    break;
            }
        }
        return rtn_status;
    }
    else
    {
        return BALANCE;
    }
}

int CmpInt(AVLNode * a,AVLNode *b)
{
    int data_a = *(int *)(a->attach);
    int data_b = *(int *)(b->attach);
    return data_b - data_a;
}

void BFSTraverse(AVLNode ** root_ptr, int maxnode)
{
    uintptr_t * dequeue = (uintptr_t *) malloc(sizeof(uintptr_t)*maxnode);
    int head=0,tail=1,len=1;
    dequeue[0] = (uintptr_t) (*root_ptr);

    while (len > 0)
    {

        AVLNode * cur = (AVLNode *) dequeue[head];
        head = (head+1)%maxnode;
        len -= 1;

        if(NULL == cur)
        {
            fprintf(stdout," ? ");
        } else
        {
            fprintf(stdout," %d ",*(int *)cur->attach);
            dequeue[tail] = (uintptr_t)(cur->lc);
            tail = (tail+1)%maxnode;
            dequeue[tail] = (uintptr_t)(cur->rc);
            tail = (tail+1)%maxnode;
            len += 2;
        }
    }
    free(dequeue);
}

void AVLDelete(AVLNode ** root_ptr, void (*back)(AVLNode*))
{
    AVLNode * cur = *root_ptr;
    if(NULL != cur->lc) AVLDelete(&cur->lc, back);
    if(NULL != cur->rc) AVLDelete(&cur->rc, back);
    if(NULL != back) back(cur);
    free(cur);
}

AVLNode * AVLQuery(AVLNode * root,AVLNode * item,int (*cmp)(AVLNode *,AVLNode *))
{
    if(NULL == root)
        return NULL;
    int res_cmp = cmp(root,item);
    if(res_cmp > 0)
    {
        return AVLQuery(root->rc,item,cmp);
    } else if(res_cmp < 0)
    {
        return AVLQuery(root->lc,item,cmp);
    } else
    {
        return root;
    }
}


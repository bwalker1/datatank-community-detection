//
//  KDTree.cpp
//  Clustering
//
//  Created by Benjamin Walker on 3/21/17.
//
//

#include "KDTree.h"
#include <math.h>
#include <algorithm>
#include <assert.h>

float distP(const pt3d &l, const pt3d &r)
{
    float ret = 0;
    for (int dim=0;dim<3;++dim)
    {
        float temp =l[dim]-r[dim];
        ret += temp*temp;
    }
    return sqrtf(ret);
}

void swap(int &l, int &r)
{
    int c = l;
    l = r;
    r = c;
}
// rvalue
float KDTree::distAccess(int id1, int id2) const
{
    assert(n*id1 + id2 < n*n);
    //turn data[n*id1 + id2];
    return C[id1][id2];
}
// lvalue
float& KDTree::distAccess(int id1, int id2)
{
    assert(n*id1 + id2 < n*n);
    //return data[n*id1 + id2];
    return C[id1][id2];
}

inline float KDTree::dist(KDNode* &l, KDNode* &r) const
{
    
    // just compute it when requested
    return distP(l->center,r->center);
}


pt3d operator+(const pt3d &l, const pt3d &r)
{
    pt3d out;
    
    for (int dim=0;dim<3;++dim)
    {
        out[dim] = l[dim];
        out[dim] += r[dim];
    }
    return out;
}

pt3d operator/=(pt3d &l, const float &r)
{
    for (int dim=0;dim<3;++dim)
    {
        l[dim]/=r;
    }
    return l;
}


KDTree::KDTree()
{
    // maybe don't do this
    exit(1);
    head = NULL;
    _size = 0;
}

// put elements into tree
// currently should be O(n log^2 n)
KDTree::KDTree(vector<pt3d> points)
{
    pts = points;
    _size = points.size();
    n = points.size();
    /*
    data = (float*) malloc(n*n*sizeof(float));*/
#ifdef checkNN
    C.resize(n);
    for (int i=0;i<n;++i)
    {
        C[i].resize(n);
    }
    for (int i=0;i<n;++i)
    {
        for (int j=0;j<n;++j)
        {
            distAccess(i, j) = distP(points[i],points[j]);
        }
    }
#endif
    head = createTree(points.begin(), points.end(), 0, NULL,false,points.begin());
}

typedef struct
{
    int dim;
    bool operator()(const pt3d &l, const pt3d &r)
    {
        //return ((float*)&l)[dim]<((float*)&r)[dim];
        return l[dim]<r[dim];
    };
} dimComp;

KDNode * createTree(vector<pt3d>::iterator begin, vector<pt3d>::iterator end, int dim, KDNode *parent, bool isLeft, vector<pt3d>::iterator absolute_begin)
{
    if (end==begin)
    {
        return NULL;
    }
    KDNode * ret = new KDNode;
    ret->parent = parent;
    ret->isLeft = isLeft;
    // base case
    if (end-begin == 1)
    {
        ret->elements = 1;
        ret->center = *begin;
        ret->id = begin->id;
        ret->dim = dim;
        ret->lchild = NULL;
        ret->rchild = NULL;
        return ret;
    }
    // choose splitting point
    // in reality we should select splitting point from fixed-size subset
    // sort along specified dimension
    dimComp comp;
    comp.dim = dim;
    sort(begin,end,comp);
    
    auto split = begin + (end-begin)/2;
    // create tree
    ret->center = *split;
    ret->elements = 1;
    ret->dim = dim;
    // add children
    ret->lchild = NULL;
    ret->rchild = NULL;
    KDNode * child1 = createTree(begin,split,(dim+1)%3,ret,true,absolute_begin);
    KDNode * child2 = createTree(split+1,end,(dim+1)%3,ret,false,absolute_begin);

    ret->lchild = child1;
    ret->rchild = child2;
    ret->id = split->id;

    return ret;
}

void KDTree::deleteNode(KDNode * node)
{
    if (!node)
        return;
    deleteNode(node->lchild);
    deleteNode(node->rchild);
    delete node;
}


// access an element. For convenience, we return a leaf (although that does make this O(logn))
KDTreeNodeAccess KDTree::getElement() const
{
    KDNode * node = head;
    //assert(node);
    while (true)
    {
        if (node->lchild)
        {
            node = node->lchild;
        } else if (node->rchild)
        {
            node = node->rchild;
        } else
        {
            break;
        }
    }
    KDTreeNodeAccess ret;
    ret.p = node;
    ret.tree = this;
    return ret;
}

KDTree::~KDTree()
{
    deleteNode(head);
    //free(data);
}

bool verifyTree(KDNode *head)
{
    bool res = true;
    assert(head);
    assert(head->id == head->center.id);

    if (head->lchild)
    {
        res = res && (head->lchild->parent == head);
        assert(res);
        res = res && (verifyTree(head->lchild));
        assert(res);
    }
    if (head->rchild)
    {
        res = res && (head->rchild->parent == head);
        assert(res);
        res = res && (verifyTree(head->rchild));
        assert(res);
    }
    return res;
}

// the meat: nearest neighbor finding algorithm
// return a pointer to the closest node found. If none are found closer than min, return null
KDNode* KDTree::findNearestNeighbor(KDNode* cur_head, KDNode* node,float &min) const
{
    assert(cur_head);
    assert(node);
    // check the head node
    float head_dist = dist(cur_head,node);

    KDNode* out = NULL;
    if (head_dist < min && cur_head != node)
    {
        out = cur_head;
        min = head_dist;
    }
    if (!cur_head->lchild && !cur_head->rchild)
    {
        return out;
    }
    KDNode *first, *second;
    // check first the side of head that node is on
    //assert(node&&cur_head);
    float splitDimDist = node->center[cur_head->dim]-cur_head->center[cur_head->dim]; // separation along dimension of split
    if (splitDimDist < 0)
    {
        first  = cur_head->lchild;
        second = cur_head->rchild;
    }
    else
    {
        first  = cur_head->rchild;
        second = cur_head->lchild;
    }
    // check the child trees (if they exist)
    if (first)
    {
        KDNode* firstNew = findNearestNeighbor(first, node, min);
        if (firstNew)
        {
            out = firstNew;
        }
    }
    if (second)
    {
        // check whether it's even worth bothering looking at the far side by computing distance between point and the hyperplane split associated with head
        float min_sep = fabsf(splitDimDist);
        
        // sometimes this is discounting valid matches - is the tree well-formed?
        if (min_sep < min)
        {
            // it is worth checking
            KDNode * secondNew = findNearestNeighbor(second, node, min);
            if (secondNew)
            {
                out = secondNew;
            }
        }
    }
    return out;
}

inline void trim(set<knnV,knncomp> &found, float &kthbest, int k)
{
    if (found.size()>k)
    {
        auto it = next(found.rbegin());
        kthbest = it->first;
        found.erase(prev(found.end()));
    }
}


// generalization of nearest-neighbor search to allow finding multiples
void KDTree::pKNN(KDNode* cur_head, KDNode* node, int k, set<knnV,knncomp> &found, float &kthbest )
{
    KDNode* first, *second;
    // check first the side of head that node is on
    //assert(node&&cur_head);
    float splitDimDist = node->center[cur_head->dim]-cur_head->center[cur_head->dim]; // separation along dimension of split
    if (splitDimDist < 0)
    {
        first  = cur_head->lchild;
        second = cur_head->rchild;
    }
    else
    {
        first  = cur_head->rchild;
        second = cur_head->lchild;
    }
    
    // check close child
    if (first)
    {
        pKNN(first, node, k,found,kthbest);
    }
    
    // check the head node
    float head_dist = dist(cur_head,node);
    
    if (head_dist < kthbest && cur_head != node)
    {
        // add to our set
        found.insert(knnV(head_dist,cur_head));
        trim(found,kthbest,k);
    }
    
    // check far child
    if (second)
    {
        // check whether it's even worth bothering looking at the far side by computing distance between point and the hyperplane split associated with head
        float min_sep = fabsf(splitDimDist);
        
        if (min_sep < kthbest)
        {
            // it is worth checking
            pKNN(second, node, k,found,kthbest);
        }
    }
}

vector<KDTreeNodeAccess> KDTree::KNN(KDTreeNodeAccess in, int k)
{
    knncomp comp;
    set<knnV,knncomp> found;
    float v = INFINITY;
    pKNN(this->head,in.p,k,found,v);
    
    // get the result into the desired form
    vector<KDTreeNodeAccess> out(k);
    auto it = found.begin();
    for (int i=0;i<k;++i)
    {
        KDTreeNodeAccess v;
        v.tree =this;
        v.p = it->second;
        out[i]=v;
        it = next(it);
    }
    return out;
}

// do the standard O(n) (shhh maybe O(n logn) here) nearest neighbor just to compare
int KDTree::naiveNearestNeighbor(KDTreeNodeAccess node, float &dist) const
{
#ifndef checkNN
    assert("Need to enable NN checking to run this\n"&&false);
#endif
    int id = node.p->id;
    float min = INFINITY;
    int min_id = -1;
    for (int i=0;i<n;++i)
    {
        if (i==id) continue;
        if (distAccess(id, i)<min)
        {
            min_id = i;
            min = distAccess(id, i);
        }
    }
    dist = min;
    return min_id;
}



KDTreeNodeAccess KDTree::findNearestNeighbor(KDTreeNodeAccess node) const
{
    KDNode* in = node.p;
    float min = INFINITY;
    KDNode* neighbor = findNearestNeighbor(head, in, min);
    assert(neighbor);
#ifdef checkNN
    float distv = -1;
    int naive = naiveNearestNeighbor(node, distv);
    assert(neighbor->id == naive);
#endif
    KDTreeNodeAccess ret;
    ret.p = neighbor;
    ret.tree = this;
    return ret;
}

KDNode* pfind(KDNode* node, int id)
{
    KDNode* res;
    if (node->id == id)
    {
        return node;
    }
    else
    {
        if (node->lchild)
        {
            res = pfind(node->lchild,id);
            if (res) return res;
        }
        if (node->rchild)
        {
            res = pfind(node->rchild,id);
            if (res) return res;
        }
    }
    return NULL;
}

KDTreeNodeAccess KDTree::find(int id)
{
    KDNode* node = pfind(head,id);
    assert(node);
    KDTreeNodeAccess out;
    out.p = node;
    out.tree = this;
    return out;
}

vector<int> idsContained(KDNode* head)
{
    //assert(head);
    vector<int> res;
    res.push_back(head->id);
    if (head->lchild)
    {
        vector<int> sub = idsContained(head->lchild);
        res.insert(res.end(),sub.begin(),sub.end());
    }
    if (head->rchild)
    {
        vector<int> sub = idsContained(head->rchild);
        res.insert(res.end(),sub.begin(),sub.end());
    }
    return res;
}

vector<int> KDTree::containedIds()
{
    vector<int> res = idsContained(this->head);
    sort(res.begin(),res.end());
    return res;
}
bool _containsID(KDNode* const head, int id)
{
    //assert(head);
    bool res = (head->id==id)||(head->lchild&&_containsID(head->lchild, id))||(head->rchild&&_containsID(head->rchild, id));
    //if (res) printf("%d\n",head->id);
    return res;
}
bool KDTree::containsID(int id) const
{
    return _containsID(head,id);
}

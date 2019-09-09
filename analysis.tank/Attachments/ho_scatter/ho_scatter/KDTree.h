//
//  KDTree.hpp
//  Clustering
//
//  Created by Benjamin Walker on 3/21/17.
//
//

#ifndef KDTree_hpp
#define KDTree_hpp

#include <stdio.h>
#include <vector>
#include <set>

// turn this on to debug and make sure it agrees with naive implementation
//#define checkNN 1

using namespace std;

class KDTree;

typedef struct
{
    //float x,y,z;
    float data[3];
    int id;
    // rvalue access
    float operator[] (int dim) const
    {
        return data[dim];
    }
    // lvalue access
    float& operator[] (int dim)
    {
        return data[dim];
    }
} pt3d;


// pt3d functions
pt3d operator+(const pt3d &l, const pt3d &r);
pt3d operator/=(pt3d &l, const float &r);
float distP(const pt3d &l, const pt3d &r);

typedef struct KDNodeStruct
{
    KDNodeStruct *lchild, *rchild,*parent;
    pt3d center;
    int elements;
    bool isLeft;
    int id;
    int dim;    // the dimension children are split along
} KDNode;

typedef vector<pair<int, double> > Cluster;

bool verifyTree(KDNode *head);

// struct used to access nodes within the tree
typedef struct KDTreeNodeAccessStruct
{
private:
    KDNode * p;
public:
    int id() {return p->id; };
    const KDTree * tree;
    friend class KDTree;
    bool operator==(const KDTreeNodeAccessStruct &other)
    {
        return this->p == other.p;
    }
} KDTreeNodeAccess;

typedef pair<float,KDNode*> knnV;

struct knncomp
{
    bool operator()(const knnV &l, const knnV &r)
    {
        return l.first < r.first;
    }
};



bool operator==(const KDTreeNodeAccess &l, const KDTreeNodeAccess &r);

KDNode * createTree(vector<pt3d>::iterator begin, vector<pt3d>::iterator end, int dim, KDNode *parent, bool isLeft, vector<pt3d>::iterator absolute_begin);
KDTreeNodeAccess merge(KDTreeNodeAccess a, KDTreeNodeAccess b);
KDNode* find(KDNode* node, int id);

class KDTree
{
protected:
    KDNode* head;
    int _size, n;
    void deleteNode(KDNode * node);
    KDNode* findNearestNeighbor(KDNode* head, KDNode* node,float &min) const;
    void pKNN(KDNode* head, KDNode* node, int k, set<knnV,knncomp> &found, float &kthbest);
    //void remove(KDNode *node);
    //void mergeTable(int a, int b);
    // node distance data (for easy complete-linkage calculation)
    vector<vector<float> > C;
    float * data;
    
    
    vector<pt3d> pts;
    
    
    // record of merges
public:
    vector<pair<int,int> > merges;
public:
    // access internal distance storage
    float& distAccess(int id1, int i2);
    float distAccess(int id1, int i2) const;
    //void nextMerge();
    //void cluster();
    bool verify() { return verifyTree(this->head); };
    int size()const { return _size; };
    //void remove(KDTreeNodeAccess node);
    KDTreeNodeAccess insert(Cluster cluster);
    KDTreeNodeAccess findNearestNeighbor(KDTreeNodeAccess node) const;
    vector<KDTreeNodeAccess> KNN(KDTreeNodeAccess node, int k);
    int naiveNearestNeighbor(KDTreeNodeAccess node, float & dist) const;
    KDTreeNodeAccess getElement() const;
    bool containsID(int id) const;
    float dist(KDNode* &l, KDNode* &r) const;
    vector<int> containedIds();
    //KDTreeNodeAccess merge(KDTreeNodeAccess &a, KDTreeNodeAccess &b);
    KDTreeNodeAccess find(int id);
    
    void balance();     // not implemented but maybe we want to rebalance the tree later on?
    
    // constructor
    KDTree();
    KDTree(vector<pt3d> points);
    
    // destructor
    ~KDTree();
};

#endif /* KDTree_hpp */

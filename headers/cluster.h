#ifndef CLUSTER
#define CLUSTER

#include <vector>
#include "image.h"
#include <set>

class cluster{
    private:
        set<image*>* clust;
    public:
        cluster();
        ~cluster();
        void insert(image*);
        int getClusterSize();
        void resetCluster();
        image* getNewCentroid(int);
        set<image*>* getClust();

};

#endif
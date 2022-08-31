#ifndef CLUSTER_ENTRY
#define CLUSTER_ENTRY

class clusterEntry{
    private:
        int clusterID;
        bool undecided;
        unsigned int clusterDistance;
    public:
        clusterEntry();
        ~clusterEntry();
        int getClusterID();
        bool isConsidered();
        unsigned int getDistanceToCluster();
        void setClusterID(int);
        void setDecidability(bool);
        void setDistance(unsigned int );
};

#endif
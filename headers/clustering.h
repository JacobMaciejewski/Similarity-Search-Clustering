#ifndef CLUSTERING
#define CLUSTERING

#include <vector>
#include "image.h"
#include "cluster.h"
#include <string>
#include "image.h"
#include "algoAux.h"
#include "hyperCube.h"
#include "LSH.h"
#include "clusterEntry.h"

using namespace std;

class clustering{
    private:
        int k; //number of clusters/centroids
        vector<cluster*> *clusters;
        vector<image*> *centroids;
        vector<image*> *imageArray;
        LSH *clusterLSH;
        hyperCube *clusterCube;
        string* method;
        vector<clusterEntry*>* entries; 
    public:
        clustering(int, int, unsigned int, unsigned int, string*, vector<image*>*);
        ~clustering();
        void nearestImagesInfo(image *targetImage, vector<image*>*, unsigned int*, int*, unsigned int (*) (vector<unsigned char>*, vector<unsigned char>*));
        void addCentroid(image*);
        void updateIthCentroid(int, image*);
        void initializeCentroids(unsigned int (*) (vector<unsigned char>*, vector<unsigned char>*));
        void Lloyd(unsigned int (*) (vector<unsigned char>*, vector<unsigned char>*));
        void updateDistanceArray(int, unsigned int*, image*, bool[], unsigned int[], unsigned int (*) (vector<unsigned char>*, vector<unsigned char>*));
        void updateProbArray(int, unsigned int, bool[], unsigned int[], double[]);
        void newCentroids();
        void getCentroidDistances(unsigned int*, unsigned int*, unsigned int (*) (vector<unsigned char>*, vector<unsigned char>*));
        void resetEntries();
        void resetClustering();
        double executeMethod(int, int, unsigned int, unsigned int, unsigned int (*) (vector<unsigned char>*, vector<unsigned char>*));
        double execute(unsigned int (*) (vector<unsigned char>*, vector<unsigned char>*));
        template<typename I> void clusterizeData(unsigned int (*) (vector<unsigned char>*, vector<unsigned char>*));
        int getNearClusterID(int, image*,  unsigned int (*) (vector<unsigned char>*, vector<unsigned char>*));
        int getMaxClusterSize();
        vector<double>* calculateSilhouettes(unsigned int (*) (vector<unsigned char>*, vector<unsigned char>*));
        double calculateIthSilhouette(int, unsigned int (*) (vector<unsigned char>*, vector<unsigned char>*), double**);
        double calculateSilhouette(unsigned int (*) (vector<unsigned char>*, vector<unsigned char>*), vector<double>**);
        template<typename I> void reverseAssignment(I, unsigned int (*) (vector<unsigned char>*, vector<unsigned char>*));
        vector<cluster*>* getClusters();
        cluster* getIthCluster(int);
        vector<image*>* getCentroids();
        image* getIthCentroid(int);
        vector<image*>* getImages();
        image* getIthImage(int);
        int getClustersNum();
        string *getMethod();
        void resetClusters();
        vector<clusterEntry*>* getEntries();
        clusterEntry* getIthEntry(int);
        hyperCube *getHyperCube();
        LSH *getLSH();
};

#endif
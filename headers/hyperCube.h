#ifndef HYPERCUBE
#define HYPERCUBE

#include <vector>
#include "hashTable.h"
#include "algoAux.h"
#include <cmath>
#include <random>
#include <algorithm>
#include <queue>
#include <set>

using namespace std;


class hyperCube{
    private:
        hashTable* cubeHashTable;
        vector<vector <unsigned char> *> *projectionFunctions;
        // max number of edges of hypercube and elements checked
        unsigned int edgeCeil;
        unsigned int checkCeil;
    public:
        hyperCube(unsigned int, unsigned int);
        ~hyperCube();
        hashTable* getHashTable();
        vector<vector <unsigned char> *>* getProjectionFuncs();
        vector<unsigned char>* getIthProjectionFunc(int);
        void insert(image *);
        void initializeHyperCube(vector<image*>*);
        unsigned int hashFunction(image *image);
        unsigned int getEdgeCeil();
        unsigned int getCheckCeil();
        vector<image*>* getNN(image*, int, unsigned int (*)(vector<unsigned char>*, vector<unsigned char>*));
        vector<image*>* radiusSearch(image*, unsigned int, unsigned int (*) (vector<unsigned char>*, vector<unsigned char>*), set<int>&);
};



// class hyperCube{
// private:
//     //the m is equal to 2 ^ 23 - 5
//     static unsigned int bigM;
//     static unsigned int smallM;
//     //array storing the result of m exponentials in h_i calculations
//     //can be reused in each h function execution, as they are static
//     static unsigned int *mExponentials;
//     static float wAspect;
//     //total number of h functions
//     static int distrVecNum;

//     unsigned int k;         //dimension of hypercube
//     static unsigned int d;         //dimension of an image
//     float r;                //radius of search
//     unsigned int hMod;  
//     vector<bucket *> *bucketArray; //bucket of an vertex
//     vector<vector <unsigned char> *> *projectionFunctions;
//     vector<vector<float> *> *disturbanceVectors; 
// public:
//     hyperCube(unsigned int k, unsigned int d, float r);
//     ~hyperCube();

//     void insert(image *);
//     unsigned int hashFunction(image *image);
//     image *NearestNeighbor(image *queryImage);
//     vector<unsigned int> *hammingNeighbors(unsigned int node);
//     void initializeHyperCube(vector<image*>*);


//     static void initMExponentials();
//     static void setMExponentials(unsigned int *);
//     unsigned int hFunction(int, image*);
//     unsigned int modularExponentiation(unsigned int, unsigned int, int);
//     unsigned int hModulo(vector<unsigned int>*, int);
    
//     static unsigned int getSmallAspectM();
//     static unsigned int getBigAspectM();
//     static unsigned int *getMExponentials();
//     static unsigned int getIthMExponential(int);
//     static float getAspectW();
//     static int getDistrVecNum();
//     static unsigned int getDimensionality();
//     vector<vector<float>*>* getDistrVectors();
//     vector<float>* getIthDistrVector(int);
// };


#endif
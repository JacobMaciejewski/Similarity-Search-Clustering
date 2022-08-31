#ifndef HASH_TABLE
#define HASH_TABLE

#define BUCKETS 460
#include "../headers/mainAux.h"
#include "../headers/image.h"
#include "../headers/bucket.h"


class hashTable{
    private:
        //the m is equal to 2 ^ 23 - 5
        static unsigned int bigM;
        static unsigned int smallM;
        //array storing the result of m exponentials in h_i calculations
        //can be reused in each h function execution, as they are static
        static unsigned int *mExponentials;
        static float wAspect;
        //total number of h functions
        static int distrVecNum;
        //number of dimensions of image vectors
        static int imgDimensions;
        // number of buckets
        static unsigned int bucketsNum;
        vector<vector<float>*> *distrVectors;
        bucket** buckets;

    public:
        hashTable();
        hashTable(unsigned int);
        ~hashTable();
        void initDistrVectors(); 
        static void initMExponentials();
        static void initMExponentials(int, int, int);
        static void setSmallM(unsigned int);
        static void setBigM(unsigned int );
        static void setWAspect(float);
        static void setDistrVecNum(unsigned int);
        static void setDimensionality(unsigned int);
        static void setMExponentials(unsigned int *);
        static void setBucketsNum(unsigned int);
        unsigned int hFunction(int, image*);
        unsigned int gFunction(image*);
        unsigned int modularExponentiation(unsigned int, unsigned int, int);
        unsigned int hModulo(vector<unsigned int>*, int);
        static unsigned int getSmallAspectM();
        static unsigned int getBigAspectM();
        static unsigned int *getMExponentials();
        static unsigned int getIthMExponential(int);
        static unsigned int getBucketsNum();
        static float getAspectW();
        static int getDistrVecNum();
        static int getDimensionality();
        unsigned int getBucket(image*);
        vector<vector<float>*>* getDistrVectors();
        vector<float>* getIthDistrVector(int);
        void pushImage(image*);
        void pushImages(vector<image*>*);
        void printBucketSizes();
        void printFirstNBucketSizes(int);
        bucket** getBuckets();
        bucket* getIthBucket(int);
};


#endif
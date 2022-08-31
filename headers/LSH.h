#ifndef LSH_
#define LSH_

#include <vector>
#include <queue>
#include "image.h"
#include "hashTable.h"
#include "imageDuplet.h"
#include "algoAux.h"
#include <set>

class LSH{
    private:
        vector<hashTable*>* hashTables;
        int hashTablesNum;
    public:
        LSH(int);
        ~LSH();
        void initializeLSH(vector<image*>*);
        void printHashBucketSizes();
        int getHashTablesNum();
        vector<hashTable*>* getHashTables();
        hashTable* getIthHashTable(int);
        vector<image*>* getNN(image*, int, unsigned int (*)(vector<unsigned char>*, vector<unsigned char>*));
        vector<image*>* radiusSearch(image*, unsigned int, unsigned int (*) (vector<unsigned char>*, vector<unsigned char>*), set<int>&);
};


#endif
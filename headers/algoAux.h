#ifndef ALGO_AUX
#define ALGO_AUX

#include "imageDuplet.h"
#include "hashTable.h"
#include <vector>
#include <queue>
#include <cmath>
#include <set>

#define MAX_CLUSTERIZATIONS 20

using namespace std;

struct compDuplet{ 
    bool operator()(imageDuplet const& p1, imageDuplet const& p2) 
    { 
        return p1.distance < p2.distance; 
    } 
};

// bool getBucketNN(image*, int, int, hashTable*, unsigned int*, unsigned int, set<int, greater<int>>&, priority_queue<imageDuplet, vector<imageDuplet>>&, unsigned int (*)(vector<unsigned char>*, vector<unsigned char>*));
unsigned int manhattanDistance(vector<unsigned char>*, vector<unsigned char>*);
void getBinaryNeighbours(unsigned int, vector<unsigned int>*);
template<typename I> I random_element(I, I);
double randomDouble(double, double);
double binarySearch(double*, int, int, double); 
#endif
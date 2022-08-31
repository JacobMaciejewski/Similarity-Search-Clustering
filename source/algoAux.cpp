#include "../headers/algoAux.h"
#include "../headers/imageDuplet.h"

template <typename I>
I random_element(I begin, I end)
{
    const unsigned long n = std::distance(begin, end);
    const unsigned long divisor = (RAND_MAX + 1) / n;

    unsigned long k;
    do { k = std::rand() / divisor; } while (k >= n);

    std::advance(begin, k);
    return begin;
}

unsigned int manhattanDistance(vector<unsigned char>* vec1, vector<unsigned char>* vec2)
{
    int vecSize = vec1->size();
    unsigned int distance = 0;

    for(int dim = 0; dim < vecSize; dim++)
    {
        distance += abs(vec1->at(dim) - vec2->at(dim));
    }
    return distance;
}

double randomDouble(double a, double b)
{
    double random = ((double) rand()) / (double) RAND_MAX;
    double diff = b - a;
    double r = random * diff;
    return a + r;
}

// returns the upper index of the area in which the x is contained
double binarySearch(double* arr, int l, int r, double x) 
{ 
    if (r >= l) { 
        int mid = l + (r - l) / 2; 
  
        // If the element is present at the middle 
        // itself 
        if (arr[mid] < x && x <= arr[mid + 1]) 
            return mid; 
  
        // If element is smaller than mid, then 
        // it can only be present in left subarray 
        if (arr[mid] > x) 
            return binarySearch(arr, l, mid - 1, x); 
  
        // Else the element can only be present 
        // in right subarray 
        return binarySearch(arr, mid + 1, r, x); 
    } 
  
    // We reach here when element is not 
    // present in array 
    return -1; 
} 
  

// producing an integer with single true bit
// shifting it in each iteration and applying xor to the current neighbour
// binary representation
// in this way we get all the neighbours with Hamming distance 1
void getBinaryNeighbours(unsigned int givenNode, vector<unsigned int> *neighbours)
{
    int neighNum = neighbours->size();
    unsigned int currentNeighbour;

    for(int i = 0; i < neighNum; i++)
    {
        currentNeighbour = givenNode ^ (1 << i);
        neighbours->at(i) = currentNeighbour;
    }
    return;
}


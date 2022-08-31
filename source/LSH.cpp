#include "../headers/LSH.h"




LSH::LSH(int hashTablesNum){
    this->hashTablesNum = hashTablesNum;
    
    this->hashTables = new vector<hashTable*>(this->hashTablesNum,NULL);
    
    for(int i = 0; i < this->hashTablesNum; i++){
        hashTables->at(i) = new hashTable();
    }
}

LSH::~LSH(){
    for(int i = 0; i < this->hashTablesNum; i++){
        delete(hashTables->at(i));
    }

    delete(hashTables);
}

void LSH::initializeLSH(vector<image*> *images){
    int totalHashTables = this->getHashTablesNum();

    for(int currHash = 0; currHash < totalHashTables; currHash++)
    {
        this->getIthHashTable(currHash)->pushImages(images);
    }

    // this->printHashBucketSizes();
    return;
}

void LSH::printHashBucketSizes()
{
    int hashesNum = this->getHashTablesNum();
    hashTable *currentHashTable;

    for(int i = 0; i < hashesNum; i++)
    {
        cout << i + 1 << "th HashTable: " << endl;
        currentHashTable = this->getIthHashTable(i);
        currentHashTable->printBucketSizes();
    }
}

int LSH::getHashTablesNum()
{
    return this->hashTablesNum;
}

vector<hashTable*>* LSH::getHashTables()
{
    return this->hashTables;
}

hashTable* LSH::getIthHashTable(int index)
{
    vector<hashTable*>* hashTablesVector = this->getHashTables();

    return hashTablesVector->at(index);
}



vector<image*>* LSH::getNN(image *query, int k, unsigned int (*distanceFunction) (vector<unsigned char>*, vector<unsigned char>*))
{
    // currently updated variables
    hashTable *currHashTable;
    bucket *currBucket;
    image *currentImage;
    // the vector representation of the given query
    vector<unsigned char>* queryVector = query->getImage();
    int queryID = query->getId();
    // the vector representation of currently checked image
    vector<unsigned char>* currentVector;
    int currentID;

    // getting the total number of hash tables in the lsh
    int hashesNum = this->getHashTablesNum();
    // number of images in the current bucket
    int bucketImagesNum;
    // the distance between current vector and the query one
    unsigned int currentDistance;
    // bucket number in current hash table
    unsigned int currBucketID;
    // priority queue used as a heap
    priority_queue<imageDuplet, vector<imageDuplet>, compDuplet> imagePool;
    // vector containing nearest images
    vector<image*> *nearestNeighbours = new vector<image*>();
    // visited vectors' ids, no need to check them again
    set<int, greater<int> > visitedIDs;
    
    // for each hashTable hash the given image and 
    // traverse its bucket searching for neighbours
    for(int currHash = 0; currHash < hashesNum; currHash++)
    {
        // getting the info about the bucket into which
        // given image is hashed in current hash table
        currHashTable = this->getIthHashTable(currHash);
        currBucketID = currHashTable->getBucket(query);
        currBucket = currHashTable->getIthBucket(currBucketID);
        bucketImagesNum = currBucket->getSize();

        for(int img = 0; img < bucketImagesNum; img++)
        {
            currentImage = currBucket->getIthImage(img);
            currentVector = currentImage->getImage();
            currentID = currentImage->getId();

            // checking for neighbours update only when new one found and it is different to the query one
            if(currentID != queryID && visitedIDs.find(currentID) == visitedIDs.end())
            {
                // add current vector's id to the searched set
                visitedIDs.insert(currentID);
                // calculate the distance between current vector and the query one
                currentDistance = distanceFunction(currentVector, queryVector); 

                // heap is full, have we found a better candidate?
                if((int) imagePool.size() == k)
                {
                    if(imagePool.top().distance > currentDistance)
                    // if so, just pop the worst previous candidate and push the new one
                    {
                        imagePool.pop();
                        imagePool.push(imageDuplet(currentDistance, currentImage));
                    }
                }else //heap not full, freely push
                {
                    imagePool.push(imageDuplet(currentDistance, currentImage));
                }  
            }
        }

    }

    // transform priority queue into vector of images
    while(!imagePool.empty())
    {
        nearestNeighbours->push_back(move(const_cast<image*>(imagePool.top().imageData)));
        imagePool.pop();
    }

    return nearestNeighbours;
}


vector<image*>* LSH::radiusSearch(image *query, unsigned int radius, unsigned int (*distanceFunction) (vector<unsigned char>*, vector<unsigned char>*), set<int>& classifiedIDs)
{
    // currently updated variables
    hashTable *currHashTable;
    bucket *currBucket;
    image *currentImage;
    // the vector representation of the given query
    vector<unsigned char>* queryVector = query->getImage();
    int queryID = query->getId();
    // the vector representation of currently checked image
    vector<unsigned char>* currentVector;
    int currentID;

    // getting the total number of hash tables in the lsh
    int hashesNum = this->getHashTablesNum();
    // number of images in the current bucket
    int bucketImagesNum;
    // the distance between current vector and the query one
    unsigned int currentDistance;
    // bucket number in current hash table
    unsigned int currBucketID;
    // priority queue used as a heap
    priority_queue<imageDuplet, vector<imageDuplet>, compDuplet> imagePool;
    // vector containing nearest images
    vector<image*> *nearestNeighbours = new vector<image*>();
    // visited vectors' ids, no need to check them again
    set<int, greater<int> > visitedIDs;
    
    // for each hashTable hash the given image and 
    // traverse its bucket searching for neighbours
    for(int currHash = 0; currHash < hashesNum; currHash++)
    {
        // getting the info about the bucket into which
        // given image is hashed in current hash table
        currHashTable = this->getIthHashTable(currHash);
        currBucketID = currHashTable->getBucket(query);
        currBucket = currHashTable->getIthBucket(currBucketID);
        bucketImagesNum = currBucket->getSize();

        for(int img = 0; img < bucketImagesNum; img++)
        {
            currentImage = currBucket->getIthImage(img);
            currentVector = currentImage->getImage();
            currentID = currentImage->getId();

            // checking for neighbours update only when new one found and it is different to the query one
            if(currentID != queryID && visitedIDs.find(currentID) == visitedIDs.end() && classifiedIDs.find(currentID) == classifiedIDs.end())
            {
                // add current vector's id to the searched set
                visitedIDs.insert(currentID);
                // calculate the distance between current vector and the query one
                currentDistance = distanceFunction(currentVector, queryVector);
                // current image within radius
                if(currentDistance <= radius)
                {
                    imagePool.push(imageDuplet(currentDistance, currentImage));
                }
            }
        } 
    }

    // transform priority queue into vector of images
    while(!imagePool.empty())
    {
        nearestNeighbours->push_back(move(const_cast<image*>(imagePool.top().imageData)));
        imagePool.pop();
    }

    return nearestNeighbours;
}
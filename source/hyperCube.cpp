#include "../headers/hyperCube.h"
#include "../headers/imageDuplet.h"

hyperCube::hyperCube(unsigned int maxChecks, unsigned int maxEdges)
{
    int projFuncsNum, bigMAspect;
    // initializing the hyper cube's hash table
    this->cubeHashTable = new hashTable();
    projFuncsNum = this->cubeHashTable->getDistrVecNum();
    bigMAspect = this->cubeHashTable->getBigAspectM();

    //creating random, uniformly distributed projection functions of 0s and 1s
    this->projectionFunctions = new vector<vector<unsigned char>*> (projFuncsNum, NULL);

    // f output will follow uniformal distribution
    random_device rd;
    mt19937 gen(rd());
    // each projection function is a vector of bits corresponding to each
    // possible output of the h function with h(max) = M
    for(int i = 0; i < projFuncsNum; i++)
    {
        vector<unsigned char>* newProjectionFunction = new vector<unsigned char>(bigMAspect);
        this->projectionFunctions->at(i) = newProjectionFunction;

        fill_n(newProjectionFunction->begin(), newProjectionFunction->size()/2, 1);
        shuffle(newProjectionFunction->begin(), newProjectionFunction->end(), gen);
    }
    // setting the hyperparameters
    this->edgeCeil = maxEdges;
    this->checkCeil = maxChecks;
}

hyperCube::~hyperCube()
{
    int projFuncsNum = this->cubeHashTable->getDistrVecNum();

    delete(this->cubeHashTable);
    for(int i = 0; i < projFuncsNum; i++){
        delete(this->projectionFunctions->at(i));
    }

    delete(this->projectionFunctions);
}

hashTable* hyperCube::getHashTable()
{
    return this->cubeHashTable;
}

vector<vector <unsigned char> *>* hyperCube::getProjectionFuncs()
{
    return this->projectionFunctions;
}

vector<unsigned char>* hyperCube::getIthProjectionFunc(int index)
{
    return this->getProjectionFuncs()->at(index);
}



//Sigoura tha thelei allages den eimai sigouros gia kapoia pramata
// hyperCube::hyperCube(unsigned int k_, unsigned int d_, float r){
//     this->k = k_;
//     hyperCube::d = d_;
//     this->hMod = pow(2, 32/this->k);
//     this->r = r;
//     hyperCube::wAspect = r * 10;

//     hyperCube::bigM = this->hMod;
//     hyperCube::smallM = pow(2, 32) - 5;

//     random_device rd;
//     mt19937 gen(rd());

//     //depends on the dimensionality of the hypercube
//     unsigned int numOfBuckets = pow(2, k);
//     this->bucketArray = new vector<bucket*> (numOfBuckets, NULL);
    
//     //initializing the buckets
//     for(unsigned int i = 0; i < numOfBuckets; i++)
//     {
//         this->bucketArray->at(i) = new bucket();
//     }

//     //creating random, uniformly distributed projection functions
//     this->projectionFunctions = new vector<vector<unsigned char>*> (k, NULL);
//     //initializing them
//     for(unsigned int i = 0; i < k; i++)
//     {
//         vector<unsigned char>* newProjectionFunction = new vector<unsigned char>(hMod);
//         this->projectionFunctions->at(i) = newProjectionFunction;

//         fill_n(newProjectionFunction->begin(), newProjectionFunction->size()/2, 1);
//         shuffle(newProjectionFunction->begin(), newProjectionFunction->end(), gen);
//     }

//     this->disturbanceVectors = initDistrVectors(k,d,wAspect);

//     this->initMExponentials();  
// }

void hyperCube::initializeHyperCube(vector<image*>* imageArray)
{
    long unsigned int totalImages = imageArray->size();
    for(long unsigned int i = 0; i < totalImages; i++)
    {
        this->insert(imageArray->at(i));
    }
}

void hyperCube::insert(image *imageToInsert)
{
    unsigned int hashResult;
    //get the vertex where we will insert
    hashResult = this->hashFunction(imageToInsert);
    this->getHashTable()->getIthBucket(hashResult)->pushImage(imageToInsert);
    return;
}


//MAY NEED TO CHANGE IT A BIT JUST FOR THE LOOKS
unsigned int hyperCube::hashFunction(image *givenImage){
    // getting the outputs of each h function
    unsigned int currentHOutput;
    unsigned char currentBit;
    unsigned int hashValue = 0;

    // hyper cube's hash table
    hashTable *cubeHash = this->getHashTable();
    // the total number of h functions that participate in g construction
    unsigned int distrVecNum = cubeHash->getDistrVecNum();

    //for every disturbance vector
    for(unsigned int i = 0; i < distrVecNum; i ++){
        // calculating the ith h function for input image
        currentHOutput = cubeHash->hFunction(i, givenImage);
        // projecting the h function value to bit space
        currentBit = this->getIthProjectionFunc(i)->at(currentHOutput);
        // simply adding current output bit to the final result
        hashValue |= currentBit;
        // making space for the next bit
        hashValue <<= 1;
    }

    // no new bit in the last iteration
    // simply shift right
    hashValue >>= 1;

    return hashValue;
}

unsigned int hyperCube::getEdgeCeil()
{
    return this->edgeCeil;
}

unsigned int hyperCube::getCheckCeil()
{
    return this->checkCeil;
}

// initializing an array of boolean representing whether ith node of hypercube has been visited
// storing node number representations in a queue, applying bfs
vector<image*>* hyperCube::getNN(image *query, int k, unsigned int (*distanceFunction) (vector<unsigned char>*, vector<unsigned char>*))
{
    auto comp = [](imageDuplet a, imageDuplet b)
    {
        return a.distance < b.distance;
    };
    // priority queue used as a heap
    priority_queue<imageDuplet, vector<imageDuplet>, decltype(comp)> imagePool(comp);
    // vector containing nearest images
    vector<image*> *nearestNeighbours = new vector<image*>();
    // total number of images checked and limit
    unsigned int checkedImages = 0;
    unsigned int imagesLimit = this->getCheckCeil();
    // total number of nodes checked and limit
    unsigned int checkedNodes = 0;
    unsigned int nodesLimit = this->getEdgeCeil();
    // total number of hypercube nodes
    int dimensions = this->getHashTable()->getDistrVecNum();
    unsigned int totalNodes = pow(2, dimensions);
    bool *visitedSet = (bool *) malloc(totalNodes * sizeof(bool));
    // queue used to store arithmetic representation of nodes in bfs
    queue<unsigned int> *toVisitSet = new queue<unsigned int>();
    // vector containing binary neighbours of current node
    vector<unsigned int> *neighbours = new vector<unsigned int>(dimensions, 0);

    // all hypercube nodes are initially unvisited
    for(unsigned int i = 0; i < totalNodes; i++)
    {
        visitedSet[i] = false;
    }
    // initial node number will be equal to query image's hash value
    unsigned int node = this->hashFunction(query);
    unsigned int neighbour;
    // query's image has to be visited, no need to check it again
    toVisitSet->push(node);
    visitedSet[node] = true;


    // variables for current bucket
    bucket *currBucket;
    vector<image*> *currImages;
    int bucketImagesNum;
    // data of the query and current image
    int currentID;
    int queryID = query->getId();
    vector<unsigned char> *currentVector;
    vector<unsigned char> *queryVector = query->getImage();
    // data of currently checked image
    image *currentImage;
    int bucketImage;
    unsigned int currentDistance;

    while(!toVisitSet->empty() && checkedNodes < nodesLimit)
    {
        // checking current node/bucket
        node = toVisitSet->front();
        toVisitSet->pop();

        currBucket = this->getHashTable()->getIthBucket(node);
        currImages = currBucket->getImages();
        bucketImagesNum = currBucket->getSize();
        // initializing the number of images checked in current bucket
        bucketImage = 0;
        // checking each new candidate till total number of images cap is reached
        // or we have checked the whole bucket
        while(bucketImage < bucketImagesNum && checkedImages < imagesLimit)
        {   //for each image in current bucket
            currentImage = currImages->at(bucketImage);
            currentVector = currentImage->getImage();
            currentID = currentImage->getId();

            // checking for neighbours update only when new one found and it is different to the query one
            if(currentID != queryID)
            {
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
                // checked another image within the bucket
            }
            // checked another image in the bucket
            checkedImages++;
            bucketImage++;
        }
        // did we reach the limit of images cap?
        if(checkedImages == imagesLimit)
        {
            break;
        }

        // getting current node's neighbours
        getBinaryNeighbours(node, neighbours);
        // checking for each neighbouring node if it has been visited
        // if not, add it to the "toVisitSet"
        int neighboursSize = neighbours->size();
        for(int i = 0; i < neighboursSize; i++)
        {
            neighbour = neighbours->at(i);
            // new unvisited node
            if(!visitedSet[neighbour])
            {
                visitedSet[neighbour] = true;
                toVisitSet->push(neighbour);
            }
        }
        // another node checked
        checkedNodes++;
    }

    // transform priority queue into vector of images
    while(!imagePool.empty())
    {
        nearestNeighbours->push_back(move(const_cast<image*>(imagePool.top().imageData)));
        imagePool.pop();
    }

    delete(toVisitSet);
    delete(neighbours);
    free(visitedSet);

    return nearestNeighbours;
}


vector<image*>* hyperCube::radiusSearch(image *query, unsigned int radius, unsigned int (*distanceFunction) (vector<unsigned char>*, vector<unsigned char>*), set<int>& classifiedIDs)
{
    auto comp = [](imageDuplet a, imageDuplet b)
    {
        return a.distance < b.distance;
    };
    // priority queue used as a heap
    priority_queue<imageDuplet, vector<imageDuplet>, decltype(comp)> imagePool(comp);
    // vector containing nearest images
    vector<image*> *nearestNeighbours = new vector<image*>();
    // total number of images checked and limit
    unsigned int checkedImages = 0;
    unsigned int imagesLimit = this->getCheckCeil();
    // total number of nodes checked and limit
    unsigned int checkedNodes = 0;
    unsigned int nodesLimit = this->getEdgeCeil();
    // total number of hypercube nodes
    int dimensions = this->getHashTable()->getDistrVecNum();
    unsigned int totalNodes = pow(2, dimensions);
    bool *visitedSet = (bool *) malloc(totalNodes * sizeof(bool));
    // queue used to store arithmetic representation of nodes in bfs
    queue<unsigned int> *toVisitSet = new queue<unsigned int>();
    // vector containing binary neighbours of current node
    vector<unsigned int> *neighbours = new vector<unsigned int>(dimensions, 0);

    // all hypercube nodes are initially unvisited
    for(unsigned int i = 0; i < totalNodes; i++)
    {
        visitedSet[i] = false;
    }
    // initial node number will be equal to query image's hash value
    unsigned int node = this->hashFunction(query);
    unsigned int neighbour;
    // query's image has to be visited, no need to check it again
    toVisitSet->push(node);
    visitedSet[node] = true;


    // variables for current bucket
    bucket *currBucket;
    vector<image*> *currImages;
    int bucketImagesNum;
    // data of the query and current image
    int currentID;
    int queryID = query->getId();
    vector<unsigned char> *currentVector;
    vector<unsigned char> *queryVector = query->getImage();
    // data of currently checked image
    image *currentImage;
    int bucketImage;
    unsigned int currentDistance;

    while(!toVisitSet->empty() && checkedNodes < nodesLimit)
    {
        // checking current node/bucket
        node = toVisitSet->front();
        toVisitSet->pop();

        currBucket = this->getHashTable()->getIthBucket(node);
        currImages = currBucket->getImages();
        bucketImagesNum = currBucket->getSize();
        // initializing the number of images checked in current bucket
        bucketImage = 0;
        // checking each new candidate till total number of images cap is reached
        // or we have checked the whole bucket
        while(bucketImage < bucketImagesNum && checkedImages < imagesLimit)
        {   //for each image in current bucket
            currentImage = currImages->at(bucketImage);
            currentVector = currentImage->getImage();
            currentID = currentImage->getId();

            // checking for neighbours update only when new one found and it is different to the query one
            if(currentID != queryID && classifiedIDs.find(currentID) == classifiedIDs.end())
            {
                // calculate the distance between current vector and the query one
                currentDistance = distanceFunction(currentVector, queryVector); 

                if(currentDistance <= radius)
                {
                    imagePool.push(imageDuplet(currentDistance, currentImage));
                }
            }
            // checked another image in the bucket
            checkedImages++;
            bucketImage++;
        }
        // did we reach the limit of images cap?
        if(checkedImages == imagesLimit)
        {
            break;
        }

        // getting current node's neighbours
        getBinaryNeighbours(node, neighbours);
        // checking for each neighbouring node if it has been visited
        // if not, add it to the "toVisitSet"
        int neighboursSize = neighbours->size();
        for(int i = 0; i < neighboursSize; i++)
        {
            neighbour = neighbours->at(i);
            // new unvisited node
            if(!visitedSet[neighbour])
            {
                visitedSet[neighbour] = true;
                toVisitSet->push(neighbour);
            }
        }
        // another node checked
        checkedNodes++;
    }

    // transform priority queue into vector of images
    while(!imagePool.empty())
    {
        nearestNeighbours->push_back(move(const_cast<image*>(imagePool.top().imageData)));
        imagePool.pop();
    }

    free(visitedSet);
    delete(toVisitSet);
    delete(neighbours);

    return nearestNeighbours;
}









//         // got to the limit of images given in current bucket
//         if(!getBucketNN(query, node, k, this->getHashTable(), &checkedImages, imagesLimit, visitedIDs, imagePool, distanceFunction))
//         {
//             break;
//         }
//         // getting current node's neighbours
//         getBinaryNeighbours(node, neighbours);
//         // checking for each neighbouring node if it has been visited
//         // if not, add it to the "toVisitSet"
//         for(int i = 0; i < neighbours->size(); i++)
//         {
//             neighbour = neighbours->at(i);
//             // new unvisited node
//             if(!visitedSet[neighbour])
//             {
//                 visitedSet[neighbour] = true;
//                 toVisitSet->push(neighbour);
//             }
//         }

//         // another node checked
//         checkedNodes++;
//     }

//     // transform priority queue into vector of images
//     while(!imagePool.empty())
//     {
//         nearestNeighbours->push_back(move(const_cast<image*>(imagePool.top().imageData)));
//         imagePool.pop();
//     }

//     return nearestNeighbours;
// }


//     bucket *currBucket = currHashTable->getIthBucket(index);
//     vector<image*> *currImages = currBucket->getImages();
//     int bucketImagesNum = currBucket->getSize();
//     // counting the number of images we checked in current bucket
//     int currImage = 0;
//     // counting the total number of images checked

//     // data of the query image
//     int queryID = query->getId();
//     vector<unsigned char> *queryVector = query->getImage();
//     // data of currently checked image
//     image *currentImage;
//     vector<unsigned char> *currentVector;
//     int currentID;

//     // checking each new candidate till total number of images cap is reached
//     // or we have checked the whole bucket
//     unsigned int currentDistance;
//     while(currImage < bucketImagesNum && *checkedImages < imagesCap)
//     {   //for each image in current bucket
//         currentImage = currImages->at(currImage);
//         currentVector = currentImage->getImage();
//         currentID = currentImage->getId();

//         // checking for neighbours update only when new one found and it is different to the query one
//         if(currentID != queryID && visitedIDs.find(currentID) == visitedIDs.end())
//         {
//             // add current vector's id to the searched set
//             visitedIDs.insert(currentID);
//             // calculate the distance between current vector and the query one
//             currentDistance = distanceFunction(currentVector, queryVector); 

//             // heap is full, have we found a better candidate?
//             if((int) imagePool.size() == k)
//             {
//                 if(imagePool.top().distance > currentDistance)
//                 // if so, just pop the worst previous candidate and push the new one
//                 {
//                     imagePool.pop();
//                     imagePool.push(imageDuplet(currentDistance, currentImage));
//                 }
//             }else //heap not full, freely push
//             {
//                 imagePool.push(imageDuplet(currentDistance, currentImage));
//             }
            
//         }
//         // checked another image within the bucket
//         *checkedImages = *checkedImages + 1;
//         currImage++;
//     }


#include "../headers/clustering.h"
#include "../headers/LSH.h"
#include "../headers/hyperCube.h"

clustering::clustering(int k, int hashTableNum, unsigned int maxChecks, unsigned int maxEdges, string* method, vector<image*> *imageArray)
{
    int totalImages;
    this->k = k;
    this->method = method;
    this->clusters = new vector<cluster*>(k, NULL);
    for(int i = 0; i < this->k; i++){
        this->clusters->at(i) = new cluster();
    }
    this->centroids = new vector<image*>(k, NULL);
    this->imageArray = imageArray;
    totalImages = this->getImages()->size();
    this->entries = new vector<clusterEntry*>(totalImages);
    //initialize the entries that will be pushed into clusters
    for(int i = 0; i < totalImages; i++)
    {
        this->entries->at(i) = new clusterEntry();
    }
    // initializing the lsh or hypercube datastructure
    // depending on the given parameters
    if(*this->getMethod() == "LSH")
    {
        this->clusterLSH = new LSH(hashTableNum);
        this->getLSH()->initializeLSH(this->imageArray);
    }else if(*this->getMethod() == "Hypercube")
    {
        this->clusterCube = new hyperCube(maxChecks, maxEdges);
        this->getHyperCube()->initializeHyperCube(this->imageArray);
    }
    cout << "Clustering class initialized" << endl;
}

clustering::~clustering(){
    if(*this->method == "LSH"){
        delete(this->clusterLSH);
    }else if(*this->method == "Hypercube"){
        delete(this->clusterCube);
    }

    int clusterNum = this->clusters->size();

    for(int i = 0; i < clusterNum; i++){
        delete(this->clusters->at(i));
    }
    delete(this->clusters);

    for(int i = 0; i < clusterNum; i++){
        if(this->centroids->at(i)->getId() == -1){
            delete(this->centroids->at(i));
        }
    }

    delete(this->centroids);

    int totalImages = this->getImages()->size();
    for(int i = 0; i < totalImages; i++){
        delete(this->entries->at(i));
    }
    delete(this->entries);
}

// returns the minimum and maximum distance between current centroids
void clustering::getCentroidDistances(unsigned int *minDistance, unsigned int *maxDistance, unsigned int (*distanceFunction) (vector<unsigned char>*, vector<unsigned char>*))
{
    int totalCentroids = this->getClustersNum();
    unsigned int currDistance;
    unsigned int currMaxDistance = 0;
    unsigned int currMinDistance = UINT_MAX;

    vector<unsigned char>* firstCentroid;
    vector<unsigned char>* secondCentroid;

    for(int i = 0; i < totalCentroids; i++)
    {
        for(int j = i + 1; j < totalCentroids; j++)
        {
            firstCentroid = this->getIthCentroid(i)->getImage();
            secondCentroid = this->getIthCentroid(j)->getImage();
            currDistance = distanceFunction(firstCentroid, secondCentroid);

            if(currDistance > currMaxDistance)
            {
                currMaxDistance = currDistance;
            }

            if(currDistance < currMinDistance)
            {
                currMinDistance = currDistance;
            }
        }
    }
    
    *minDistance = currMinDistance;
    *maxDistance = currMaxDistance;
    return;
}


template <typename DataType>
void clustering::clusterizeData(unsigned int (*distanceFunction) (vector<unsigned char>*, vector<unsigned char>*))
{
    LSH *LSH_hasher;
    hyperCube *cube_hasher;
    // getting initial centroids used to construct the first clusters
    this->initializeCentroids(distanceFunction);

    if(*this->getMethod() == "LSH")
    {
        LSH_hasher = this->getLSH();
    }else if(*this->getMethod() == "Hypercube")
    {
        cube_hasher = this->getHyperCube();
    }

    // for each iteration clusterize according to requested technique
    for(int iteration = 0; iteration < MAX_CLUSTERIZATIONS; iteration++)
    {
        if(*this->getMethod() == "Classic")
        {
            this->Lloyd(distanceFunction);
        }else if(*this->getMethod() == "LSH")
        {
            this->reverseAssignment(LSH_hasher, distanceFunction);
        }else
        {
            this->reverseAssignment(cube_hasher, distanceFunction);
        }
        // clusters should not be dropped in the final iteration
        if(iteration != MAX_CLUSTERIZATIONS - 1)
        {
            this->newCentroids();
            this->resetClusters();
        }
        cout << iteration + 1 << " phase" << endl;
    }
    cout << endl;
    return;
}





template <typename DataType>
void clustering::reverseAssignment(DataType hasher, unsigned int (*distanceFunction) (vector<unsigned char>*, vector<unsigned char>*))
{
    int totalImages = this->getImages()->size();
    // different distances
    unsigned int maxCenDis, minCenDis;
    unsigned int currDistance;
    int totalCentroids = this->getClustersNum();
    // getting the maximum and minimum distance between current centroids
    this->getCentroidDistances(&minCenDis, &maxCenDis, distanceFunction);
    unsigned int initialRadius = minCenDis / 2;
    unsigned int currRadius = initialRadius;
    // data concerning the images that will be clusterized in current execution
    vector<image*> *currCandidates;
    image *currCandidate;
    int currCandidateID;
    // current image and current centroid
    image *centroid;
    clusterEntry *currEntry;
    // containing the ids of images that have been assigned to a specific cluster
    set<int> clusteredIDs;
    while(currRadius < maxCenDis)
    {
        // for each centroid, do radius search and decide clusterization for found images
        for(int currCentroid = 0; currCentroid < totalCentroids; currCentroid++)
        {
            centroid = this->getIthCentroid(currCentroid);
            // getting candidates for current centroid within current radius
            currCandidates = hasher->radiusSearch(centroid, currRadius, distanceFunction, clusteredIDs);
            while(!currCandidates->empty())
            {
                currCandidate = currCandidates->front();
                currCandidates->erase(currCandidates->begin());
                currCandidateID = currCandidate->getId();
                // found an image that has not been assigned
                if(clusteredIDs.find(currCandidateID) == clusteredIDs.end())
                {
                    // has current image been found within the radius
                    // of another centroid in current iteration?
                    currEntry = this->getIthEntry(currCandidateID);
                    // calculate the distance of current image with current centroid
                    currDistance = distanceFunction(currCandidate->getImage(), centroid->getImage());
                    if(currEntry->isConsidered())
                    {
                        // checking if we found a centroid that is more near
                        // the current image
                        if(currDistance < currEntry->getDistanceToCluster())
                        {
                            currEntry->setClusterID(currCentroid);
                            currEntry->setDistance(currDistance);
                        }
                    }else
                    {
                        currEntry->setClusterID(currCentroid);
                        // current image is considered in current iteration
                        currEntry->setDecidability(true);
                        currEntry->setDistance(currDistance);
                    }
                    
                }
            } 
            delete(currCandidates);  
        }
        // having considered all the centroids
        // the debatable images can now be assigned to actual
        // nearest centroid

        for(int entry = 0; entry < totalImages; entry++)
        {
            currEntry = this->getIthEntry(entry);
            // is image considered in current iteration and was not previous appointed?
            if(clusteredIDs.find(entry) == clusteredIDs.end() && currEntry->isConsidered())
            {
                clusteredIDs.insert(entry);
            }
        }
        currRadius = currRadius * 2;
    }

    int clusterIndex;
    unsigned int centroidDistance;
    // assigning found images into their clusters
    for(int entry = 0; entry < totalImages; entry++)
    {
        currEntry = this->getIthEntry(entry);
        currCandidate = this->getIthImage(entry);
        // cluster was found for specific image
        if(clusteredIDs.find(entry) != clusteredIDs.end())
        {
            clusterIndex = currEntry->getClusterID();
            this->getIthCluster(clusterIndex)->insert(currCandidate);
        }else
        {//no cluster for current image, adding it to the unidentified cluster images
            nearestImagesInfo(currCandidate, this->getCentroids(), &centroidDistance, &clusterIndex, distanceFunction);
            this->getIthCluster(clusterIndex)->insert(currCandidate);
        }   
    }
    // all entries initially not considered in the next range search call
    this->resetEntries();
    return;
}

void clustering::resetClustering()
{
    int totalClusters = this->getClustersNum();
    image *currentCentroid;
    this->resetClusters();
    // delete centroids that are not real images from previous method
    for(int i = 0; i < totalClusters; i++)
    {
        currentCentroid = this->getIthCentroid(i);

        if(currentCentroid->getId() == -1)
        {
            delete(currentCentroid);
        }
    }
}

// executing one of the 3 possible clusterization methods
// based on the input index
double clustering::executeMethod(int index, int hashTableNum, unsigned int maxChecksHyperCube, unsigned int probesHyperCube, unsigned int (*distanceFunction) (vector<unsigned char>*, vector<unsigned char>*))
{
    LSH *LSH_hasher;
    hyperCube *cube_hasher;

    // clearing data from previous method
    if(index != 0)
    {
        this->resetClustering();
    }

    if(index == 1)
    {
        LSH_hasher = new LSH(hashTableNum);
        LSH_hasher->initializeLSH(this->getImages());
    }else if(index == 2)
    {
        cube_hasher = new hyperCube(maxChecksHyperCube, probesHyperCube);
        cube_hasher->initializeHyperCube(this->getImages());
    }

    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    // getting initial centroids used to construct the first clusters
    this->initializeCentroids(distanceFunction);

    // for each iteration clusterize according to requested technique
    for(int iteration = 0; iteration < MAX_CLUSTERIZATIONS; iteration++)
    {
        if(index == 0)//Lloyd
        {
            this->Lloyd(distanceFunction);
        }else if(index == 1)//LSH
        {
            this->reverseAssignment(LSH_hasher, distanceFunction);
        }else//hyperCube
        {
            this->reverseAssignment(cube_hasher, distanceFunction);
        }
        // clusters should not be dropped in the final iteration
        if(iteration != MAX_CLUSTERIZATIONS - 1)
        {
            this->newCentroids();
            this->resetClusters();
        }
        cout << iteration + 1 << " phase" << endl;
    }
    cout << endl;

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    // calculate the time passed
    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);

    if(index == 1)
    {
        delete(LSH_hasher);
    }else if(index == 2)
    {
        delete(cube_hasher);
    }

    return time_span.count();
}


void clustering::resetEntries()
{
    int totalImages = this->getImages()->size();

    for(int i = 0; i < totalImages; i++)
    {
        this->getIthEntry(i)->setDecidability(false);
    }
    return;
}

void clustering::Lloyd(unsigned int (*distanceFunction) (vector<unsigned char>*, vector<unsigned char>*)){
    unsigned int imageArraySize = imageArray->size();
    image* currImage;
    unsigned int distance;
    // the final cluster into which current image will be put into
    int wantedIndex;
    cluster* wantedCluster;

    for(unsigned int i = 0; i < imageArraySize; i++){
        currImage = imageArray->at(i);
        nearestImagesInfo(currImage, this->getCentroids(), &distance, &wantedIndex, distanceFunction);
        wantedCluster = this->getIthCluster(wantedIndex);
        wantedCluster->insert(currImage);
    }
}


void clustering::addCentroid(image *newCentroid)
{
    this->getCentroids()->push_back(newCentroid);
}

void clustering::updateIthCentroid(int index, image *newCentroid)
{
    this->getCentroids()->at(index) = newCentroid;
}

void clustering::nearestImagesInfo(image *targetImage, vector<image*> *neighbours, unsigned int *distance, int *neighbourIndex, unsigned int (*distanceFunction) (vector<unsigned char>*, vector<unsigned char>*))
{
    int totalNeighbours = neighbours->size();
    int minIndex = 0;
    // setting the minimum distance as the distance with the first neighbour
    unsigned int minDistance = distanceFunction(targetImage->getImage(), neighbours->at(0)->getImage());
    unsigned int currDistance;


    for(int i = 1; i < totalNeighbours; i++)
    {
        currDistance = distanceFunction(targetImage->getImage(), neighbours->at(i)->getImage());

        if(currDistance < minDistance)
        {
            minDistance = currDistance;
            minIndex = i;
        }
    }

    *neighbourIndex = minIndex;
    *distance = minDistance;
    return;
}

void clustering::updateDistanceArray(int totalImages, unsigned int *maxDistance, image *currentCentroid, bool* isCentroid, unsigned int* distanceArray, unsigned int (*distanceFunction) (vector<unsigned char>*, vector<unsigned char>*))
{
    image *currentImage;
    unsigned int currMaxDistance = 0;
    unsigned int currentDistance;

    for(int i = 0; i < totalImages; i++)
    {
        currentImage = this->getIthImage(i);

        if(!isCentroid[i])
        {
            currentDistance = distanceFunction(currentImage->getImage(), currentCentroid->getImage());

            if(currentDistance < distanceArray[i])
            {
                distanceArray[i] = currentDistance;
            }

            if(currentDistance > currMaxDistance)
            {
                currMaxDistance = currentDistance;
            }
        }else
        {
            distanceArray[i] = 0;
        }
    }

    *maxDistance = currMaxDistance;
    return;
}

void clustering::updateProbArray(int totalImages, unsigned int maxDistance, bool* isCentroid,unsigned int* distanceArray, double* probArray)
{
    int probSize = totalImages + 1;
    double sumIncrease;

    for(int i = 1; i < probSize; i++)
    {
        // if current element is a centroid, distance between previous sum is zero
        if(isCentroid[i - 1])
        {
            sumIncrease = 0.0;
        }else
        {
            sumIncrease = double( distanceArray[i - 1]/ double(maxDistance));
            sumIncrease = sumIncrease * sumIncrease;
        }
        probArray[i] = probArray[i - 1] + sumIncrease;
    }
    return;
}

void clustering::initializeCentroids(unsigned int (*distanceFunction) (vector<unsigned char>*, vector<unsigned char>*))
{
    int totalCentroids = this->getClustersNum();
    // total number of images
    int imagesNum = this->getImages()->size();
    // used to check if we want to check current point (not a centroid)
    bool isCentroid[imagesNum];
    fill_n(isCentroid, imagesNum, false);
    // array containing the increasing sum, used for probability distribution
    double probArray[imagesNum + 1];
    fill_n(probArray, imagesNum + 1, 0.0);
    double randomNum;
    // contains the distances of each image from nearest centroid
    unsigned int distanceArray[imagesNum];
    fill_n(distanceArray, imagesNum, UINT_MAX);
    unsigned int maxDistance;
    // new centroid and its index
    image *currCentroid;
    int newCentroidIndex;

    currCentroid = this->getIthImage(rand() % imagesNum);
    this->updateIthCentroid(0, currCentroid);
    isCentroid[currCentroid->getId()] = true;

    // initialize random distribution of doubles
    default_random_engine re;

    for(int i = 1; i < totalCentroids; i++)
    {
        // distance and probability function array update
        updateDistanceArray(imagesNum, &maxDistance, currCentroid, isCentroid, distanceArray, distanceFunction);
        updateProbArray(imagesNum, maxDistance, isCentroid, distanceArray, probArray);
        uniform_real_distribution<double> unif(probArray[0], probArray[imagesNum]);
        // get a random float between the min and max probability number
        // randomNum = randomDouble(probArray[0], probArray[imagesNum]);
        randomNum = unif(re);
        // get the index of the new centroid
        newCentroidIndex = binarySearch(probArray, 0, imagesNum, randomNum);
        // adding the new centroid into the centroid pool
        isCentroid[newCentroidIndex] = true;
        currCentroid = getIthImage(newCentroidIndex);
        this->updateIthCentroid(i, currCentroid);
        isCentroid[currCentroid->getId()] = true;
    }
    return;
}

// this function simply updates the centroids
// new centroids are the median vectors of previous clusters
void clustering::newCentroids()
{
    int totalCentroids = this->getClustersNum();
    int imageSize = this->getIthImage(0)->getImage()->size();
    image* newCentroid;

    for(int i = 0; i < totalCentroids; i++)
    {
        newCentroid = this->getIthCluster(i)->getNewCentroid(imageSize);
        if(this->centroids->at(i)->getId() == -1){
            delete(this->centroids->at(i));
        }
        this->centroids->at(i) = newCentroid;
    }
}


vector<cluster*>* clustering::getClusters()
{
    return this->clusters;
}

cluster* clustering::getIthCluster(int index)
{
    return this->clusters->at(index);
}

vector<image*>* clustering::getCentroids()
{
    return this->centroids;
}

image* clustering::getIthCentroid(int index)
{
    return this->centroids->at(index);
}

vector<image*>* clustering::getImages()
{
    return this->imageArray;
}

image* clustering::getIthImage(int index)
{
    return this->imageArray->at(index);
}

int clustering::getClustersNum()
{
    return this->k;
}

string* clustering::getMethod()
{
    return this->method;
}

vector<clusterEntry*>* clustering::getEntries()
{
    return this->entries;
}

clusterEntry* clustering::getIthEntry(int index)
{
    return this->entries->at(index);
}

hyperCube* clustering::getHyperCube()
{
    return this->clusterCube;
}

LSH* clustering::getLSH()
{
    return this->clusterLSH;
}

void clustering::resetClusters(){
    int sizeOfClusters = clusters->size();

    for(int i = 0; i < sizeOfClusters; i++){
        this->getIthCluster(i)->resetCluster();
    }
}

// executes consecutive clusterizations based on given technique
double clustering::execute(unsigned int (*distanceFunction) (vector<unsigned char>*, vector<unsigned char>*))
{
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    if(*this->getMethod() == "LSH")
    {
        this->clusterizeData<LSH*>(distanceFunction);
    }else
    {
        this->clusterizeData<hyperCube*>(distanceFunction);
    }
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    // calculate the time passed
    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
    return time_span.count();
}

// returns the id of the second nearest cluster to given image
int clustering::getNearClusterID(int currentCluster, image *currentImage, unsigned int (*distanceFunction) (vector<unsigned char>*, vector<unsigned char>*))
{
    int totalClusters = this->getClustersNum();
    int neighbourClusterID;
    unsigned int smallestDistance = UINT_MAX;
    unsigned int currentDistance;
    image* compareCentroid;

    // for each centroid not belonging to image's cluster
    // check for the nearest one
    for(int i = 0; i < totalClusters; i++)
    {
        if(currentCluster != i)
        {
            compareCentroid = this->getIthCentroid(i);
            currentDistance = distanceFunction(compareCentroid->getImage(), currentImage->getImage());
            if(currentDistance < smallestDistance)
            {
                smallestDistance = currentDistance;
                neighbourClusterID = i;
            }
        }
    }

    return neighbourClusterID;
}

// returns the size of the biggest cluster
int clustering::getMaxClusterSize()
{
    int currSize;
    int maxSize = 0;
    int totalClusters = this->getClustersNum();

    for(int i = 0; i < totalClusters; i++)
    {
        currSize = this->getIthCluster(i)->getClusterSize();

        if(currSize > maxSize)
        {
            maxSize = currSize;
        }
    }
    return maxSize;
}

double clustering::calculateSilhouette(unsigned int (*distanceFunction) (vector<unsigned char>*, vector<unsigned char>*), vector<double>** silhouettes)
{
    *silhouettes = this->calculateSilhouettes(distanceFunction);
    int totalClusters = this->getClustersNum();
    // contains the final silhouette of all clusters
    double silhouette = 0.0;

    for(int i = 0; i < totalClusters; i++)
    {
        silhouette += (*silhouettes)->at(i);
    }

    silhouette /= (double)(totalClusters);

    return silhouette;
}


// returns a vector containing the silhouette of each cluster
vector<double>* clustering::calculateSilhouettes(unsigned int (*distanceFunction) (vector<unsigned char>*, vector<unsigned char>*))
{
    int totalCentroids = this->getClustersNum();
    int maxClusterSize = this->getMaxClusterSize();
    vector<double>* silhouettes = new vector<double>(totalCentroids);
    double i_silhouette;

    double** distanceMatrix = (double**)malloc(sizeof(double*)*maxClusterSize);
    for(int i = 0; i < maxClusterSize; i++){
        distanceMatrix[i] = (double*)malloc(sizeof(double)*maxClusterSize);
    }
    


    for(int i = 0; i < totalCentroids; i++)
    {
        i_silhouette = calculateIthSilhouette(i, distanceFunction, distanceMatrix);
        silhouettes->at(i) = i_silhouette;
    }

    for(int i = 0; i < maxClusterSize; i++){
        free(distanceMatrix[i]);
    }
    free(distanceMatrix);

    return silhouettes;
}

double clustering::calculateIthSilhouette(int index, unsigned int (*distanceFunction) (vector<unsigned char>*, vector<unsigned char>*), double **distanceMatrix)
{
    set<image*>* currentCluster = this->getIthCluster(index)->getClust();
    set<image*>* otherCluster;
    int neighbourClusterID;
    int currClusterSize, otherClusterSize;
    // to avoid distance calculation with the same items
    int outerIter = 0;
    int innerIter = 0;
    double a_i = 0.0;
    double b_i = 0.0;
    double partial_add;
    // silhouette of current cluster
    double silhouette = 0.0;

    // number of elements in current cluster
    currClusterSize = this->getIthCluster(index)->getClusterSize() - 1;

    for(image* elem : *currentCluster)
    {
        neighbourClusterID = this->getNearClusterID(index, elem, distanceFunction);
        otherCluster = this->getIthCluster(neighbourClusterID)->getClust();
        // getting the size of the second nearest cluster to current image
        otherClusterSize = this->getIthCluster(neighbourClusterID)->getClusterSize();

        //calculating 
        innerIter = 0;
        a_i = 0.0;
        for(image* elem1 : *currentCluster)
        {
            if(innerIter != outerIter)
            {
                if(innerIter >= outerIter)
                {
                    partial_add = (double) ((double) (distanceFunction(elem->getImage(), elem1->getImage())) / (double) (currClusterSize));
                    distanceMatrix[outerIter][innerIter] = partial_add;
                    a_i = a_i + partial_add;
                }else
                {
                    // distance already calculated, just add it
                    a_i = a_i + distanceMatrix[innerIter][outerIter];
                }
            }
            innerIter++;
        }

        b_i = 0.0;
        for(image* elem2 : *otherCluster)
        {
            partial_add = (double) ((double) (distanceFunction(elem->getImage(), elem2->getImage())) / (double) (otherClusterSize));
            b_i = b_i + partial_add;
        }

        outerIter++;

        if(a_i == b_i)
        {
            silhouette += 0;
        }else if(a_i < b_i)
        {
            silhouette += 1.0 - a_i / b_i;
        }else
        {
            silhouette += b_i / a_i - 1.0;
        }
    }

    // Current cluster is not empty 
    if(currClusterSize != -1 && currClusterSize != 0)
        silhouette = silhouette / (double)(currClusterSize);
    
    return silhouette;
}

template void clustering::reverseAssignment<LSH*>(LSH*, unsigned int (*) (vector<unsigned char>*, vector<unsigned char>*));
template void clustering::reverseAssignment<hyperCube*>(hyperCube*, unsigned int (*) (vector<unsigned char>*, vector<unsigned char>*));
template void clustering::clusterizeData<LSH*>(unsigned int (*) (vector<unsigned char>*, vector<unsigned char>*));
template void clustering::clusterizeData<hyperCube*>(unsigned int (*) (vector<unsigned char>*, vector<unsigned char>*));

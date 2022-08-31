#include "../headers/cluster.h"
#include <bits/stdc++.h> 

cluster::cluster(){
    this->clust = new set<image*>;
}

cluster::~cluster(){
    delete(this->clust);
}

int cluster::getClusterSize(){
    return this->clust->size();
}

void cluster::insert(image* newImage){
    this->clust->insert(newImage);
}

image* cluster::getNewCentroid(int imageSize){

    int clusterSize = this->getClusterSize();         //getting the size of this cluster
    unsigned char xDimensionOfImages[clusterSize];              //creating an array which holds the value of a dimension for each image
    //int sizeOfArray = sizeof(xDimensionOfImages)/sizeof(xDimensionOfImages[0]);
    vector<unsigned char>* newCoordinates = new vector<unsigned char>(imageSize);
    
    for(int i = 0; i < imageSize; i++){
        int j = 0;
        for(image* currImage : *this->clust){
            xDimensionOfImages[j] = currImage->getImage()->at(i); //getting values from the images of the cluster
            j++;
        }
        sort(xDimensionOfImages, xDimensionOfImages + clusterSize); //sorting the array of one column of coordinates
        int medianIndex = clusterSize/2;                            // getting the median value
        unsigned char medianCoordinate = xDimensionOfImages[medianIndex];   
        
        newCoordinates->at(i) = medianCoordinate; //building the coordinates of the new image/centroid
    }

    //creating a new image
    image* newCentroid = new image(-1, newCoordinates);
    
    //updating the new centroid
    return newCentroid;
    
}

void cluster::resetCluster(){
    delete(this->clust);
    this->clust = new set<image*>;
}

set<image*>* cluster::getClust()
{
    return this->clust;
}
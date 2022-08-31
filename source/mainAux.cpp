#include "../headers/mainAux.h"


void getArgumentsLSH(int argc, char **argv, string *inputFile, string *queryFile, string *outputFile, int *hFuncs, int *hashNum, int *neighbours, double *radius)
{
    *hFuncs = 4;
    *hashNum = 5;
    *neighbours = 1;
    *radius = 10000.0;

    for(int i = 0; i < argc; i++)
    {
        if(!strcmp(argv[i], "-d"))
        {
            *inputFile = string(argv[i+1]);
        }
        else if(!strcmp(argv[i], "-q"))
        {
            *queryFile = string(argv[i+1]);
        }
        else if(!strcmp(argv[i], "-k"))
        {
          *hFuncs = atoi(argv[i+1]);
        }
        else if(!strcmp(argv[i], "-L"))
        {
          *hashNum = atoi(argv[i+1]);
        }
        else if(!strcmp(argv[i], "-o"))
        {
          *outputFile = string(argv[i+1]);
        }
        else if(!strcmp(argv[i], "-N"))
        {
          *neighbours = atoi(argv[i+1]);
        }
        else if(!strcmp(argv[i], "-R"))
        {
          *radius = atof(argv[i+1]);
        }
    }
}

void getArgumentsCube(int argc, char **argv, string *inputFile, string *queryFile, string *outputFile, int *dimensions, int *maxPointsCheck,int *probes, int *neighbours, double *radius)
{
    *dimensions = 14;
    *maxPointsCheck = 10;
    *probes = 2;
    *neighbours = 1;
    *radius = 10000.0;

    for(int i = 0; i < argc; i++)
    {
        if(!strcmp(argv[i], "-d"))
        {
            *inputFile = string(argv[i+1]);
        }
        else if(!strcmp(argv[i], "-q"))
        {
            *queryFile = string(argv[i+1]);
        }
        else if(!strcmp(argv[i], "-k"))
        {
          *dimensions = atoi(argv[i+1]);
        }
        else if(!strcmp(argv[i], "-M"))
        {
          *maxPointsCheck = atoi(argv[i+1]);
        }
        else if(!strcmp(argv[i], "-probes"))
        {
          *probes = atoi(argv[i+1]);
        }
        else if(!strcmp(argv[i], "-o"))
        {
          *outputFile = string(argv[i+1]);
        }
        else if(!strcmp(argv[i], "-N"))
        {
          *neighbours = atoi(argv[i+1]);
        }
        else if(!strcmp(argv[i], "-R"))
        {
          *radius = atof(argv[i+1]);
        }
    }
}

void getArgumentsCluster(int argc, char ** argv, string* inputFile, string * configFile, string * outputFile, bool* complete, string* method){
    
  *complete = false;
  *method = "All";  
  
  for(int i = 0; i < argc; i++)
  {
      if(!strcmp(argv[i], "-i"))
      {
        *inputFile = string(argv[i+1]);
      }
      else if(!strcmp(argv[i], "-c"))
      {
        *configFile = string(argv[i+1]);
      }
      else if(!strcmp(argv[i], "-o"))
      {
        *outputFile = string(argv[i+1]);
      }
      else if(!strcmp(argv[i], "-complete"))
      {
        *complete = true;
      }
      else if(!strcmp(argv[i], "-m"))
      {
        *method = string(argv[i+1]);
      }
  }
}


void readMetadata(fstream *imageFile, int *magicNumber, int *numImages, int *rows, int *columns)
{
  readInteger(imageFile, magicNumber);
  readInteger(imageFile, numImages);
  readInteger(imageFile, rows);
  readInteger(imageFile, columns);
  return;
}

void readInteger(fstream* imageFile, int* num){
  imageFile->read((char*)num, 4);
  *num = htonl(*num);
  return;
}

vector<unsigned char>* getImage(fstream* file, int size){
  unsigned char tempImage[size];
  file->read((char*)tempImage, size);
  vector<unsigned char> *v = new vector<unsigned char>(tempImage, tempImage + size);
  return v;
}



//returns a vector containing images in form of vectors
void getImages(int numImages, int pixels, vector<image*> **imagesArray, fstream *imageFile)
{
  *imagesArray = new vector<image*>(numImages, NULL);
  //current image
  image *currImage;

  for(int img = 0; img < numImages; img++)
  {
    currImage = new image(img, getImage(imageFile, pixels));
    (*imagesArray)->at(img) = currImage;
  }
  return;
}


void initializeData(int *pixels, int *numImages, vector<image*> **imagesArray, fstream **imageFile, string inputFile){
  int magicNumber, rows, columns;

  *imageFile = new fstream(inputFile, ios::in | ios::binary);
  //getting the metadata of the input file
  readMetadata(*imageFile, &magicNumber, numImages, &rows, &columns);
  //total number of pixels for each image
  *pixels = rows * columns;
  //getting a vector containing all the images
  getImages(*numImages, *pixels, imagesArray, *imageFile);
  
  (*imageFile)->close();
  delete(*imageFile);
  return;
}

void initHashTableMetadata(int hFuncs, int pixels, float radius)
{
  unsigned int smallM, bigM;
  float wAspect = 10.0;

  bigM = pow(2, 32/hFuncs);
  smallM = pow(2, 32) - 5;

  hashTable::setSmallM(smallM);
  hashTable::setBigM(bigM);
  hashTable::setWAspect(wAspect);
  hashTable::setDistrVecNum(hFuncs);
  hashTable::setDimensionality(pixels);
  hashTable::initMExponentials();

  return;
}

void initiHashTableMetada(int hFuncsNum, int dimensionsHyperCube, int pixels, int numImages, int index)
{
  unsigned int smallM, bigM;
  float wAspect = 10.0;
  int LSHBuckets = (numImages / 128);

  if(LSHBuckets == 0)
  {
    LSHBuckets = LSHBuckets + 1;
  }

  smallM = pow(2, 32) - 5;

  hashTable::setSmallM(smallM);
  hashTable::setWAspect(wAspect);
  hashTable::setDimensionality(pixels);

  // classic or lsh
  if(index == 0 || index == 1)
  {
      bigM = pow(2, 32/hFuncsNum);
      hashTable::setBigM(bigM);
      hashTable::setDistrVecNum(hFuncsNum);
      hashTable::setBucketsNum(LSHBuckets);
  }else
  {
      bigM = pow(2, 32/dimensionsHyperCube);
      hashTable::setBigM(bigM);
      hashTable::setDistrVecNum(dimensionsHyperCube);
      hashTable::setBucketsNum(pow(2, dimensionsHyperCube));
  }
  return;
}

void initHashTableMetadata(int k, int pixels, int totalImages, float radius, string &method)
{
  unsigned int smallM, bigM;
  float wAspect = 10.0;
  int LSHBuckets = (totalImages / 128);

  bigM = pow(2, 32/k);
  smallM = pow(2, 32) - 5;

  if(LSHBuckets == 0)
  {
    LSHBuckets = LSHBuckets + 1;
  }

  hashTable::setSmallM(smallM);
  hashTable::setBigM(bigM);
  hashTable::setWAspect(wAspect);
  hashTable::setDistrVecNum(k);
  hashTable::setDimensionality(pixels);
  hashTable::initMExponentials();

  if(method == "LSH")
  {
    hashTable::setBucketsNum(LSHBuckets);
  }else
  {
    hashTable::setBucketsNum(pow(2, k));
  }
  
}

void readClusterConfig(string path, int* numClusters, int* numHashtables, int* numHfunctions, int* MhyperCube, int* hyperCubeDimension, int* hyperCubeProbes){
  ifstream input(path);
  //initializing values to check later if there is an initialization value from the config file
  *numClusters = -1;
  *numHashtables = 3;
  *numHfunctions = 4;
  *MhyperCube = 10;
  *hyperCubeDimension = 3;
  *hyperCubeProbes = 2;

  for(string line; getline(input,line); ){
    if(line.find("number_of_clusters:") != string::npos){
      *numClusters = extractInt(&line);
    }else if(line.find("number_of_vector_hash_tables:") != string::npos){
      *numHashtables = extractInt(&line);
    }else if(line.find("number_of_vector_hash_functions:") != string::npos){
      *numHfunctions = extractInt(&line);
    }else if(line.find("max_number_M_hypercube:") != string::npos){
      *MhyperCube = extractInt(&line);
    }else if(line.find("number_of_hypercube_dimensions:") != string::npos){
      *hyperCubeDimension = extractInt(&line);
    }else if(line.find("number_of_probes:") != string::npos){
      *hyperCubeProbes = extractInt(&line);
    }
  }

  if(*numClusters == -1){
    cout << "ERROR ~ Number of Clusters is uninitialized - Check config file!" << endl;
  }

  return;
}

int extractInt(string* line){
  stringstream ss;     
  
  ss << *line; 
  
  string temp; 
  int found; 
  while (!ss.eof()) { 
    ss >> temp; 

    if (stringstream(temp) >> found) 
      return found;
  } 
  return -1;
}

vector<image*>* getTrueNN(int N, image* queryImage, vector<image*>* imageArray, unsigned int (*distanceFunction) (vector<unsigned char>*, vector<unsigned char>*)){
  priority_queue<imageDuplet, vector<imageDuplet>, compDuplet> nearestNeighbours;
  unsigned int currDistance;
  image* currImage;
  int imageArraySize = imageArray->size();

  for(int i = 0; i < imageArraySize; i++){
    currImage = imageArray->at(i);
    currDistance = distanceFunction(queryImage->getImage(), currImage->getImage());
    
    if((int)nearestNeighbours.size() != N){
      nearestNeighbours.push(imageDuplet(currDistance, currImage));
    }else{
      if(currDistance < nearestNeighbours.top().distance){
        nearestNeighbours.pop();
        nearestNeighbours.push(imageDuplet(currDistance, currImage));
      }
    }
  }
  
  vector<image*>* wantedVector = new vector<image*>(N,NULL);
  for(int i = 0; i < N; i++){
    wantedVector->at(i) = nearestNeighbours.top().imageData;
    nearestNeighbours.pop();
  }

  return wantedVector;
}

void writeOutputLSH(int imageNum, image* queryImage, vector<image*>* nearestNeighbours, vector<image*>* nearestNeighboursTrue, duration<double> timeElapsed, duration<double> timeElapsedTrue, vector<image*>* nearestRadiusNeighbours, ofstream* outfile, unsigned int (*distanceFunction) (vector<unsigned char>*, vector<unsigned char>*)){
  int Nminus1 = nearestNeighbours->size()-1;
  int N = Nminus1 + 1;
  *outfile << "Query: " << imageNum << endl;
  
  for(int i = Nminus1; i >= 0; i--){
    *outfile << "Nearest neighbour-" << N-i << ": " << nearestNeighbours->at(i)->getId() << endl;
    *outfile << "distanceLSH: " << manhattanDistance(nearestNeighbours->at(i)->getImage(), queryImage->getImage()) << endl;
    *outfile << "distanceTrue: " << manhattanDistance(nearestNeighboursTrue->at(i)->getImage(), queryImage->getImage()) << endl;
  }

  *outfile << "tLSH:" << timeElapsed.count() << endl;
  *outfile << "tTrue:" << timeElapsedTrue.count() << endl;
  *outfile << "R-nearest neighbours:" << endl;

  int radiusNNsize = nearestRadiusNeighbours->size()-1;
  for(int i = radiusNNsize; i >= 0; i--){
    *outfile << nearestRadiusNeighbours->at(i)->getId() << endl;
  }
}

void writeOutputHyperCube(int imageNum, image* queryImage, vector<image*>* nearestNeighbours, vector<image*>* nearestNeighboursTrue, duration<double> timeElapsed, duration<double> timeElapsedTrue, vector<image*>* nearestRadiusNeighbours, ofstream* outfile, unsigned int (*distanceFunction) (vector<unsigned char>*, vector<unsigned char>*)){
  int Nminus1 = nearestNeighbours->size()-1;
  int N = Nminus1 + 1;
  *outfile << "Query: " << imageNum << endl;
  
  for(int i = Nminus1; i >= 0; i--){
    *outfile << "Nearest neighbour-" << N-i << ": " << nearestNeighbours->at(i)->getId() << endl;
    *outfile << "distanceHyperCube: " << manhattanDistance(nearestNeighbours->at(i)->getImage(), queryImage->getImage()) << endl;
    *outfile << "distanceTrue: " << manhattanDistance(nearestNeighboursTrue->at(i)->getImage(), queryImage->getImage()) << endl;
  }

  *outfile << "tHyperCube:" << timeElapsed.count() << endl;
  *outfile << "tTrue:" << timeElapsedTrue.count() << endl;
  *outfile << "R-nearest neighbours:" << endl;

  int radiusNNsize = nearestRadiusNeighbours->size()-1;
  for(int i = radiusNNsize; i >= 0; i--){
    *outfile << nearestRadiusNeighbours->at(i)->getId() << endl;
  }
}

void writeOutputCluster(string* method, clustering* myClustering, double elapsedTime, double silhouette, vector<double>* silhouettes, bool complete, ofstream* output, int index){
  // total number of clusters
  int totalClusters = myClustering->getClustersNum();
  cluster* currCluster;
  image* currCentroid;
  int imageSize, silhouettesSize;

  vector<cluster*>* clusters = myClustering->getClusters();
  vector<image*>* centroids = myClustering->getCentroids();
  if(index == 0){
    *method = "Classic";
  }else if(index == 1){
    *method = "LSH";
  }else if(index == 2){
    *method = "Hypercube";
  }

  *output << "Algorithm: " << *method << endl;
  
  for(int i = 0; i < totalClusters; i++){
    currCluster = clusters->at(i);
    currCentroid = centroids->at(i);
    imageSize = currCentroid->getImage()->size();
    
    *output << "CLUSTER-" << i+1 << " {size: " << currCluster->getClusterSize() << ", centroid: [";
    
    for(int j = 0; j < imageSize-1; j++){
      *output << (int) currCentroid->getImage()->at(j) << ", ";
    }
    *output << (int) currCentroid->getImage()->at(imageSize-1) << "]}" << endl;
  }

  *output << "clustering_time: " << elapsedTime << endl;

  *output << "Silhouette: [";
  silhouettesSize = silhouettes->size();
  for(int i = 0; i < silhouettesSize; i++){
    *output << silhouettes->at(i) << ", ";
  }
  *output << silhouette << "]" << endl;



  if(complete)
  {
    for(int cluster = 0; cluster < totalClusters; cluster++)
    {
      
      currCentroid = centroids->at(cluster);
      imageSize = currCentroid->getImage()->size();
      *output << "CLUSTER-" << cluster+1 << " {"; 
      //centroid
      *output << "[";
      for(int i = 0; i < imageSize-1; i++){
        *output << (int) currCentroid->getImage()->at(i) << ", ";
      }
      *output << (int) currCentroid->getImage()->at(imageSize-1) << "]";

      set<image*>* clusterSet;
      clusterSet = clusters->at(cluster)->getClust();

      for(image* imageFromSet : *clusterSet){
        *output << ", " << imageFromSet->getId();
      }
      *output << "}" << endl;

    }
  }    

  // freeing the memory of the silhouette vector
  delete silhouettes;
}

void deleteforLSH(LSH** myLSH, vector<image*>** images, int numImages){
    delete(*myLSH);
    free(hashTable::getMExponentials());
    for(int i = 0; i < numImages; i++){
        delete((*images)->at(i));
    }
    delete(*images);
}

void deleteforHyperCube(hyperCube** myHyperCube, vector<image*>** images, int numImages){
    delete(*myHyperCube);
    free(hashTable::getMExponentials());
    for(int i = 0; i < numImages; i++){
        delete((*images)->at(i));
    }
    delete(*images);

}


void deleteforCluster(clustering** myCluster, vector<image*>** images, int numImages){
    delete(*myCluster);
    free(hashTable::getMExponentials());
    for(int i = 0; i < numImages; i++){
        delete((*images)->at(i));
    }
    delete(*images);
}

void deleteQueries(vector<image*>** queries){
  int queriesSize = (*queries)->size();

  for(int i = 0; i < queriesSize; i++){
    delete((*queries)->at(i));
  }
  delete(*queries);
}
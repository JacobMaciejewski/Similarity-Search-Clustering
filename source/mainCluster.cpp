#include "../headers/mainAux.h"


using namespace std;
using namespace std::chrono;

int main(int argc, char** argv){
    string inputFile, outputFile, configFile, method;
    bool complete;
    int clustersNum, hashTablesNum, hfuncsNum, maxChecksHyperCube, dimensionsHyperCube, probesHyperCube;
    int pixels;
    int numImages;
    double elapsedTime, silhouetteTotal;
    vector<double>* silhouettes;

    //array containing references to images pixels values
    vector<image*> *imagesArray;
    fstream *imageFile;
    
    //getting arguments from terminal
    getArgumentsCluster(argc, argv, &inputFile, &configFile, &outputFile, &complete, &method);
    //reading the config file and initializing values
    readClusterConfig(configFile, &clustersNum, &hashTablesNum, &hfuncsNum, &maxChecksHyperCube, &dimensionsHyperCube, &probesHyperCube);
    // initializing images from train file
    initializeData(&pixels, &numImages, &imagesArray, &imageFile, inputFile);
    // initializing images from query file
    //creating output file
    ofstream output(outputFile);
    clustering *myClustering;
    // initializing random number generator
    srand(time(NULL));
    if(method == "LSH")
    {
        initHashTableMetadata(hfuncsNum, pixels, numImages, 1.0, method);
        // initializing a new cluster
        myClustering = new clustering(clustersNum, hashTablesNum, (unsigned int) maxChecksHyperCube, (unsigned int)probesHyperCube, &method, imagesArray);
        // executing clustering and returning the elapsed time
        elapsedTime = myClustering->execute(manhattanDistance);
        // computing the silhouettes
        silhouetteTotal = myClustering->calculateSilhouette(manhattanDistance, &silhouettes);
        // printing output into file
        writeOutputCluster(&method, myClustering, elapsedTime, silhouetteTotal, silhouettes, complete, &output, -1);
        // delete clusters and the left out memory
        deleteforCluster(&myClustering, &imagesArray, numImages);
    }else if(method == "Hypercube" || method == "Classic")
    {
        initHashTableMetadata(dimensionsHyperCube, pixels, numImages, 1.0, method);
        // initializing a new cluster
        myClustering = new clustering(clustersNum, hashTablesNum, (unsigned int) maxChecksHyperCube, (unsigned int)probesHyperCube, &method, imagesArray);
        // executing clustering and returning the elapsed time
        elapsedTime = myClustering->execute(manhattanDistance);
        // computing the silhouettes
        silhouetteTotal = myClustering->calculateSilhouette(manhattanDistance, &silhouettes);
        // printing output into file
        writeOutputCluster(&method, myClustering, elapsedTime, silhouetteTotal, silhouettes, complete, &output, -1);
        // delete clusters and the left out memory
        deleteforCluster(&myClustering, &imagesArray, numImages);
    }else if(method == "All")
    {
        clustering *myClustering = new clustering(clustersNum, hashTablesNum, (unsigned int) maxChecksHyperCube, (unsigned int)probesHyperCube, &method, imagesArray);
        // initializing m modulo exponentials, keen for all methods
        for(int i = 0; i < 3; i++)
        {
            initiHashTableMetada(hfuncsNum, dimensionsHyperCube, pixels, numImages, i);
            hashTable::initMExponentials(i, hfuncsNum, dimensionsHyperCube);
            elapsedTime = myClustering->executeMethod(i, hashTablesNum, (unsigned int) maxChecksHyperCube, (unsigned int)probesHyperCube, manhattanDistance);
            // computing the silhouettes
            silhouetteTotal = myClustering->calculateSilhouette(manhattanDistance, &silhouettes);
            writeOutputCluster(&method, myClustering, elapsedTime, silhouetteTotal, silhouettes, complete, &output, i);
        }
        // delete clusters and the left out memory
        deleteforCluster(&myClustering, &imagesArray, numImages);
    }
    return 0;   

}
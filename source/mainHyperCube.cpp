#include "../headers/mainAux.h"
#include "../headers/hyperCube.h"

using namespace std;
using namespace std::chrono;

int main(int argc, char** argv){
    string inputFile, queryFile, outputFile;
    int dimensions, maxChecks, neighbours, probes;
    int pixels;
    int numImages, numImagesQuery;
    double radius;
    string method = "Hypercube";
    //array containing references to images pixels values
    vector<image*> *imagesArray, *queryArray;
    fstream *imageFile, *queryFileStream;  
    vector<image*> *nearestNeighbours, *nearestRadiusNeighbours, *nearestNeighboursTrue;
    image* currImage;

    srand(time(NULL));
    //getting arguments from terminal
    getArgumentsCube(argc, argv, &inputFile, &queryFile, &outputFile,  &dimensions, &maxChecks, &probes, &neighbours, &radius);
    //initializing images from train file
    initializeData(&pixels, &numImages, &imagesArray, &imageFile, inputFile);
    // initialize hashtable metadata
    initHashTableMetadata(dimensions, pixels, numImages, (float)radius, method);
    hyperCube* myHyperCube = new hyperCube(maxChecks, probes);
    //initializing the Hyper Cube with the train images
    myHyperCube->initializeHyperCube(imagesArray);
    //this is the output file
    ofstream outfile(outputFile);
    do{
        //initializing images from query file
        initializeData(&pixels, &numImagesQuery, &queryArray, &queryFileStream, queryFile);

        for(int i = 0; i < numImagesQuery; i++){
            currImage = queryArray->at(i);

            high_resolution_clock::time_point t1 = high_resolution_clock::now();
            //getting the N nearest neighbours
            nearestNeighbours = myHyperCube->getNN(currImage, neighbours, manhattanDistance);
            high_resolution_clock::time_point t2 = high_resolution_clock::now();
            //calculating elapsed time
            duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
            //getting the nearest neighbour
            //nearestNeighbour = nearestNeighbours->at(neighbours-1);
            //getting the nearest neighbours within radius R
            set<int> checkedImages;
            nearestRadiusNeighbours = myHyperCube->radiusSearch(currImage, radius, manhattanDistance, checkedImages);

            //caclulate TRUE distance
            t1 = high_resolution_clock::now();
            nearestNeighboursTrue = getTrueNN(neighbours, currImage, imagesArray, manhattanDistance);
            t2 = high_resolution_clock::now();
            //calculating elapsed time
            duration<double> time_spanTrue = duration_cast<duration<double>>(t2 - t1);
            writeOutputHyperCube(i, currImage ,nearestNeighbours, nearestNeighboursTrue, time_span, time_spanTrue, nearestRadiusNeighbours, &outfile, manhattanDistance);
            delete(nearestNeighbours);
            delete(nearestNeighboursTrue);
            delete(nearestRadiusNeighbours);
        }

        deleteQueries(&queryArray);

        cout << "Insert new query file(path) or exit" << endl;
        cin >> queryFile;
        
    }while(queryFile != "exit");
    
    deleteforHyperCube(&myHyperCube, &imagesArray, numImages);

    return 0;
}
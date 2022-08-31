#include "../headers/mainAux.h"
#include "../headers/LSH.h"


using namespace std;
using namespace std::chrono;

int main(int argc, char** argv){
    string inputFile, queryFile, outputFile;
    int hFuncs, hashTablesNum, neighbours;
    int pixels;
    int numImages, numImagesQuery;
    double radius;
    //array containing references to images pixels values
    vector<image*> *imagesArray, *queryArray;
    fstream *imageFile, *queryFileStream;
    string method = "LSH";
    vector<image*> *nearestNeighbours, *nearestRadiusNeighbours ,*nearestNeighboursTrue;
    image* currImage;
    srand(time(NULL));

    //getting arguments from terminal
    getArgumentsLSH(argc, argv, &inputFile, &queryFile, &outputFile, &hFuncs, &hashTablesNum, &neighbours, &radius);
    //initializing images from train file
    initializeData(&pixels, &numImages, &imagesArray, &imageFile, inputFile);
    // initialize hashtable metadata
    initHashTableMetadata(hFuncs, pixels, numImages, (float) radius, method);
    //this is the output file
    ofstream outfile(outputFile);
    LSH* myLSH = new LSH(hashTablesNum);
    //initializing the LSH with the train images
    myLSH->initializeLSH(imagesArray);
    
    do{  
        //initializing images from query file
        initializeData(&pixels, &numImagesQuery, &queryArray, &queryFileStream, queryFile);

        for(int i = 0; i < numImagesQuery; i++){

            currImage = queryArray->at(i);

            high_resolution_clock::time_point t1 = high_resolution_clock::now();
            //getting the N nearest neighbours
            nearestNeighbours = myLSH->getNN(currImage, neighbours, manhattanDistance);
            high_resolution_clock::time_point t2 = high_resolution_clock::now();
            //calculating elapsed time
            duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
            //getting the nearest neighbours within radius R
            set<int> checkedImages;
            nearestRadiusNeighbours = myLSH->radiusSearch(currImage, radius, manhattanDistance, checkedImages);

            //caclulate TRUE distance
            t1 = high_resolution_clock::now();
            nearestNeighboursTrue = getTrueNN(neighbours, currImage, imagesArray, manhattanDistance);
            t2 = high_resolution_clock::now();
            //calculating elapsed time
            duration<double> time_spanTrue = duration_cast<duration<double>>(t2 - t1);
            writeOutputLSH(i, currImage ,nearestNeighbours, nearestNeighboursTrue, time_span, time_spanTrue, nearestRadiusNeighbours, &outfile, manhattanDistance);
            delete(nearestNeighbours);
            delete(nearestNeighboursTrue);
            delete(nearestRadiusNeighbours);
        }

        deleteQueries(&queryArray);
        cout << "Insert new query file(path) or exit" << endl;
        cin >> queryFile;
    }while(queryFile != "exit");

    deleteforLSH(&myLSH, &imagesArray, numImages);

    
    return 0;
}
class hashTable;
class clustering;
class LSH;
class hyperCube;

#ifndef MAIN_AUX
#define MAIN_AUX
#include <iostream>
#include <string>
#include <fstream>
#include <iterator>
#include <string.h>
#include <vector>
#include <array>
#include <bits/stdc++.h>
#include <arpa/inet.h>
#include <random>
#include <cstdlib>
#include "cmath"
#include "hashTable.h"
#include "image.h"
#include "bucket.h"
#include "clustering.h"
#include "algoAux.h"
#include <chrono>
#include <ratio>
#include <ctime>
#include <iostream>
#include <fstream>


using namespace std;
using namespace std::chrono;
void getArgumentsLSH(int, char **, string *, string *, string *, int *, int *, int *, double *);
void getArgumentsCube(int , char **, string *, string *, string *, int *, int *, int*, int *, double *);
void getArgumentsCluster(int, char **, string*, string *, string *, bool*, string*);
void initializeData(int*, int*, vector<image*>**, fstream**, string);;
void readMetadata(fstream*, int*, int*, int*, int*);
void readInteger(fstream*, int*);
vector<unsigned char>* getImage(fstream*, int);
void getImages(int, int, vector<image*>**, fstream*);
void initHashTableMetadata(int, int, float);
void initHashTableMetadata(int, int, int, float, string&);
void readClusterConfig(string, int*, int*, int*, int*, int*, int*);
vector<unsigned char>* getImage(fstream*, int);
int extractInt(string*);

vector<image*>* getTrueNN(int, image*, vector<image*>*, unsigned int (*distanceFunction) (vector<unsigned char>*, vector<unsigned char>*));
void writeOutputLSH(int, image*, vector<image*>*, vector<image*>*, duration<double>, duration<double>, vector<image*>*, ofstream*, unsigned int (*distanceFunction) (vector<unsigned char>*, vector<unsigned char>*));
void writeOutputHyperCube(int, image*, vector<image*>*, vector<image*>*, duration<double>, duration<double>, vector<image*>*, ofstream*, unsigned int (*distanceFunction) (vector<unsigned char>*, vector<unsigned char>*));
void writeOutputCluster(string*, clustering*, double, double, vector<double>*,bool, ofstream*, int);
void deleteforLSH(LSH**, vector<image*>**, int);
void deleteforHyperCube(hyperCube**, vector<image*>**, int);
void deleteforCluster(clustering**, vector<image*>**, int);
void deleteQueries(vector<image*>**);
void initiHashTableMetada(int, int, int, int, int);
#endif
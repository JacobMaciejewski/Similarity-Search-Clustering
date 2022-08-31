
#include "../headers/hashTable.h"

unsigned int hashTable::bigM;
unsigned int hashTable::smallM;
unsigned int * hashTable::mExponentials;
unsigned int hashTable::bucketsNum;
float hashTable::wAspect;
int hashTable::distrVecNum;
int hashTable::imgDimensions;

hashTable::hashTable()
{
    unsigned int totalBuckets = this->getBucketsNum();
    //allocating memory for the array of disturbance vectors
    this->distrVectors = new vector<vector<float>*>(this->getDistrVecNum(), NULL);
    //all disturbance vectors initially have zero magnitude
    for(int i = 0; i < this->distrVecNum; i++)
    {
        this->distrVectors->at(i) = new vector<float>(this->getDimensionality(), 0.0);
    }

    //initializing disturbance vectors
    this->initDistrVectors();

    this->buckets = (bucket**)malloc(sizeof(bucket*)*totalBuckets);
    for(unsigned int i = 0; i < totalBuckets; i++){
        this->buckets[i] = new bucket();
    }
}

hashTable::hashTable(unsigned int cubeDimensionality)
{
    unsigned int totalCorners = pow(2, cubeDimensionality);
    //allocating memory for the array of disturbance vectors
    this->distrVectors = new vector<vector<float>*>(this->getDistrVecNum(), NULL);
    //all disturbance vectors initially have zero magnitude
    for(int i = 0; i < this->distrVecNum; i++)
    {
        this->distrVectors->at(i) = new vector<float>(this->getDimensionality(), 0.0);
    }

    //initializing disturbance vectors
    this->initDistrVectors();

    this->buckets = (bucket**)malloc(sizeof(bucket*)*totalCorners);
    for(unsigned int i = 0; i < totalCorners; i++){
        this->buckets[i] = new bucket();
    }
}

hashTable::~hashTable()
{

    for(int i = 0; i < this->distrVecNum; i++)
    {
        delete(this->distrVectors->at(i));
    }

    delete(this->distrVectors);

    for(unsigned int i = 0; i < this->getBucketsNum(); i++){
        delete(this->buckets[i]);
    }

    free(this->buckets);
}

//initializing k real, uniformly distributed disturbance vectors
void hashTable::initDistrVectors()
{
    int totalVectors = this->getDistrVecNum();
    int totalDimensions = this->getDimensionality();
    //initialize uniform random distribution
    random_device rd;
    mt19937 e2(rd());
    //getting random float numbers from 0 to w aspect
    uniform_real_distribution<> dist(0, this->getAspectW());

    //initialize each disturbance vector with real values
    for(int index = 0; index < totalVectors; index++)
    {
        //updating each coordinate of current disturbance vector
        for(int dimension = 0; dimension < totalDimensions; dimension++)
        {
            this->getIthDistrVector(index)->at(dimension) = dist(e2);
        }
    }

    return;
}

void hashTable::initMExponentials()
{
    unsigned int m = hashTable::getSmallAspectM();
    unsigned int M = hashTable::getBigAspectM();
    //used to update exponential modulo of m in each iteration
    //is equal to m mod M
    unsigned int mModulo = m % M;
    unsigned int dims = hashTable::getDimensionality();
    unsigned int *expModulos;

    expModulos = (unsigned int *) malloc(sizeof(unsigned int) * dims);
    hashTable::setMExponentials(expModulos);
    //m ^ 0 is equal to 1, the coefficient of a_0
    expModulos[0] = 1;

    //initializing the modulos of exponential m factors using formulas described in h modulo function
    for(unsigned int i = 1; i < dims; i++)
    {
        expModulos[i] = (expModulos[i - 1] * mModulo) % M;
    }
    return;
}

void hashTable::initMExponentials(int index, int hfuncsNum, int dimensionsHyperCube)
{
    unsigned int m = hashTable::getSmallAspectM();
    unsigned int M = hashTable::getBigAspectM();
    //used to update exponential modulo of m in each iteration
    //is equal to m mod M
    unsigned int mModulo = m % M;
    unsigned int dims = hashTable::getDimensionality();
    unsigned int *expModulos;
    // should update the exponentials array
    // if the k is different in LSH and HyperCube method
    if(index != 0 && hfuncsNum != dimensionsHyperCube)
    {
        free(hashTable::mExponentials);
        expModulos = (unsigned int *) malloc(sizeof(unsigned int) * dims);
    }else if(index == 0) //first method, memory allocation is necessary
    {
        expModulos = (unsigned int *) malloc(sizeof(unsigned int) * dims);
    }
    
    hashTable::setMExponentials(expModulos);
    //m ^ 0 is equal to 1, the coefficient of a_0
    expModulos[0] = 1;

    //initializing the modulos of exponential m factors using formulas described in h modulo function
    for(unsigned int i = 1; i < dims; i++)
    {
        expModulos[i] = (expModulos[i - 1] * mModulo) % M;
    }
    return;
}

void hashTable::setSmallM(unsigned int m)
{
    hashTable::smallM = m;
}

void hashTable::setBigM(unsigned int M)
{
    hashTable::bigM = M;
}

void hashTable::setWAspect(float w)
{
    hashTable::wAspect = w;
}

void hashTable::setDistrVecNum(unsigned int hFuncs)
{
    hashTable::distrVecNum = hFuncs;
}
void hashTable::setDimensionality(unsigned int dim)
{
    hashTable::imgDimensions = dim;
}

void hashTable::setMExponentials(unsigned int *array)
{
    hashTable::mExponentials = array;
}

void hashTable::setBucketsNum(unsigned int buckets)
{
    hashTable::bucketsNum = buckets;
}

//return the output value for input vector for h_i(x)
unsigned int hashTable::hFunction(int functionIndex, image *inputVector)
{
    unsigned int hFunctionValue;
    //total number of dimensions of image vectors
    int totalPixels = this->getDimensionality();
    //parts of a_i formula a_ij = x_i - s_ij / w
    float sub, div, currentPixel, currDisturbance;
    //input vector will be substracted by disturbance vector and normalized by w aspect
    vector<unsigned int>* normalizedInputVector = new vector<unsigned int>(totalPixels,0);

    //constructing the normalized vector following x_i - s_ij / w
    //for each of the pixels
    for(int dimension = 0; dimension < totalPixels; dimension++)
    {
        currentPixel = (float)inputVector->getIthPixel(dimension);
        currDisturbance = this->getIthDistrVector(functionIndex)->at(dimension);

        if(currentPixel > currDisturbance){
            sub = currentPixel - currDisturbance;
        }else{
            sub = currDisturbance - currentPixel;
        }

        div = float(sub / this->getAspectW());
        normalizedInputVector->at(dimension) = floor(div);
    }

    //following the recursive modulo formula in order to avoid overflow
    hFunctionValue = this->hModulo(normalizedInputVector, this->getDimensionality());
    delete normalizedInputVector;
    return hFunctionValue;
}

//We are trying to calculate sum{x_i} % M, where x_i = a_i * m ^ {d - 1 - i}
//a_i being natural and d being the dimension of the given vector, avoiding overflow

//It can be proven by induction that this is equal to incrementally adding the modulo
//of the previous x_{i-1} to the modulo of x_{i} and modulling the whole result

//x_i modulo calculation can be simplified by employing the modular multiplication formula
//where we have to calculate the modulo of the exponential coefficient m ^ d - 1 - i

//Using the modular multiplication formula, we can store the result of the previous 
//exponential coefficient modulo, multiply it by the modulo of m (which will be calculated one time as an constant)
//and modulo the result of the multiplication. In this way we employ only two mathematical operations in each iteration.

//As a result, we get a polynomial coefficient modulation algorithm with O(d) complexity (d being the image dimensionality)

unsigned int hashTable::hModulo(vector<unsigned int> *polyCoeffs, int totalCoeffs)
{
    unsigned int expModulo, aModulo, xModulo, result;
    unsigned int M = this->getBigAspectM();

    //modulo of the exponential m ^ d - i - 1
    expModulo = getIthMExponential(0);
    //modulo of current a_i coefficient
    aModulo = polyCoeffs->at(0) % M;
    //following the described formula for the first term x_0 
    xModulo = (expModulo * aModulo) % M;

    result = xModulo;


    for(int dimension = 1; dimension < totalCoeffs; dimension++)
    {
        aModulo = polyCoeffs->at(dimension) % M;
        //updating the exponential m modulo using the previous value and the set mModulo
        expModulo = getIthMExponential(dimension);
        //calculating the multiplication modulo of current iteration
        xModulo = (aModulo * expModulo) % M;
        //updated result is equal to the previous plus the current one, modulo
        result = (result + xModulo) % M;
    }

    return result;
}

unsigned int* hashTable::getMExponentials()
{
    return hashTable::mExponentials;
}
unsigned int hashTable::getIthMExponential(int index)
{
    return hashTable::getMExponentials()[index];
}

//returns a value that moded by the number of buckers, returns the
//corresponding bucket to given vector
unsigned int hashTable::gFunction(image *inputVector)
{   
    int totalhFuncs = this->getDistrVecNum();
    int shiftSize = totalhFuncs/sizeof(unsigned int); 
    unsigned int finalValue = 0;
    
    finalValue = this->hFunction(0, inputVector);

    for(int i = 1; i < totalhFuncs; i++){
        finalValue <<= shiftSize;
        finalValue |= this->hFunction(i,inputVector);
    }

    return finalValue;
}

unsigned int hashTable::getBucket(image *inputVector)
{
    return this->gFunction(inputVector) % this->getBucketsNum();
}

unsigned int hashTable::getSmallAspectM()
{
    return hashTable::smallM;
}

unsigned int hashTable::getBigAspectM()
{
    return hashTable::bigM;
}

unsigned int hashTable::getBucketsNum()
{
    return hashTable::bucketsNum;
}


float hashTable::getAspectW()
{
    return hashTable::wAspect;
}

int hashTable::getDistrVecNum()
{
    return hashTable::distrVecNum;
}

int hashTable::getDimensionality()
{
    return hashTable::imgDimensions;
}


vector<vector<float>*>* hashTable::getDistrVectors()
{
    return this->distrVectors;
}


vector<float>* hashTable::getIthDistrVector(int index)
{
    return this->getDistrVectors()->at(index);
}

void hashTable::pushImage(image* newImage){
    int wantedBucket = this->getBucket(newImage);
    this->buckets[wantedBucket]->pushImage(newImage);
}

void hashTable::pushImages(vector<image*>* images){
    int totalImages = images->size();
    for(int i = 0; i < totalImages; i++){
        pushImage(images->at(i));
    }
}

void hashTable::printBucketSizes()
{
    unsigned int totalBuckets = this->getBucketsNum();
    for(unsigned int currBucket = 0; currBucket < totalBuckets; currBucket++)
    {
        cout << "[" << currBucket + 1 << "]: " << this->buckets[currBucket]->getSize() << " entries" << endl;
    }
    cout << endl;
}

void hashTable::printFirstNBucketSizes(int buckets)
{
    for(int currBucket = 0; currBucket < buckets; currBucket++)
    {
        cout << "[" << currBucket + 1 << "]: " << this->buckets[currBucket]->getSize() << " entries" << endl;
    }
    cout << endl;
}

bucket** hashTable::getBuckets(){
    return this->buckets;
}

bucket* hashTable::getIthBucket(int index)
{
    return this->buckets[index];
}
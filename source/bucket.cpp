#include "../headers/bucket.h"
//TO MAKEFILE
bucket::bucket(){
    this->images = new vector<image*>;
}

bucket::~bucket(){
    delete(this->images);
}

vector<image*>* bucket::getImages(){
    return this->images;
}

image* bucket::getIthImage(int index){
    return this->images->at(index);
}

void bucket::pushImage(image* newImage){
    this->images->push_back(newImage);
}

int bucket::getSize(){
    return this->getImages()->size();
}



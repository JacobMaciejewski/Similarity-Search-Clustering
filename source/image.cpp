#include "../headers/image.h"

image::image(int id, vector<unsigned char>* image){
    this->imageData = image;
    this->id = id;
}

image::~image(){
    delete(this->imageData);
}

vector<unsigned char>* image::getImage(){
    return this->imageData;
}

int image::getId(){
    return this->id;
}

void image::setId(unsigned int newId){
    this->id = newId;
}

void image::setImage(vector<unsigned char>* newImage){
    this->imageData = newImage;
}

unsigned char image::getIthPixel(int index){
    return this->getImage()->at(index);
}

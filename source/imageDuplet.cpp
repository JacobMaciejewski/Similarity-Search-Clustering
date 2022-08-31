#include "../headers/imageDuplet.h"

// struct compDuplet{ 
//     bool operator()(imageDuplet & p1, imageDuplet & p2) 
//     { 
//         return p1.distance < p2.distance; 
//     } 
// };

imageDuplet::imageDuplet(int givenDistance, image *givenImage)
{
    this->distance = givenDistance;
    this->imageData = givenImage;
}

imageDuplet::~imageDuplet()
{

}

vector<unsigned char>* imageDuplet::getImageVector()
{
    image *dupletImage = this->getImage();
    return dupletImage->getImage();
}

int imageDuplet::getImageId()
{
    image *dupletImage = this->getImage();
    return dupletImage->getId();

}

image* imageDuplet::getImage()
{
    return this->imageData;
}

unsigned int imageDuplet::getDistance()
{
    return this->distance;
} 
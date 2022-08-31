#ifndef ALGO_IMG_DUPLET
#define ALGO_IMG_DUPLET
#include "../headers/image.h"

// duplet class containing image's data and its distance to queue
class imageDuplet{
    private:
    public:
        unsigned int distance;
        image *imageData;
        imageDuplet(int, image*);
        ~imageDuplet();
        vector<unsigned char>* getImageVector();
        int getImageId();
        image* getImage();
        unsigned int getDistance();      
};



#endif
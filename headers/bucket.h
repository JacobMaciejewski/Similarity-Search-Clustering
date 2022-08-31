#ifndef BUCKET
#define BUCKET

#include <vector>
#include "image.h"

class bucket{
    private:
        vector<image*> *images;
    public:
        bucket();
        ~bucket();
        vector<image*> *getImages();
        image* getIthImage(int);
        int getSize();
        void pushImage(image*);
};

#endif
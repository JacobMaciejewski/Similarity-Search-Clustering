#ifndef IMAGE
#define IMAGE
#include <vector>
using namespace std;

class image{
    private:
        vector<unsigned char>* imageData;
        int id;
    public:
        image(int, vector<unsigned char>*);
        ~image();
        vector<unsigned char> *getImage();
        int getId();
        void setId(unsigned int);
        void setImage(vector<unsigned char>*);
        unsigned char getIthPixel(int);
};

#endif
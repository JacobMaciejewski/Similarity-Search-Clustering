 #include "../headers/clusterEntry.h"
 
clusterEntry::clusterEntry()
{
    this->setClusterID(-1);
    this->setDecidability(false);
    this->setDistance(0);
}


clusterEntry::~clusterEntry()
{

}

int clusterEntry::getClusterID()
{
    return this->clusterID;
}

bool clusterEntry::isConsidered()
{
    return this->undecided;
}

unsigned int clusterEntry::getDistanceToCluster()
{
    return this->clusterDistance;
}


void clusterEntry::setClusterID(int givenID)
{
    this->clusterID = givenID;
}

void clusterEntry::setDecidability(bool state)
{
    this->undecided = state;
}

void clusterEntry::setDistance(unsigned int givenDistance)
{
    this->clusterDistance = givenDistance;
}
/**********************************************************************************************
************************************ O(6) MODEL MONTE CODE ************************************
***********************************************************************************************
* Lauren E. Hayward Sierens
***********************************************************************************************
* File:   Lattice.h (Abstract Class)
**********************************************************************************************/

#ifndef LATTICE_H
#define LATTICE_H

#include <cmath>

class Lattice 
{ 
  public:
    typedef unsigned int  uint;
  
  protected: 
    uint   N_;          //total number of lattice sites
    uint   z_;          //number of nearest neighbouring sites for each site
    uint** neighbours_; //coordinates of each vertex's nearest neighbours
    
  public:
    Lattice();
    virtual ~Lattice();
    
    //pure virtual methods (to be implemented by all child classes):
    virtual uint  getNeighbour(uint i, uint j) = 0;
    virtual void  printParams() = 0;
    virtual void  printNeighbours() = 0;
    
    //getter methods:
    uint getN();
    uint getZ();
};

#endif  // LATTICE_H

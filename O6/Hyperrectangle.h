/**********************************************************************************************
************************************ O(6) MODEL MONTE CODE ************************************
***********************************************************************************************
* Lauren E. Hayward Sierens
***********************************************************************************************
* File:   Hyperrectangle.h
* Description: Note that L_[0] == Lx = Lattice length along x
*                        L_[1] == Ly = Lattice length along y
*                          .
*                          .
*                          .
*
*              The vertices are labelled along the x-direction, followed by the 
*              y-direction, etc. 
*              So for example, the first Lx*Ly (in the x-y plane in D>=2) are:
*                (Lx*Ly - Lx)    (Lx*Ly - Lx + 1)    . . .    (Lx*Ly - 1)
*                     .                 .                          .
*                     .                 .                          .
*                     .                 .                          .
*                     Lx               Lx+1          . . .      2Lx - 1
*                     0                 1            . . .       Lx-1
**********************************************************************************************/

#ifndef HYPERRECTANGLE_H 
#define HYPERRECTANGLE_H

#include <fstream>
#include "Lattice.h"

class Hyperrectangle : public Lattice
{ 
  public:
    typedef unsigned int  uint;
    
  private:
    uint  D_; //dimension
    uint* L_; //length in each dimension
    
    void   initNeighbours();
    int    round(double num);
    //void   trimWhiteSpace(std::string* word);
    uint   nextInDir(uint dir); //used by initNeighbours method
    uint   uintPower(uint base, uint exp);
    
  public:
    Hyperrectangle(std::ifstream* fin, std::string fileName);
    virtual ~Hyperrectangle();
    
    virtual uint getNeighbour(uint i, uint j);
    virtual void printParams();
    virtual void printNeighbours();
    
    //getter methods:
    uint  getD();
    uint* getL();
};

#endif  // HYPERRECTANGLE_H

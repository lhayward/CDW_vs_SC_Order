/**********************************************************************************************
******************************** CLASSICAL REPLICA MONTE CODE *********************************
***********************************************************************************************
* Lauren Hayward
***********************************************************************************************
* File:   Hypercube.h
* Description: The vertices of a are labelled along the x-direction, followed by the 
*              y-direction, etc. 
*              So for example, the first L^2 vertices (in the x-y plane in D>=2) are:
*                   (L^2-L)  (L^2-L+1)  . . .  (L^2-1)
*                      .         .                .
*                      .         .                .
*                      .         .                .
*                      L        L+1     . . .    2L-1
*                      0         1      . . .    L-1
**********************************************************************************************/

#ifndef HYPERCUBE_H 
#define HYPERCUBE_H

#include <fstream>
#include "Lattice.h"

class Hypercube : public Lattice
{ 
  public:
    typedef unsigned int  uint;
    
  private:
    uint   D_; //dimension
    
    double fillCylinder(bool* regionA, double fracA);
    void   fill_KP_A(bool* regionA);
    void   fill_KP_B(bool* regionA);
    void   fill_KP_C(bool* regionA);
    void   fill_LW_a(bool* regionA);
    void   fill_LW_b(bool* regionA);
    void   fill_LW_c(bool* regionA);
    void   fill_LW_d(bool* regionA);
    void   fillRect2D(bool* regionA, uint xMin, uint xMax, uint yMin, uint yMax);
    void   init_N_and_neighbours();
    int    round(double num);
    void   trimWhiteSpace(std::string* word);
    uint   uintPower(uint base, uint exp);
    
  public:
    Hypercube(uint L, uint D);
    Hypercube(std::ifstream* fin, std::string fileName);
    virtual ~Hypercube();
    
    virtual uint getNeighbour(uint i, uint j);
    //virtual bool* getRegionA(double fracA);
    std::pair<std::string,bool*> getRegionA(std::string regionAInputStr);
    virtual void printParams();
    virtual void printNeighbours();
    
    //getter methods:
    uint getD();
};

#endif  // HYPERCUBE_H

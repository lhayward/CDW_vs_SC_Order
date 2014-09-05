/*********************************************************************************************
***************************************** CDW vs. SC *****************************************
**********************************************************************************************
* Lauren Hayward
**********************************************************************************************
* File:   VecND.h
*********************************************************************************************/

#ifndef VECND_H
#define VECND_H

#include "MersenneTwister.h"

using namespace std;

class VecND 
{
  private:
  public:
    typedef unsigned int  uint;
    
    static int numVecND;
    
    uint    N_;  //dimensionality of the vector
    double* v_;  //components of the ND vector
  
    VecND(uint N, VecND* oldVec); //copy constructor
    VecND(uint N, double val);
    VecND(uint N, int val);
    VecND(uint N, MTRand* randomGen);
    virtual ~VecND();
  
    void   add(VecND* vec2);
    void   clear();
    double dot(VecND* vec2);
    double dotForRange(VecND* vec2, uint start, uint end);
    VecND* getAbsComponents();
    VecND* getMultiple(double c);
    VecND* getReflection(VecND* r);
    double getSquare();
    double getSquareForRange(uint start, uint end);
    VecND* getSqComponents();
    void   multiply(double c);
    void   normalize();
    void   print();
    void   reflectAndNormalize(VecND* r);
    void   subtract(VecND* vec2);
};

#endif  /* VECND_H */


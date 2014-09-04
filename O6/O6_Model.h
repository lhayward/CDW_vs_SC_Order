/**********************************************************************************************
******************************** CLASSICAL REPLICA MONTE CODE *********************************
***********************************************************************************************
* Lauren Hayward
***********************************************************************************************
* File:   O6_Model.h
**********************************************************************************************/

#ifndef O6_MODEL
#define O6_MODEL

#include <string>
#include "Hypercube.h"
//#include "IsingSpins.h"
#include "Model.h"

class O6_Model : public Model
{ 
  public:
    typedef unsigned int  uint;
    
  private:
    double  lambda_;  //helicity modulus for spatial variations of the CDW order
    double  g_;       //anisotropy in the energy between the CDW and SC directions
    double  gPrime_;  //quartic anisotropy
    double  w_;       //controls the symmetry of the CDW order (stripe for w<0, checkerboard
                      //when w>0)
    uint    D_;       //dimension
    uint    L_;       //hypercube linear length
    uint    N_;       //number of spins living on the hypercube   
    
    Hypercube*  hcube_; //the hypercubic lattice on which the d.o.f. live
    //IsingSpins* spins_; //the degrees of freedom (d.o.f.) for the model
    
    void   localUpdate       (MTRand* randomGen);
    uint   uintPower(uint base, uint exp);
    
  public:
    O6_Model(std::ifstream* fin, std::string outFileName, Lattice* lattice);
    virtual ~O6_Model();
    
    virtual void makeMeasurement();
    virtual void printParams    ();
    virtual void printSpins     ();
    virtual void randomize      (MTRand* randomGen);
    virtual void setT           (double newT);
    virtual void sweep          (MTRand* randomGen);
    virtual void updateEnergy   ();
    virtual void writeBin       (int binNum, int numMeas);
};  

#endif  // O6_MODEL

/**********************************************************************************************
************************************ O(6) MODEL MONTE CODE ************************************
***********************************************************************************************
* Lauren E. Hayward Sierens
***********************************************************************************************
* File:   O6_Model.h
**********************************************************************************************/

#ifndef O6_MODEL
#define O6_MODEL

#include <string>
#include "Hypercube.h"
#include "Model.h"
#include "Vector_NDim.h"
#include "VectorSpins.h"

class O6_Model : public Model
{ 
  public:
    typedef unsigned int  uint;
    
  private:
    static const int VECTOR_SPIN_DIM = 6; //dimensionality of the spins on the lattice
    
    //parameters in the O(6) model:
    double  lambda_;  //helicity modulus for spatial variations of the CDW order
    double  g_;       //anisotropy in the energy between the CDW and SC directions
    double  gPrime_;  //quartic anisotropy
    double  w_;       //controls the symmetry of the CDW order (stripe for w<0, checkerboard
                      //when w>0)
                      
    //parameters of the lattice:
    uint    D_;       //dimension
    uint    N_;       //number of spins living on the hypercube
    uint*   L_;       //lattice's linear length in each dimension
    
    Hypercube*   hcube_; //the hypercubic lattice on which the d.o.f. live
    VectorSpins* spins_; //the degrees of freedom (d.o.f.) for the model
    Vector_NDim* mag_;   //total magnetization of the spins
    
    void   localUpdate(MTRand* randomGen);
    uint   uintPower  (uint base, uint exp);
    
  public:
    O6_Model(std::ifstream* fin, std::string outFileName, Lattice* lattice, MTRand* randomGen);
    virtual ~O6_Model();
    
    virtual void makeMeasurement    ();
    virtual void printParams        ();
    virtual void printSpins         ();
    virtual void randomizeLattice   (MTRand* randomGen);
    virtual void setT               (double newT);
    virtual void sweep              (MTRand* randomGen);
    virtual void updateEnergy       ();
    virtual void updateMagnetization();
    virtual void writeBin           (int binNum, int numMeas);
};  

#endif  // O6_MODEL

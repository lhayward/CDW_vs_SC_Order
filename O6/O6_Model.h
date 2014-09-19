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
#include "Hyperrectangle.h"
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
    double  Vz_;      //interlayer coupling for CDW order (only applicable when D_=3)
    double  VzPrime_; //interlayer coupling for SC order (only applicable when D_=3)
                      
    //parameters of the lattice:
    uint    D_;       //dimension
    uint    N_;       //number of spins living on the hyperrectangle
    uint*   L_;       //lattice's linear length in each dimension
    
    Hyperrectangle* hrect_; //the hyperrectangular lattice on which the d.o.f. live
    VectorSpins*    spins_; //the degrees of freedom (d.o.f.) for the model
    
    double       getHelicityModulus(int dir);
    double       getIsingOrder     ();
    double       getEnergy         ();
    Vector_NDim* getMagnetization  ();
    uint         uintPower         (uint base, uint exp);
    void         updateObservables ();
    
  public:
    O6_Model(std::ifstream* fin, std::string outFileName, Lattice* lattice, MTRand* randomGen);
    virtual ~O6_Model();
    
    virtual void localUpdate        (MTRand* randomGen);
    virtual void makeMeasurement    ();
    virtual void printParams        ();
    virtual void printSpins         ();
    virtual void randomizeLattice   (MTRand* randomGen);
    virtual void setT               (double newT);
    virtual void sweep              (MTRand* randomGen);
    virtual void writeBin           (int binNum, int numMeas);
};  

#endif  // O6_MODEL

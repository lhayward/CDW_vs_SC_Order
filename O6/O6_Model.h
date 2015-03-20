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
    static const uint VECTOR_SPIN_DIM = 6; //dimensionality of the spins on the lattice
    
    //parameters in the O(6) model:
    double  lambda_;  //helicity modulus for spatial variations of the CDW order
    double  g_;       //anisotropy in the energy between the CDW and SC directions
    double  gPrime_;  //quartic anisotropy
    double  w_;       //controls the symmetry of the CDW order (stripe for w<0, checkerboard
                      //when w>0)
    double r_;        //coefficient of the additional quartic term
    double  Vz_;      //interlayer coupling for CDW order (only applicable when D_=3)
    double  VzPrime_; //interlayer coupling for SC order (only applicable when D_=3)
    double  sigma_;   //related to the standard deviation of the disorder term 
                      //(std. dev. = sqrt(2)*sigma)
                      
    //parameters of the lattice:
    uint    D_;       //dimension
    uint    N_;       //number of spins living on the hyperrectangle
    uint    Nxy_;     //number of spins in each layer along the xy-plane
    uint    Lz_;      //system length along the z direction (equals 1 when D=2)
    uint*   L_;       //lattice's linear length in each dimension
    
    Hyperrectangle* hrect_; //the hyperrectangular lattice on which the d.o.f. live
    VectorSpins*    spins_; //the degrees of freedom (d.o.f.) for the model
    
    double** h_; //random field (for the disorder term in the Hamiltonian)
    
    bool*           inCluster_;  //boolean array for wolffUpdate indicating whether or not each
                                 //spin is in the cluster (redundant information to cluster 
                                 //vec, but stored for efficiency purposes)
                                 
    //acceptance rates for the two kinds of updates:
    int numAccept_local_;
    int numAccept_clust_;
    
    //variables related to storing the cluster sizes to generate a histogram:
    static const bool writeClusts_ = false;  //should we write data for the cluster histograms?
    uint*             clustSizes_;
    uint*             clustSizes_accepted_;
    uint*             clustSizes_rejected_;
    
    void         clearCluster          (std::vector<uint>* cluster);
    void         flipCluster           (std::vector<uint>* cluster, Vector_NDim* r);
    double       getClusterOnSiteEnergy(std::vector<uint>* cluster);
    double       getEnergy             ();
    double       getHelicityModulus    (int dir);
    double       getIsingOrder         ();
    Vector_NDim* getMagnetization      ();
    uint         uintPower             (uint base, uint exp);
    void         updateObservables     ();
    void         wolffUpdate           (MTRand &randomGen, uint start, uint end, bool pr);
    
  public:
    O6_Model(std::ifstream* fin, std::string outFileName, Lattice* lattice, MTRand &randomGen);
    virtual ~O6_Model();
    
    
    virtual void changeT            (double newT);
    virtual void localUpdate        (MTRand &randomGen);
    virtual void makeMeasurement    ();
    virtual void markWarmupDone     ();
    virtual void printParams        ();
    virtual void printSpins         ();
    virtual void randomizeLattice   (MTRand &randomGen);
    virtual void sweep              (MTRand &randomGen, bool pr);
    virtual void writeBin           (int binNum, int numMeas, int sweepsPerMeas);
    virtual void writeClustHistoData(std::string fileName);
    virtual void zeroMeasurements   ();
};  

#endif  // O6_MODEL

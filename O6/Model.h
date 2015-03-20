/**********************************************************************************************
************************************ O(6) MODEL MONTE CODE ************************************
***********************************************************************************************
* Lauren E. Hayward Sierens
***********************************************************************************************
* File:   Model.h (Abstract Class)
**********************************************************************************************/

#ifndef MODEL_H
#define MODEL_H

#include <stdlib.h>
#include <string>
#include "Measure.h"
#include "MersenneTwister.h"

class Model 
{ 
  public:
    typedef unsigned int  uint;
    
  protected:
    double        J_;       //coupling 
    double        T_;       //current temperature
    Measure       measures; //observables to record
    bool          warmupDone; //has the system finished warming up?
    std::ofstream fout;
    
  public:
    Model(std::ifstream* fin, std::string outFileName);
    virtual ~Model();
    
    //methods implemented in Model class:
    //void zeroMeasurements();
    
    //methods that can be overwritten by child classes:
    virtual void markWarmupDone  ();
    virtual void printParams     ();
    virtual void changeT         (double newT);
    virtual void zeroMeasurements();
    
    //pure virtual methods (to be implemented by all child classes):
    virtual void localUpdate        (MTRand &randomGen) = 0;
    virtual void makeMeasurement    () = 0;
    virtual void printSpins         () = 0;
    virtual void randomizeLattice   (MTRand &randomGen) = 0;
    virtual void sweep              (MTRand &randomGen, bool pr) = 0;
    virtual void writeBin           (int binNum, int numMeas, int sweepsPerMeas) = 0;
    virtual void writeClustHistoData(std::string fileName) = 0;
};

#endif  // MODEL_H

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
    double        energy_;  //current energy
    Measure       measures; //observables to record
    std::ofstream fout;
    
  public:
    Model(std::ifstream* fin, std::string outFileName);
    virtual ~Model();
    
    //methods implemented in Model class:
    double getEnergy();
    void   zeroMeasurements();
    
    //methods that can be overwritten by child classes:
    virtual void printParams();
    virtual void setT       (double newT);
    
    //pure virtual methods (to be implemented by all child classes):
    virtual void makeMeasurement    () = 0;
    virtual void printSpins         () = 0;
    virtual void randomizeLattice   (MTRand* randomGen) = 0;
    virtual void sweep              (MTRand* randomGen) = 0;
    virtual void updateEnergy       () = 0;
    virtual void updateMagnetization() = 0;
    virtual void writeBin           (int binNum, int numMeas) = 0;
};

#endif  // MODEL_H

/**********************************************************************************************
************************************ O(6) MODEL MONTE CODE ************************************
***********************************************************************************************
* Lauren E. Hayward Sierens
***********************************************************************************************
* File: Main.cpp
* Note: MersenneTwister.h was taken from external sources and is used as the random number
*       generator
**********************************************************************************************/

//#include <cmath>
#include <ctime>
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <typeinfo>
#include "FileReading.h"
#include "Hyperrectangle.h"
#include "MersenneTwister.h"
#include "Lattice.h"
#include "Model.h"
#include "O6_Model.h"
#include "SimParameters.h"

typedef unsigned long ulong;
typedef unsigned int  uint;

std::string getFileSuffix(int argc, char** argv);
Lattice* readLattice(std::string latticeName, std::string fileName, std::string startStr);
Model* readModel(std::string modelName, std::string inFileName, std::string startStr, 
                 std::string outFileName, Lattice* lattice, MTRand* randomGen);

/**********************************************************************************************
******************************************** main *********************************************
**********************************************************************************************/
int main(int argc, char** argv) 
{
  MTRand         randomGen;
  SimParameters* params;
  Lattice*       lattice;
  Model*         model;
  double         T; //current temperature
  time_t         sec1, sec2;  //for timing
  
  //variables related to input/output data from/to files:
  std::string fileSuffix      = getFileSuffix( argc, argv );
  std::string paramFileName   = "params" + fileSuffix + ".txt";
  std::string simParamStr     = "SIMULATION PARAMETERS";
  std::string latticeParamStr = "LATTICE PARAMETERS";
  std::string modelParamStr   = "MODEL PARAMETERS";
  std::string outFileName     = "bins" + fileSuffix + ".txt";
  
  std::cout.precision(8);
  std::cout << "\nParameter File: " << paramFileName << "\n" << std::endl;
  
  params = new SimParameters( paramFileName, simParamStr );
  params->print();
  
  lattice = readLattice( params->latticeType_, paramFileName, latticeParamStr );
  lattice->printParams();

  model = readModel( params->modelName_, paramFileName, modelParamStr, outFileName, lattice,
                     params->randomGen_ );
  model->printParams();
  
  std::cout << "\n*** STARTING SIMULATION ***\n" << std::endl;
  sec1 = time (NULL);
  //loop over the different temperatures:
  for( uint TIndex=0; TIndex<(params->TList_->size()); TIndex++)
  {
    T = (*(params->TList_))[TIndex];
    std::cout << "******** T = " << T << " (Temperature #" << (TIndex+1) << ") ********"
              << std::endl;
    model->setT(T);
    model->randomizeLattice( params->randomGen_ );
    
    //equilibrate:
    for( uint i=0; i<params->numWarmUpSweeps_; i++ )
    { model->sweep( params->randomGen_ ); }
    
    //loop over Monte Carlo bins:
    for( uint i=0; i<params->numBins_; i++ )
    {
      model->zeroMeasurements();
      //perform the measurements for one bin:
      for( uint j=0; j<params->measPerBin_; j++ )
      {
        //perform the sweeps for one measurement:
        for( uint k=0; k<params->sweepsPerMeas_; k++ )
        { model->sweep( params->randomGen_ ); }
        model->makeMeasurement();
      } //loop over measurements
      model->writeBin((i+1), params->measPerBin_);
      
      if( (i+1)%100==0 )
      { std::cout << (i+1) << " Bins Complete" << std::endl; }
    } //loop over bins
    std::cout << std::endl;
  } //temperature loop 
  
  sec2 = time(NULL);
  std::cout << "Time: " << (sec2 - sec1) << " seconds" << std::endl;
  std::cout << "\n*** END OF SIMULATION ***\n" << std::endl;
  return 0;
} //closes main

/**************************** getFileSuffix(int argc, char** argv) ***************************/
std::string getFileSuffix(int argc, char** argv)
{
  std::string result = "";
  
  if( argc > 1 )
  { result = "_" + std::string(argv[1]); }
  
  return result;
}

/****** readLattice(std::string latticeName, std::string fileName, std::string startStr) *****/
Lattice* readLattice(std::string latticeName, std::string fileName, std::string startStr)
{
  std::ifstream fin;
  fin.open(fileName.c_str());
  
  if( fin.is_open() )
  { FileReading::readUntilFound(&fin, startStr); }
  
  return new Hyperrectangle(&fin, fileName);
}

/**************** readModel(std::string modelName, std::string inFileName, ... ****************
*****************           std::string startStr, std::string outFileName,     ****************
*****************           Lattice* lattice, MTRand* randomGen)               ***************/
Model* readModel(std::string modelName, std::string inFileName, std::string startStr, 
                 std::string outFileName, Lattice* lattice, MTRand* randomGen )
{
  Model* result=NULL;
  std::ifstream fin;
  fin.open(inFileName.c_str());
  
  if( fin.is_open() )
  { FileReading::readUntilFound(&fin, startStr); }
  
  if( modelName == "o6" )
  { 
    result = new O6_Model(&fin, outFileName, lattice, randomGen); 
  }
  else
  {
    std::cout << "ERROR in readModel(std::string modelName, std::string inFileName, "
              << "std::string startStr, Lattice* lattice): the model name is not valid, so a "
              << "NULL Model object will be returned\n" << std::endl;
  }
  return result;
}
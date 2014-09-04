/**********************************************************************************************
************************************ O(6) MODEL MONTE CODE ************************************
***********************************************************************************************
* Lauren E. Hayward Sierens
***********************************************************************************************
* File:   O6_Model.cpp
**********************************************************************************************/

#include <fstream>
#include <iostream>
#include <string>
#include <typeinfo>
#include "FileReading.h"
#include "O6_Model.h"

//typdef needed because uint is a return types:
typedef O6_Model::uint uint;

/********** O6_Model(std::ifstream* fin, std::string outFileName, Lattice* lattice) ***********
**************************************** (constructor) ***************************************/
O6_Model::O6_Model(std::ifstream* fin, std::string outFileName, Lattice* lattice)
  : Model(fin, outFileName)
{
  const char EQUALS_CHAR = '=';
  
  std::cout.precision(15);
  
  if( fin!=NULL && fin->is_open() )
  {
    lambda_ = FileReading::readDouble(fin, EQUALS_CHAR);
    g_      = FileReading::readDouble(fin, EQUALS_CHAR);
    gPrime_ = FileReading::readDouble(fin, EQUALS_CHAR);
    w_      = FileReading::readDouble(fin, EQUALS_CHAR);
  
    if( lattice != NULL )
    {
      hcube_ = dynamic_cast<Hypercube *>(lattice);
      if(hcube_)
      {
        D_  = hcube_->getD();
        L_  = hcube_->getL();
        N_  = hcube_->getN();
      
        //spins_ = new IsingSpins(alpha_, N1_);
      }
      else
      {
        std::cout << "ERROR in O6_Model constructor:\n" 
                  << "  A lattice of type Hypercube is required.\n"
                  << "  A lattice of type " << typeid(*lattice).name() << " was given.\n" 
                  << std::endl;
      }
    }
    else
    {
      std::cout << "ERROR in O6_Model constructor: The passed Lattice object is not "
                << "valid\n" << std::endl; 
    }
  }
  else
  { 
    std::cout << "ERROR in Model constructor: could not read from input file\n" << std::endl; 
  }
}

/******************************* ~O6_Model() (destructor) *******************************/
O6_Model::~O6_Model()
{
  //delete the IsingSpins object:
  /*if( spins_ != NULL )
  { delete spins_; }
  spins_ = NULL;*/
}

/******************************* localUpdate(MTRand* randomGen) ******************************/
void O6_Model::localUpdate(MTRand* randomGen)
{
  uint latticeSite; //randomly selected spin location
  
  latticeSite = randomGen->randInt(N_-1);
}

/************************************* makeMeasurement() *************************************/
void O6_Model::makeMeasurement()
{
  double energyPerSpin = energy_/(1.0*N_);
  
  measures.accumulate( "E",   energyPerSpin ) ;
  measures.accumulate( "ESq", pow(energyPerSpin,2) );
}

/*************************************** printParams() ***************************************/
void O6_Model::printParams()
{
  std::cout << "O(6) Model Parameters:\n"
            << "---------------------------" << std::endl;
  Model::printParams();
  std::cout << "  lambda = " << lambda_ << "\n"
            << "       g = " << g_      << "\n"
            << "      g' = " << gPrime_ << "\n"
            << "       w = " << w_      << "\n"
            << std::endl;
}

/**************************************** printSpins() ***************************************/
void O6_Model::printSpins()
{ 
  //spins_->print(); 
}

/******************************** randomize(MTRand* randomGen) *******************************/
void O6_Model::randomize(MTRand* randomGen)
{
  //spins_->randomize(randomGen, regionA_);
}

/************************************* setT(double newT) *************************************/
void O6_Model::setT(double newT)
{ 
  Model::setT(newT); 
}

/********************************** sweep(MTRand* randomGen) *********************************/
void O6_Model::sweep(MTRand* randomGen)
{
  for( uint i=0; i<N_; i++ )
  { localUpdate(randomGen); }
}

/******************************** uintPower(int base, int exp) *******************************/
uint O6_Model::uintPower(uint base, uint exp)
{
  uint result = 1;
  for(uint i=1; i<=exp; i++)
  { result *= base; } 
  
  return result;
} //uintPower method

/*************************************** updateEnergy() **************************************/
void O6_Model::updateEnergy()
{
  energy_=0;
  energy_ *= -J_;
}

/***************************** writeBin(int binNum, int numMeas) *****************************/
void O6_Model::writeBin(int binNum, int numMeas)
{
  fout << L_ << '\t' << T_ << '\t' << binNum;
  measures.writeAverages(&fout, numMeas);
  fout << std::endl;
}
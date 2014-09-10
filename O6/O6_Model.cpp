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

/***************** O6_Model(std::ifstream* fin, std::string outFileName, ... ******************
******************          Lattice* lattice, MTRand* randomGen)             ******************
******************                       (constructor)                       *****************/
O6_Model::O6_Model(std::ifstream* fin, std::string outFileName, Lattice* lattice,
                   MTRand* randomGen)
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
      hrect_ = dynamic_cast<Hyperrectangle *>(lattice);
      if(hrect_)
      {
        D_  = hrect_->getD();
        
        if( D_==2 || D_==3 )
        {
          L_  = hrect_->getL();
          N_  = hrect_->getN();
      
          spins_ = new VectorSpins(N_, VECTOR_SPIN_DIM);
          randomizeLattice(randomGen);
          
          updateEnergy();
          mag_ = new Vector_NDim(VECTOR_SPIN_DIM,0);
          updateMagnetization();
        }
        else
        {
          std::cout << "ERROR in O6_Model constructor:\n" 
                    << "  The Hyperrectangle lattice must have dimension 2 or 3." 
                    << std::endl;
        }
      }
      else
      {
        std::cout << "ERROR in O6_Model constructor:\n" 
                  << "  A lattice of type Hyperrectangle is required.\n"
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
  //delete the VectorSpins object:
  if( spins_ != NULL )
  { delete spins_; }
  spins_ = NULL;
  
  if( mag_ != NULL )
  { delete mag_; }
  mag_ = NULL;
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
            << "---------------------" << std::endl;
  Model::printParams();
  std::cout << "  lambda = " << lambda_ << "\n"
            << "       g = " << g_      << "\n"
            << "      g' = " << gPrime_ << "\n"
            << "       w = " << w_      << "\n"
            << std::endl;
}

/**************************************** printSpins() ***************************************/
void O6_Model::printSpins()
{ spins_->print(); }

/**************************** randomizeLattice(MTRand* randomGen) ****************************/
void O6_Model::randomizeLattice(MTRand* randomGen)
{ spins_->randomize(randomGen); }

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
  Vector_NDim* currSpin;
  Vector_NDim* neighbour_x;
  Vector_NDim* neighbour_y;
  double       sumCDWSqs    = 0;
  double       energy1      = 0;
  double       energyLambda = 0;
  double       energyg      = 0;
  double       energygPrime = 0;
  double       energyw      = 0;
  
  energy_=0;
  
  for( uint i=0; i<N_; i++ )
  { 
    currSpin    = spins_->getSpin(i);
    neighbour_x = spins_->getSpin( hrect_->getNeighbour(i,0) );
    neighbour_y = spins_->getSpin( hrect_->getNeighbour(i,1) );
    
    energy1 += currSpin->dotForRange( neighbour_x, 0, 1 ) 
               + currSpin->dotForRange( neighbour_y, 0, 1 ); 
    
    energyLambda += currSpin->dotForRange( neighbour_x, 2, VECTOR_SPIN_DIM-1 ) 
                    + currSpin->dotForRange( neighbour_y, 2, VECTOR_SPIN_DIM-1 );
  }
  energy1 *= -1;
  energyLambda *= -1*lambda_;
  
  for( uint i=0; i<N_; i++ )
  { 
    sumCDWSqs = spins_->getSpin(i)->getSquareForRange(2,VECTOR_SPIN_DIM-1); 
    energyg += sumCDWSqs;
    energygPrime += pow(sumCDWSqs,2.0);
  }
  energyg *= (g_ + 4.0*(lambda_-1.0))/2.0;
  energygPrime *= gPrime_/2.0;
  
  for( uint i=0; i<N_; i++ )
  { energyw += pow( spins_->getSpin(i)->getSquareForRange(2,3), 2.0 ) 
               + pow( spins_->getSpin(i)->getSquareForRange(4,5), 2.0 ); }
  energyw *= w_/2.0;
  
  energy_ = J_*(energy1 + energyLambda + energyg + energygPrime + energyw);
}

/*********************************** updateMagnetization() ***********************************/
void O6_Model::updateMagnetization()
{
  mag_->clear();
  
  for( uint i=0; i<N_; i++ )
  { mag_->add( spins_->getSpin(i) ); }
}

/***************************** writeBin(int binNum, int numMeas) *****************************/
void O6_Model::writeBin(int binNum, int numMeas)
{
  fout << L_[0] << '\t' << T_ << '\t' << binNum;
  measures.writeAverages(&fout, numMeas);
  fout << std::endl;
}
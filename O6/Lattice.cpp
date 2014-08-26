/**********************************************************************************************
******************************** CLASSICAL REPLICA MONTE CODE *********************************
***********************************************************************************************
* Lauren Hayward
***********************************************************************************************
* File:   Lattice.cpp (Abstract Class)
**********************************************************************************************/

#include <iostream>
#include "FileReading.h"
#include "Lattice.h"

//typdef needed because uint is a return types:
typedef Lattice::uint uint;

/******************************** Lattice(int L) (constructor) *******************************/
Lattice::Lattice(uint L)
{
  L_ = L;
  N_ = 0;
  neighbours_ = NULL;
  isValid_ = false; //must be a member of one of the child classes and have the neighbours
                    //array not NULL to be a valid Lattice object
}

/************** Lattice(std::ifstream* fin, std::string fileName) (constructor) **************/
Lattice::Lattice(std::ifstream* fin, std::string fileName)
{
  const char EQUALS_CHAR = '=';
  
  if( fin!=NULL && fin->is_open() )
  { L_ = FileReading::readUint(fin, EQUALS_CHAR); }
  else
  {
    std::cout << "ERROR in Lattice constructor: could not read from file \"" << fileName 
              << "\"\n" << std::endl;
    L_=0;
  }
  
  N_ = 0;
  neighbours_ = NULL;
  isValid_ = false; //must be a member of one of the child classes and have the neighbours
                    //array not NULL to be a valid Lattice object
}

/********************************** ~Lattice() (destructor) **********************************/
Lattice::~Lattice()
{
  if( neighbours_ != NULL )
  {
    for(uint i=0; i<N_; i++)
    { delete[] neighbours_[i]; }
    delete[] neighbours_; 
  }
  neighbours_=NULL;
}

/*********************************** Public Getter Methods: **********************************/
bool Lattice::isValid(){ return isValid_; }
uint Lattice::getL()   { return L_; }
uint Lattice::getN()   { return N_; }
uint Lattice::getZ()   { return z_; }
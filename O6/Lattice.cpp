/**********************************************************************************************
************************************ O(6) MODEL MONTE CODE ************************************
***********************************************************************************************
* Lauren E. Hayward Sierens
***********************************************************************************************
* File:   Lattice.cpp (Abstract Class)
**********************************************************************************************/

#include <iostream>
#include "FileReading.h"
#include "Lattice.h"

//typdef needed because uint is a return types:
typedef Lattice::uint uint;

/********************************** Lattice() (constructor) **********************************/
Lattice::Lattice()
{
  N_ = 0;
  neighbours_ = NULL;
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
uint Lattice::getN()   { return N_; }
uint Lattice::getZ()   { return z_; }
/**********************************************************************************************
************************************ O(6) MODEL MONTE CODE ************************************
***********************************************************************************************
* Lauren E. Hayward Sierens
***********************************************************************************************
* File:   Hypercube.cpp
**********************************************************************************************/

#include <iostream>
#include <sstream>
#include <stdlib.h>
#include "FileReading.h"
#include "Hypercube.h"

//typdef needed because uint is a return types:
typedef Hypercube::uint uint;

/************* Hypercube(std::ifstream* fin, std::string fileName) (constructor) *************/
Hypercube::Hypercube(std::ifstream* fin, std::string fileName)
  : Lattice()
{
  const char EQUALS_CHAR = '=';
  const char LIST_START_CHAR = '[';
  const char LIST_END_CHAR = ']';
  
  //read in D_ from the file:
  if( fin!=NULL && fin->is_open() )
  { 
    D_ = FileReading::readUint     (fin, EQUALS_CHAR); 
    L_ = FileReading::readUintArray(fin, D_, 1, EQUALS_CHAR, LIST_START_CHAR, LIST_END_CHAR);
    
    N_ = 1;
    for( uint i=0; i<D_; i++ )
    { N_ *= L_[i]; }
    
    z_ = 2*D_;
    initNeighbours();
  }
  else
  {
    std::cout << "ERROR in Hypercube constructor: could not read from file \"" << fileName 
              << "\"\n" << std::endl;
    D_=1;
  }
}

/********************************* ~Hypercube() (destructor) *********************************/
Hypercube::~Hypercube()
{
  //delete the neighbours_ array:
  for(uint i=0; i<N_; i++)
  { 
    if( neighbours_[i] != NULL )
    { delete[] neighbours_[i]; }
    neighbours_[i] = NULL; 
  }
  if( neighbours_ != NULL )
  { delete[] neighbours_; }
  neighbours_ = NULL;
} // ~Hypercube

/******************************** getNeighbour(uint i, uint j) *******************************/
uint Hypercube::getNeighbour(uint i, uint j)
{
  uint result = 0;
  
  if( (neighbours_ != NULL) && i<N_ && j<(2*D_) )
  { result = neighbours_[i][j]; }
  else
  { 
    std::cout << "ERROR in Hypercube::getNeighbour(uint i, uint j): NULL neighbours_ array or "
              << "index out of bounds" << std::endl; 
  }
  
  return result;
}

/************************************** initNeighbours() *************************************/
void Hypercube::initNeighbours()
{
  uint next;     //used to naively estimate neighbour (before boundary correction) 
  uint next_mod; //used to determine if we are on boundary

  //initialize the neighbours_ array (note periodic boundary conditions):
  neighbours_ = new uint*[N_];
  for( uint i=0; i<N_; i++ )
  { neighbours_[i] = new uint[2*D_]; }
  for( uint i=0; i<N_; i++ )
  { 
    for( uint j=0; j<D_; j++ ) 
    {
      next     = nextInDir(j);
      next_mod = nextInDir(j+1);
      
      neighbours_[i][j] = i + next;
      //fix at the boundaries:
      if( neighbours_[i][j]%next_mod < next )
      { neighbours_[i][j] -= next_mod; }
      
      //initialize the corresponding neighbour (note that this information is redundant for a
      //hypercube, but we include it since some Model classes won't assume a hypercubic 
      //lattice)
      neighbours_[ neighbours_[i][j] ] [j + D_] = i;
    } //j
  } //i
}

/************************************ nextInDir(uint dir) ************************************/
uint Hypercube::nextInDir(uint dir)
{
  uint result = 1;
  for(uint i=0; i<dir; i++)
  { result *= L_[i]; } 
  
  return result;
} //nextInDir method

/*************************************** printParams() ***************************************/
void Hypercube::printParams()
{
  std::cout << "Hypercube Parameters:\n"
            << "--------------------\n"
            << "                Dimension D = " << D_ << "\n";
  if( L_[0] == L_[1] )
  { std::cout << "     Lattice Length Lx = Ly = " << L_[0] << "\n"; }
  else
  {
    std::cout << "          Lattice Length Lx = " << L_[0] << "\n"
              << "          Lattice Length Ly = " << L_[1] << "\n";
  }
  
  //if D=3, also print the length along the z direction:
  if( D_==3 )
  { std::cout << "          Lattice Length Lz = " << L_[2] << "\n"; }
  
  std::cout << "  Number of Lattice Sites N = " << N_ << "\n"
            << std::endl;
  
} //printParams method

/****************************************** printNeighbours() ******************************************/
void Hypercube::printNeighbours()
{
  std::cout << "Hypercube Neighbours list:" << std::endl;

  //print the neighbours_ array:
  for( uint i=0; i<N_; i++ )
  {
    std::cout.width(4);
    std::cout << "    " << i << ": ";
    for( uint j=0; j<(2*D_); j++ )
    {
      std::cout.width(4);
      std::cout << neighbours_[i][j] << " ";
    } //j
    std::cout << std::endl;
  } //i 
  std::cout << std::endl;
} //printNeighbours method

/******************************************* round ******************************************/
int Hypercube::round(double num)
{
    int result;
    
    if( num < 0.0 )
    { result = (int)ceil(num - 0.5); }
    else
    { result = (int)floor(num + 0.5); }
    
    return result;
}

/***************************** trimWhiteSpace(std::string* word) *****************************/
/*void Hypercube::trimWhiteSpace(std::string* word)
{
  std::size_t index;
  
  //trim the leading whitespace:
  index = word->find_first_not_of(" \t\n");
  *word = word->substr(index);

  //trim the trailing whitespace:
  index = word->find_last_not_of(" \t\n");
  *word = word->substr(0,index+1);
}*/

/******************************** uintPower(int base, int exp) *******************************/
uint Hypercube::uintPower(uint base, uint exp)
{
  uint result = 1;
  for(uint i=1; i<=exp; i++)
  { result *= base; } 
  
  return result;
} //uintPower method

/*********************************** Public Getter Methods: **********************************/
uint  Hypercube::getD(){ return D_; }
uint* Hypercube::getL(){ return L_; }

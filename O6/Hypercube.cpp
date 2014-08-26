/**********************************************************************************************
******************************** CLASSICAL REPLICA MONTE CODE *********************************
***********************************************************************************************
* Lauren Hayward
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

/*************************** Hypercube(int L, int D) (constructor) ***************************/
Hypercube::Hypercube(uint L, uint D)
  : Lattice(L)
{
  isValid_=true;
  D_ = D;
  z_ = 2*D;
  init_N_and_neighbours();
}

/************* Hypercube(std::ifstream* fin, std::string fileName) (constructor) *************/
Hypercube::Hypercube(std::ifstream* fin, std::string fileName)
  : Lattice(fin, fileName)
{
  const char EQUALS_CHAR = '=';
  
  isValid_ = true;
  
  //read in D from the file:
  if( fin!=NULL && fin->is_open() )
  { D_ = FileReading::readUint(fin, EQUALS_CHAR); }
  else
  {
    std::cout << "ERROR in Hypercube constructor: could not read from file \"" << fileName 
              << "\"\n" << std::endl;
    D_=1;
    isValid_ = false;
  }
  
  z_ = 2*D_;
  init_N_and_neighbours();
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

/*********************** fillCylinder(bool* regionA, double inputFracA) **********************/
double Hypercube::fillCylinder(bool* regionA, double inputFracA)
{
  uint NInLastDim    = (uint)round(inputFracA*L_);  //number of spins from the L spins in the
                                                    //Dth dimension that will be in region A
  uint NInA          = NInLastDim*uintPower(L_,D_-1);  //total number of spins in region A
  double outputFracA = (double)NInLastDim/L_;
  
  for( uint i=0; i<NInA; i++ )
  { regionA[i] = 1; }
  
  return outputFracA;
}

/**********************************************************************************************
* NOTE: The Kitaev-Preskill regions A,B,C can be illustrated as:
* 
*                   L               
*         <------------------->     
*          ___________________      
*         |                   |     
*         |                   |     
*         |                   |     
*         |_________          |     
*         | A  |    |         |     
*         |____|  C |         |     
*         | B  |    |         |     
*         |____|____|_________|     
*         <--------->               
*             L/2                   
*                                   
**********************************************************************************************/

/********************************** fill_KP_A(bool* regionA) *********************************/
void Hypercube::fill_KP_A(bool* regionA)
{ fillRect2D(regionA, 0, L_/4, L_/4, L_/2); }

/********************************** fill_KP_B(bool* regionA) *********************************/
void Hypercube::fill_KP_B(bool* regionA)
{ fillRect2D(regionA, 0, L_/4, 0, L_/4); }

/********************************** fill_KP_C(bool* regionA) *********************************/
void Hypercube::fill_KP_C(bool* regionA)
{ fillRect2D(regionA, L_/4, L_/2, 0, L_/2); }

/**********************************************************************************************
* NOTE: The Levin-Wen regions a,b,c,d can be illustrated as:
* 
*                             L                           
*         <--------------------------------------->       
*          _______________________________________        
*         |                                       |       
*         |                                       |       
*         |                                       |       
*         |                                       |       
*         |                                       |       
*         |                                       |       
*         |                                       |       
*         |     ______________                    |       
*         |    |      b       |                   |       
*         |    |______________|                   |       
*         |    | c  |    | d  |                   |       
*         |    |____|____|____|                   |       
*         |    |      a       |                   |       
*         |    |______________|                   |       
*         |                                       |       
*         |_______________________________________|       
*              <-------------->                           
*                    3L/8                                 
*                                                         
**********************************************************************************************/

/********************************** fill_LW_a(bool* regionA) *********************************/
void Hypercube::fill_LW_a(bool* regionA)
{ fillRect2D(regionA, L_/8, L_/2, L_/8, L_/4); }

/********************************** fill_LW_b(bool* regionA) *********************************/
void Hypercube::fill_LW_b(bool* regionA)
{ fillRect2D(regionA, L_/8, L_/2, 3*L_/8, L_/2); }

/********************************** fill_LW_c(bool* regionA) *********************************/
void Hypercube::fill_LW_c(bool* regionA)
{ fillRect2D(regionA, L_/8, L_/4, L_/4, 3*L_/8); }

/********************************** fill_LW_d(bool* regionA) *********************************/
void Hypercube::fill_LW_d(bool* regionA)
{ fillRect2D(regionA, 3*L_/8, L_/2, L_/4, 3*L_/8); }

/*********** fillRect2D(bool* regionA, uint xMin, uint xMax, uint yMin, uint yMax) ***********/
void Hypercube::fillRect2D(bool* regionA, uint xMin, uint xMax, uint yMin, uint yMax)
{
  for( uint y=yMin; y<yMax; y++ )
  {
    for( uint x=xMin; x<xMax; x++ )
    { regionA[y*L_ + x] = 1; }
  }
}

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

/********************************** getRegionA(double fracA) *********************************/
/*bool* Hypercube::getRegionA(double fracA)
{ 
  bool* regionA = new bool[N_];
  
  uint NInLastDim = (uint)round(fracA*L_);  //number of spins from the L spins in the Dth
                                            //dimension that will be in region A
  uint NInA       = NInLastDim*uintPower(L_,D_-1);  //total number of spins in region A
  
  for( uint i=0; i<NInA; i++ )
  { regionA[i] = 1; }
  for( uint i=NInA; i<N_; i++ )
  { regionA[i] = 0; }
  
  return regionA;
}*/

/************************** getRegionA(std::string regionAInputStr) **************************/
std::pair<std::string,bool*> Hypercube::getRegionA(std::string regionAInputStr)
{
  bool*             regionA = new bool[N_];
  double            inputFracA; //used when the region type is a cylinder
  double            outputFracA;
  std::size_t       commaIndex;
  std::stringstream ssout;
  std::string       regionType_1;
  std::string       regionType_2;  //more specific than regionType_1
  std::string       regionAOutputString;
  
  for( uint i=0; i<N_; i++ )
  { regionA[i] = 0; }
  
  commaIndex = regionAInputStr.find_last_of(",");
  regionType_1 = regionAInputStr.substr( 0, std::min(commaIndex,regionAInputStr.length()) );
  
  if( regionType_1 == "cylinder" && commaIndex<(regionAInputStr.length()-1) )
  {
    inputFracA = strtod( (regionAInputStr.substr(commaIndex+1)).c_str(), NULL);
    
    if( inputFracA > 0 && inputFracA <= 1 )
    {
      outputFracA = fillCylinder(regionA, inputFracA);
      ssout << "cylinder with fraction " << outputFracA << " of the spins";
      regionAOutputString = ssout.str();
    }
    else
    { regionAOutputString = "none"; }
  } // if for cylinder case
  
  //the Kitaev-Preskill regions {A,B,C,AB,AC,BC,ABC} can only be implemented on a D=2 lattice 
  //with L a multiple of 4:
  else if(regionType_1 == "KP" && commaIndex<(regionAInputStr.length()-1) && D_==2 && L_%4==0)
  {
    regionType_2 = regionAInputStr.substr( commaIndex + 1 );
    trimWhiteSpace(&regionType_2);
    regionAOutputString = "Kitaev-Preskill Region ";
    
    //Check for which KP regions {A,B,C} we want to include to make the complete region:
    if( (regionType_2.find("A") < regionType_2.length()) )
    { 
      fill_KP_A(regionA);
      regionAOutputString.append("A");
    }
    if( (regionType_2.find("B") < regionType_2.length()) )
    { 
      fill_KP_B(regionA);
      regionAOutputString.append("B");
    }
    if( (regionType_2.find("C") < regionType_2.length()) )
    { 
      fill_KP_C(regionA);
      regionAOutputString.append("C");
    }
    
    //if none of the regions {A,B,C} were specified:
    if( regionAOutputString == "Kitaev-Preskill Region " )
    { regionAOutputString = "none"; }
  } //Kitaev-Preskill regions
  
  //the Levin-Wen regions {A1,A2,A3,A4} can only be implemented on a D=2 lattice with L a 
  //multiple of 8:
  else if(regionType_1 == "LW" && commaIndex<(regionAInputStr.length()-1) && D_==2 && L_%8==0)
  {
    regionType_2 = regionAInputStr.substr( commaIndex + 1 );
    trimWhiteSpace(&regionType_2);
    
    if( regionType_2 == "A1" )
    {
      fill_LW_a(regionA);
      fill_LW_b(regionA);
      fill_LW_c(regionA);
      fill_LW_d(regionA);
      regionAOutputString = "Levin-Wen Region A1 (donut)";
    }
    else if( regionType_2 == "A2" )
    {
      fill_LW_a(regionA);
      fill_LW_b(regionA);
      fill_LW_c(regionA);
      regionAOutputString = "Levin-Wen Region A2 (C shape)";
    }
    else if( regionType_2 == "A3" )
    {
      fill_LW_a(regionA);
      fill_LW_b(regionA);
      fill_LW_d(regionA);
      regionAOutputString = "Levin-Wen Region A3 (backwards C shape)";
    }
    else if( regionType_2 == "A4" )
    {
      fill_LW_a(regionA);
      fill_LW_b(regionA);
      regionAOutputString = "Levin-Wen Region A4 (bars)";
    }
    else
    { regionAOutputString = "none"; }
  } //Levin-Wen regions
  
  else
  { regionAOutputString = "none"; }
  
  return std::pair<std::string,bool*>(regionAOutputString, regionA);
}

/********************************** init_N_and_neighbours() **********************************/
void Hypercube::init_N_and_neighbours()
{
  if( isValid_ )
  {
    N_ = uintPower(L_,D_);
  
    //initialize the neighbours_ array (note periodic boundary conditions):
    neighbours_ = new uint*[N_];
    for( uint i=0; i<N_; i++ )
    { neighbours_[i] = new uint[2*D_]; }
    for( uint i=0; i<N_; i++ )
    { 
      for( uint j=0; j<D_; j++ ) 
      {
        neighbours_[i][j] = i + uintPower(L_,j);
        //fix at the boundaries:
        if( neighbours_[i][j]%uintPower(L_,(j+1)) < uintPower(L_,j) )
        { neighbours_[i][j] -= uintPower(L_,(j+1)); }
        
        //initialize the corresponding neighbour (note that this information is redundant for a
        //hypercube, but we include it since some Model classes won't assume a hypercubic 
        //lattice)
        neighbours_[ neighbours_[i][j] ] [j + D_] = i;
      } //j
    } //i
  }
  else
  {
    std::cout << "ERROR in Hypercube::init_N_and_neighbours(): the Hypercube object is not "
              << "valid\n" << std::endl;
  }
}

/*************************************** printParams() ***************************************/
void Hypercube::printParams()
{
  if( isValid_ )
  {
    std::cout << "Hypercube Parameters:\n"
              << "--------------------\n"
              << "           Lattice Length L = " << L_ << "\n"
              << "                Dimension D = " << D_ << "\n"
              << "  Number of Lattice Sites N = " << N_ << "\n"
              << std::endl;
  }
  else
  {
    std::cout << "ERROR in Hypercube::printParams(): the Hypercube object is not valid\n" 
             << std::endl;
  }
  
} //printParams method

/****************************************** printNeighbours() ******************************************/
void Hypercube::printNeighbours()
{
  if( isValid_ )
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
  }
  else
  {
    std::cout << "ERROR in Hypercube::printNeighbours(): the Hypercube object is not valid\n" 
              << std::endl;
  }
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
void Hypercube::trimWhiteSpace(std::string* word)
{
  std::size_t index;
  
  //trim the leading whitespace:
  index = word->find_first_not_of(" \t\n");
  *word = word->substr(index);

  //trim the trailing whitespace:
  index = word->find_last_not_of(" \t\n");
  *word = word->substr(0,index+1);
}

/******************************** uintPower(int base, int exp) *******************************/
uint Hypercube::uintPower(uint base, uint exp)
{
  uint result = 1;
  for(uint i=1; i<=exp; i++)
  { result *= base; } 
  
  return result;
} //uintPower method

/*********************************** Public Getter Methods: **********************************/
uint Hypercube::getD(){ return D_; }

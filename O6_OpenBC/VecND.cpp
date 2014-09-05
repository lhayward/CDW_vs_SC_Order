/*********************************************************************************************
***************************************** CDW vs. SC *****************************************
**********************************************************************************************
* Lauren Hayward
**********************************************************************************************
* File:   VecND.cpp
*********************************************************************************************/
#include <cmath>
#include <cstdlib>
#include <iostream>
#include "VecND.h"

using namespace std;

/************************************ VecND (constructor) ***********************************/
VecND::VecND(uint N, VecND* oldVec)
{
  N_ = N;
  v_ = new double[N];
  
  for( uint i=0; i<N_; i++)
  { v_[i] = oldVec->v_[i]; }
  numVecND++;
}

/************************************ VecND (constructor) ***********************************/  
VecND::VecND(uint N, double val) 
{
  N_ = N;
  v_ = new double[N];
  
  for( uint i=0; i<N_; i++ )
  { v_[i] = val; }
  numVecND++;
}

/************************************ VecND (constructor) ***********************************/  
VecND::VecND(uint N, int val) 
{
  N_ = N;
  v_ = new double[N];
  
  for( uint i=0; i<N_; i++ )
  { v_[i] = val; }
  numVecND++;
}

/************************************ VecND (constructor) *************************************
* This is a method for generating a random point on an N-dimensional unit hypersphere. For 
* efficiency reasons, different techniques are used for different values of N:
* 
* For N<=3:
* The method generates a random point insides an N-dimensional hypercube of volumes 2^N. If 
* this point is outside the inscribed unit hypersphere, then is it is rejected. Otherwise, the 
* point is accepted and it is projected onto the surface of the unit hypersphere.
*
* For N>=4:
* The method generates N random coordinates according to a normal distribution with mean 0 and 
* variance 1. Each coordinate is then divided by the norm of the resulting vector so that the 
* final point is on the surface of the N-dimensional unit hypersphere.
*
* Note: It is also possible to use the method proposed by Marsaglia in 1972 for the N=4 case.
**********************************************************************************************/

VecND::VecND(uint N, MTRand* randomGen)
{
  N_ = N;
  v_ = new double[N];
  double S=0;
  
  //hypercube method is more efficient for N_<=3
  if(N_<=3)
  {
    for( uint i=0; i<N_; i++ )
    {
      v_[i] = -1.0 + 2*(randomGen->randDblExc());
      S += pow(v_[i],2);
    } //for
  
    while( S >= 1 )
    {
      S=0;
      for( uint i=0; i<N_; i++ )
      {
        v_[i] = -1.0 + 2*(randomGen->randDblExc());
        S += pow(v_[i],2);
      } //for
    } //while
  } //if
  
  //normal distribution method is more efficient for N_>=4:
  else
  {
    for( uint i=0; i<N_; i++ )
    {
      v_[i] = randomGen->randNorm(0,1);
      S += pow(v_[i],2);
    }
  } //else
  
  //normalize:
  S = pow(S,0.5);
  
  for( uint i=0; i<N_; i++ )
  { v_[i] = v_[i]/S; }
  numVecND++;
}

/************************************ ~VecND (destructor) ***********************************/ 
VecND::~VecND()
{
  delete[] v_;
  numVecND--;
}

/******************************************** add *******************************************/ 
void VecND::add(VecND* vec2)
{
  for( uint i=0; i<N_; i++ )
  { v_[i] += vec2->v_[i]; }
}

/******************************************* clear *******************************************/ 
void VecND::clear()
{
  for( uint i=0; i<N_; i++ )
  { v_[i] = 0; }
}

/******************************************** dot *******************************************/
double VecND::dot(VecND* vec2)
{
  double result=0;
  
  for( uint i=0; i<N_; i++ )
  { result += v_[i]*vec2->v_[i]; }
  
  return result;
}

/**************************************** dotForRange ****************************************/
double VecND::dotForRange(VecND* vec2, uint start, uint end)
{
  double result=0;
  uint index1=start;
  uint index2=end;
  
  if(index1 >= N_)
  { index1 = (N_-1); }
  if(index2 >= N_)
  { index2 = (N_-1); }
  
  for( uint i=index1; i<=index2; i++ )
  { result += v_[i]*vec2->v_[i]; } 
  
  return result;
}

/**************************************** getMultiple ****************************************/
VecND* VecND::getMultiple(double c)
{
  VecND* result = new VecND(N_,0);
  
  for( uint i=0; i<N_; i++ )
  { result->v_[i] = c*v_[i]; }
  
  return result;
}

/*************************************** getReflection ***************************************/
VecND* VecND::getReflection(VecND* r)
{
  VecND* result = new VecND(N_, this);
  result->reflectAndNormalize(r);
  return result;
}

/***************************************** getSquare ****************************************/
double VecND::getSquare()
{
  double result = 0;
  
  for( uint i=0; i<N_; i++ )
  { result += pow(v_[i],2); }
  
  return result;
}

/************************************* getSquareForRange *************************************/
double VecND::getSquareForRange(uint start, uint end)
{
  double result = 0;
  uint index1=start;
  uint index2=end;
  
  if(index1 >= N_)
  { index1 = (N_-1); }
  if(index2 >= N_)
  { index2 = (N_-1); }
  
  for( uint i=index1; i<=index2; i++ )
  { result += pow(v_[i],2); }
  
  return result;
}

/************************************** getSqComponents **************************************/
VecND* VecND::getSqComponents()
{
  VecND* result = new VecND(N_,0);
  
  for( uint i=0; i<N_; i++ )
  { result->v_[i] = pow(v_[i],2); }
  
  return result;
}

/***************************************** multiply *****************************************/
void VecND::multiply(double c)
{
  for( uint i=0; i<N_; i++ )
  { v_[i] = c*v_[i]; }
}

/***************************************** normalize *****************************************/
void VecND::normalize()
{
  double len = pow(getSquare(),0.5);
  
  for( uint i=0; i<N_; i++ )
  { v_[i] = v_[i]/len; }
}

/******************************************* print ******************************************/
void VecND::print()
{
  cout << "[ "; 
  for( uint i=0; i<(N_-1); i++ )
  { cout << v_[i] << " , "; }
  cout << v_[N_-1] << " ]" << endl;
}

/************************************ reflectAndNormalize ************************************/
void VecND::reflectAndNormalize(VecND* r)
{
  double dotProd = dot(r); 
  
  for( uint i=0; i<N_; i++ )
  { v_[i] = v_[i] - 2*dotProd*r->v_[i]; }
  
  //ensure that the norm is still one (to avoid rounding errors building up):
  normalize();
}

/***************************************** subtract *****************************************/
void VecND::subtract(VecND* vec2)
{
  for( uint i=0; i<N_; i++ )
  { v_[i] -= vec2->v_[i]; }
}

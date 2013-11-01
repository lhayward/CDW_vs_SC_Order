/*********************************************************************************************
***************************************** CDW vs. SC *****************************************
**********************************************************************************************
* Lauren Hayward
**********************************************************************************************
* File:   Simulation.h
*********************************************************************************************/
#ifndef SIMULATION_H
#define SIMULATION_H

#include <vector>
#include "MersenneTwister.h"
#include "VecND.h"

//using namespace std;

class Simulation 
{
public:
  typedef unsigned int  uint;
  
private:
  double               J,sigmaBar,T;
  double               energy;  //total energy
  VecND*               mag;     //total magnetization
  uint                 spinDim;
  int                  L,N;
  int                  maxZ;  //max. coordinate number
  int                  numWarmUpSweeps, sweepsPerMeas, measPerBin, numBins;
  VecND**              spins;
  int*                 coordNums;
  int**                crossProds;  //2D array, where crossProds[i][a] stores TWICE the cross 
                                    //product of r_i with bond vector "a" ( a=0 => \hat{x},
                                    //a=1 => -\hat{x}, a=2 => \hat{y}, a=3 => -\hat{y} )
  int**                neighbours;
  const char*          outputFileName;
  MTRand*              randomGen;
  //vector<int>*         cluster; //vector storing the sites of spins in the cluster
  //bool*                inCluster;  //boolean array indicating whether or not each spin
                                   //is in the cluster (redundant information to cluster vec,
                                   //but stored for efficiency purposes)
  //vector<int>*         buffer;
  vector<double>*      TList;
  
  void   calculateEnergy();
  void   calculateIsingOrder();
  void   calculateMagnetization();
  //void   clearCluster();
  //void   flipCluster(VecND* r);
  ////void   flipSpin(int site, Vec4D* r);  //for original Wolff cluster algorithm
  ////bool   isInCluster(int site);
  //double getClusterOnSiteEnergy();
  double getCorrelation(int i, int j);
  double getDiamag();
  double getHelicityModulus(int dir);
  double getSF();
  double getSumPsiSq();
  void   metropolisStep();
  //void   printCluster();
  void   printLattice();
  void   printNeighbours();
  void   randomizeLattice();
  void   setUpNeighbours();
  void   sweep();
  //void   wolffStep();
  
public:
  Simulation(double J, double sigmaBar, vector<double>* TList, int L,
             MTRand* randomGen, int numWarmUpSweeps, int sweepsPerMeas, int measPerBin, 
             int numBins, const char* outputFileName);
  virtual ~Simulation();
  void    runSim();  
};

#endif  /* SIMULATION_H */


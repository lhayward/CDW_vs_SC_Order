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
  double               isingOrderParam;
  VecND*               mag;     //total magnetization
  uint                 spinDim;
  int                  L,N;
  int                  numWarmUpSweeps, sweepsPerMeas, measPerBin, numBins;
  VecND**              spins;
  int**                neighbours;
  const char*          outputFileName;
  MTRand*              randomGen;
  vector<int>*         cluster; //vector storing the sites of spins in the cluster
  bool*                inCluster;  //boolean array indicating whether or not each spin
                                   //is in the cluster (redundant information to cluster vec,
                                   //but stored for efficiency purposes)
  vector<int>*         buffer;
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
  double getCPhi(int i, int j);
  double getHelicityModulus(int dir);
  double getSF();
  double getSFPhi();
  void   metropolisStep();
  //void   printCluster();
  void   printLattice();
  void   randomizeLattice();
  void   setUpNeighbours();
  void   sweep();
  //void   wolffStep();
  
  //methods for generating a random point on an N-dimensional unit hypersphere (proposed by 
  //Marsaglia for N=4):
  //VecND* getRandomVecND_3();
public:
  Simulation(double J, double sigmaBar, vector<double>* TList, int L,
             MTRand* randomGen, int numWarmUpSweeps, int sweepsPerMeas, int measPerBin, 
             int numBins, const char* outputFileName);
  virtual ~Simulation();
  void    runSim();  
};

#endif  /* SIMULATION_H */


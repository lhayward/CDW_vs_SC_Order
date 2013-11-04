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
  double               J,lambda,g,w,T;
  double               energy;  //total energy
  double               isingOrderParam;
  VecND*               mag;     //total magnetization
  uint                 spinDim, extraSpinLoc;
  int                  L,N;
  int                  maxZ;  //max. coordinate number
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
  void   clearCluster();
  void   flipCluster(VecND* r);
  //void   flipSpin(int site, Vec4D* r);  //for original Wolff cluster algorithm
  double getClusterOnSiteEnergy();
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
  void   wolffStep();
  
public:
  Simulation(double J, double lambda, double g, double w, vector<double>* TList, int L,
             MTRand* randomGen, int numWarmUpSweeps, int sweepsPerMeas, int measPerBin, 
             int numBins, const char* outputFileName);
  virtual ~Simulation();
  void    runSim();  
};

#endif  /* SIMULATION_H */


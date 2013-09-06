/*********************************************************************************************
***************************************** CDW vs. SC *****************************************
**********************************************************************************************
* Lauren Hayward
**********************************************************************************************
* File:   Simulation.cpp
*********************************************************************************************/

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include "MersenneTwister.h"
#include "Simulation.h"

//using namespace std;

/********************************* Simulation (constructor) *********************************/ 
Simulation::Simulation(double J, double lambda, double g, double w, vector<double>* TList,
                       int L, MTRand* randomGen, int numWarmUpSweeps, int sweepsPerMeas, 
                       int measPerBin, int numBins, const char* outputFileName)
{
  spinDim               = 6;
  this->J               = J;
  this->lambda          = lambda;
  this->g               = g;
  this->w               = w;
  this->TList           = TList;
  this->L               = L;
  this->N               = L*L;
  this->randomGen       = randomGen;
  this->numWarmUpSweeps = numWarmUpSweeps;
  this->sweepsPerMeas   = sweepsPerMeas;
  this->measPerBin      = measPerBin;
  this->numBins         = numBins;
  this->outputFileName  = outputFileName;
  
  spins = new VecND*[N];
  //for( int i=0; i<N; i++ )
  //{ spins[i] = getRandomVecND_1(); }
  randomizeLattice();
  
  neighbours = new int*[N];
  for( int i=0; i<N; i++ )
  { neighbours[i] = new int[4]; }
  
  setUpNeighbours();
  
  calculateEnergy();
  mag = new VecND(spinDim,0);
  calculateMagnetization();
  
  cluster  = new vector<int>;
  buffer   = new vector<int>;
}

/********************************* ~Simulation (destructor) *********************************/
Simulation::~Simulation()
{
  for( int i=0; i<N; i++ )
  { delete spins[i]; }
  delete[] spins;
  
  for(int i=0; i<N; i++)
  { delete[] neighbours[i]; }
  delete[] neighbours; 
  
  if(mag!=NULL)
  { delete mag; }
  mag = NULL;
  
  delete cluster;
  delete buffer;
}

/****************************************** runSim ******************************************/ 
void Simulation::runSim()
{
  unsigned int      TIndex;
  ofstream          outFile;
  string            spinFileName;
  //double            currM;
  double            currPsiSq;
  VecND*            currn;
  VecND*            currnSq;
  double            aveE, aveESq;
  double            aveHelicityX, aveHelicityY;
  double            aveM, aveMAbs, aveMSq, aveM4;
  double            avePsiSq, avePsi4;
  //double            aveSF, aveSFPhi;
  VecND* aven   = new VecND(spinDim,0);
  VecND* avenSq = new VecND(spinDim,0);
  
  
  outFile.open(outputFileName);
  outFile.precision(20);
  
  //loop over all temperatures:
  for(TIndex=0; TIndex<(TList->size()); TIndex++)
  {
    //T = TMax - TIndex*(TMax-TMin)/( (double)max(1,(numTSteps-1)) );
    T = TList->at(TIndex);
    std::cout << "\n******** L = " << L << ", T = " << T << " (Temperature #" << (TIndex+1) << ") ********" << std::endl;
    
    //take out the following two lines if you want the system to use its previous state from 
    //the last temperature (i.e. cooling the system):
    randomizeLattice();
    
    //do warm-up sweeps:
    for( int i=0; i<numWarmUpSweeps; i++ )
    { sweep(); }
    
    for( int i=0; i<numBins; i++ )
    {
      aveE=0;
      aveESq=0;
      aveHelicityX=0;
      aveHelicityY=0;
      aveM=0;
      aveMAbs=0;
      aveMSq=0;
      aveM4=0;
      avePsiSq=0;
      avePsi4=0;
      aven->clear();
      avenSq->clear();
      //aveSF = 0;
      //aveSFPhi = 0;
  
      for( int j=0; j<measPerBin; j++ )
      { 
        //perform one measurement:
        for( int k=0; k<sweepsPerMeas; k++ )
        { sweep(); }
        
        //put this section in to avoid rounding errors (rather than trusting the energy and
        //magnetization we have been updating ourselves):
        calculateEnergy();
        calculateMagnetization();
        calculateIsingOrder();
        
        currn = mag->getMultiple(1.0/N);
        currnSq = mag->getSqComponents();
        currnSq->multiply(1.0/(N*N));
        //currM = currnSq->v_[2] + currnSq->v_[3] - currnSq->v_[4] - currnSq->v_[5];
        currPsiSq = currnSq->v_[0] + currnSq->v_[1];
        
        aveE         += energy/N;
        aveESq       += pow( (energy/N), 2);
        //aveSF      += getSF();
        //aveSFPhi   += getSFPhi();
        aveHelicityX += getHelicityModulus(0);
        aveHelicityY += getHelicityModulus(1);
        aveM         += isingOrderParam;
        aveMAbs      += abs(isingOrderParam);
        aveMSq       += pow(isingOrderParam, 2);
        aveM4        += pow(isingOrderParam, 4);
        avePsiSq     += currPsiSq;
        avePsi4      += pow(currPsiSq,2);
        aven->add(currn);
        avenSq->add(currnSq);
        
        if(currn!=NULL)
        { delete currn; }
        currn = NULL;
        
        if(currnSq!=NULL)
        { delete currnSq; }
        currnSq=NULL;
        
      } //j (measPerBin)
      //calculate average for the bin:
      aveE         = aveE/measPerBin;
      aveESq       = aveESq/measPerBin;
      aveHelicityX = aveHelicityX/measPerBin;
      aveHelicityY = aveHelicityY/measPerBin;
      //aveSF      = aveSF/measPerBin;
      //aveSFPhi   = aveSFPhi/measPerBin;
      aveM         = aveM/measPerBin;
      aveMAbs      = aveMAbs/measPerBin;
      aveMSq       = aveMSq/measPerBin;
      aveM4        = aveM4/measPerBin;
      avePsiSq     = avePsiSq/measPerBin;
      avePsi4      = avePsi4/measPerBin;
      aven->multiply(1.0/measPerBin);
      avenSq->multiply(1.0/measPerBin);
      
      outFile << L << '\t' << T << '\t' << (i+1) << '\t' << aveE << '\t' << aveESq << '\t'
              << aveHelicityX << '\t' << aveHelicityY << '\t'
              << aveM << '\t' << aveMAbs << '\t' << aveMSq << '\t' << aveM4 << '\t'
              << avePsiSq << '\t' << avePsi4 << '\t';
      for( uint j=0; j<spinDim; j++ )
      { outFile << aven->v_[j] << '\t'; }
      for( uint j=0; j<spinDim; j++ )
      { outFile << avenSq->v_[j] << '\t'; }
      outFile << std::endl;
      
      std::cout << (i+1) << " Bins Complete" << std::endl; 
    } //i (bins)
  }  //closes T loop
  
  outFile.close();
}

/************************************** calculateEnergy *************************************/  
void Simulation::calculateEnergy()
{
  int    i;
  energy=0;
  double energy1=0;
  double energyLambda=0;
  double energyg=0;
  double energyw=0;
  
  for( i=0; i<N; i++ )
  { energy1 += spins[i]->dotForRange(spins[neighbours[i][0]],0,1) 
               + spins[i]->dotForRange(spins[neighbours[i][2]],0,1); }
  energy1 *= -1;
  
  for( i=0; i<N; i++ )
  { energyLambda += spins[i]->dotForRange(spins[neighbours[i][0]],2,spinDim-1) 
                    + spins[i]->dotForRange(spins[neighbours[i][2]],2,spinDim-1); }
  energyLambda *= -1*lambda;
  
  for( i=0; i<N; i++ )
  { energyg += spins[i]->getSquareForRange(2,spinDim-1); }
  energyg *= (g + 4.0*(lambda-1.0))/2.0;
  
  for( i=0; i<N; i++ )
  { energyw += pow(spins[i]->getSquareForRange(2,3),2.0) 
               + pow(spins[i]->getSquareForRange(4,5),2.0); }
  energyw *= w/2.0;
  
  energy = J*(energy1 + energyLambda + energyg + energyw);
}

/************************************ calculateIsingOrder ************************************/ 
void Simulation::calculateIsingOrder()
{
  isingOrderParam=0;
  VecND* orderParamComponents = new VecND(spinDim,0);
  VecND* sqN;  //squared components of spins[i] for a given i
  
  for( int i=0; i<N; i++ )
  { 
    sqN = spins[i]->getSqComponents();
    sqN->multiply(1.0/(N*N));
    orderParamComponents->add(sqN); 
    
    if( sqN!=NULL )
    { delete sqN; }
    sqN=NULL;
  }

  isingOrderParam = orderParamComponents->v_[2] + orderParamComponents->v_[3]
                    - orderParamComponents->v_[4] - orderParamComponents->v_[5];
  delete orderParamComponents;
}

/********************************** calculateMagnetization **********************************/ 
void Simulation::calculateMagnetization()
{
  int i;
  mag->clear();
  
  for( i=0; i<N; i++ )
  { mag->add(spins[i]); }

}

/*************************************** clearCluster ***************************************/ 
void Simulation::clearCluster()
{
  while( !cluster->empty() )
  { cluster->pop_back(); }
}

/**************************************** flipCluster *****************************************
* This function flips all spins in the cluster by using the passed vector r.
**********************************************************************************************/
void Simulation::flipCluster(VecND* r)
{
  int clustSize = (int)cluster->size();
  
  for( int i=0; i<clustSize; i++ )
  { spins[cluster->at(i)]->reflectAndNormalize(r); } //for loop
}

/***************************************** flipSpin ******************************************
* The function reflects the spin at the passed site about the hyperplane orthogonal to the 
* passed vector r.
*********************************************************************************************/
/*void Simulation::flipSpin(int site, Vec4D* r)
{
  int numNeighbours = 4;
  double deltaE;
  Vec4D* nnSum = new Vec4D(spinDim,0);
  
  //loop to calculate the nearest neighbour sum:
  for( int i=0; i<numNeighbours; i++ )
  { nnSum->add(spins[neighbours[site][i]]); }
  
  deltaE = 2*J*(r->dot(spins[site]))*(r->dot(nnSum));
  mag->subtract(spins[site]);
  spins[site]->reflect(r);
  mag->add(spins[site]);
  
  energy += deltaE;
  
  cluster->push_back(site);
  buffer->push_back(site);
  
  //delete the vector storing the nearest neighbour sum:
  if(nnSum!=NULL)
  { delete nnSum; }
  nnSum = NULL;
}*/

/*********************************** getClusterOnSiteEnergy ***********************************
* This function calculates the onsite (non-interacting) part of the energy corresponding to 
* the spins in the cluster.
**********************************************************************************************/
double Simulation::getClusterOnSiteEnergy()
{
  int clustSize = (int)cluster->size();
  int site;
  double energyg=0;
  double energyw=0;
  
  for( int i=0; i<clustSize; i++ )
  {
    site = cluster->at(i);
    energyg += spins[site]->getSquareForRange(2,spinDim-1);
    energyw += pow(spins[site]->getSquareForRange(2,3),2.0) 
                 + pow(spins[site]->getSquareForRange(4,5),2.0);
  } //for loop
  
  energyg *= (g + 4.0*(lambda-1.0))/2.0;
  energyw *= w/2.0;
  return J*(energyg + energyw);
}

/************************************** getCorrelation ***************************************/
double Simulation::getCorrelation(int i, int j)
{
  return spins[i]->dot(spins[j]);
}

/****************************************** getCPhi ******************************************/
double Simulation::getCPhi(int i, int j)
{
  return spins[i]->dotForRange(spins[j],2,3);
}

/************************************* getHelicityModulus *************************************
* Calculates and returns the helicity modulus for the desired lattice directions.
* Note: dir=0 corresponds to the x-direction
*       dir=1 corresponds to the y-direction
**********************************************************************************************/
double Simulation::getHelicityModulus(int dir)
{
  double helicityMod = 0;
  int neighDir;  //direction (for neighbour's array) corresponding to the specified dir
  VecND* neigh; //current neighbour
  double sum1;  //first sum in formula for helicity modulus
  double sum2;  //second sum in formula for helicity modulus
  
  //compute the helicity modulus if a valid direction was passed (otherwise just return 0):
  if(dir==0 || dir==1)
  {
     neighDir = 2*dir;  //0 if x-direction (dir=0), 2 if y-direction (dir=1)
     sum1=0;
     sum2=0;
     
     //loop over all spins:
     for(int i=0; i<N; i++)
     {
      neigh = spins[neighbours[i][neighDir]];
      
      sum1 += spins[i]->dotForRange(neigh,0,1);
      sum2 += spins[i]->v_[0]*neigh->v_[1] - spins[i]->v_[1]*neigh->v_[0];
     }
     helicityMod = sum1/(1.0*N) - J/(T*N)*pow(sum2,2); 
  } //end if
  
  return helicityMod;
}

/************************************** getRandomVecND_3 **************************************
* This is a method for generating a random point on a 4D unit sphere. This method was proposed by
* Marsaglia in 1972. This method turns out to be the most efficient in 4D when compared to the
* other two methods above for generating a random point on the surface of a 4D sphere.
**********************************************************************************************/
/*VecND* Simulation::getRandomVecND_3()
{
  double x1, x2, x3, x4;
  double S1, S2;
  double C;  //factor we need to multiply by at the end
  
  x1 = -1.0 + 2*(randomGen->randDblExc());
  x2 = -1.0 + 2*(randomGen->randDblExc());
  S1 = pow(x1,2) + pow(x2,2);
  
  while( S1 >= 1 )
  {
    x1 = -1.0 + 2*(randomGen->randDblExc());
    x2 = -1.0 + 2*(randomGen->randDblExc());
    S1 = pow(x1,2) + pow(x2,2);
  }
  
  x3 = -1.0 + 2*(randomGen->randDblExc());
  x4 = -1.0 + 2*(randomGen->randDblExc());
  S2 = pow(x3,2) + pow(x4,2);
  
  while( S2 >= 1 )
  {
    x3 = -1.0 + 2*(randomGen->randDblExc());
    x4 = -1.0 + 2*(randomGen->randDblExc());
    S2 = pow(x3,2) + pow(x4,2);
  }
  
  C = pow(abs(1.0-S1)/S2,0.5);
  
  return new VecND(x1,x2,x3*C,x4*C);
}*/

/******************************************* getSF *******************************************/
double Simulation::getSF()
{
  double SF = 0;
  
  ////sum over all pairs of spins (ensure j<i for no double counting):
  //sum over all pairs of spins:
  for( int i=0; i<N; i++ )
  {
    for( int j=0; j<N; j++ )
    { SF += getCorrelation(i,j); }
  }
  SF /= N;
  
  return SF;
}

/****************************************** getSFPhi *****************************************/
double Simulation::getSFPhi()
{
  double SFPhi = 0;
  
  ////sum over all pairs of spins (ensure j<i for no double counting):
  //sum over all pairs of spins:
  for( int i=0; i<N; i++ )
  {
    for( int j=0; j<N; j++ )
    { SFPhi += getCPhi(i,j); }
  }
  SFPhi /= N;
  
  return SFPhi;
}

/**************************************** isInCluster ****************************************
* This function checks whether or not spin at the passed site is in the cluster.
*********************************************************************************************/  
bool Simulation::isInCluster(int site)
{
  int i;
  int clustSize = (int)cluster->size();
  bool found = false;
  
  i=0;
  while( !found && i<clustSize )
  {
    if( cluster->at(i)==site )
    { found = true; }
    i++;
  }  //closes while loop
  
  return found;
}

/************************************** metropolisStep *************************************/  
void Simulation::metropolisStep()
{
  int numNeighbours = 4;
  int site;
  double deltaE;
  VecND* sNew = new VecND(spinDim, randomGen);
  VecND* nnSum = new VecND(spinDim,0);
  
  site = randomGen->randInt(N-1);
  
  //loop to calculate the nearest neighbour sum:
  for( int i=0; i<numNeighbours; i++ )
  { nnSum->add(spins[neighbours[site][i]]); }
  
  deltaE = J*( -1*(nnSum->dotForRange(sNew,0,1) - nnSum->dotForRange(spins[site],0,1)) 
               - lambda*(nnSum->dotForRange(sNew,2,spinDim-1) - nnSum->dotForRange(spins[site],2,spinDim-1)) 
               + (g + 4*(lambda-1.0))/2.0*(sNew->getSquareForRange(2,spinDim-1) 
                                          - spins[site]->getSquareForRange(2,spinDim-1)) 
               + w/2.0*(pow(sNew->getSquareForRange(2,3),2.0) + pow(sNew->getSquareForRange(4,5),2.0) 
                        - pow(spins[site]->getSquareForRange(2,3),2.0) - pow(spins[site]->getSquareForRange(4,5),2.0)) );
  
  if( deltaE<=0 || randomGen->randDblExc() < exp(-deltaE/T) )
  {
    //Calculate energy and mag before writing to file, not now:
    //energy += deltaE;
    //mag->subtract(spins[site]);
    //mag->add(sNew);
    
    //delete the vector storing the old state of the spin:
    if(spins[site]!=NULL)
    { delete spins[site]; }
    spins[site] = NULL;
    
    spins[site] = sNew;
  }
  else
  {
    //delete the vector storing the rejected new state of the spin:
    if(sNew!=NULL)
    { delete sNew; }
    sNew = NULL;  
  }
  
  //delete the vector storing the nearest neighbour sum:
  if(nnSum!=NULL)
  { delete nnSum; }
  nnSum = NULL;
}

/*************************************** printCluster ***************************************/
/*void Simulation::printCluster()
{
  int clustSize = (int)cluster->size();
  
  for( int i=0; i<clustSize; i++ )
  {
    std::cout << cluster->at(i) << "   ";
  }
  std::cout << std::endl;
}*/

/**************************************** printLattice ***************************************/
void Simulation::printLattice()
{
  for( int i=0; i<N; i++ )
  {
    std::cout << "Spin " << i << ": ";
    spins[i]->print();
  }
  std::cout << std::endl;
}

/************************************** randomizeLattice *************************************/
void Simulation::randomizeLattice()
{
  for( int i=0; i<N; i++ )
  { spins[i] = new VecND(spinDim,randomGen); }
}

/************************************** setUpNeighbours **************************************
* This functions sets up the "neighbours" array for a square lattice with periodic boundary
* conditions. 
*********************************************************************************************/
void Simulation::setUpNeighbours()
{
  int i,x,y;
  
  //loop to assign the neighbours (with periodic boundary conditions):
  for( i=0; i<N; i++ )
  {
    //assign the x-direction neighbours:
    x = i+1;
    if( (i+1)%L==0 )
    { x = x - L; }
    neighbours[i][0] = x; //+x neighbour of i is x
    neighbours[x][1] = i; //-x neighbour of x is i
    
    //assign the y-direction neighbours:
    y = i+L;
        if( y%(L*L) < L )
        { y = y - L*L; }
    neighbours[i][2] = y;  //+y neighbour of i is y
    neighbours[y][3] = i;  //-y neighbour of y is i  
  }  //closes for loop
}

/******************************************* sweep ******************************************/  
void Simulation::sweep()
{
  int i;
  int N1 = N/2; //number of Metropolis steps before Wolff step
  int N2 = N - N1;  //number of Metropolis steps after Wolff step
  
  for( i=0; i<N1; i++ )
  { metropolisStep(); }
  
  wolffStep();
  
  for( i=0; i<N2; i++ )
  { metropolisStep(); }
}

/***************************************** wolffStep ****************************************/
void Simulation::wolffStep()
{
  int numNeighbours = 4;
  int i, site;
  double PAdd;
  double exponent;  //exponent for PAdd
  VecND* r = new VecND(spinDim,randomGen);
  VecND* reflectedSpin;
  double onsiteEnergy_initial;
  double onsiteEnergy_final;
  double onsiteEnergy_diff;
  double PAcceptCluster;
  
  site = randomGen->randInt(N-1);
  //flipSpin(site, r);
  cluster->push_back(site);
  buffer->push_back(site);
  
  while( !buffer->empty() )
  {
    site = buffer->back();
    buffer->pop_back();
    for( i=0; i<numNeighbours; i++ )
    {
      //Note: spins[site] is not flipped yet so we have to consider the energy difference that
      //      would result if it were already flipped:
      reflectedSpin = spins[site]->getReflection(r);
      exponent = (2.0*J/T)*( r->dot(reflectedSpin) ) * ( r->dot(spins[neighbours[site][i]]) );
      if (exponent < 0 )
      {
        PAdd = 1.0 - exp(exponent);
        if( (randomGen->randDblExc() < PAdd) && !(isInCluster(neighbours[site][i])) )
        { 
          //flipSpin(neighbours[site][i],r); 
          cluster->push_back( neighbours[site][i] );
          buffer ->push_back( neighbours[site][i] );
        }
      }
      
      if(reflectedSpin!=NULL)
      { delete reflectedSpin; }
      reflectedSpin = NULL;
    }  //for loop over neighbours
  } //while loop for buffer
  
  //flip the cluster if it is accepted:
  onsiteEnergy_initial = getClusterOnSiteEnergy();
  flipCluster(r);
  onsiteEnergy_final = getClusterOnSiteEnergy();
  onsiteEnergy_diff = onsiteEnergy_final - onsiteEnergy_initial;
  
  //If the onsite energy diff. is negative, accept the cluster move. If the onsite energy 
  //diff. is positive, accept the cluster move with probability exp^(-beta*onsiteEnery_diff).
  //Note: cluster is already flipped currently, so this is a check to see if we need to flip 
  //      it back:
  if( onsiteEnergy_diff > 0 )
  {
    PAcceptCluster = exp(-1.0*onsiteEnergy_diff/T);
    //check if we need to flip the cluster back to its original orientation:
    if( randomGen->randDblExc() > PAcceptCluster )
    { flipCluster(r); }
  }
  
  clearCluster();
  
  if(r!=NULL)
  { delete r; }
  r = NULL;
}

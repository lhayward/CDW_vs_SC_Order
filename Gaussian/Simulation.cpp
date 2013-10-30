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
Simulation::Simulation(double J, double sigmaBar, vector<double>* TList, int L, 
                       MTRand* randomGen, int numWarmUpSweeps, int sweepsPerMeas, 
                       int measPerBin, int numBins, const char* outputFileName)
{
  spinDim               = 2;
  maxZ                  = 4;
  this->J               = J;
  this->sigmaBar        = sigmaBar;
  this->TList           = TList;
  this->T               = TList->at(0);
  this->L               = L;
  this->N               = L*L;
  this->randomGen       = randomGen;
  this->numWarmUpSweeps = numWarmUpSweeps;
  this->sweepsPerMeas   = sweepsPerMeas;
  this->measPerBin      = measPerBin;
  this->numBins         = numBins;
  this->outputFileName  = outputFileName;
  
  spins = new VecND*[N+1];
  spins[N] = new VecND(spinDim,0);  //The (N+1)st spin is an effective neighbour for spins
                                    //on the lattice boundary. The value of this spin is zero
                                    //so that it does not affect values of observables.
  
  neighbours = new int*[N];
  for( int i=0; i<N; i++ )
  { neighbours[i] = new int[maxZ]; }
  crossProds = new int*[N];
  for( int i=0; i<N; i++ )
  { crossProds[i] = new int[maxZ]; }
  coordNums = new int[N];
  setUpNeighbours();
  
  randomizeLattice();
  
  calculateEnergy();
  mag = new VecND(spinDim,0);
  calculateMagnetization();
  
  /*inCluster = new bool[N];
  for( int i=0; i<N; i++ )
  { inCluster[i] = 0; }
  
  cluster  = new vector<int>;
  buffer   = new vector<int>;*/
}

/********************************* ~Simulation (destructor) *********************************/
Simulation::~Simulation()
{
  for( int i=0; i<(N+1); i++ )
  { delete spins[i]; }
  delete[] spins;
  
  for(int i=0; i<N; i++)
  { delete[] neighbours[i]; }
  delete[] neighbours; 
  
  for(int i=0; i<N; i++)
  { delete[] crossProds[i]; }
  delete[] crossProds; 
  
  if( coordNums!=NULL )
  { delete[] coordNums; }
  
  if(mag!=NULL)
  { delete mag; }
  mag = NULL;
  
  /*delete[] inCluster;
  delete cluster;
  delete buffer;*/
}

/****************************************** runSim ******************************************/ 
void Simulation::runSim()
{
  unsigned int      TIndex;
  ofstream          outFile;
  double            currPsiSq;
  VecND*            currn;
  VecND*            currnSq;
  double            aveE, aveESq;
  double            aveDiamag;
  double            aveHelicityX, aveHelicityY;
  double            avePsiSq, avePsi4;
  //double            aveSF;
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
      aveDiamag=0;
      aveHelicityX=0;
      aveHelicityY=0;
      avePsiSq=0;
      avePsi4=0;
      aven->clear();
      avenSq->clear();
      //aveSF = 0;
  
      for( int j=0; j<measPerBin; j++ )
      { 
        //perform one measurement:
        for( int k=0; k<sweepsPerMeas; k++ )
        { sweep(); }
        
        //put this section in to avoid rounding errors (rather than trusting the energy and
        //magnetization we have been updating ourselves):
        calculateEnergy();
        calculateMagnetization();
        
        currn = mag->getMultiple(1.0/N);
        currnSq = mag->getSqComponents();
        currnSq->multiply(1.0/(N*N));
        currPsiSq = currnSq->v_[0] + currnSq->v_[1];
        
        aveE         += energy/N;
        aveESq       += pow( (energy/N), 2);
        //aveSF      += getSF();
        aveDiamag    += getDiamag()/(1.0*T);
        aveHelicityX += getHelicityModulus(0);
        aveHelicityY += getHelicityModulus(1);
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
      aveDiamag    = aveDiamag/measPerBin;
      aveHelicityX = aveHelicityX/measPerBin;
      aveHelicityY = aveHelicityY/measPerBin;
      //aveSF      = aveSF/measPerBin;
      avePsiSq     = avePsiSq/measPerBin;
      avePsi4      = avePsi4/measPerBin;
      aven->multiply(1.0/measPerBin);
      avenSq->multiply(1.0/measPerBin);
      
      outFile << L << '\t' << T << '\t' << (i+1) << '\t' << aveE << '\t' << aveESq << '\t'
              << aveDiamag << '\t' << aveHelicityX << '\t' << aveHelicityY << '\t'
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
  energy=0;
  double energyOnSite = 0;
  double energyNeigh = 0;
  
  for( int i=0; i<N; i++ )
  { energyOnSite += (coordNums[i] + sigmaBar)*spins[i]->getSquare(); }
  energyOnSite *= 1.0/2.0;
  
  for( int i=0; i<N; i++ )
  { energyNeigh += spins[i]->dot( spins[neighbours[i][0]] )
                   + spins[i]->dot( spins[neighbours[i][2]] ); }
  energyNeigh *= -1;
  
  energy = J*(energyOnSite + energyNeigh);
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
/*void Simulation::clearCluster()
{
  int site;
  
  while( !cluster->empty() )
  { 
    site = cluster->back();
    inCluster[site]=0;
    cluster->pop_back(); 
  }
  
  for( int i=0; i<N; i++ )
  {
    if( inCluster[i]==1 )
    { cout << "*** 1 ***" << endl; }
  }
}*/

/**************************************** flipCluster *****************************************
* This function flips all spins in the cluster by using the passed vector r.
**********************************************************************************************/
/*void Simulation::flipCluster(VecND* r)
{
  int clustSize = (int)cluster->size();
  
  for( int i=0; i<clustSize; i++ )
  { spins[cluster->at(i)]->reflectAndNormalize(r); } //for loop
}*/

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
/*double Simulation::getClusterOnSiteEnergy()
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
}*/

/***************************************** getDiamag *****************************************/
double Simulation::getDiamag()
{
  double diamag = 0;  // M/B
  double term1 = 0; //First term in the expression for M/B
  double term2 = 0; //Second term in the expression for M/B
  
  for( int i=0; i<N; i++ )
  {
    //sum over indices corresponding to a=+\hat{x} and a=+\hat{y}:
    for( int j=0; j<maxZ; j=j+2 )
    { term1 += crossProds[i][j]*crossProds[i][j]*spins[i]->dot( spins[ neighbours[i][j] ] ); }
  }
  term1 *= -1*J/(16.0*N);
  
  for( int i=0; i<N; i++ )
  {
    //sum over indices corresponding to a=+\hat{x} and a=+\hat{y}:
    for( int j=0; j<maxZ; j=j+2 )
    { term2 += crossProds[i][j]*( spins[i]->v_[0]*spins[neighbours[i][j]]->v_[1] ) 
                                  - spins[i]->v_[1]*spins[neighbours[i][j]]->v_[0]; }
  }
  term2 = J*J*term2*term2/(16.0*T*N);
  
  diamag = term1 + term2;
  return diamag;
}

/************************************** getCorrelation ***************************************/
double Simulation::getCorrelation(int i, int j)
{
  return spins[i]->dot(spins[j]);
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
     helicityMod = sum1/(1.0*N) - J/(1.0*T*N)*pow(sum2,2); 
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

/**************************************** isInCluster ****************************************
* This function checks whether or not spin at the passed site is in the cluster.
*********************************************************************************************/  
/*bool Simulation::isInCluster(int site)
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
}*/

/************************************** metropolisStep *************************************/  
void Simulation::metropolisStep()
{
  int site;
  double deltaE;
  double mean = 0;
  double stddev;
  VecND* sNew; // = new VecND(spinDim, randomGen);
  VecND* nnSum = new VecND(spinDim,0);
  
  //Generate the new spin using a Gaussian distribution based on the on-site term in the 
  //Hamiltonian:
  site = randomGen->randInt(N-1);
  stddev = sqrt( T/(J*1.0*(coordNums[site] + sigmaBar)) );;
  sNew = new VecND( spinDim, randomGen, mean, stddev );
  
  //loop to calculate the nearest neighbour sum:
  for( int i=0; i<maxZ; i++ )
  { nnSum->add(spins[neighbours[site][i]]); }
  
  //The on-site contribution is taken care of by the normal distribution used to generate the
  //random new spin. Therefore, delta(E) is from nearest neighbour interactions only:
  deltaE = -1*J*( nnSum->dot(sNew) - nnSum->dot(spins[site]) );
  
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

void Simulation::printNeighbours()
{
  for( int i=0; i<N; i++ )
  {
    std::cout << "Spin " << i << ": " << std::endl;
    std::cout << "  Co-ordination Number = " << coordNums[i] << std::endl;
    
    std::cout << "  Neighbours: [  ";
    for( int j=0; j<maxZ; j++ )
    { std::cout << neighbours[i][j] << "  "; }
    std::cout << "]" << std::endl;
    
    std::cout << "  Cross Products: [  ";
    for( int j=0; j<maxZ; j++ )
    { std::cout << crossProds[i][j] << "  "; }
    std::cout << "]" << std::endl;
    
  }
  std::cout << std::endl;
}

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
  double mean;
  double stddev;
  
  for( int i=0; i<N; i++ )
  {
    mean = 0;
    stddev = sqrt( T/(J*1.0*(coordNums[i] + sigmaBar)) );
    spins[i] = new VecND(spinDim,randomGen, mean, stddev); 
  }
}

/************************************** setUpNeighbours **************************************
* This functions sets up the "neighbours" and "coordNums" arrays for a square lattice with 
* open boundary conditions. 
*********************************************************************************************/
void Simulation::setUpNeighbours()
{
  int ix,  iy;  //coordinates of spin i
  int inx, iny; //coordinates of spin i's neighbour
  
  //loop to assign the neighbours (with open boundary conditions):
  for( int i=0; i<N; i++ )
  {
    coordNums[i] = 0;
    
    ix = i%L;
    iy = (i-ix)/L;
    
    //-----------------------------------------------------------------------------------------
    //assign the neighbour to the right (+x direction):
    inx=ix+1;
    iny=iy;
    //if spin i is not on the right-hand boundary:
    if( inx < L )
    {
      coordNums[i]++;
      neighbours[i][0] = iny*L + inx;
      crossProds[i][0] = -1*(2*iy - L + 1);
      
    }
    //if spin i is on the right-hand boundary:
    else
    {
      neighbours[i][0] = N; //neighbour is the extra "spin" at the end of the spin array
                            //(with value 0)
      crossProds[i][0] = 0;
    }
    
    //-----------------------------------------------------------------------------------------
    //assign the neighbour to the left (-x direction):
    inx=ix-1;
    iny=iy;
    //if spin i is not on the left-hand boundary:
    if( inx >= 0 )
    {
      coordNums[i]++;
      neighbours[i][1] = iny*L + inx;
      crossProds[i][1] = 2*iy - L + 1;
      
    }
    //if spin i is on the left-hand boundary:
    else
    {
      neighbours[i][1] = N; //neighbour is the extra "spin" at the end of the spin array
                            //(with value 0)
      crossProds[i][1] = 0;
    }
    
    //-----------------------------------------------------------------------------------------
    //assign the neighbour to the top (+y direction):
    inx=ix;
    iny=iy+1;
    //if spin i is not on the top boundary:
    if( iny < L)
    {
      coordNums[i]++;
      neighbours[i][2] = iny*L + inx;
      crossProds[i][2] = 2*ix - L + 1;
      
    }
    //if spin i is on the top boundary:
    else
    {
      neighbours[i][2] = N; //neighbour is the extra "spin" at the end of the spin array
                            //(with value 0)
      crossProds[i][2] = 0;
    }
    
    //-----------------------------------------------------------------------------------------
    //assign the neighbour to the bottom (-y direction):
    inx=ix;
    iny=iy-1;
    //if spin i is not on the bottom boundary:
    if( iny >= 0 )
    {
      coordNums[i]++;
      neighbours[i][3] = iny*L + inx;
      crossProds[i][3] = -1*(2*ix - L + 1);
      
    }
    //if spin i is on the bottom boundary:
    else
    {
      neighbours[i][3] = N; //neighbour is the extra "spin" at the end of the spin array
                            //(with value 0)
      crossProds[i][3] = 0;
    }
    
  }  //closes for loop
}

/******************************************* sweep ******************************************/  
void Simulation::sweep()
{ 
  /*int i;
  int N1 = N/2; //number of Metropolis steps before Wolff step
  int N2 = N - N1;  //number of Metropolis steps after Wolff step
  
  for( i=0; i<N1; i++ )
  { metropolisStep(); }
  
  if( lambda==1 )
  { wolffStep(); }
  
  for( i=0; i<N2; i++ )
  { metropolisStep(); }*/
  
  for( int i=0; i<N; i++ )
  { metropolisStep(); }
}

/***************************************** wolffStep ****************************************/
/*void Simulation::wolffStep()
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
  inCluster[site] = 1;
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
        //if( (randomGen->randDblExc() < PAdd) && !(isInCluster(neighbours[site][i])) )
        if( !( inCluster[ neighbours[site][i] ] ) && (randomGen->randDblExc() < PAdd) )
        { 
          //flipSpin(neighbours[site][i],r); 
          inCluster[ neighbours[site][i] ] = 1;
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
  flipCluster(r); //flip in order to calculate the final energy
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
}*/

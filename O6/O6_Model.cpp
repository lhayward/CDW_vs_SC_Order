/**********************************************************************************************
************************************ O(6) MODEL MONTE CODE ************************************
***********************************************************************************************
* Lauren E. Hayward Sierens
***********************************************************************************************
* File:   O6_Model.cpp
**********************************************************************************************/

#include <fstream>
#include <iostream>
#include <string>
#include <typeinfo>
#include "FileReading.h"
#include "O6_Model.h"

//typdef needed because uint is a return types:
typedef O6_Model::uint uint;

/***************** O6_Model(std::ifstream* fin, std::string outFileName, ... ******************
******************          Lattice* lattice, MTRand* randomGen)             ******************
******************                       (constructor)                       *****************/
O6_Model::O6_Model(std::ifstream* fin, std::string outFileName, Lattice* lattice,
                   MTRand* randomGen)
  : Model(fin, outFileName)
{
  const char EQUALS_CHAR = '=';
  
  std::cout.precision(15);
  
  if( fin!=NULL && fin->is_open() )
  {
    lambda_ = FileReading::readDouble(fin, EQUALS_CHAR);
    g_      = FileReading::readDouble(fin, EQUALS_CHAR);
    gPrime_ = FileReading::readDouble(fin, EQUALS_CHAR);
    w_      = FileReading::readDouble(fin, EQUALS_CHAR);
  
    if( lattice != NULL )
    {
      hrect_ = dynamic_cast<Hyperrectangle *>(lattice);
      if(hrect_)
      {
        D_  = hrect_->getD();
        
        //This model only works in two or three dimension:
        if( D_==2 || D_==3 )
        {
          L_  = hrect_->getL();
          N_  = hrect_->getN();
          
          Vz_      = 0;
          VzPrime_ = 0;
          
          //if D_=3, then the input file should also specify the interlayer coupling strength:
          if( D_==3 )
          {
            Vz_      = FileReading::readDouble(fin, EQUALS_CHAR);
            VzPrime_ = FileReading::readDouble(fin, EQUALS_CHAR);
          } //if for D_==3
          
          mag_ = new Vector_NDim(VECTOR_SPIN_DIM,0);
          
          spins_ = new VectorSpins(N_, VECTOR_SPIN_DIM);
          randomizeLattice(randomGen);
        } //if for dimension
        else
        {
          std::cout << "ERROR in O6_Model constructor:\n" 
                    << "  The Hyperrectangle lattice must have dimension 2 or 3." 
                    << std::endl;
        }
      }
      else
      {
        std::cout << "ERROR in O6_Model constructor:\n" 
                  << "  A lattice of type Hyperrectangle is required.\n"
                  << "  A lattice of type " << typeid(*lattice).name() << " was given.\n" 
                  << std::endl;
      }
    }
    else
    {
      std::cout << "ERROR in O6_Model constructor: The passed Lattice object is not "
                << "valid\n" << std::endl; 
    }
  }
  else
  { 
    std::cout << "ERROR in Model constructor: could not read from input file\n" << std::endl; 
  }
}

/******************************* ~O6_Model() (destructor) *******************************/
O6_Model::~O6_Model()
{
  //delete the VectorSpins object:
  if( spins_ != NULL )
  { delete spins_; }
  spins_ = NULL;
  
  if( mag_ != NULL )
  { delete mag_; }
  mag_ = NULL;
}

/************************************* getHelicityModulus *************************************
* Calculates and returns the helicity modulus for the desired lattice directions.
* Note: dir=0 corresponds to the x-direction
*       dir=1 corresponds to the y-direction
*       dir=2 corresponds to the z-direction (only applicable when D_=3)
**********************************************************************************************/
double O6_Model::getHelicityModulus(int dir)
{
  double       helicityMod = 0; //resulting helicity modulus
  //int neighDir;  //direction (for neighbour's array) corresponding to the specified dir
  Vector_NDim* currSpin;        //current spin
  Vector_NDim* neigh;           //current neighbour
  double       sum1;            //first sum in formula for helicity modulus
  double       sum2;            //second sum in formula for helicity modulus
  
  //compute the helicity modulus if a valid direction was passed (otherwise just return 0):
  if( dir==0 || dir==1 || (dir==2 && D_==3) )
  {
     //neighDir = 2*dir;  //0 if x-direction (dir=0), 2 if y-direction (dir=1)
     sum1=0;
     sum2=0;
     
     //loop over all spins:
     for(int i=0; i<N_; i++)
     {
      currSpin = spins_->getSpin(i);
      neigh    = spins_->getSpin( hrect_->getNeighbour(i,dir) ); //nearest neighbour along +x
      //neigh = spins[neighbours[i][neighDir]];
      
      sum1 += currSpin->dotForRange(neigh,0,1);
      sum2 += ( currSpin->v_[0] * neigh->v_[1] ) - ( currSpin->v_[1] * neigh->v_[0] );
     }
     helicityMod = sum1/(1.0*N_) - J_/(T_*N_)*pow(sum2,2); 
  } //end if
  
  return helicityMod;
}

/******************************* localUpdate(MTRand* randomGen) ******************************/
void O6_Model::localUpdate(MTRand* randomGen)
{
  const uint   dir_x = 0;
  const uint   dir_y = 1;
  const uint   dir_z = 2;
  
  uint         latticeSite; //randomly selected spin location
  double       deltaE;
  double       oldCDWSumSqs;
  double       newCDWSumSqs;
  Vector_NDim* spin_old;  //previous state of spin (at randomly selected lattice site)
  Vector_NDim* spin_new;  //new proposed value for spin
  Vector_NDim* nnSum_xy;
  Vector_NDim* nnSum_z;
  
  //randomly generate a new spin:
  spin_new = new Vector_NDim(VECTOR_SPIN_DIM, randomGen);
  
  //randomly select a spin on the lattice:
  latticeSite = randomGen->randInt(N_-1);
  spin_old    = spins_->getSpin(latticeSite);
  
  //loop to calculate the nearest neighbour sum for the xy-plane:
  nnSum_xy = new Vector_NDim(VECTOR_SPIN_DIM, 0);
  for( int i=dir_x; i<=dir_y; i++ )
  { 
    nnSum_xy->add( spins_->getSpin( hrect_->getNeighbour( latticeSite, i    ) ) ); 
    nnSum_xy->add( spins_->getSpin( hrect_->getNeighbour( latticeSite, i+D_ ) ) );
  }
  
  //calculate the nearest neighbour sum for the z-direction (if D_=3):
  nnSum_z  = new Vector_NDim(VECTOR_SPIN_DIM, 0);
  if( D_==3 )
  {
    nnSum_z->add( spins_->getSpin( hrect_->getNeighbour( latticeSite, dir_z    ) ) ); 
    nnSum_z->add( spins_->getSpin( hrect_->getNeighbour( latticeSite, dir_z+D_ ) ) );
  }
  
  oldCDWSumSqs = spin_old->getSquareForRange( 2, VECTOR_SPIN_DIM-1 );
  newCDWSumSqs = spin_new->getSquareForRange( 2,VECTOR_SPIN_DIM-1 );
  
  //calculate the energy change for the proposed move:
  deltaE = J_*( -1*( nnSum_xy->dotForRange(spin_new,0,1) 
                     - nnSum_xy->dotForRange(spin_old,0,1) )
               - lambda_*( nnSum_xy->dotForRange(spin_new,2,VECTOR_SPIN_DIM-1) 
                           - nnSum_xy->dotForRange(spin_old,2,VECTOR_SPIN_DIM-1) ) 
               + (g_ + 4*(lambda_-1.0))/2.0*( newCDWSumSqs - oldCDWSumSqs )
               + gPrime_/2.0*( pow(newCDWSumSqs,2.0) - pow(oldCDWSumSqs,2.0) )
               + w_/2.0*( pow(spin_new->getSquareForRange(2,3),2.0) 
                         + pow(spin_new->getSquareForRange(4,5),2.0) 
                         - pow(spin_old->getSquareForRange(2,3),2.0) 
                         - pow(spin_old->getSquareForRange(4,5),2.0) ) );
  //if D_=3 then also include the contribution to deltaE from interlayer coupling:
  if( D_==3 )
  {
    deltaE += J_*( -Vz_*( nnSum_z->dotForRange(spin_new,2,VECTOR_SPIN_DIM-1) 
                          - nnSum_z->dotForRange(spin_old,2,VECTOR_SPIN_DIM-1) ) 
                   - VzPrime_*( nnSum_z->dotForRange(spin_new,0,1) 
                                - nnSum_z->dotForRange(spin_old,0,1) ) );
  }
  
  //if the move is accepted:
  if( deltaE<=0 || randomGen->randDblExc() < exp(-deltaE/T_) )
  { 
    //delete the vector storing the old state of the spin:
    if(spin_old!=NULL)
    { delete spin_old; }
    spin_old = NULL;
    
    spins_->setSpin( latticeSite, spin_new );
    energy_ += deltaE;
  }
  //otherwise, the move is rejected:
  else
  {
    //delete the vector storing the rejected new state of the spin:
    if(spin_new!=NULL)
    { delete spin_new; }
    spin_new = NULL;  
  }
  
  //delete the vectors storing the nearest neighbour sums:
  if(nnSum_xy!=NULL)
  { delete nnSum_xy; }
  nnSum_xy = NULL;
  
  if(nnSum_z!=NULL)
  { delete nnSum_z; }
  nnSum_z = NULL;
}

/************************************* makeMeasurement() *************************************/
void O6_Model::makeMeasurement()
{
  double energyPerSpin = energy_/(1.0*N_);
  
  measures.accumulate( "E",   energyPerSpin ) ;
  measures.accumulate( "ESq", pow(energyPerSpin,2) );
}

/*************************************** printParams() ***************************************/
void O6_Model::printParams()
{
  std::cout << "O(6) Model Parameters:\n"
            << "---------------------" << std::endl;
  Model::printParams();
  std::cout << "  lambda = " << lambda_ << "\n"
            << "       g = " << g_      << "\n"
            << "      g' = " << gPrime_ << "\n"
            << "       w = " << w_      << "\n";
  //if D_==3, then also print the interlayer coupling:
  if( D_==3 )
  {
    std::cout << "      Vz = " << Vz_      << "\n"
              << "     Vz' = " << VzPrime_ << "\n";
  }
  std::cout << std::endl;
}

/**************************************** printSpins() ***************************************/
void O6_Model::printSpins()
{ spins_->print(); }

/**************************** randomizeLattice(MTRand* randomGen) ****************************/
void O6_Model::randomizeLattice(MTRand* randomGen)
{ 
  spins_->randomize(randomGen);
  updateEnergy();
  updateMagnetization();
}

/************************************* setT(double newT) *************************************/
void O6_Model::setT(double newT)
{ 
  Model::setT(newT); 
}

/********************************** sweep(MTRand* randomGen) *********************************/
void O6_Model::sweep(MTRand* randomGen)
{
  for( uint i=0; i<N_; i++ )
  { localUpdate(randomGen); }
  updateEnergy();
  updateMagnetization();
}

/******************************** uintPower(int base, int exp) *******************************/
uint O6_Model::uintPower(uint base, uint exp)
{
  uint result = 1;
  for(uint i=1; i<=exp; i++)
  { result *= base; } 
  
  return result;
} //uintPower method

/*************************************** updateEnergy() **************************************/
void O6_Model::updateEnergy()
{
  Vector_NDim* currSpin;
  Vector_NDim* neighbour_x;
  Vector_NDim* neighbour_y;
  Vector_NDim* neighbour_z;
  double       sumCDWSqs     = 0;
  double       energy1       = 0;
  double       energyLambda  = 0;
  double       energyg       = 0;
  double       energygPrime  = 0;
  double       energyw       = 0;
  double       energyVz      = 0;
  double       energyVzPrime = 0;
  
  energy_=0;
  
  for( uint i=0; i<N_; i++ )
  { 
    currSpin    = spins_->getSpin(i);
    neighbour_x = spins_->getSpin( hrect_->getNeighbour(i,0) ); //nearest neighbour along +x
    neighbour_y = spins_->getSpin( hrect_->getNeighbour(i,1) ); //nearest neighbour along +y
    
    energy1 += currSpin->dotForRange( neighbour_x, 0, 1 ) 
               + currSpin->dotForRange( neighbour_y, 0, 1 ); 
    
    energyLambda += currSpin->dotForRange( neighbour_x, 2, VECTOR_SPIN_DIM-1 ) 
                    + currSpin->dotForRange( neighbour_y, 2, VECTOR_SPIN_DIM-1 );
  }
  energy1 *= -1;
  energyLambda *= -1*lambda_;
  
  for( uint i=0; i<N_; i++ )
  { 
    sumCDWSqs = spins_->getSpin(i)->getSquareForRange(2,VECTOR_SPIN_DIM-1); 
    energyg += sumCDWSqs;
    energygPrime += pow(sumCDWSqs,2.0);
  }
  energyg *= (g_ + 4.0*(lambda_-1.0))/2.0;
  energygPrime *= gPrime_/2.0;
  
  for( uint i=0; i<N_; i++ )
  { energyw += pow( spins_->getSpin(i)->getSquareForRange(2,3), 2.0 ) 
               + pow( spins_->getSpin(i)->getSquareForRange(4,5), 2.0 ); }
  energyw *= w_/2.0;
  
  //in three dimensions, then also add the energy due to interlayer coupling:
  if( D_==3 )
  {
    for( uint i=0; i<N_; i++ )
    { 
      currSpin    = spins_->getSpin(i);
      neighbour_z = spins_->getSpin( hrect_->getNeighbour(i,2) ); //nearest neighbour along +z
      
      energyVz      += currSpin->dotForRange( neighbour_z, 2, VECTOR_SPIN_DIM-1 ); 
      energyVzPrime += currSpin->dotForRange( neighbour_z, 0, 1 );
    }
    energyVz      *= -1*Vz_;
    energyVzPrime *= -1*VzPrime_;
  }//if for D_==3
  
  energy_ = J_*(energy1 + energyLambda + energyg + energygPrime + energyw 
                + energyVz + energyVzPrime);
}

/*********************************** updateMagnetization() ***********************************/
void O6_Model::updateMagnetization()
{
  if( mag_ != NULL )
  { mag_->clear(); }
  else
  { mag_ = new Vector_NDim(VECTOR_SPIN_DIM,0); }
  
  for( uint i=0; i<N_; i++ )
  { mag_->add( spins_->getSpin(i) ); }
}

/***************************** writeBin(int binNum, int numMeas) *****************************/
void O6_Model::writeBin(int binNum, int numMeas)
{
  //if this is the first bin being written to file, then also write a line of text explaining
  //each column:
  if( binNum == 1)
  {
    fout << "# L \t T \t binNum";
    measures.writeMeasNames(&fout);
    fout << std::endl;
  }
  fout << L_[0] << '\t' << T_ << '\t' << binNum;
  measures.writeAverages(&fout, numMeas);
  fout << std::endl;
}
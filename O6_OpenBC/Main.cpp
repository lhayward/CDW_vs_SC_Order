/*********************************************************************************************
***************************************** CDW vs. SC *****************************************
**********************************************************************************************
* Lauren Hayward
**********************************************************************************************
* File: Main.cpp
* Note: MersenneTwister.h was taken from external sources and is used as the random number
*       generator
*********************************************************************************************/

#include <cmath>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include "MersenneTwister.h"
#include "Simulation.h"

//using namespace std;

typedef unsigned long ulong;

void   printDoubleVec(vector<double>* vec);
void   readInput  (const char* inFileName, double* J, double* lambda, double* g, double* w, 
                   vector<double>* TList, int* LMin, int* LMax, int* LSteps, ulong* randomSeed, 
                   int* numWarmUpSweeps,  int* sweepsPerMeas, int* measPerBin, int* numBins);
double readDouble    (ifstream* in, char delim);
void   readDoubleList(vector<double>* list, ifstream* in, char delim, char startChar, char endChar);
int    readInt       (ifstream* in, char delim); 
ulong  readLongInt   (ifstream* in, char delim);

/*********************************************************************************************
******************************************** main ********************************************
*********************************************************************************************/
int main(int argc, char** argv) 
{
  //variables:
  int           i, L; 
  MTRand*       randomGen;
  Simulation*   sim;
  ostringstream converter;
  string        inFileName = "input";
  string        outFileName;
  time_t        sec1, sec2;
  
  std::cout.precision(20);
  
  //variables to be read from file (values here are default values): 
  double          J=2.0, lambda=1, g=0, w=0;
  int             LMax=20, LMin=10, LSteps=2;
  vector<double>* TList = new vector<double>;
  ulong           randomSeed = 12345678;
  int             numWarmUpSweeps = 100, sweepsPerMeas = 1, measPerBin=1000, numBins = 10;
  
  if( argc >1 )
  { inFileName = inFileName + '_' + argv[1]; }
  inFileName = inFileName + ".txt";
  std::cout << inFileName << std::endl;
  
  std::cout << "\n***STARTING SIMULATION***\n" << std::endl;
  
  readInput(inFileName.c_str(), &J, &lambda, &g, &w, TList, &LMin, &LMax, &LSteps, &randomSeed, 
            &numWarmUpSweeps, &sweepsPerMeas, &measPerBin, &numBins);
  randomGen = new MTRand(randomSeed);
  
  std::cout << "                               J = " << J << std::endl
            << "                          lambda = " << lambda << std::endl
            << "                               g = " << g << std::endl
            << "                               w = " << w << std::endl;
  std::cout << "                Temperature List = ";
  printDoubleVec(TList);
  
  std::cout << "                       Minimum L = " << LMin << std::endl
            << "                       Maximum L = " << LMax << std::endl
            << "               Number of L Steps = " << LSteps << std::endl
            << "           Random Generator Seed = " << randomSeed << std::endl
            << "         Number of Warmup Sweeps = " << numWarmUpSweeps << std::endl
            << "Number of Sweeps per Measurement = " << sweepsPerMeas << std::endl
            << "  Number of Measurements per Bin = " << measPerBin << std::endl 
            << "                  Number of Bins = " << numBins << std::endl << std::endl;
  
  for(i=0; i<LSteps; i++)
  {
    L = LMin + i*(LMax-LMin)/max((LSteps-1),1);
    
    converter.str("");
    if( argc > 1 )
    { converter << "bins_" << argv[1] << ".txt"; }
    else
    { converter << "bins_L" << L << ".txt"; }
    outFileName = converter.str();
    std::cout << "OUTFILE = " << outFileName<< std::endl;
    sim = new Simulation(J, lambda, g, w, TList, L, randomGen, numWarmUpSweeps, sweepsPerMeas, 
                         measPerBin, numBins, outFileName.c_str());
    sec1 = time (NULL);
    sim->runSim();
    sec2 = time(NULL);
    std::cout << "Time for L = " << L << ": " << (sec2 - sec1) << " seconds" << std::endl;
    
    if(sim!=NULL)
    { delete sim; }
    sim = NULL;
  }

  std::cout << "\n***END OF SIMULATION***\n" << std::endl;
  return 0;
} //closes main

/************************************** printDoubleVec **************************************/
void printDoubleVec(vector<double>* vec)
{
  size_t i;
  
  std::cout << "[ ";
  for( i=0; i<vec->size() ;i++ )
  { 
    std::cout << vec->at(i); 
    if( i<(vec->size() - 1) )
    { std::cout << ","; }
    std::cout << " ";
  }
  std::cout << "]" << std::endl;
}


/***************************************** readInput ****************************************/
void readInput(const char* inFileName, double* J, double* lambda, double* g, double* w,
               vector<double>* TList, int* LMin, int* LMax, int* LSteps, ulong* randomSeed, 
               int* numWarmUpSweeps, int* sweepsPerMeas, int* measPerBin, int* numBins)
{
  const char EQUALS_CHAR = '=';
  const char LIST_START_CHAR = '[';
  const char LIST_END_CHAR = ']';
  ifstream   inFile;
  
  inFile.open(inFileName);
  
  if( inFile.is_open() )
  {
    *J               = readDouble (&inFile, EQUALS_CHAR);
    *lambda          = readDouble (&inFile, EQUALS_CHAR);
    *g               = readDouble (&inFile, EQUALS_CHAR);
    *w               = readDouble (&inFile, EQUALS_CHAR);
    readDoubleList(TList, &inFile, EQUALS_CHAR, LIST_START_CHAR, LIST_END_CHAR);
    *LMin            = readInt    (&inFile, EQUALS_CHAR);
    *LMax            = readInt    (&inFile, EQUALS_CHAR);
    *LSteps          = readInt    (&inFile, EQUALS_CHAR);
    *randomSeed      = readLongInt(&inFile, EQUALS_CHAR);
    *numWarmUpSweeps = readInt    (&inFile, EQUALS_CHAR);
    *sweepsPerMeas   = readInt    (&inFile, EQUALS_CHAR);
    *measPerBin      = readInt    (&inFile, EQUALS_CHAR);
    *numBins         = readInt    (&inFile, EQUALS_CHAR);
  }
  else{ std::cout << "Could not find file \"" << inFileName << "\"" << std::endl; }
  
  inFile.close();
}

/**************************************** readDouble ****************************************/
double readDouble(ifstream* in, char delim)
{
  string currLine;
  int index;
  
  getline(*in, currLine);
  index = currLine.find_last_of(delim);
  return strtod( (currLine.substr(index+1)).c_str(), NULL);
}

/************************************** readDoubleList **************************************/
void readDoubleList(vector<double>* list, ifstream* in, char delim, char startChar, char endChar)
{
  string currLine;
  size_t delimIndex, startIndex, endIndex, commaIndex;
  
  getline(*in, currLine);
  delimIndex = currLine.find_last_of(delim);
  currLine = currLine.substr(delimIndex+1);
  
  startIndex = currLine.find_first_of(startChar);
  endIndex = currLine.find_last_of(endChar);
  currLine = currLine.substr(startIndex+1, (endIndex - startIndex - 1) );
  
  commaIndex = currLine.find_first_of(",");
  while( commaIndex>=0 && commaIndex<currLine.size() )
  {
    list->push_back( strtod( (currLine.substr(0, commaIndex)).c_str(), NULL) );
    currLine = currLine.substr(commaIndex+1);
    commaIndex = currLine.find_first_of(",");
  }
  list->push_back( strtod( currLine.c_str(), NULL) );
}

/****************************************** readInt *****************************************/
int readInt(ifstream* in, char delim)
{
  string currLine;
  int index;
  
  getline(*in, currLine);
  index = currLine.find_last_of(delim);
  return strtol( (currLine.substr(index+1)).c_str(), NULL, 0);
}

/**************************************** readLongInt ***************************************/
ulong readLongInt(ifstream* in, char delim)
{
  string currLine;
  int index;
  
  getline(*in, currLine);
  index = currLine.find_last_of(delim);
  return strtoul( (currLine.substr(index+1)).c_str(), NULL, 0);
}


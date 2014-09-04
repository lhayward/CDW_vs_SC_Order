/**********************************************************************************************
******************************** CLASSICAL REPLICA MONTE CODE *********************************
***********************************************************************************************
* Lauren Hayward
***********************************************************************************************
* File:   Model.cpp (Abstract Class)
**********************************************************************************************/

#include <fstream>
#include <iostream>
#include <string>
#include "FileReading.h"
#include "Model.h"

/************** Model(std::ifstream* fin, std::string outFileName) (constructor) *************/
Model::Model(std::ifstream* fin, std::string outFileName)
{ 
  const char EQUALS_CHAR = '=';
  
  if( fin!=NULL && fin->is_open() )
  {
    J_ = FileReading::readDouble(fin, EQUALS_CHAR);
  }
  else
  { 
    std::cout << "ERROR in Model constructor: could not read from input file\n" << std::endl; 
  }
  
  fout.open(outFileName.c_str());
  fout.precision(15);
  
  energy_=0;
  //initialize the temperature (should be changed by user to desired temperature before
  //starting the simulation):
  T_ = 1.0;
  
  //Add measurement names to Measure object:
  measures.insert("E");
  measures.insert("ESq");
}

/*********************************** ~Model() (destructor) ***********************************/
Model::~Model()
{ fout.close(); }

/*************************************** printParams() ***************************************/
void Model::printParams()
{
  std::cout << "       J = " << J_ << "\n";
}

/************************************* setT(double newT) *************************************/
void Model::setT(double newT)
{ T_ = newT; }

/************************************* zeroMeasurements() ************************************/
void Model::zeroMeasurements()
{ measures.zero(); }

/*********************************** Public Getter Methods: **********************************/
double Model::getEnergy(){ return energy_; }
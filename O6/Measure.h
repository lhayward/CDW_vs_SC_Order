/**********************************************************************************************
******************************** CLASSICAL REPLICA MONTE CODE *********************************
***********************************************************************************************
* Lauren Hayward
***********************************************************************************************
* File:   Measure.h
**********************************************************************************************/

#ifndef MEASURE_H
#define MEASURE_H

#include <map>
#include <string>
#include <vector>

class Measure 
{ 
  public:
    typedef unsigned int  uint;
    
  private:
    std::map<std::string,double> measurements;
    std::vector<std::string> measStrings;
    
  public:
    Measure();
    virtual ~Measure();
    
    void accumulate   (std::string label, double newMeas);
    void insert       (std::string label);
    void print        ();
    void writeAverages(std::ofstream* fout, uint numMeas);
    void zero         ();
};

#endif  // MEASURE_H

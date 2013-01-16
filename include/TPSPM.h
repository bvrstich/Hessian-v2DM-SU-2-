#ifndef TPSPM_H
#define TPSPM_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

class TPTPM;

#include "RecMat.h"

/**
 * @author Brecht Verstichel
 * @date 14-01-2013\n\n
 * This class TPSPM is a class written for the singly-traced TPTPM matrix object,
 * being a rectangular matrix on TP for the rows, and SP for the columns, it inherits alle the function from its mother 
 * RecMat, a rectangular matrix, some special member functions and two lists 
 */
class TPSPM : public RecMat {

   public:
      
      //constructor
      TPSPM();

      //copy constructor
      TPSPM(const TPSPM &);

      //destructor
      virtual ~TPSPM();

      using RecMat::operator=;

      using RecMat::operator();

      //easy to access the numbers, in sp mode
      double operator()(int S,int a,int b,int c,int d,int e,int z) const;

      void dpt(double,const TPM &);

      void dpt(double,const PHM &);

      void dpt3(double,const TPTPM &);

   private:

};

#endif

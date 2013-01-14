#ifndef SPSPM_H
#define SPSPM_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

#include "Matrix.h"

/**
 * @author Brecht Verstichel
 * @date 23-11-2012\n\n
 * This class SPSPM is a class written for matrices of two particle matrices, it inherits alle the function from its mother 
 * Matrix, some special member functions and two lists that give the relationship between the sp and the tp 
 * basis.
 */

class SPSPM : public Matrix {

   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param spmm_p the SPSPM you want to print
    */
   friend ostream &operator<<(ostream &output,const SPSPM &spmm_p);

   public:
      
      //constructor
      SPSPM();

      //copy constructor
      SPSPM(const SPSPM &);

      //destructor
      virtual ~SPSPM();

      using Matrix::operator=;

      using Matrix::operator();

      const double &operator()(int,int,int,int) const;

      static int gn();

      static int gspmm2s(int,int);

      static int gs2spmm(int,int);

      static void init();

      static void clear();

   private:

      //!list relating the single-particle space to the SPSPM basis
      static vector< vector<int> > spmm2s;

      //!list relating the single-particle space to the SPSPM basis
      static int **s2spmm;

};

#endif

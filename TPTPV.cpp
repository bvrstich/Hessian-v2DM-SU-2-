#include <iostream>
#include <cmath>
#include <fstream>

using std::ostream;
using std::ofstream;
using std::ifstream;
using std::cout;
using std::endl;
using std::ios;

#include "include.h"

double *TPTPV::norm;
int TPTPV::n;

/**
 * initialize the static lists
 */
void TPTPV::init(){

   n = TPTPM::gn();

   norm = new double [n];

   int hess = 0;

   for(int i = 0;i < TPM::gn();++i)
      for(int j = i;j < TPM::gn();++j){

         if(i == j)
            norm[hess] = 0.5;
         else
            norm[hess] = std::sqrt(0.5);

         ++hess;

      }

}

/**
 * deallocate the static lists
 */
void TPTPV::clear(){

   delete [] norm;

}

/**
 * standard constructor:
 */
TPTPV::TPTPV() {

   tpvv = new double [n];
   
}

/**
 * copy constructor
 * @param tpvv_c object that will be copied into this.
 */
TPTPV::TPTPV(const TPTPV &tpvv_c){

   tpvv = new double [n];

   int inc = 1;

   dcopy_(&n,tpvv_c.gpointer(),&inc,tpvv,&inc);
   
}

/**
 * destructor
 */
TPTPV::~TPTPV(){ 

   delete [] tpvv;
   
}

/**
 * access to the norm from outside of the class
 * @param i hess index, if a == b norm = 1.0/sqrt(2.0)
 */
double TPTPV::gnorm(int i){

   return norm[i];

}

/**
 * access to the norm from outside of the class, in tp mode
 * @param I row index
 * @param J column index
 */
double TPTPV::gnorm(int I,int J){

   return norm[TPTPM::gt2tpmm(I,J)];

}

/**
 * @return the pointer to the vector, for blas and lapack stuff
 */
const double *TPTPV::gpointer() const{

   return tpvv;

}

/** 
 * access to the numbers in the vector, const
 */
const double &TPTPV::operator[](int i) const{

   return tpvv[i];

}

/** 
 * access to the numbers in the vector, no const
 */
double &TPTPV::operator[](int i){

   return tpvv[i];

}

/**
 * convert a 2DM/TPM to vector form
 * @param tpm the 2DM in question
 */
void TPTPV::convert(const TPM &tpm){

   int I,J;

   for(int i = 0;i < n;++i){

      I = TPTPM::gtpmm2t(i,0);
      J = TPTPM::gtpmm2t(i,1);

      tpvv[i] = 2.0 * norm[i] * tpm(I,J);

   }

}

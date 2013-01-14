#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>

using std::ostream;
using std::ofstream;
using std::ifstream;
using std::vector;
using std::cout;
using std::endl;
using std::ios;

#include "include.h"

/**
 * standard constructor:
 */
TPSPM::TPSPM() : RecMat(TPTPM::gn(),SPSPM::gn()) { }

/**
 * copy constructor: constructs RecMat object
 * @param tpspm_c object that will be copied into this.
 */
TPSPM::TPSPM(const TPSPM &tpspm_c) : RecMat(tpspm_c){ }

/**
 * destructor
 */
TPSPM::~TPSPM(){ }

/**
 * access the elements of the matrix in sp mode, antisymmetry is automatically accounted for:\n\n
 * @param S spinindex of the tp part of the matrix
 * @param a first sp index that forms the tp row index i together with b, part of the TPSPM row index
 * @param b second sp index that forms the tp row index i together with a, part of the TPSPM row index
 * @param c first sp index that forms the tp column index j together with d, part of the TPSPM row index
 * @param d second sp index that forms the tp column index j together with c, part of the TPSPM row index
 * @param e first sp index that forms part of the TPSPM colum index together with b
 * @param z second sp index that forms part of the TPSPM colum index together with a
 * @return the number on place TPSPM(i,j) with the right phase.
 */
double TPSPM::operator()(int S,int a,int b,int c,int d,int e,int z) const{

   if(S == 0){

      int I = TPM::gs2t(0,a,b);
      int J = TPM::gs2t(0,c,d);

      return (*this)(TPTPM::gt2tpmm(S,I,J),SPSPM::gs2spmm(e,z));

   }
   else{

      if( (a == b) || (c == d) )
         return 0;
      else{

         int I = TPM::gs2t(1,a,b);
         int J = TPM::gs2t(1,c,d);

         int phase = 1;

         if(a > b)
            phase *= -1;
         if(c > d)
            phase *= -1;

         return phase * (*this)(TPTPM::gt2tpmm(S,I,J),SPSPM::gs2spmm(e,z));

      }

   }

}

/**
 * construct the singly-traced direct product of two TPM objects
 * @param tpm input TPM object
 */
void TPSPM::dpt(double scale,const TPM &tpm){

   int S,I,J;

   int a,b,c,d;

   int e,z;

   for(int i = 0;i < TPTPM::gn();++i){

      S = TPTPM::gtpmm2t(i,0);
      I = TPTPM::gtpmm2t(i,1);
      J = TPTPM::gtpmm2t(i,2);

      a = TPM::gt2s(S,I,0);
      b = TPM::gt2s(S,I,1);
      c = TPM::gt2s(S,J,0);
      d = TPM::gt2s(S,J,1);

      for(int j = 0;j < SPSPM::gn();++j){

         e = SPSPM::gspmm2s(j,0);
         z = SPSPM::gspmm2s(j,1);

         (*this)(i,j) = 0.0;

         for(int l = 0;l < Tools::gM();++l)
            (*this)(i,j) += (tpm(S,a,b,e,l) * tpm(S,c,d,z,l) + tpm(S,a,b,z,l) * tpm(S,c,d,e,l))/ (TPM::gnorm(e,l) * TPM::gnorm(z,l));
            
         (*this)(i,j) *= scale;

      }

   }

}

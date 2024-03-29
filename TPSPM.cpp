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

/**
 * construct the antisymmetrized singly-traced direct product of two PHM objects
 * @param phm input PHM object
 */
void TPSPM::dpt(double scale,const PHM &phm){

   int M = Tools::gM();
   int M2 = M*M;
   int M3 = M2*M;
   int M4 = M3*M;

   double *phmarray = new double [2 * M4];

   phm.convert(phmarray);

   int S,I,J;

   int sign;

   int a,b,c,d;

   int e,t;

   for(int i = 0;i < TPTPM::gn();++i){

      S = TPTPM::gtpmm2t(i,0);

      sign = 1 - 2*S;

      I = TPTPM::gtpmm2t(i,1);
      J = TPTPM::gtpmm2t(i,2);

      a = TPM::gt2s(S,I,0);
      b = TPM::gt2s(S,I,1);
      c = TPM::gt2s(S,J,0);
      d = TPM::gt2s(S,J,1);

      for(int j = 0;j < SPSPM::gn();++j){

         e = SPSPM::gspmm2s(j,0);
         t = SPSPM::gspmm2s(j,1);

         (*this)(i,j) = 0.0;

         for(int Z = 0;Z < 2;++Z){

            for(int l = 0;l < M;++l){

               double ward = phmarray[a + d*M + e*M2 + l*M3 + Z*M4] * phmarray[c + b*M + t*M2 + l*M3 + Z*M4] 

                  + phmarray[a + d*M + t*M2 + l*M3 + Z*M4] * phmarray[c + b*M + e*M2 + l*M3 + Z*M4] 

                  + sign * ( phmarray[b + d*M + e*M2 + l*M3 + Z*M4] * phmarray[c + a*M + t*M2 + l*M3 + Z*M4] 

                        + phmarray[b + d*M + t*M2 + l*M3 + Z*M4] * phmarray[c + a*M + e*M2 + l*M3 + Z*M4] )

                  + sign * ( phmarray[a + c*M + e*M2 + l*M3 + Z*M4] * phmarray[d + b*M + t*M2 + l*M3 + Z*M4] 

                        + phmarray[a + c*M + t*M2 + l*M3 + Z*M4] * phmarray[d + b*M + e*M2 + l*M3 + Z*M4] )

                  + phmarray[b + c*M + e*M2 + l*M3 + Z*M4] * phmarray[d + a*M + t*M2 + l*M3 + Z*M4] 

                  + phmarray[b + c*M + t*M2 + l*M3 + Z*M4] * phmarray[d + a*M + e*M2 + l*M3 + Z*M4];

               (*this)(i,j) += (2*Z + 1.0) * Tools::g6j(0,0,Z,S) * ward;

            }
         }

         (*this)(i,j) *= scale;

      }
   }

   delete [] phmarray;

}

/**
 * construct the antisymmetrized triply-traced direct product of two DPM objects
 * @param tpmm TPTPM
 */
void TPSPM::dpt3(double scale,const TPTPM &tpmm){

   int M = Tools::gM();

   int S,I,J;

   int e,t;

   for(int i = 0;i < TPTPM::gn();++i){

      S = TPTPM::gtpmm2t(i,0);

      I = TPTPM::gtpmm2t(i,1);
      J = TPTPM::gtpmm2t(i,2);

      for(int j = 0;j < SPSPM::gn();++j){

         e = SPSPM::gspmm2s(j,0);
         t = SPSPM::gspmm2s(j,1);

         (*this)(i,j) = 0.0;

         //first S_el = 0
         double ward = 0.0;

         for(int l = 0;l < M;++l){

            int K = TPM::gs2t(0,e,l);
            int L = TPM::gs2t(0,t,l);

            ward += tpmm(S,I,J,0,K,L)/ ( TPM::gnorm(e,l) * TPM::gnorm(t,l) );

         }

         (*this)(i,j) += ward;

         ward = 0.0;
         
         //then S_el = 1
         for(int l = 0;l < e;++l){

            int K = TPM::gs2t(1,e,l);
            int L = TPM::gs2t(1,t,l);

            ward += tpmm(S,I,J,1,K,L);

         }

         for(int l = e + 1;l < t;++l){

            int K = TPM::gs2t(1,e,l);
            int L = TPM::gs2t(1,t,l);

            ward -= tpmm(S,I,J,1,K,L);

         }

         for(int l = t + 1;l < M;++l){

            int K = TPM::gs2t(1,e,l);
            int L = TPM::gs2t(1,t,l);

            ward += tpmm(S,I,J,1,K,L);

         }

         (*this)(i,j) += 3.0 * ward;


         (*this)(i,j) *= 0.5 * scale;

      }
   }

}

/**
 * construct the antisymmetrized triply- skew traced direct product of two PPHM objects
 */
void TPSPM::dpw3(double scale,const PPHM &pphm){

   int M = Tools::gM();
   int M2 = M*M;
   int M3 = M2*M;
   int M4 = M3*M;
   int M5 = M4*M;
   int M6 = M5*M;

   double **ppharray = new double * [2];

   ppharray[0] = new double [4*M6];
   ppharray[1] = new double [M6];

   pphm.convert_st(ppharray);

   int S,I,J_index;

   int sign;

   int a,b,c,d;

   int e,t;

   for(int i = 0;i < TPTPM::gn();++i){

      S = TPTPM::gtpmm2t(i,0);

      sign = 1 - 2*S;

      I = TPTPM::gtpmm2t(i,1);
      J_index = TPTPM::gtpmm2t(i,2);

      a = TPM::gt2s(S,I,0);
      b = TPM::gt2s(S,I,1);

      c = TPM::gt2s(S,J_index,0);
      d = TPM::gt2s(S,J_index,1);

      for(int j = 0;j < SPSPM::gn();++j){

         e = SPSPM::gspmm2s(j,0);
         t = SPSPM::gspmm2s(j,1);

         (*this)(i,j) = 0.0;

         double ward;

         //first S = 1/2
         for(int J = 0;J < 2;++J){

            ward = 0.0;

            for(int k = 0;k < M;++k)
               for(int S_mn = 0;S_mn < 2;++S_mn)
                  for(int m = 0;m < M;++m)
                     for(int n = m + S_mn;n < M;++n){

                        ward += ppharray[0][m + n*M + e*M2 + k*M3 + d*M4 + a*M5 + S_mn*M6 + 2*J*M6] * ppharray[0][m + n*M + t*M2 + k*M3 + b*M4 + c*M5 + S_mn*M6 + 2*J*M6]

                           + ppharray[0][m + n*M + e*M2 + k*M3 + b*M4 + c*M5 + S_mn*M6 + 2*J*M6] * ppharray[0][m + n*M + t*M2 + k*M3 + d*M4 + a*M5 + S_mn*M6 + 2*J*M6]

                           + sign * ( ppharray[0][m + n*M + e*M2 + k*M3 + d*M4 + b*M5 + S_mn*M6 + 2*J*M6] * ppharray[0][m + n*M + t*M2 + k*M3 + a*M4 + c*M5 + S_mn*M6 + 2*J*M6]

                                 + ppharray[0][m + n*M + e*M2 + k*M3 + a*M4 + c*M5 + S_mn*M6 + 2*J*M6] * ppharray[0][m + n*M + t*M2 + k*M3 + d*M4 + b*M5 + S_mn*M6 + 2*J*M6] )

                           + sign * ( ppharray[0][m + n*M + e*M2 + k*M3 + c*M4 + a*M5 + S_mn*M6 + 2*J*M6] * ppharray[0][m + n*M + t*M2 + k*M3 + b*M4 + d*M5 + S_mn*M6 + 2*J*M6]

                                 + ppharray[0][m + n*M + e*M2 + k*M3 + b*M4 + d*M5 + S_mn*M6 + 2*J*M6] * ppharray[0][m + n*M + t*M2 + k*M3 + c*M4 + a*M5 + S_mn*M6 + 2*J*M6] )

                           + ppharray[0][m + n*M + e*M2 + k*M3 + c*M4 + b*M5 + S_mn*M6 + 2*J*M6] * ppharray[0][m + n*M + t*M2 + k*M3 + a*M4 + d*M5 + S_mn*M6 + 2*J*M6]

                           + ppharray[0][m + n*M + e*M2 + k*M3 + a*M4 + d*M5 + S_mn*M6 + 2*J*M6] * ppharray[0][m + n*M + t*M2 + k*M3 + c*M4 + b*M5 + S_mn*M6 + 2*J*M6];


                     }

            (*this)(i,j) += (2*J + 1.0) * Tools::g6j(0,0,S,J) * ward;

         }

         ward = 0.0;

         //then S = 3/2
         for(int k = 0;k < M;++k)
            for(int m = 0;m < M;++m)
               for(int n = m + 1;n < M;++n){

                  ward += ppharray[1][m + n*M + e*M2 + k*M3 + d*M4 + a*M5] * ppharray[1][m + n*M + t*M2 + k*M3 + b*M4 + c*M5]

                     + ppharray[1][m + n*M + e*M2 + k*M3 + b*M4 + c*M5] * ppharray[1][m + n*M + t*M2 + k*M3 + d*M4 + a*M5]

                     + sign * ( ppharray[1][m + n*M + e*M2 + k*M3 + d*M4 + b*M5] * ppharray[1][m + n*M + t*M2 + k*M3 + a*M4 + c*M5]

                           + ppharray[1][m + n*M + e*M2 + k*M3 + a*M4 + c*M5] * ppharray[1][m + n*M + t*M2 + k*M3 + d*M4 + b*M5] )

                     + sign * ( ppharray[1][m + n*M + e*M2 + k*M3 + c*M4 + a*M5] * ppharray[1][m + n*M + t*M2 + k*M3 + b*M4 + d*M5]

                           + ppharray[1][m + n*M + e*M2 + k*M3 + b*M4 + d*M5] * ppharray[1][m + n*M + t*M2 + k*M3 + c*M4 + a*M5] )

                     + ppharray[1][m + n*M + e*M2 + k*M3 + c*M4 + b*M5] * ppharray[1][m + n*M + t*M2 + k*M3 + a*M4 + d*M5]

                     + ppharray[1][m + n*M + e*M2 + k*M3 + a*M4 + d*M5] * ppharray[1][m + n*M + t*M2 + k*M3 + c*M4 + b*M5];

               }

         (*this)(i,j) += 6.0 * Tools::g6j(0,0,S,1) * ward;

         (*this)(i,j) *= scale;

      }
   }

   delete [] ppharray[0];
   delete [] ppharray[1];

   delete [] ppharray;


}

/**
 * construct the doubly-skew traced and once regular traced,  direct product of two PPHM objects
 * @param pphm input PPHM
 */
void TPSPM::dptw2(double scale,const PPHM &pphm){

   int M = Tools::gM();
   int M2 = M*M;
   int M3 = M2*M;
   int M4 = M3*M;
   int M5 = M4*M;
   int M6 = M5*M;

   double **ppharray = new double * [2];

   ppharray[0] = new double [4*M6];
   ppharray[1] = new double [M6];

   pphm.convert(ppharray);

   int S,I,J_index;

   int a,b,c,d;

   int e,t;

   for(int i = 0;i < TPTPM::gn();++i){

      S = TPTPM::gtpmm2t(i,0);

      I = TPTPM::gtpmm2t(i,1);
      J_index = TPTPM::gtpmm2t(i,2);

      a = TPM::gt2s(S,I,0);
      b = TPM::gt2s(S,I,1);

      c = TPM::gt2s(S,J_index,0);
      d = TPM::gt2s(S,J_index,1);

      for(int j = 0;j < SPSPM::gn();++j){

         e = SPSPM::gspmm2s(j,0);
         t = SPSPM::gspmm2s(j,1);

         (*this)(i,j) = 0.0;

         double ward = 0.0;

         //first S = 1/2
         for(int s = 0;s < M;++s)
            for(int S_kl = 0;S_kl < 2;++S_kl)
               for(int k = 0;k < M;++k)
                  for(int l = k + S_kl;l < M;++l){

                     ward += ppharray[0][a + b*M + s*M2 + k*M3 + l*M4 + e*M5 + S*M6 + 2*S_kl*M6] * ppharray[0][c + d*M + s*M2 + k*M3 + l*M4 + t*M5 + S*M6 + 2*S_kl*M6]

                        + ppharray[0][a + b*M + s*M2 + k*M3 + l*M4 + t*M5 + S*M6 + 2*S_kl*M6] * ppharray[0][c + d*M + s*M2 + k*M3 + l*M4 + e*M5 + S*M6 + 2*S_kl*M6];

                  }

         (*this)(i,j) += ward;

         
         //S = 3/2, only if:
         if(S == 1){

            ward = 0.0;

            for(int s = 0;s < M;++s)
               for(int k = 0;k < M;++k)
                  for(int l = k + 1;l < M;++l){

                     ward += ppharray[1][a + b*M + s*M2 + k*M3 + l*M4 + e*M5] * ppharray[1][c + d*M + s*M2 + k*M3 + l*M4 + t*M5]

                        + ppharray[1][a + b*M + s*M2 + k*M3 + l*M4 + t*M5] * ppharray[1][c + d*M + s*M2 + k*M3 + l*M4 + e*M5];

                  }

            (*this)(i,j) += 2.0 * ward;

         }

         (*this)(i,j) *= scale/(2.0*S + 1.0);

      }
   }

   delete [] ppharray[0];
   delete [] ppharray[1];

   delete [] ppharray;


}

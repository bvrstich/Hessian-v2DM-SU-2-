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

int ***TPTPM::t2tpmm;
vector< vector<int> > TPTPM::tpmm2t;

/**
 * initialize the static lists
 */
void TPTPM::init(){

   t2tpmm = new int ** [2];

   for(int S = 0;S < 2;++S){

      t2tpmm[S] = new int * [TPM::gdim(S)];

      for(int i = 0;i < TPM::gdim(S);++i)
      t2tpmm[S][i] = new int [TPM::gdim(S)];

   }

   vector<int> v(3);

   int tpmm = 0;

   for(int S = 0;S < 2;++S){

      for(int i = 0;i < TPM::gdim(S);++i)
         for(int j = i;j < TPM::gdim(S);++j){

            v[0] = S;
            v[1] = i;
            v[2] = j;

            tpmm2t.push_back(v);

            t2tpmm[S][i][j] = tpmm;
            t2tpmm[S][j][i] = tpmm;

            ++tpmm;

         }

   }

}

/**
 * deallocate the static lists
 */
void TPTPM::clear(){

   for(int S = 0;S < 2;++S){

      for(int i = 0;i < TPM::gdim(S);++i)
         delete [] t2tpmm[S][i];

      delete [] t2tpmm[S];

   }

   delete [] t2tpmm;

}

/**
 * standard constructor:
 */
TPTPM::TPTPM() : Matrix(tpmm2t.size()) { }

/**
 * copy constructor: constructs Matrix object of dimension M*(M - 1)/2 and fills it with the content of matrix tpmm_c
 * @param tpmm_c object that will be copied into this.
 */
TPTPM::TPTPM(const TPTPM &tpmm_c) : Matrix(tpmm_c){ }

/**
 * destructor
 */
TPTPM::~TPTPM(){ }

/**
 * access the elements of the matrix in tp mode
 * @param S spin index of the first two indices
 * @param I first tp index that forms the tpmm row index i together with J
 * @param J second tp index that forms the tpmm row index i together with I
 * @param S_ spin index of the second two indices
 * @param K first tp index that forms the tpmm column index j together with L
 * @param L second tp index that forms the tpmm column index j together with K
 * @return the number on place TPTPM(i,j)
 */
double TPTPM::operator()(int S,int I,int J,int S_,int K,int L) const{

   int i = t2tpmm[S][I][J];
   int j = t2tpmm[S_][K][L];

   return (*this)(i,j);

}

ostream &operator<<(ostream &output,const TPTPM &tpmm_p){

   int S,I,J,S_,K,L;

   int a,b,c,d;
   int e,z,t,h;

   for(int i = 0;i < TPTPM::gn();++i){

      S = tpmm_p.tpmm2t[i][0];
      I = tpmm_p.tpmm2t[i][1];
      J = tpmm_p.tpmm2t[i][2];

      a = TPM::gt2s(S,I,0);
      b = TPM::gt2s(S,I,1);

      c = TPM::gt2s(S,J,0);
      d = TPM::gt2s(S,J,1);

      for(int j = i;j < TPTPM::gn();++j){

         S_ = tpmm_p.tpmm2t[j][0]; 
         K = tpmm_p.tpmm2t[j][1];
         L = tpmm_p.tpmm2t[j][2];

         e = TPM::gt2s(S_,K,0);
         z = TPM::gt2s(S_,K,1);

         t = TPM::gt2s(S_,L,0);
         h = TPM::gt2s(S_,L,1);

         output << i << "\t" << j << "\t|\t(" << S << ")\t" << I << "\t" << J << "\t(" << S_ << ")\t" << K << "\t" << L << "\t|\t" << 

            "(" << a << "," << b << "," << c << "," << d << ")\t(" << e << "," << z << "," << t << "," << h << ")\t|\t" << tpmm_p(i,j) << endl;

      }

   }

   return output;

}

/**
 * @return the dimension of a TPTPM matrix
 */
int TPTPM::gn(){

   return tpmm2t.size();

}

/**
 * access to the lists from outside the class
 */
int TPTPM::gt2tpmm(int S,int I,int J){

   return t2tpmm[S][I][J];

}

/**
 * access to the lists from outside the class
 * @param option == 0 return S, == 1 return a, == 2 return b
 */
int TPTPM::gtpmm2t(int i,int option){

   return tpmm2t[i][option];

}

/**
 * fill a TPTPM object by antisymmetrizing a direct product of two PHM's
 */
void TPTPM::dp(const PHM &phm){

   int M = Tools::gM();
   int M2 = M*M;
   int M3 = M2*M;
   int M4 = M3*M;

   double *phmarray = new double [2 * M4];

   phm.convert(phmarray);

   int S,S_;

   int sign,sign_;

   int I,J,K,L;

   int a,b,c,d;
   int e,z,t,h;

   for(int i = 0;i < gn();++i){

      S = tpmm2t[i][0];

      sign = 1 - 2*S;

      I = tpmm2t[i][1];
      J = tpmm2t[i][2];

      a = TPM::gt2s(S,I,0);
      b = TPM::gt2s(S,I,1);

      c = TPM::gt2s(S,J,0);
      d = TPM::gt2s(S,J,1);

      for(int j = i;j < gn();++j){

         S_ = tpmm2t[j][0];

         sign_ = 1 - 2*S_;

         K = tpmm2t[j][1];
         L = tpmm2t[j][2];

         e = TPM::gt2s(S_,K,0);
         z = TPM::gt2s(S_,K,1);

         t = TPM::gt2s(S_,L,0);
         h = TPM::gt2s(S_,L,1);

         (*this)(i,j) = 0.0;

         for(int Z = 0;Z < 2;++Z){

            double ward = phmarray[a + d*M + e*M2 + h*M3 + Z*M4] * phmarray[c + b*M + t*M2 + z*M3 + Z*M4]

               + phmarray[a + d*M + t*M2 + z*M3 + Z*M4] * phmarray[c + b*M + e*M2 + h*M3 + Z*M4]

               + sign_ * (phmarray[a + d*M + z*M2 + h*M3 + Z*M4] * phmarray[c + b*M + t*M2 + e*M3 + Z*M4]

                     + phmarray[a + d*M + t*M2 + e*M3 + Z*M4] * phmarray[c + b*M + z*M2 + h*M3 + Z*M4])

               + sign_ * (phmarray[a + d*M + e*M2 + t*M3 + Z*M4] * phmarray[c + b*M + h*M2 + z*M3 + Z*M4]

                     + phmarray[a + d*M + h*M2 + z*M3 + Z*M4] * phmarray[c + b*M + e*M2 + t*M3 + Z*M4])

               + phmarray[a + d*M + z*M2 + t*M3 + Z*M4] * phmarray[c + b*M + h*M2 + e*M3 + Z*M4] 

               + phmarray[a + d*M + h*M2 + e*M3 + Z*M4] * phmarray[c + b*M + z*M2 + t*M3 + Z*M4]

               + sign * ( phmarray[b + d*M + e*M2 + h*M3 + Z*M4] * phmarray[c + a*M + t*M2 + z*M3 + Z*M4]

                     + phmarray[b + d*M + t*M2 + z*M3 + Z*M4] * phmarray[c + a*M + e*M2 + h*M3 + Z*M4]

                     + sign_ * (phmarray[b + d*M + z*M2 + h*M3 + Z*M4] * phmarray[c + a*M + t*M2 + e*M3 + Z*M4]

                        + phmarray[b + d*M + t*M2 + e*M3 + Z*M4] * phmarray[c + a*M + z*M2 + h*M3 + Z*M4])

                     + sign_ * (phmarray[b + d*M + e*M2 + t*M3 + Z*M4] * phmarray[c + a*M + h*M2 + z*M3 + Z*M4]

                        + phmarray[b + d*M + h*M2 + z*M3 + Z*M4] * phmarray[c + a*M + e*M2 + t*M3 + Z*M4])

                     + phmarray[b + d*M + z*M2 + t*M3 + Z*M4] * phmarray[c + a*M + h*M2 + e*M3 + Z*M4] 

                     + phmarray[b + d*M + h*M2 + e*M3 + Z*M4] * phmarray[c + a*M + z*M2 + t*M3 + Z*M4] )

               + sign * ( phmarray[a + c*M + e*M2 + h*M3 + Z*M4] * phmarray[d + b*M + t*M2 + z*M3 + Z*M4]

               + phmarray[a + c*M + t*M2 + z*M3 + Z*M4] * phmarray[d + b*M + e*M2 + h*M3 + Z*M4]

               + sign_ * (phmarray[a + c*M + z*M2 + h*M3 + Z*M4] * phmarray[d + b*M + t*M2 + e*M3 + Z*M4]

                     + phmarray[a + c*M + t*M2 + e*M3 + Z*M4] * phmarray[d + b*M + z*M2 + h*M3 + Z*M4])

               + sign_ * (phmarray[a + c*M + e*M2 + t*M3 + Z*M4] * phmarray[d + b*M + h*M2 + z*M3 + Z*M4]

                     + phmarray[a + c*M + h*M2 + z*M3 + Z*M4] * phmarray[d + b*M + e*M2 + t*M3 + Z*M4])

               + phmarray[a + c*M + z*M2 + t*M3 + Z*M4] * phmarray[d + b*M + h*M2 + e*M3 + Z*M4] 

               + phmarray[a + c*M + h*M2 + e*M3 + Z*M4] * phmarray[d + b*M + z*M2 + t*M3 + Z*M4] )

               + phmarray[b + c*M + e*M2 + h*M3 + Z*M4] * phmarray[d + a*M + t*M2 + z*M3 + Z*M4]

                     + phmarray[b + c*M + t*M2 + z*M3 + Z*M4] * phmarray[d + a*M + e*M2 + h*M3 + Z*M4]

                     + sign_ * (phmarray[b + c*M + z*M2 + h*M3 + Z*M4] * phmarray[d + a*M + t*M2 + e*M3 + Z*M4]

                        + phmarray[b + c*M + t*M2 + e*M3 + Z*M4] * phmarray[d + a*M + z*M2 + h*M3 + Z*M4])

                     + sign_ * (phmarray[b + c*M + e*M2 + t*M3 + Z*M4] * phmarray[d + a*M + h*M2 + z*M3 + Z*M4]

                        + phmarray[b + c*M + h*M2 + z*M3 + Z*M4] * phmarray[d + a*M + e*M2 + t*M3 + Z*M4])

                     + phmarray[b + c*M + z*M2 + t*M3 + Z*M4] * phmarray[d + a*M + h*M2 + e*M3 + Z*M4] 

                     + phmarray[b + c*M + h*M2 + e*M3 + Z*M4] * phmarray[d + a*M + z*M2 + t*M3 + Z*M4];

            (*this)(i,j) += (2*Z + 1.0) * Tools::g6j(0,0,Z,S) * Tools::g6j(0,0,Z,S_) * ward;

         }//end over Z loop

      }
   }

   delete [] phmarray;

}

/**
 * construct a TPTPM by doubly-tracing the direct product of a DPM with itsself
 */
void TPTPM::dpt2(const DPM &dpm){

   int M = Tools::gM();
   int M2 = M*M;
   int M3 = M2*M;
   int M4 = M3*M;
   int M5 = M4*M;
   int M6 = M5*M;

   double **dpmarray = new double * [2];
   
   dpmarray[0] = new double [4*M6];
   dpmarray[1] = new double [M6];

   dpm.convert(dpmarray);

   int S,S_;

   int I,J,K,L;

   int a,b,c,d;
   int e,z,t,h;

   for(int i = 0;i < gn();++i){

      S = tpmm2t[i][0];

      I = tpmm2t[i][1];
      J = tpmm2t[i][2];

      a = TPM::gt2s(S,I,0);
      b = TPM::gt2s(S,I,1);

      c = TPM::gt2s(S,J,0);
      d = TPM::gt2s(S,J,1);

      for(int j = i;j < gn();++j){

         S_ = tpmm2t[j][0];

         K = tpmm2t[j][1];
         L = tpmm2t[j][2];

         e = TPM::gt2s(S_,K,0);
         z = TPM::gt2s(S_,K,1);

         t = TPM::gt2s(S_,L,0);
         h = TPM::gt2s(S_,L,1);

         double ward = 0.0;

         //first S = 1/2
         for(int r = 0;r < M;++r)
            for(int s = 0;s < M;++s){

               ward += dpmarray[0][e + z*M + r*M2 + a*M3 + b*M4 + s*M5 + S_*M6 + 2*S*M6] * dpmarray[0][t + h*M + r*M2 + c*M3 + d*M4 + s*M5 + S_*M6 + 2*S*M6]

                  + dpmarray[0][e + z*M + r*M2 + c*M3 + d*M4 + s*M5 + S_*M6 + 2*S*M6] * dpmarray[0][t + h*M + r*M2 + a*M3 + b*M4 + s*M5 + S_*M6 + 2*S*M6];

            }

         (*this)(i,j) = 2.0/( (2*S_ + 1.0) * (2*S + 1.0) ) * ward;

         //then S = 3/2: only contributes if both S and S_ == 1 !
         if(S == 1 && S_ == 1){

            ward = 0.0;

            for(int r = 0;r < M;++r)
               for(int s = 0;s < M;++s){

                  ward += dpmarray[1][e + z*M + r*M2 + a*M3 + b*M4 + s*M5] * dpmarray[1][t + h*M + r*M2 + c*M3 + d*M4 + s*M5]

                     + dpmarray[1][e + z*M + r*M2 + c*M3 + d*M4 + s*M5] * dpmarray[1][t + h*M + r*M2 + a*M3 + b*M4 + s*M5];

               }

            (*this)(i,j) += 4.0 / 9.0 * ward;

         }

      }
   }

   this->symmetrize();

   delete [] dpmarray[0];
   delete [] dpmarray[1];

   delete [] dpmarray;

}

/**
 * construct a TPTPM by twice skew-tracing a direct product of two PPHM's
 */
void TPTPM::dpw2(const PPHM &pphm){

   int M = Tools::gM();
   int M2 = M*M;
   int M3 = M2*M;
   int M4 = M3*M;
   int M5 = M4*M;
   int M6 = M5*M;

   double **ppharray = new double * [2];

   ppharray[0] = new double [4*M6];
   ppharray[1] = new double [M6];

   pphm.convert_st2(ppharray);

   int S,S_;

   int sign,sign_;

   int I,J_index,K,L;

   int a,b,c,d;
   int e,z,t,h;

   for(int i = 0;i < gn();++i){

      S = tpmm2t[i][0];

      sign = 1 - 2*S;

      I = tpmm2t[i][1];
      J_index = tpmm2t[i][2];

      a = TPM::gt2s(S,I,0);
      b = TPM::gt2s(S,I,1);

      c = TPM::gt2s(S,J_index,0);
      d = TPM::gt2s(S,J_index,1);

      for(int j = i;j < gn();++j){

         S_ = tpmm2t[j][0];

         sign_ = 1 - 2*S_;

         K = tpmm2t[j][1];
         L = tpmm2t[j][2];

         e = TPM::gt2s(S_,K,0);
         z = TPM::gt2s(S_,K,1);

         t = TPM::gt2s(S_,L,0);
         h = TPM::gt2s(S_,L,1);

         double ward = 0.0;

         (*this)(i,j) = 0.0;

         //first S = 1/2
         for(int J = 0;J < 2;++J)
            for(int J_ = 0;J_ < 2;++J_){

               ward = 0.0;

               for(int k = 0;k < M;++k)
                  for(int m = 0;m < M;++m){

                     ward += ( ppharray[0][k + d*M + a*M2 + m*M3 + h*M4 + e*M5 + J*M6 + 2*J_*M6] * ppharray[0][k + b*M + c*M2 + m*M3 + z*M4 + t*M5 + J*M6 + 2*J_*M6]

                        + ppharray[0][k + d*M + a*M2 + m*M3 + z*M4 + t*M5 + J*M6 + 2*J_*M6] * ppharray[0][k + b*M + c*M2 + m*M3 + h*M4 + e*M5 + J*M6 + 2*J_*M6]

                        + sign_ * ( ppharray[0][k + d*M + a*M2 + m*M3 + h*M4 + z*M5 + J*M6 + 2*J_*M6] * ppharray[0][k + b*M + c*M2 + m*M3 + e*M4 + t*M5 + J*M6 + 2*J_*M6]

                        + ppharray[0][k + d*M + a*M2 + m*M3 + e*M4 + t*M5 + J*M6 + 2*J_*M6] * ppharray[0][k + b*M + c*M2 + m*M3 + h*M4 + z*M5 + J*M6 + 2*J_*M6] )

                        + sign_ * ( ppharray[0][k + d*M + a*M2 + m*M3 + t*M4 + e*M5 + J*M6 + 2*J_*M6] * ppharray[0][k + b*M + c*M2 + m*M3 + z*M4 + h*M5 + J*M6 + 2*J_*M6]

                        + ppharray[0][k + d*M + a*M2 + m*M3 + z*M4 + h*M5 + J*M6 + 2*J_*M6] * ppharray[0][k + b*M + c*M2 + m*M3 + t*M4 + e*M5 + J*M6 + 2*J_*M6] )

                        + ppharray[0][k + d*M + a*M2 + m*M3 + t*M4 + z*M5 + J*M6 + 2*J_*M6] * ppharray[0][k + b*M + c*M2 + m*M3 + e*M4 + h*M5 + J*M6 + 2*J_*M6]

                        + ppharray[0][k + d*M + a*M2 + m*M3 + e*M4 + h*M5 + J*M6 + 2*J_*M6] * ppharray[0][k + b*M + c*M2 + m*M3 + t*M4 + z*M5 + J*M6 + 2*J_*M6] 

                        /***********/

                        + sign * ( ppharray[0][k + d*M + b*M2 + m*M3 + h*M4 + e*M5 + J*M6 + 2*J_*M6] * ppharray[0][k + a*M + c*M2 + m*M3 + z*M4 + t*M5 + J*M6 + 2*J_*M6]

                        + ppharray[0][k + d*M + b*M2 + m*M3 + z*M4 + t*M5 + J*M6 + 2*J_*M6] * ppharray[0][k + a*M + c*M2 + m*M3 + h*M4 + e*M5 + J*M6 + 2*J_*M6]

                        + sign_ * ( ppharray[0][k + d*M + b*M2 + m*M3 + h*M4 + z*M5 + J*M6 + 2*J_*M6] * ppharray[0][k + a*M + c*M2 + m*M3 + e*M4 + t*M5 + J*M6 + 2*J_*M6]

                        + ppharray[0][k + d*M + b*M2 + m*M3 + e*M4 + t*M5 + J*M6 + 2*J_*M6] * ppharray[0][k + a*M + c*M2 + m*M3 + h*M4 + z*M5 + J*M6 + 2*J_*M6] )

                        + sign_ * ( ppharray[0][k + d*M + b*M2 + m*M3 + t*M4 + e*M5 + J*M6 + 2*J_*M6] * ppharray[0][k + a*M + c*M2 + m*M3 + z*M4 + h*M5 + J*M6 + 2*J_*M6]

                        + ppharray[0][k + d*M + b*M2 + m*M3 + z*M4 + h*M5 + J*M6 + 2*J_*M6] * ppharray[0][k + a*M + c*M2 + m*M3 + t*M4 + e*M5 + J*M6 + 2*J_*M6] )

                        + ppharray[0][k + d*M + b*M2 + m*M3 + t*M4 + z*M5 + J*M6 + 2*J_*M6] * ppharray[0][k + a*M + c*M2 + m*M3 + e*M4 + h*M5 + J*M6 + 2*J_*M6]

                        + ppharray[0][k + d*M + b*M2 + m*M3 + e*M4 + h*M5 + J*M6 + 2*J_*M6] * ppharray[0][k + a*M + c*M2 + m*M3 + t*M4 + z*M5 + J*M6 + 2*J_*M6] )

                        /************/

                        + sign * ( ppharray[0][k + c*M + a*M2 + m*M3 + h*M4 + e*M5 + J*M6 + 2*J_*M6] * ppharray[0][k + b*M + d*M2 + m*M3 + z*M4 + t*M5 + J*M6 + 2*J_*M6]

                        + ppharray[0][k + c*M + a*M2 + m*M3 + z*M4 + t*M5 + J*M6 + 2*J_*M6] * ppharray[0][k + b*M + d*M2 + m*M3 + h*M4 + e*M5 + J*M6 + 2*J_*M6]

                        + sign_ * ( ppharray[0][k + c*M + a*M2 + m*M3 + h*M4 + z*M5 + J*M6 + 2*J_*M6] * ppharray[0][k + b*M + d*M2 + m*M3 + e*M4 + t*M5 + J*M6 + 2*J_*M6]

                        + ppharray[0][k + c*M + a*M2 + m*M3 + e*M4 + t*M5 + J*M6 + 2*J_*M6] * ppharray[0][k + b*M + d*M2 + m*M3 + h*M4 + z*M5 + J*M6 + 2*J_*M6] )

                        + sign_ * ( ppharray[0][k + c*M + a*M2 + m*M3 + t*M4 + e*M5 + J*M6 + 2*J_*M6] * ppharray[0][k + b*M + d*M2 + m*M3 + z*M4 + h*M5 + J*M6 + 2*J_*M6]

                        + ppharray[0][k + c*M + a*M2 + m*M3 + z*M4 + h*M5 + J*M6 + 2*J_*M6] * ppharray[0][k + b*M + d*M2 + m*M3 + t*M4 + e*M5 + J*M6 + 2*J_*M6] )

                        + ppharray[0][k + c*M + a*M2 + m*M3 + t*M4 + z*M5 + J*M6 + 2*J_*M6] * ppharray[0][k + b*M + d*M2 + m*M3 + e*M4 + h*M5 + J*M6 + 2*J_*M6]

                        + ppharray[0][k + c*M + a*M2 + m*M3 + e*M4 + h*M5 + J*M6 + 2*J_*M6] * ppharray[0][k + b*M + d*M2 + m*M3 + t*M4 + z*M5 + J*M6 + 2*J_*M6] )

                        /*********/
                        
                        +  ppharray[0][k + c*M + b*M2 + m*M3 + h*M4 + e*M5 + J*M6 + 2*J_*M6] * ppharray[0][k + a*M + d*M2 + m*M3 + z*M4 + t*M5 + J*M6 + 2*J_*M6]

                        + ppharray[0][k + c*M + b*M2 + m*M3 + z*M4 + t*M5 + J*M6 + 2*J_*M6] * ppharray[0][k + a*M + d*M2 + m*M3 + h*M4 + e*M5 + J*M6 + 2*J_*M6]

                        + sign_ * ( ppharray[0][k + c*M + b*M2 + m*M3 + h*M4 + z*M5 + J*M6 + 2*J_*M6] * ppharray[0][k + a*M + d*M2 + m*M3 + e*M4 + t*M5 + J*M6 + 2*J_*M6]

                        + ppharray[0][k + c*M + b*M2 + m*M3 + e*M4 + t*M5 + J*M6 + 2*J_*M6] * ppharray[0][k + a*M + d*M2 + m*M3 + h*M4 + z*M5 + J*M6 + 2*J_*M6] )

                        + sign_ * ( ppharray[0][k + c*M + b*M2 + m*M3 + t*M4 + e*M5 + J*M6 + 2*J_*M6] * ppharray[0][k + a*M + d*M2 + m*M3 + z*M4 + h*M5 + J*M6 + 2*J_*M6]

                        + ppharray[0][k + c*M + b*M2 + m*M3 + z*M4 + h*M5 + J*M6 + 2*J_*M6] * ppharray[0][k + a*M + d*M2 + m*M3 + t*M4 + e*M5 + J*M6 + 2*J_*M6] )

                        + ppharray[0][k + c*M + b*M2 + m*M3 + t*M4 + z*M5 + J*M6 + 2*J_*M6] * ppharray[0][k + a*M + d*M2 + m*M3 + e*M4 + h*M5 + J*M6 + 2*J_*M6]

                        + ppharray[0][k + c*M + b*M2 + m*M3 + e*M4 + h*M5 + J*M6 + 2*J_*M6] * ppharray[0][k + a*M + d*M2 + m*M3 + t*M4 + z*M5 + J*M6 + 2*J_*M6] );

                  }

               (*this)(i,j) += 2.0 * (2*J_ + 1.0) * (2*J + 1.0) * Tools::g6j(0,0,S,J) * Tools::g6j(0,0,S_,J_) * ward;

            }

         //then the S = 3/2 part
         ward = 0.0;

         for(int k = 0;k < M;++k)
            for(int m = 0;m < M;++m){

               ward += ppharray[1][k + d*M + a*M2 + m*M3 + h*M4 + e*M5] * ppharray[1][k + b*M + c*M2 + m*M3 + z*M4 + t*M5]

                  + ppharray[1][k + d*M + a*M2 + m*M3 + z*M4 + t*M5] * ppharray[1][k + b*M + c*M2 + m*M3 + h*M4 + e*M5]

                  + sign_ * ( ppharray[1][k + d*M + a*M2 + m*M3 + h*M4 + z*M5] * ppharray[1][k + b*M + c*M2 + m*M3 + e*M4 + t*M5]

                        + ppharray[1][k + d*M + a*M2 + m*M3 + e*M4 + t*M5] * ppharray[1][k + b*M + c*M2 + m*M3 + h*M4 + z*M5] )

                  + sign_ * ( ppharray[1][k + d*M + a*M2 + m*M3 + t*M4 + e*M5] * ppharray[1][k + b*M + c*M2 + m*M3 + z*M4 + h*M5]

                        + ppharray[1][k + d*M + a*M2 + m*M3 + z*M4 + h*M5] * ppharray[1][k + b*M + c*M2 + m*M3 + t*M4 + e*M5] )

                  + ppharray[1][k + d*M + a*M2 + m*M3 + t*M4 + z*M5] * ppharray[1][k + b*M + c*M2 + m*M3 + e*M4 + h*M5]

                  + ppharray[1][k + d*M + a*M2 + m*M3 + e*M4 + h*M5] * ppharray[1][k + b*M + c*M2 + m*M3 + t*M4 + z*M5] 

                  /***********/

                  + sign * ( ppharray[1][k + d*M + b*M2 + m*M3 + h*M4 + e*M5] * ppharray[1][k + a*M + c*M2 + m*M3 + z*M4 + t*M5]

                        + ppharray[1][k + d*M + b*M2 + m*M3 + z*M4 + t*M5] * ppharray[1][k + a*M + c*M2 + m*M3 + h*M4 + e*M5]

                        + sign_ * ( ppharray[1][k + d*M + b*M2 + m*M3 + h*M4 + z*M5] * ppharray[1][k + a*M + c*M2 + m*M3 + e*M4 + t*M5]

                           + ppharray[1][k + d*M + b*M2 + m*M3 + e*M4 + t*M5] * ppharray[1][k + a*M + c*M2 + m*M3 + h*M4 + z*M5] )

                        + sign_ * ( ppharray[1][k + d*M + b*M2 + m*M3 + t*M4 + e*M5] * ppharray[1][k + a*M + c*M2 + m*M3 + z*M4 + h*M5]

                           + ppharray[1][k + d*M + b*M2 + m*M3 + z*M4 + h*M5] * ppharray[1][k + a*M + c*M2 + m*M3 + t*M4 + e*M5] )

                        + ppharray[1][k + d*M + b*M2 + m*M3 + t*M4 + z*M5] * ppharray[1][k + a*M + c*M2 + m*M3 + e*M4 + h*M5]

                        + ppharray[1][k + d*M + b*M2 + m*M3 + e*M4 + h*M5] * ppharray[1][k + a*M + c*M2 + m*M3 + t*M4 + z*M5] )

                  /************/

                  + sign * ( ppharray[1][k + c*M + a*M2 + m*M3 + h*M4 + e*M5] * ppharray[1][k + b*M + d*M2 + m*M3 + z*M4 + t*M5]

                        + ppharray[1][k + c*M + a*M2 + m*M3 + z*M4 + t*M5] * ppharray[1][k + b*M + d*M2 + m*M3 + h*M4 + e*M5]

                        + sign_ * ( ppharray[1][k + c*M + a*M2 + m*M3 + h*M4 + z*M5] * ppharray[1][k + b*M + d*M2 + m*M3 + e*M4 + t*M5]

                           + ppharray[1][k + c*M + a*M2 + m*M3 + e*M4 + t*M5] * ppharray[1][k + b*M + d*M2 + m*M3 + h*M4 + z*M5] )

                        + sign_ * ( ppharray[1][k + c*M + a*M2 + m*M3 + t*M4 + e*M5] * ppharray[1][k + b*M + d*M2 + m*M3 + z*M4 + h*M5]

                           + ppharray[1][k + c*M + a*M2 + m*M3 + z*M4 + h*M5] * ppharray[1][k + b*M + d*M2 + m*M3 + t*M4 + e*M5] )

                        + ppharray[1][k + c*M + a*M2 + m*M3 + t*M4 + z*M5] * ppharray[1][k + b*M + d*M2 + m*M3 + e*M4 + h*M5]

                        + ppharray[1][k + c*M + a*M2 + m*M3 + e*M4 + h*M5] * ppharray[1][k + b*M + d*M2 + m*M3 + t*M4 + z*M5] )

                  /*********/

                  +  ppharray[1][k + c*M + b*M2 + m*M3 + h*M4 + e*M5] * ppharray[1][k + a*M + d*M2 + m*M3 + z*M4 + t*M5]

                  + ppharray[1][k + c*M + b*M2 + m*M3 + z*M4 + t*M5] * ppharray[1][k + a*M + d*M2 + m*M3 + h*M4 + e*M5]

                  + sign_ * ( ppharray[1][k + c*M + b*M2 + m*M3 + h*M4 + z*M5] * ppharray[1][k + a*M + d*M2 + m*M3 + e*M4 + t*M5]

                        + ppharray[1][k + c*M + b*M2 + m*M3 + e*M4 + t*M5] * ppharray[1][k + a*M + d*M2 + m*M3 + h*M4 + z*M5] )

                  + sign_ * ( ppharray[1][k + c*M + b*M2 + m*M3 + t*M4 + e*M5] * ppharray[1][k + a*M + d*M2 + m*M3 + z*M4 + h*M5]

                        + ppharray[1][k + c*M + b*M2 + m*M3 + z*M4 + h*M5] * ppharray[1][k + a*M + d*M2 + m*M3 + t*M4 + e*M5] )

                  + ppharray[1][k + c*M + b*M2 + m*M3 + t*M4 + z*M5] * ppharray[1][k + a*M + d*M2 + m*M3 + e*M4 + h*M5]

                  + ppharray[1][k + c*M + b*M2 + m*M3 + e*M4 + h*M5] * ppharray[1][k + a*M + d*M2 + m*M3 + t*M4 + z*M5];

            }

         (*this)(i,j) += 4.0 * 9.0 * Tools::g6j(0,0,S,1) * Tools::g6j(0,0,S_,1) * ward;

      }
   }

   this->symmetrize();

   delete [] ppharray[0];
   delete [] ppharray[1];

   delete [] ppharray;

}

/**
 * construct a TPTPM by skew-tracing and regular tracing a direct product of two PPHM's: NOT SYMMETRIC!
 */
void TPTPM::dptw(const PPHM &pphm){

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

   int S,S_;

   int sign;

   int I,J_index,K,L;

   int a,b,c,d;
   int e,z,t,h;

   for(int i = 0;i < gn();++i){

      S = tpmm2t[i][0];

      sign = 1 - 2*S;

      I = tpmm2t[i][1];
      J_index = tpmm2t[i][2];

      a = TPM::gt2s(S,I,0);
      b = TPM::gt2s(S,I,1);

      c = TPM::gt2s(S,J_index,0);
      d = TPM::gt2s(S,J_index,1);

      for(int j = 0;j < gn();++j){//not SYMMETRIC!

         S_ = tpmm2t[j][0];

         K = tpmm2t[j][1];
         L = tpmm2t[j][2];

         e = TPM::gt2s(S_,K,0);
         z = TPM::gt2s(S_,K,1);

         t = TPM::gt2s(S_,L,0);
         h = TPM::gt2s(S_,L,1);

         double ward = 0.0;

         (*this)(i,j) = 0.0;

         //first S = 1/2
         for(int J = 0;J < 2;++J){

            ward = 0.0;

            for(int k = 0;k < M;++k)
               for(int s = 0;s < M;++s){

                  ward += ppharray[0][e + z*M + s*M2 + k*M3 + d*M4 + a*M5 + S_*M6 + 2*J*M6] *  ppharray[0][t + h*M + s*M2 + k*M3 + b*M4 + c*M5 + S_*M6 + 2*J*M6]

                     + ppharray[0][e + z*M + s*M2 + k*M3 + b*M4 + c*M5 + S_*M6 + 2*J*M6] *  ppharray[0][t + h*M + s*M2 + k*M3 + d*M4 + a*M5 + S_*M6 + 2*J*M6]

                     + sign * (ppharray[0][e + z*M + s*M2 + k*M3 + d*M4 + b*M5 + S_*M6 + 2*J*M6] *  ppharray[0][t + h*M + s*M2 + k*M3 + a*M4 + c*M5 + S_*M6 + 2*J*M6]

                           + ppharray[0][e + z*M + s*M2 + k*M3 + a*M4 + c*M5 + S_*M6 + 2*J*M6] *  ppharray[0][t + h*M + s*M2 + k*M3 + d*M4 + b*M5 + S_*M6 + 2*J*M6] )

                     + sign * ( ppharray[0][e + z*M + s*M2 + k*M3 + c*M4 + a*M5 + S_*M6 + 2*J*M6] *  ppharray[0][t + h*M + s*M2 + k*M3 + b*M4 + d*M5 + S_*M6 + 2*J*M6]

                           + ppharray[0][e + z*M + s*M2 + k*M3 + b*M4 + d*M5 + S_*M6 + 2*J*M6] *  ppharray[0][t + h*M + s*M2 + k*M3 + c*M4 + a*M5 + S_*M6 + 2*J*M6] )

                     + ppharray[0][e + z*M + s*M2 + k*M3 + c*M4 + b*M5 + S_*M6 + 2*J*M6] *  ppharray[0][t + h*M + s*M2 + k*M3 + a*M4 + d*M5 + S_*M6 + 2*J*M6]

                     + ppharray[0][e + z*M + s*M2 + k*M3 + a*M4 + d*M5 + S_*M6 + 2*J*M6] *  ppharray[0][t + h*M + s*M2 + k*M3 + c*M4 + b*M5 + S_*M6 + 2*J*M6];

               }

            (*this)(i,j) += 2.0 / (2*S_ + 1.0) * (2*J + 1.0) * Tools::g6j(0,0,S,J) * ward;

         }

         //only if S_ == 1 contribution from S = 3/2
         if(S_ == 1){

            ward = 0.0;

            //then S = 3/2
            for(int k = 0;k < M;++k)
               for(int s = 0;s < M;++s){

                  ward += ppharray[1][e + z*M + s*M2 + k*M3 + d*M4 + a*M5] *  ppharray[1][t + h*M + s*M2 + k*M3 + b*M4 + c*M5]

                     + ppharray[1][e + z*M + s*M2 + k*M3 + b*M4 + c*M5] *  ppharray[1][t + h*M + s*M2 + k*M3 + d*M4 + a*M5]

                     + sign * (ppharray[1][e + z*M + s*M2 + k*M3 + d*M4 + b*M5] *  ppharray[1][t + h*M + s*M2 + k*M3 + a*M4 + c*M5]

                           + ppharray[1][e + z*M + s*M2 + k*M3 + a*M4 + c*M5] *  ppharray[1][t + h*M + s*M2 + k*M3 + d*M4 + b*M5] )

                     + sign * ( ppharray[1][e + z*M + s*M2 + k*M3 + c*M4 + a*M5] *  ppharray[1][t + h*M + s*M2 + k*M3 + b*M4 + d*M5]

                           + ppharray[1][e + z*M + s*M2 + k*M3 + b*M4 + d*M5] *  ppharray[1][t + h*M + s*M2 + k*M3 + c*M4 + a*M5] )

                     + ppharray[1][e + z*M + s*M2 + k*M3 + c*M4 + b*M5] *  ppharray[1][t + h*M + s*M2 + k*M3 + a*M4 + d*M5]

                     + ppharray[1][e + z*M + s*M2 + k*M3 + a*M4 + d*M5] *  ppharray[1][t + h*M + s*M2 + k*M3 + c*M4 + b*M5];

               }

            (*this)(i,j) += 4.0 * Tools::g6j(0,0,S,1) * ward;

         }

      }
   }

   delete [] ppharray[0];
   delete [] ppharray[1];

   delete [] ppharray;

}

/**
 * construct a TPTPM by skew-tracing and regular tracing a direct product of two PPHM's: NOT SYMMETRIC!
 */
void TPTPM::dpt2(const PPHM &pphm){

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

   int S,S_;

   int I,J_index,K,L;

   int a,b,c,d;
   int e,z,t,h;

   for(int i = 0;i < gn();++i){

      S = tpmm2t[i][0];

      I = tpmm2t[i][1];
      J_index = tpmm2t[i][2];

      a = TPM::gt2s(S,I,0);
      b = TPM::gt2s(S,I,1);

      c = TPM::gt2s(S,J_index,0);
      d = TPM::gt2s(S,J_index,1);

      for(int j = i;j < gn();++j){

         S_ = tpmm2t[j][0];

         K = tpmm2t[j][1];
         L = tpmm2t[j][2];

         e = TPM::gt2s(S_,K,0);
         z = TPM::gt2s(S_,K,1);

         t = TPM::gt2s(S_,L,0);
         h = TPM::gt2s(S_,L,1);

         double ward = 0.0;

         (*this)(i,j) = 0.0;

         //first S = 1/2
         for(int s = 0;s < M;++s)
            for(int r = 0;r < M;++r){

               ward += ppharray[0][a + b*M + s*M2 + e*M3 + z*M4 + r*M5 + S*M6 + 2*S_*M6] * ppharray[0][c + d*M + s*M2 + t*M3 + h*M4 + r*M5 + S*M6 + 2*S_*M6]

                  + ppharray[0][a + b*M + s*M2 + t*M3 + h*M4 + r*M5 + S*M6 + 2*S_*M6] * ppharray[0][c + d*M + s*M2 + e*M3 + z*M4 + r*M5 + S*M6 + 2*S_*M6];

            }

         (*this)(i,j) += 2.0 / ( (2.0*S + 1.0) * (2.0*S_ + 1.0) ) * ward;

         if(S == 1 && S_ == 1){//only then S = 3/2 contribution

            ward = 0.0;

            for(int s = 0;s < M;++s)
               for(int r = 0;r < M;++r){

                  ward += ppharray[1][a + b*M + s*M2 + e*M3 + z*M4 + r*M5] * ppharray[1][c + d*M + s*M2 + t*M3 + h*M4 + r*M5]

                     + ppharray[1][a + b*M + s*M2 + t*M3 + h*M4 + r*M5] * ppharray[1][c + d*M + s*M2 + e*M3 + z*M4 + r*M5];

               }

            (*this)(i,j) += 4.0/9.0 * ward;

         }

      }
   }

   this->symmetrize();

   delete [] ppharray[0];
   delete [] ppharray[1];

   delete [] ppharray;

}

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

int **SPSPM::s2spmm;
vector< vector<int> > SPSPM::spmm2s;

/**
 * initialize the static lists
 */
void SPSPM::init(){

   int M = Tools::gM();

   s2spmm = new int * [M];

   for(int a = 0;a < M;++a)
      s2spmm[a] = new int [M];

   vector<int> v(2);

   int spmm = 0;

   for(int a = 0;a < M;++a)
      for(int b = a;b < M;++b){

         v[0] = a;
         v[1] = b;

         spmm2s.push_back(v);

         s2spmm[a][b] = spmm;
         s2spmm[b][a] = spmm;

         ++spmm;

      }

}

/**
 * deallocate the static lists
 */
void SPSPM::clear(){

   int M = Tools::gM();

   for(int a = 0;a < M;++a)
      delete [] s2spmm[a];

   delete [] s2spmm;

}

/**
 * standard constructor:
 */
SPSPM::SPSPM() : Matrix(spmm2s.size()) { }

/**
 * copy constructor: constructs Matrix object of dimension M*(M - 1)/2 and fills it with the content of matrix tpmm_c
 * @param spmm_c object that will be copied into this.
 */
SPSPM::SPSPM(const SPSPM &spmm_c) : Matrix(spmm_c){ }

/**
 * destructor
 */
SPSPM::~SPSPM(){ }

/**
 * access the elements of the matrix in sp mode
 * @param a first tp index that forms the tpmm row index i together with J
 * @param b second tp index that forms the tpmm row index i together with I
 * @param c first tp index that forms the tpmm column index j together with L
 * @param d second tp index that forms the tpmm column index j together with K
 * @return the number on place SPSPM(i,j)
 */
const double &SPSPM::operator()(int a,int b,int c,int d) const{

   int i = s2spmm[a][b];
   int j = s2spmm[c][d];

   return (*this)(i,j);

}

ostream &operator<<(ostream &output,const SPSPM &spmm_p){

   int a,b,c,d;

   for(int i = 0;i < SPSPM::gn();++i){

      a = SPSPM::gspmm2s(i,0);
      b = SPSPM::gspmm2s(i,1);

      for(int j = i;j < SPSPM::gn();++j){

         c = SPSPM::gspmm2s(j,0);
         d = SPSPM::gspmm2s(j,1);

         output << i << "\t" << j << "\t|\t" << a << "\t" << b << "\t" << c << "\t" << d << "\t|\t" << spmm_p(i,j) << endl;

      }

   }

   return output;

}

/**
 * @return the dimension of a SPSPM matrix
 */
int SPSPM::gn(){

   return spmm2s.size();

}

/**
 * access to the lists from outside the class
 */
int SPSPM::gs2spmm(int a,int b){

   return s2spmm[a][b];

}

/**
 * access to the lists from outside the class
 * @param if == 0 return a, == 1 return b
 */
int SPSPM::gspmm2s(int i,int option){

   return spmm2s[i][option];

}

/**
 * construct the doubly traced direct product of two TPM's
 */
void SPSPM::dpt2(double scale,const TPSPM &tpspm){

   int a,c,e,t;

   for(int i = 0;i < gn();++i){

      a = spmm2s[i][0];
      c = spmm2s[i][1];

      for(int j = i;j < gn();++j){

         e = spmm2s[j][0];
         t = spmm2s[j][1];

         (*this)(i,j) = 0.0;

         for(int l = 0;l < Tools::gM();++l){

            //S = 0 part
            (*this)(i,j) += tpspm(0,a,l,c,l,e,t)/( TPM::gnorm(a,l) * TPM::gnorm(c,l) );

            //S = 1 part
            (*this)(i,j) += 3.0 * tpspm(0,a,l,c,l,e,t);

         }

         (*this)(i,j) *= 0.5 * scale;

      }

   }

   this->symmetrize();

}

/**
 * construct the doubly traced direct product of two TPM's
 */
void SPSPM::dpt2(double scale,const TPM &tpm){

   int a,c,e,t;

   for(int i = 0;i < gn();++i){

      a = spmm2s[i][0];
      c = spmm2s[i][1];

      for(int j = i;j < gn();++j){

         e = spmm2s[j][0];
         t = spmm2s[j][1];

         (*this)(i,j) = 0.0;

         for(int l = 0;l < Tools::gM();++l)
            for(int k = 0;k < Tools::gM();++k){

               //S = 0 part
               (*this)(i,j) += (tpm(0,a,k,e,l) * tpm(0,c,k,t,l) + tpm(0,a,k,t,l) * tpm(0,c,k,e,l)) / ( TPM::gnorm(a,k) * TPM::gnorm(c,k) * TPM::gnorm(e,l) * TPM::gnorm(t,l) );

               //S = 1 part
               (*this)(i,j) += 3.0 * (tpm(1,a,k,e,l) * tpm(1,c,k,t,l) + tpm(1,a,k,t,l) * tpm(1,c,k,e,l));

            }

         (*this)(i,j) *= 0.5 * scale;

      }

   }

   this->symmetrize();

}

/**
 * construct the doubly traced direct product of two PHM's
 */
void SPSPM::dpt2(double scale,const PHM &phm){

   int M = Tools::gM();
   int M2 = M*M;
   int M3 = M2*M;
   int M4 = M3*M;

   double *phmarray = new double [2 * M4];

   phm.convert(phmarray);

   int a,c,e,t;

   for(int i = 0;i < gn();++i){

      a = spmm2s[i][0];
      c = spmm2s[i][1];

      for(int j = i;j < gn();++j){

         e = spmm2s[j][0];
         t = spmm2s[j][1];

         (*this)(i,j) = 0.0;

         for(int l = 0;l < Tools::gM();++l)
            for(int k = 0;k < Tools::gM();++k){

               //S = 0 part
               (*this)(i,j) += phmarray[a + k*M + e*M2 + l*M3] * phmarray[c + k*M + t*M2 + l*M3] + phmarray[a + k*M + t*M2 + l*M3] * phmarray[c + k*M + e*M2 + l*M3];

               //S = 1 part
               (*this)(i,j) += 3.0 * ( phmarray[a + k*M + e*M2 + l*M3 + M4] * phmarray[c + k*M + t*M2 + l*M3 + M4]
               
                  + phmarray[a + k*M + t*M2 + l*M3 + M4] * phmarray[c + k*M + e*M2 + l*M3 + M4] );

            }

         (*this)(i,j) *= 0.5 * scale;

      }

   }

   delete [] phmarray;

   this->symmetrize();

}

/**
 * construct the doubly traced direct product of two PHM's
 */
void SPSPM::dpt4(double scale,const TPSPM &tpspm){

   int a,c,e,t;

   for(int i = 0;i < gn();++i){

      a = spmm2s[i][0];
      c = spmm2s[i][1];

      for(int j = i;j < gn();++j){

         e = spmm2s[j][0];
         t = spmm2s[j][1];

         (*this)(i,j) = 0.0;

         for(int l = 0;l < Tools::gM();++l){

            //S = 0 part
            (*this)(i,j) += tpspm(0,a,l,c,l,e,t) / (TPM::gnorm(a,l) * TPM::gnorm(c,l)) ;

            //S = 1 part
            (*this)(i,j) += 3.0 * tpspm(1,a,l,c,l,e,t);

         }

         (*this)(i,j) *= 0.25 * scale;

      }
   }

   this->symmetrize();

}

/**
 * construct the quartly skew-traced direct product of two PPHM's
 */
void SPSPM::dpw4(double scale,const PPHM &pphm){

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

   int a,c,e,t;

   for(int i = 0;i < gn();++i){

      a = spmm2s[i][0];
      c = spmm2s[i][1];

      for(int j = i;j < gn();++j){

         e = spmm2s[j][0];
         t = spmm2s[j][1];

         (*this)(i,j) = 0.0;

         double ward = 0.0;

         //first S = 1/2
         for(int S_kl = 0;S_kl < 2;++S_kl)
            for(int k = 0;k < M;++k)
               for(int l = k + S_kl;l < M;++l)
                  for(int S_mn = 0;S_mn < 2;++S_mn)
                     for(int m = 0;m < M;++m)
                        for(int n = m + S_mn;n < M;++n){

                           ward += ppharray[0][k + l*M + a*M2 + m*M3 + n*M4 + e*M5 + S_kl*M6 + 2*S_mn*M6] *  ppharray[0][k + l*M + c*M2 + m*M3 + n*M4 + t*M5 + S_kl*M6 + 2*S_mn*M6]

                              + ppharray[0][k + l*M + a*M2 + m*M3 + n*M4 + t*M5 + S_kl*M6 + 2*S_mn*M6] *  ppharray[0][k + l*M + c*M2 + m*M3 + n*M4 + e*M5 + S_kl*M6 + 2*S_mn*M6];

                        }

         (*this)(i,j) += ward;

         ward = 0.0;

         //then S = 3/2
         for(int k = 0;k < M;++k)
            for(int l = k + 1;l < M;++l)
               for(int m = 0;m < M;++m)
                  for(int n = m + 1;n < M;++n){

                     ward += ppharray[1][k + l*M + a*M2 + m*M3 + n*M4 + e*M5] *  ppharray[1][k + l*M + c*M2 + m*M3 + n*M4 + t*M5]

                        + ppharray[1][k + l*M + a*M2 + m*M3 + n*M4 + t*M5] *  ppharray[1][k + l*M + c*M2 + m*M3 + n*M4 + e*M5];


                  }

         (*this)(i,j) += 2.0 * ward;

         (*this)(i,j) *= scale;

      }
   }

   delete [] ppharray[0];
   delete [] ppharray[1];

   delete [] ppharray;

   this->symmetrize();

}

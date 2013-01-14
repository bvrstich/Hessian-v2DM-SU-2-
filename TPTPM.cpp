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

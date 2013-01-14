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

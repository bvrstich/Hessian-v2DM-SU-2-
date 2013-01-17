#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

using std::ostream;
using std::ofstream;
using std::endl;

#include "include.h"

vector< vector<int> > *PPHM::pph2s;
int *****PPHM::s2pph;

/**
 * initialize the static lists
 */
void PPHM::init(){

   int M = Tools::gM();

   //first allocation
   pph2s = new vector< vector<int> > [2];//two total spinblocks

   s2pph = new int **** [2];//two spinblocks

   s2pph[0] = new int *** [2];//for the S = 1/2, we have that S_ab can be 0 or 1

   for(int S_ab = 0;S_ab < 2;++S_ab){//loop and allocate

      s2pph[0][S_ab] = new int ** [M];

      for(int a = 0;a < M;++a){

         s2pph[0][S_ab][a] = new int * [M];

         for(int b = 0;b < M;++b)
            s2pph[0][S_ab][a][b] = new int [M];

      }

   }

   s2pph[1] = new int *** [1];//for the S = 3/2, we have that S_ab can be only 1

   s2pph[1][0] = new int ** [M];

   for(int a = 0;a < M;++a){//loop and allocate

      s2pph[1][0][a] = new int * [M];

      for(int b = 0;b < M;++b)
         s2pph[1][0][a][b] = new int [M];

   }

   //initialize the lists
   int pph = 0;

   vector<int> v(4);

   //S = 1/2 S_ab = 0: a <= b, c
   for(int a = 0;a < M;++a)
      for(int b = a;b < M;++b)
         for(int c = 0;c < M;++c){

            v[0] = 0;//S_ab

            v[1] = a;
            v[2] = b;
            v[3] = c;

            pph2s[0].push_back(v);

            s2pph[0][0][a][b][c] = pph;
            s2pph[0][0][b][a][c] = pph;

            ++pph;

         }

   //S = 1/2, S_ab = 1, a < b ,c
   for(int a = 0;a < M;++a)
      for(int b = a + 1;b < M;++b)
         for(int c = 0;c < M;++c){

            v[0] = 1;//S_ab

            v[1] = a;
            v[2] = b;
            v[3] = c;

            pph2s[0].push_back(v);

            s2pph[0][1][a][b][c] = pph;
            s2pph[0][1][b][a][c] = pph;

            ++pph;

         }

   //re-init teller for block S = 3/2
   pph = 0;

   //S = 3/2, S_ab = 1, a < b, c
   for(int a = 0;a < M;++a)
      for(int b = a + 1;b < M;++b)
         for(int c = 0;c < M;++c){

            v[0] = 1;//S_ab

            v[1] = a;
            v[2] = b;
            v[3] = c;

            pph2s[1].push_back(v);

            s2pph[1][0][a][b][c] = pph;
            s2pph[1][0][b][a][c] = pph;

            ++pph;

         }

}

/**
 * deallocate the static lists
 */
void PPHM::clear(){

   int M = Tools::gM();

   //first delete S = 1/2 part
   for(int S_ab = 0;S_ab < 2;++S_ab){

      for(int a = 0;a < M;++a){

         for(int b = 0;b < M;++b)
            delete [] s2pph[0][S_ab][a][b];

         delete [] s2pph[0][S_ab][a];

      }

      delete [] s2pph[0][S_ab];

   }

   //then the S = 3/2 part
   for(int a = 0;a < M;++a){

      for(int b = 0;b < M;++b)
         delete [] s2pph[1][0][a][b];

      delete [] s2pph[1][0][a];

   }

   delete [] s2pph[1][0];

   for(int S = 0;S < 2;++S)
      delete [] s2pph[S];

   delete [] s2pph;

   delete [] pph2s;

}

/**
 * standard constructor: constructs BlockMatrix object with 2 blocks, for S = 1/2 and 3/2.
 */
PPHM::PPHM() : BlockMatrix(2) {

   //set the dimension and the degeneracies of the blocks
   this->setMatrixDim(0,pph2s[0].size(),2);//S=1/2 block
   this->setMatrixDim(1,pph2s[1].size(),4);//S=3/2 block

}

/**
 * copy constructor: constructs BlockMatrix object with two blocks, on for S=1/2 and one for S=3/2, and copies the content of the pphm_c blocks into it,
 * @param pphm_c PPHM object to be copied into (*this)
 */
PPHM::PPHM(const PPHM &pphm_c) : BlockMatrix(pphm_c) { }

/**
 * Destructor
 */
PPHM::~PPHM(){ }


/**
 * access the elements of the matrix in sp mode, special symmetry and antisymmetry relations are automatically accounted for:\n\n
 * @param S The block index, when == 0 then access the block S = 1/2, for block == 1 we access the S = 3/2.
 * @param S_ab The intermediate spinquantumnumber of a and b.
 * @param a first sp index that forms the pph row index i together with b, c and S_ab in block S
 * @param b second sp index that forms the pph row index i together with a, c and S_ab in block S
 * @param c third sp index that forms the pph row index i together with a, b and S_ab in block S
 * @param S_de The intermediate spinquantumnumber of d and e.
 * @param d first sp index that forms the pph column index j together with e, z and S_de in block S
 * @param e second sp index that forms the pph column index j together with d, z and S_de in block S
 * @param z third sp index that forms the pph column index j together with d, e and S_de in block S
 * @return the number on place PPHM(S,i,j) with the right phase.
 */
double PPHM::operator()(int S,int S_ab,int a,int b,int c,int S_de,int d,int e,int z) const {

   int i,j;

   int phase_i = get_inco(S,S_ab,a,b,c,i);

   if(phase_i == 0)
      return 0;

   int phase_j = get_inco(S,S_de,d,e,z,j);

   if(phase_j == 0)
      return 0;

   return phase_i*phase_j* (*this)(S,i,j);

}

/** 
 * Static member function that gets the pph-index and phase corresponding to the sp indices S,S_ab,a,b,c.
 * @param S block index of the state, 0 -> S = 1/2, 1 -> S = 3/2
 * @param S_ab intermediate spincoupling of a and b. = 0 or 1
 * @param a first sp orbital
 * @param b second sp orbital
 * @param c third sp orbital
 * @param i the corresponding pph index will be stored in this int after calling the function
 * @return the phase needed to get to a normal ordering of indices that corresponds to a pph index i
 */
int PPHM::get_inco(int S,int S_ab,int a,int b,int c,int &i){

   if(S == 0){//S = 1/2

      if(S_ab == 0){//symmetric in spatial sp's

         if(a <= b)
            i = s2pph[0][0][a][b][c];
         else
            i = s2pph[0][0][b][a][c];

         return 1;

      }
      else{//antisymmetric in spatial sp's

         if(a == b)
            return 0;

         if(a < b){

            i = s2pph[0][1][a][b][c];

            return 1;

         }
         else{

            i = s2pph[0][1][b][a][c];

            return -1;

         }

      }

   }
   else{//S = 3/2

      if(S_ab == 0)
         return 0;

      if(a == b)
         return 0;

      if(a < b){

         i = s2pph[1][0][a][b][c];

         return 1;

      }
      else{

         i = s2pph[1][0][b][a][c];

         return -1;

      }

   }

}

/**
 * The spincoupled T2 map, maps a TPM onto a PPHM object. See notes for more info
 * @param tpm Input TPM matrix
 */
void PPHM::T(const TPM &tpm){

   SPM spm;
   spm.bar(1.0/(Tools::gN() - 1.0),tpm);

   int a,b,c,d,e,z;
   int S_ab,S_de;

   double norm_ab,norm_de;
   int sign_ab,sign_de;

   //first S = 1/2
   for(int i = 0;i < gdim(0);++i){

      S_ab = pph2s[0][i][0];

      a = pph2s[0][i][1];
      b = pph2s[0][i][2];
      c = pph2s[0][i][3];

      sign_ab = 1 - 2*S_ab;

      norm_ab = 1.0;

      if(a == b)
         norm_ab /= std::sqrt(2.0);

      for(int j = i;j < gdim(0);++j){

         S_de = pph2s[0][j][0];

         d = pph2s[0][j][1];
         e = pph2s[0][j][2];
         z = pph2s[0][j][3];

         sign_de = 1 - 2*S_de;

         norm_de = 1.0;

         if(d == e)
            norm_de /= std::sqrt(2.0);


         //start the map:
         (*this)(0,i,j) = 0.0;

         //tp(1)
         if(c == z)
            if(S_ab == S_de)
               (*this)(0,i,j) += tpm(S_ab,a,b,d,e);

         if(a == d){

            //sp(1) first term
            if(b == e)
               if(S_ab == S_de)
                  (*this)(0,i,j) += norm_ab * norm_de * spm(c,z);

            //tp(2)
            double ward = 0.0;

            for(int Z = 0;Z < 2;++Z)
               ward += (2*Z + 1.0) * Tools::g9j(0,Z,S_ab,S_de) * tpm(Z,c,e,z,b);

            ward *= norm_ab * norm_de * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_de + 1.0) );

            if(c == e)
               ward *= std::sqrt(2.0);

            if(z == b)
               ward *= std::sqrt(2.0);

            (*this)(0,i,j) -= ward;

         }

         if(b == d){

            //sp(1) second term
            if(a == e)
               if(S_ab == S_de)
                  (*this)(0,i,j) += sign_ab * norm_ab * norm_de * spm(c,z);

            //tp(3)
            double ward = 0.0;

            for(int Z = 0;Z < 2;++Z)
               ward += (2*Z + 1.0) * Tools::g9j(0,Z,S_ab,S_de) * tpm(Z,c,e,z,a);

            ward *= norm_ab * norm_de * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_de + 1.0) );

            if(c == e)
               ward *= std::sqrt(2.0);

            if(z == a)
               ward *= std::sqrt(2.0);

            (*this)(0,i,j) -= sign_ab * ward;

         }

         //tp(4)
         if(a == e){

            double ward = 0.0;

            for(int Z = 0;Z < 2;++Z)
               ward += (2*Z + 1.0) * Tools::g9j(0,Z,S_ab,S_de) * tpm(Z,c,d,z,b);

            ward *= norm_ab * norm_de * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_de + 1.0) );

            if(c == d)
               ward *= std::sqrt(2.0);

            if(z == b)
               ward *= std::sqrt(2.0);

            (*this)(0,i,j) -= sign_de * ward;

         }

         //tp(5)
         if(b == e){

            double ward = 0.0;

            for(int Z = 0;Z < 2;++Z)
               ward += (2*Z + 1.0) * Tools::g9j(0,Z,S_ab,S_de) * tpm(Z,c,d,z,a);

            ward *= norm_ab * norm_de * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_de + 1.0) );

            if(c == d)
               ward *= std::sqrt(2.0);

            if(z == a)
               ward *= std::sqrt(2.0);

            (*this)(0,i,j) -= sign_ab * sign_de * ward;

         }

      }

   }

   //the easier S = 3/2 part:
   for(int i = 0;i < this->gdim(1);++i){

      a = pph2s[1][i][1];
      b = pph2s[1][i][2];
      c = pph2s[1][i][3];

      for(int j = i;j < this->gdim(1);++j){

         d = pph2s[1][j][1];
         e = pph2s[1][j][2];
         z = pph2s[1][j][3];

         //init
         (*this)(1,i,j) = 0.0;

         //tp(1)
         if(c == z)
            (*this)(1,i,j) += tpm(1,a,b,d,e);

         if(a == d){

            //sp(1)
            if(b == e)
               (*this)(1,i,j) += spm(c,z);

            //tp(2)
            double ward = 0.0;

            for(int Z = 0;Z < 2;++Z)
               ward += (2*Z + 1.0) * Tools::g6j(0,0,1,Z) * tpm(Z,c,e,z,b);

            if(c == e)
               ward *= std::sqrt(2.0);

            if(z == b)
               ward *= std::sqrt(2.0);

            (*this)(1,i,j) -= ward;

         }

         //tp(3)
         if(b == d){

            double ward = 0.0;

            for(int Z = 0;Z < 2;++Z)
               ward += (2*Z + 1.0) * Tools::g6j(0,0,1,Z) * tpm(Z,c,e,z,a);

            if(c == e)
               ward *= std::sqrt(2.0);

            if(z == a)
               ward *= std::sqrt(2.0);

            (*this)(1,i,j) += ward;

         }

         //tp(5)
         if(b == e){

            double ward = 0.0;

            for(int Z = 0;Z < 2;++Z)
               ward += (2*Z + 1.0) * Tools::g6j(0,0,1,Z) * tpm(Z,c,d,z,a);

            if(c == d)
               ward *= std::sqrt(2.0);

            if(z == a)
               ward *= std::sqrt(2.0);

            (*this)(1,i,j) -= ward;

         }

      }

   }

   this->symmetrize();

}

ostream &operator<<(ostream &output,const PPHM &pphm_p){

   for(int S = 0;S < pphm_p.gnr();++S){

      output << S << "\t" << pphm_p.gdim(S) << "\t" << pphm_p.gdeg(S) << std::endl;
      output << std::endl;

      for(int i = 0;i < pphm_p.gdim(S);++i)
         for(int j = 0;j < pphm_p.gdim(S);++j){

            output << S << "\t" << i << "\t" << j << "\t|\t" << 

               pphm_p.pph2s[S][i][0] << "\t" << pphm_p.pph2s[S][i][1] << "\t" << pphm_p.pph2s[S][i][2] << "\t" << pphm_p.pph2s[S][i][3] << 

               "\t" << pphm_p.pph2s[S][j][0] << "\t" << pphm_p.pph2s[S][j][1] << "\t" << pphm_p.pph2s[S][j][2] << "\t" << pphm_p.pph2s[S][j][3] 

               << "\t" << pphm_p(S,i,j) << endl;

         }

      output << endl;

   }

   return output;

}

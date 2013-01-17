#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <cmath>

using std::ostream;
using std::ofstream;
using std::cout;
using std::endl;
using std::vector;

#include "include.h"

vector< vector<int> > PHM::ph2s;
int **PHM::s2ph;

/**
 * initialize the statics
 */
void PHM::init(){

   int M = Tools::gM();

   s2ph = new int * [M];
   s2ph[0] = new int [M*M];

   for(int i = 1;i < M;++i)
      s2ph[i] = s2ph[i - 1] + M;

   //initialisation of the two arrays
   int ph = 0;

   vector<int> v(2);

   for(int a = 0;a < M;++a)
      for(int b = 0;b < M;++b){

         v[0] = a;
         v[1] = b;

         ph2s.push_back(v);

         s2ph[a][b] = ph;

         ++ph;

      }
}

/** 
 * deallocate the statics
 */
void PHM::clear(){

   delete [] s2ph[0];
   delete [] s2ph;

}

/**
 * standard constructor: constructs BlockMatrix object with 2 blocks, for S = 0 and 1 of dimension M*M/4.
 */
PHM::PHM() : BlockMatrix(2) {

   //set the dimension of the blocks
   this->setMatrixDim(0,ph2s.size(),1);
   this->setMatrixDim(1,ph2s.size(),3);

}

/**
 * copy constructor: constructs BlockMatrix object with two blocks of dimension M*M/4 and copies the content of phm_c into it,
 * @param phm_c PHM to be copied into (*this)
 */
PHM::PHM(const PHM &phm_c) : BlockMatrix(phm_c){ }

/**
 * destructor: if counter == 1 the memory for the static lists ph2s en s2ph will be deleted.
 */
PHM::~PHM(){ }

/**
 * access the elements of the matrix in sp mode, 
 * @param S The spin of the block you want to access
 * @param a first sp index that forms the ph row index i in block S together with b
 * @param b second sp index that forms the ph row index i in block S together with a
 * @param c first sp index that forms the ph column index j in block S together with d
 * @param d second sp index that forms the ph column index j in block S together with c
 * @return the number on place PHM(i,j)
 */
double &PHM::operator()(int S,int a,int b,int c,int d){

   int i = s2ph[a][b];
   int j = s2ph[c][d];

   return (*this)(S,i,j);

}

/**
 * access the elements of the matrix in sp mode, const version
 * @param S The spin of the block you want to access
 * @param a first sp index that forms the ph row index i in block S together with b
 * @param b second sp index that forms the ph row index i in block S together with a
 * @param c first sp index that forms the ph column index j in block S together with d
 * @param d second sp index that forms the ph column index j in block S together with c
 * @return the number on place PHM(i,j)
 */
double PHM::operator()(int S,int a,int b,int c,int d) const{

   int i = s2ph[a][b];
   int j = s2ph[c][d];

   return (*this)(S,i,j);

}

ostream &operator<<(ostream &output,const PHM &phm_p){

   for(int S = 0;S < phm_p.gnr();++S){

      output << S << "\t" << phm_p.gdim(S) << "\t" << phm_p.gdeg(S) << std::endl;
      output << std::endl;

      for(int i = 0;i < phm_p.gdim(S);++i)
         for(int j = 0;j < phm_p.gdim(S);++j){

            output << S << "\t" << i << "\t" << j << "\t|\t" << phm_p.ph2s[i][0] << "\t" << phm_p.ph2s[i][1]

               << "\t" << phm_p.ph2s[j][0] << "\t" << phm_p.ph2s[j][1] << "\t" << phm_p(S,i,j) << endl;

         }

      output << endl;

   }

   return output;

}

/**
 * De G map, maps a TPM object on a PHM object.
 * @param tpm input TPM
 */
void PHM::G(const TPM &tpm){

   //construct the SPM corresponding to the TPM
   SPM spm;
   spm.bar(1.0/(Tools::gN() - 1.0),tpm);

   int a,b,c,d;

   for(int S = 0;S < 2;++S){

      for(int i = 0;i < gdim(S);++i){

         a = ph2s[i][0];
         b = ph2s[i][1];

         for(int j = i;j < gdim(S);++j){

            c = ph2s[j][0];
            d = ph2s[j][1];

            //tp part
            (*this)(S,i,j) = - Tools::g6j(0,0,S,0)*tpm(0,a,d,c,b)/(TPM::gnorm(a,d)*TPM::gnorm(c,b)) - 3.0*Tools::g6j(0,0,S,1)*tpm(1,a,d,c,b);

            //sp part
            if(b == d)
               (*this)(S,i,j) += spm(a,c);

         }
      }

   }

   this->symmetrize();

}

/**
 * access to the lists from outside of the class, same for S = 0/1
 */
int PHM::gs2ph(int a,int b){

   return s2ph[a][b];

}

/**
 * access to the lists from outside of the class, same for S = 0/1
 */
int PHM::gph2s(int i,int option){

   return s2ph[i][option];

}

/**
 * convert a PHM to an array, for fast access
 */
void PHM::convert(double *array) const {

   int M = Tools::gM();
   int M2 = M * M;
   int M3 = M2 * M;
   int M4 = M3 * M;

   for(int S = 0;S < 2;++S)
      for(int a = 0;a < M;++a)
         for(int b = 0;b < M;++b)
            for(int c = 0;c < M;++c)
               for(int d = 0;d < M;++d)
                  array[a + b*M + c*M2 + d*M3 + S*M4] = (*this)(S,a,b,c,d);

}

/**
 * The bar function that maps a PPHM object onto a PHM object by tracing away the first pair of incdices of the PPHM
 * @param scale factor to scale the result with
 * @param pphm Input PPHM object
 */
void PHM::bar(double scale,const PPHM &pphm){

   int a,b,c,d;

   double ward,hard;

   for(int S = 0;S < 2;++S){//loop over spinblocks of PHM

      for(int i = 0;i < this->gdim(S);++i){

         a = ph2s[i][0];
         b = ph2s[i][1];

         for(int j = i;j < this->gdim(S);++j){

            c = ph2s[j][0];
            d = ph2s[j][1];

            (*this)(S,i,j) = 0.0;

            //first the S = 1/2 block of the PPHM matrix
            for(int S_ab = 0;S_ab < 2;++S_ab)
               for(int S_de = 0;S_de < 2;++S_de){

                  ward = 2.0 * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_de + 1.0) ) * Tools::g6j(0,0,S,S_ab) * Tools::g6j(0,0,S,S_de);

                  for(int l = 0;l < Tools::gM();++l){

                     hard = ward * pphm(0,S_ab,l,a,b,S_de,l,c,d);

                     //norms
                     if(l == a)
                        hard *= std::sqrt(2.0);

                     if(l == c)
                        hard *= std::sqrt(2.0);

                     (*this)(S,i,j) += hard;

                  }

               }

            //then the S = 3/2 block
            if(S == 1)
               for(int l = 0;l < Tools::gM();++l)
                  (*this)(S,i,j) += 4.0/3.0 * pphm(1,1,l,a,b,1,l,c,d);

            (*this)(S,i,j) *= scale;

         }
      }

   }

   this->symmetrize();

}

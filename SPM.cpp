#include <iostream>
#include <fstream>

using std::ostream;
using std::endl;

#include "include.h"

/**
 * constructor, makes matrix of dimension M
 */
SPM::SPM() : Matrix(Tools::gM()) { }

/**
 * copy constructor
 * @param spm_copy content of this matrix will be copied into the constructed matrix
 */
SPM::SPM(const SPM &spm_copy) : Matrix(spm_copy) { }

/**
 * destructor
 */
SPM::~SPM(){ }

ostream &operator<<(ostream &output,const SPM &spm_p){

   for(int i = 0;i < Tools::gM();++i)
      for(int j = 0;j < Tools::gM();++j)
         output << i << "\t" << j << "\t" << spm_p(i,j) << endl;

   return output;

}

/**
 * Trace out a set of indices to create the "bar" matrix of a TPM
 * @param scale the factor u want the SPM to be scaled with (1/N-1 for normal sp density matrix)
 * @param tpm the TPM out of which the SPM will be filled
 */
void SPM::bar(double scale,const TPM &tpm){

   for(int a = 0;a < Tools::gM();++a)
      for(int c = a;c < Tools::gM();++c){

         (*this)(a,c) = 0.0;

         for(int b = 0;b < Tools::gM();++b){

            //S = 0 stuk
            (*this)(a,c) += tpm(0,a,b,c,b)/(TPM::gnorm(a,b)*TPM::gnorm(c,b));

            //S = 1 stuk: hier kan nooit a = b en c = d wegens antisymmetrie
            (*this)(a,c) += 3.0*tpm(1,a,b,c,b);

         }

         //nog schalen
         (*this)(a,c) *= 0.5*scale;

      }

   this->symmetrize();

}

/**
 * Trace out a set of indices to create the "bar" matrix of a PHM, slight difference from the bar(TPM) function (normalization of the tp basisset).
 * @param scale the factor u want the SPM to be scaled with
 * @param phm the PHM out of which the SPM will be filled
 */
void SPM::bar(double scale,const PHM &phm){

   int M = Tools::gM();

   for(int a = 0;a < M;++a)
      for(int c = a;c < M;++c){

         (*this)(a,c) = 0.0;

         for(int b = 0;b < M;++b){

            //S = 0 stuk
            for(int S = 0;S < 2;++S)
               (*this)(a,c) += phm.gdeg(S)*phm(S,a,b,c,b);

         }

         //nog schalen
         (*this)(a,c) *= 0.5*scale;

      }

   this->symmetrize();

}

/** 
 * This bar function maps a PPHM object directly onto a SPM object, scaling it with a factor scale
 * @param scale the scalefactor
 * @param pphm Input PPHM object
 */
void SPM::bar(double scale,const PPHM &pphm){

   int M = Tools::gM();

   for(int a = 0;a < M;++a)
      for(int c = a;c < M;++c){

         (*this)(a,c) = 0.0;

         //first S = 1/2 part
         for(int S_lk = 0;S_lk < 2;++S_lk){

            for(int l = 0;l < M;++l){

               for(int k = 0;k < l;++k)//k < l
                  (*this)(a,c) += pphm(0,S_lk,l,k,a,S_lk,l,k,c);

               //k == l norm correction
               (*this)(a,c) += 2.0 * pphm(0,S_lk,l,l,a,S_lk,l,l,c);

               for(int k = l + 1;k < M;++k)//k > l
                  (*this)(a,c) += pphm(0,S_lk,l,k,a,S_lk,l,k,c);

            }
         }
         
         //then S = 3/2 part:
         for(int l = 0;l < M;++l)
            for(int k = 0;k < M;++k)
               (*this)(a,c) += 2.0 * pphm(1,1,l,k,a,1,l,k,c);

         //scaling
         (*this)(a,c) *= scale;

      }

   this->symmetrize();

}

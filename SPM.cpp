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

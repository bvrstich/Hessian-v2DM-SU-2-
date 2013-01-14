#include <iostream>
#include <cmath>
#include <fstream>

using std::ostream;
using std::ofstream;
using std::ifstream;
using std::cout;
using std::endl;
using std::ios;

#include "include.h"

/**
 * standard constructor:
 */
Gradient::Gradient() {

   gradient = new double [TPTPM::gn() + 1];
   
}

/**
 * copy constructor
 * @param gradient_c object that will be copied into this.
 */
Gradient::Gradient(const Gradient &gradient_c){

   int n = TPTPM::gn() + 1;

   gradient = new double [n];

   int inc = 1;

   dcopy_(&n,gradient_c.gpointer(),&inc,gradient,&inc);
   
}

/**
 * destructor
 */
Gradient::~Gradient(){ 

   delete [] gradient;
   
}

/**
 * @return the pointer to the vector, for blas and lapack stuff, const version
 */
const double *Gradient::gpointer() const{

   return gradient;

}

/**
 * @return the pointer to the vector, for blas and lapack stuff
 */
double *Gradient::gpointer(){

   return gradient;

}



/** 
 * access to the numbers in the vector, const
 */
const double &Gradient::operator[](int i) const{

   return gradient[i];

}

/** 
 * access to the numbers in the vector, no const
 */
double &Gradient::operator[](int i){

   return gradient[i];

}

/**
 * construct 'minus' the gradient vector
 * @param t potential scaling parameter
 * @param ham hamiltonian object
 * @param P object containing the inverse of the matrix constraints p and q.
 */
void Gradient::construct(double t,const TPM &ham,const SUP &P){

   int I,J;

#ifdef __Q_CON
   TPM QQ;
   QQ.Q(1,P.gQ());
#endif

#ifdef __G_CON
   TPM GG;
   GG.G(1,P.gG());
#endif

#ifdef __T1_CON
   TPM TT1;
   TT1.T(1,P.gT1());
#endif

#ifdef __T2_CON
   TPM TT2;
   TT2.T(P.gT2());
#endif

   for(int i = 0;i < TPTPM::gn();++i){

      I = TPTPM::gtpmm2t(i,0);
      J = TPTPM::gtpmm2t(i,1);

      gradient[i] = t * P.gI()(I,J) - ham(I,J);

#ifdef __Q_CON
      gradient[i] += t * QQ(I,J);
#endif

#ifdef __G_CON
      gradient[i] += t * GG(I,J);
#endif

#ifdef __T1_CON
      gradient[i] += t * TT1(I,J);
#endif

#ifdef __T2_CON
      gradient[i] += t * TT2(I,J);
#endif

      gradient[i] *= 2.0 * TPTPV::gnorm(i);

   }

   //last part of right-hand side (lagrange multiplier)
   gradient[TPTPM::gn()] = 0.0;

}

/**
 * convert a 2DM/TPM to vector form
 * @param tpm the 2DM in question
 */
void Gradient::convert(const TPM &tpm){

   int I,J;

   for(int i = 0;i < TPTPM::gn();++i){

      I = TPTPM::gtpmm2t(i,0);
      J = TPTPM::gtpmm2t(i,1);

      gradient[i] = 2.0 * TPTPV::gnorm(i) * tpm(I,J);

   }

}

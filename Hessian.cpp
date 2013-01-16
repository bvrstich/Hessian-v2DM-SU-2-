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
Hessian::Hessian() : Matrix(TPTPM::gn() + 1) { }

/**
 * copy constructor: constructs Matrix object of dimension M*(M - 1)/2 and fills it with the content of matrix hess_c
 * @param hess_c object that will be copied into this.
 */
Hessian::Hessian(const Hessian &hess_c) : Matrix(hess_c){ }

/**
 * destructor
 */
Hessian::~Hessian(){ }

ostream &operator<<(ostream &output,const Hessian &hess_p){

   int S,I,J,S_,K,L;

   int a,b,c,d;
   int e,z,t,h;

   for(int i = 0;i < TPTPM::gn();++i){

      S = TPTPM::gtpmm2t(i,0);
      I = TPTPM::gtpmm2t(i,1);
      J = TPTPM::gtpmm2t(i,2);

      a = TPM::gt2s(S,I,0);
      b = TPM::gt2s(S,I,1);

      c = TPM::gt2s(S,J,0);
      d = TPM::gt2s(S,J,1);

      for(int j = i;j < TPTPM::gn();++j){

         S_ = TPTPM::gtpmm2t(j,0); 
         K = TPTPM::gtpmm2t(j,1);
         L = TPTPM::gtpmm2t(j,2);

         e = TPM::gt2s(S_,K,0);
         z = TPM::gt2s(S_,K,1);

         t = TPM::gt2s(S_,L,0);
         h = TPM::gt2s(S_,L,1);

         output << i << "\t" << j << "\t|\t(" << S << ")\t" << I << "\t" << J << "\t(" << S_ << ")\t" << K << "\t" << L << "\t|\t" << 

            "(" << a << "," << b << "," << c << "," << d << ")\t(" << e << "," << z << "," << t << "," << h << ")\t|\t" << hess_p(i,j) << endl;

      }

   }

   return output;

}

/**
 * construct the I part of the hessian matrix
 */
void Hessian::I(const TPM &tpm){

   int S,I,J,S_,K,L;

   for(int i = 0;i < TPTPM::gn();++i){

      S = TPTPM::gtpmm2t(i,0);
      I = TPTPM::gtpmm2t(i,1);
      J = TPTPM::gtpmm2t(i,2);

      for(int j = i;j < TPTPM::gn();++j){

         S_ = TPTPM::gtpmm2t(j,0);
         K = TPTPM::gtpmm2t(j,1);
         L = TPTPM::gtpmm2t(j,2);

         if(S == S_)
            (*this)(i,j) = 2.0 * (2.0*S + 1.0) * Gradient::gnorm(i) * Gradient::gnorm(j) * ( tpm(S,I,K) * tpm(S,J,L) + tpm(S,I,L) * tpm(S,J,K) );
         else
            (*this)(i,j) = 0.0;

      }

   }

}

/**
 * construct the lagrange multiplier part of the Hessian
 */
void Hessian::lagr(){

   int S,I,J;

   for(int i = 0;i < TPTPM::gn();++i){

      S = TPTPM::gtpmm2t(i,0);
      I = TPTPM::gtpmm2t(i,1);
      J = TPTPM::gtpmm2t(i,2);

      if(I == J)
         (*this)(i,TPTPM::gn()) = std::sqrt(2.0*S + 1.0);
      else
         (*this)(i,TPTPM::gn()) = 0.0;

   }

}

/**
 * construct the Q part of the hessian
 */
void Hessian::Q(const TPM &Q){

   int N = Tools::gN();

   TPM Q2;
   Q2.squaresym(Q);

   SPM Q2bar;
   Q2bar.bar(8.0/(N*(N - 1.0)*(N - 1.0)),Q2);

   double Q2trace = 16 * Q2.trace()/ (N*N*(N - 1.0)*(N - 1.0));

   TPSPM dpt;
   dpt.dpt(1.0/(N - 1.0),Q);

   SPSPM dpt2;
   dpt2.dpt2(1.0/((N - 1.0)*(N - 1.0)),Q);

   int S,I,J,S_,K,L;

   //first store everything in ward, then multiply with norms and add to (*this)!
   double ward;

   int a,b,c,d;
   int e,z,t,h;

   int sign,sign_;

   for(int i = 0;i < TPTPM::gn();++i){

      S = TPTPM::gtpmm2t(i,0);

      sign = 1 - 2*S;

      I = TPTPM::gtpmm2t(i,1);
      J = TPTPM::gtpmm2t(i,2);

      a = TPM::gt2s(S,I,0);
      b = TPM::gt2s(S,I,1);
      c = TPM::gt2s(S,J,0);
      d = TPM::gt2s(S,J,1);

      for(int j = i;j < TPTPM::gn();++j){

         S_ = TPTPM::gtpmm2t(j,0);

         sign_ = 1 - 2*S_;

         K = TPTPM::gtpmm2t(j,1);
         L = TPTPM::gtpmm2t(j,2);

         e = TPM::gt2s(S_,K,0);
         z = TPM::gt2s(S_,K,1);
         t = TPM::gt2s(S_,L,0);
         h = TPM::gt2s(S_,L,1);

         ward = 0.0;

         if(S == S_)
            ward += 2.0/(2.0*S + 1.0) * ( Q(S,I,K) * Q(S,J,L) + Q(S,I,L) * Q(S,J,K) );

         if(I == J){

            if(K == L)
               ward += Q2trace;

            ward += 8.0/(N*(N - 1.0)) * Q2(S_,K,L);

            if(z == h)
               ward -= TPM::gnorm(e,z) * TPM::gnorm(t,h) * Q2bar(e,t);

            if(e == h)
               ward -= sign_ * TPM::gnorm(e,z) * TPM::gnorm(t,h) * Q2bar(z,t);

            if(z == t)
               ward -= sign_ * TPM::gnorm(e,z) * TPM::gnorm(t,h) * Q2bar(e,h);

            if(e == t)
               ward -= TPM::gnorm(e,z) * TPM::gnorm(t,h) * Q2bar(z,h);

         }

         if(K == L){

            ward +=  8.0/(N*(N - 1.0)) * Q2(S,I,J);

            if(b == d)
               ward -= TPM::gnorm(a,b) * TPM::gnorm(c,d) * Q2bar(a,c);

            if(a == d)
               ward -= sign * TPM::gnorm(a,b) * TPM::gnorm(c,d) * Q2bar(b,c);

            if(b == c)
               ward -= sign * TPM::gnorm(a,b) * TPM::gnorm(c,d) * Q2bar(a,d);

            if(a == c)
               ward -= TPM::gnorm(a,b) * TPM::gnorm(c,d) * Q2bar(b,d);

         }

         if(b == d){

            ward -= TPM::gnorm(a,b) * TPM::gnorm(c,d) * dpt(S_,e,z,t,h,a,c);

            if(z == h)
               ward += TPM::gnorm(a,b) * TPM::gnorm(c,d) * TPM::gnorm(e,z) * TPM::gnorm(t,h) * dpt2(a,c,e,t);

            if(e == h)
               ward += sign_ * TPM::gnorm(a,b) * TPM::gnorm(c,d) * TPM::gnorm(e,z) * TPM::gnorm(t,h) * dpt2(a,c,z,t);

            if(z == t)
               ward += sign_ * TPM::gnorm(a,b) * TPM::gnorm(c,d) * TPM::gnorm(e,z) * TPM::gnorm(t,h) * dpt2(a,c,e,h);

            if(e == t)
               ward += TPM::gnorm(a,b) * TPM::gnorm(c,d) * TPM::gnorm(e,z) * TPM::gnorm(t,h) * dpt2(a,c,z,h);

         }

         if(a == d){

            ward -= sign * TPM::gnorm(a,b) * TPM::gnorm(c,d) * dpt(S_,e,z,t,h,b,c);

            if(z == h)
               ward += sign * TPM::gnorm(a,b) * TPM::gnorm(c,d) * TPM::gnorm(e,z) * TPM::gnorm(t,h) * dpt2(b,c,e,t);

            if(e == h)
               ward += sign * sign_ * TPM::gnorm(a,b) * TPM::gnorm(c,d) * TPM::gnorm(e,z) * TPM::gnorm(t,h) * dpt2(b,c,z,t);

            if(z == t)
               ward += sign * sign_ * TPM::gnorm(a,b) * TPM::gnorm(c,d) * TPM::gnorm(e,z) * TPM::gnorm(t,h) * dpt2(b,c,e,h);

            if(e == t)
               ward += sign * TPM::gnorm(a,b) * TPM::gnorm(c,d) * TPM::gnorm(e,z) * TPM::gnorm(t,h) * dpt2(b,c,z,h);

         }

         if(b == c){

            ward -= sign * TPM::gnorm(a,b) * TPM::gnorm(c,d) * dpt(S_,e,z,t,h,a,d);

            if(z == h)
               ward += sign * TPM::gnorm(a,b) * TPM::gnorm(c,d) * TPM::gnorm(e,z) * TPM::gnorm(t,h) * dpt2(a,d,e,t);

            if(e == h)
               ward += sign * sign_ * TPM::gnorm(a,b) * TPM::gnorm(c,d) * TPM::gnorm(e,z) * TPM::gnorm(t,h) * dpt2(a,d,z,t);

            if(z == t)
               ward += sign * sign_ * TPM::gnorm(a,b) * TPM::gnorm(c,d) * TPM::gnorm(e,z) * TPM::gnorm(t,h) * dpt2(a,d,e,h);

            if(e == t)
               ward += sign * TPM::gnorm(a,b) * TPM::gnorm(c,d) * TPM::gnorm(e,z) * TPM::gnorm(t,h) * dpt2(a,d,z,h);

         }

         if(a == c){

            ward -= TPM::gnorm(a,b) * TPM::gnorm(c,d) * dpt(S_,e,z,t,h,b,d);

            if(z == h)
               ward += TPM::gnorm(a,b) * TPM::gnorm(c,d) * TPM::gnorm(e,z) * TPM::gnorm(t,h) * dpt2(b,d,e,t);

            if(e == h)
               ward += sign_ * TPM::gnorm(a,b) * TPM::gnorm(c,d) * TPM::gnorm(e,z) * TPM::gnorm(t,h) * dpt2(b,d,z,t);

            if(z == t)
               ward += sign_ * TPM::gnorm(a,b) * TPM::gnorm(c,d) * TPM::gnorm(e,z) * TPM::gnorm(t,h) * dpt2(b,d,e,h);

            if(e == t)
               ward += TPM::gnorm(a,b) * TPM::gnorm(c,d) * TPM::gnorm(e,z) * TPM::gnorm(t,h) * dpt2(b,d,z,h);

         }

         if(z == h)
            ward -= TPM::gnorm(e,z) * TPM::gnorm(t,h) * dpt(S,a,b,c,d,e,t);

         if(e == h)
            ward -= sign_ * TPM::gnorm(e,z) * TPM::gnorm(t,h) * dpt(S,a,b,c,d,z,t);

         if(z == t)
            ward -= sign_ * TPM::gnorm(e,z) * TPM::gnorm(t,h) * dpt(S,a,b,c,d,e,h);

         if(e == t)
            ward -= TPM::gnorm(e,z) * TPM::gnorm(t,h) * dpt(S,a,b,c,d,z,h);

         //finally
         (*this)(i,j) += ward * Gradient::gnorm(i) * Gradient::gnorm(j) * (2.0*S + 1.0) * (2.0*S_ + 1.0);

      }
   }

}

/**
 * construct the G part of the Hessian
 */
void Hessian::G(const PHM &G){

}

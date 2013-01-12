#include <iostream>
#include <fstream>
#include <cmath>

using std::cout;
using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;

#include "include.h"

/**
 * standard constructor
 */
Newton::Newton(){

   H = new Hessian();

   gradient = new Gradient();

}

/**
 * copy constructor
 */
Newton::Newton(const Newton &newton_c){

   H = new Hessian(newton_c.gH());

   gradient = new Gradient(newton_c.gGradient());

}

/**
 * destructor
 */
Newton::~Newton(){

   delete H;

   delete gradient;

}

/**
 * @return the full Hessian matrix, read only
 */
const Hessian &Newton::gH() const {

   return *H;

}

/**
 * @return the full Gradient, read only
 */
const Gradient &Newton::gGradient() const {

   return *gradient;

}

/**
 * construct the different parts of the linear system
 * @param t potential scaling parameter
 * @param ham hamiltonian object
 * @param P object containing the inverse of the matrix constraints p and q.
 */
void Newton::construct(double t,const TPM &ham,const SUP &P){

   //first construct the gradient
   gradient(t,ham,P);

   //construct the p part of the hessian
   H->I(P.gI());

#ifdef __Q_CON
   H->Q(P.gQ());
#endif

#ifdef __G_CON
   H->G(P.gG());
#endif

#ifdef __T1_CON
   H->T(P.gT1());
#endif

#ifdef __T2_CON
   H->T(P.gT2());
#endif

   H->dscal(t);

   //the constraint/lagrange multiplier part of the Hessian
   H->lagr();

   //last element zero!
   (*H)(TPTPM::gn(),TPTPM::gn()) = 0.0;
   
   H->symmetrize();

   //and last but not least, solve the system
   H->solve_sy(x);

}

/**
 * construct 'minus' the gradient vector
 * @param t potential scaling parameter
 * @param ham hamiltonian object
 * @param P object containing the inverse of the matrix constraints p and q.
 */
void Newton::gradient(double t,const TPM &ham,const SUP &P){

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

      x[i] = t * P.gI()(I,J) - ham(I,J);

#ifdef __Q_CON
      x[i] += t * QQ(I,J);
#endif

#ifdef __G_CON
      x[i] += t * GG(I,J);
#endif

#ifdef __T1_CON
      x[i] += t * TT1(I,J);
#endif

#ifdef __T2_CON
      x[i] += t * TT2(I,J);
#endif

      x[i] *= 2.0 * TPTPV::gnorm(i);

   }

   //last part of right-hand side (lagrange multiplier)
   x[TPTPM::gn()] = 0.0;

}

/**
 * access to the x-elements, the solution, as it were.
 * elements are already transformed to matrix mode
 */
double Newton::gx(int I,int J) const {

   int i = TPTPM::gt2tpmm(I,J);

   return x[i];

}

/**
 * access to the x-elements, the solution, as it were in sp mode
 * elements are already transformed to matrix mode
 */
double Newton::gx(int a,int b,int c,int d) const {

   if( (a == b) || (c == d) )
      return 0.0;

   int phase = 1;

   if(a > b)
      phase *= -1;

   if(c > d)
      phase *= -1;

   int I = TPM::gs2t(a,b);
   int J = TPM::gs2t(c,d);

   int i = TPTPM::gt2tpmm(I,J);

   return phase * x[i];

}

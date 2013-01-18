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
 * @param D metric matrix
 */
void Newton::construct(const SUP &D){

   //construct the p part of the hessian
   H->I(D.gI());

#ifdef __Q_CON
   H->Q(D.gQ());
#endif

#ifdef __G_CON
   H->G(D.gG());
#endif

#ifdef __T1_CON
   H->T(D.gT1());
#endif

#ifdef __T2_CON
   H->T(D.gT2());
#endif

   //the constraint/lagrange multiplier part of the Hessian
   H->lagr();

   //last element zero!
   (*H)(TPTPM::gn(),TPTPM::gn()) = 0.0;
   
   H->symmetrize();

   //and last but not least, solve the system
   H->solve_sy(gradient->gpointer());

}


/**
 *  set the right-hand side of the Newton equations
 */
void Newton::set_rhs(const TPM &tpm){

   gradient->convert(tpm);

   //last part of right-hand side (lagrange multiplier)
   (*gradient)[TPTPM::gn()] = 0.0;

}

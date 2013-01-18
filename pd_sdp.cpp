/**
 * @mainpage 
 * This is an implementation of the dual only, potential reduction interior point method
 * for optimizing the second order density matrix using the P, Q, G and T_1 N-representability conditions.
 * Compiling can be done with the options PQ, PQG, PQGT1, PQGT2 and PQGT (for all conditions active) with logical consequences for the program.
 * @author Brecht Verstichel, Ward Poelmans
 * @date 22-02-2010
 */

#include <iostream>
#include <fstream>
#include <cmath>

using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;

//includes all important headers and defines which conditions are
//going to be used:
#include "include.h"

/**
 * In the main the actual program is run.\n 
 * We start from the unity density matrix normed on the particle number and minimize the 
 * ojective function:\n\n
 * Tr (Gamma H) - t * ln(det P(Gamma)) \n\n
 * Once the minimum is found the parameter t is reduced and a new search is initiated,
 * this goes on until convergence is reached.\n
 * The potential is minimized using the Newton-Raphson method and the resulting linear system
 * is solved via the linear conjugate gradient method.
 */
int main(void) {

   cout.precision(10);

   const int M = 6;//dim sp hilbert space
   const int N = 5;//nr of particles

   Tools::init(M,N);

   TPM::init();
   PHM::init();
   DPM::init();
   PPHM::init();

   TPTPM::init();
   SPSPM::init();

   Gradient::init();

   SUP::init();
   EIG::init();

   Newton newton;

   TPM ham;
   ham.hubbard(1.0);

   SUP Z;
   Z.init_Z();

   SUP X;
   X.init_X(1000.0,ham,Z);

   int dim = SUP::gdim();

   //eerste primal dual gap:
   double pd_gap = Z.ddot(X);
   double energy = (Z.gI()).ddot(ham);

   double center_dev = Z.center_dev(X);

   //eerst centering
   double gamma = 1.0;

   double tolerance = 1.0e-4;

   //flag == 0 : initiele centering run (tot op tolerance)
   //flag == 1 : doe een stap met gamma = 0
   //flag == 2 : doe een stap met gamma = 1
   //flag == 3 : game over man
   int flag = 0;

   double a;//stapgrootte

   int iter = 0;

   while(flag != 3){

      cout << (Z.gI()).trace() << "\t" << pd_gap << "\t" << center_dev << "\t" << energy << endl;

      //matrix D aanmaken voor de hessiaan van het duale stelsel
      SUP D;
      D.D(Z,X);

      //D inverteren voor de hessiaan van het primale stelsel
      SUP D_inv(D);
      D_inv.invert();

      //rechterlid maken van stelsel dat moet worden opgelost:
      SUP B(Z);

      //invert B
      B.invert();

      //schalen met 
      B.dscal(gamma*pd_gap/dim);

      B -= X;

      //collaps B onto b to construct the right hand side of the primal Newton equation
      TPM b;

      b.collaps(0,B);

      newton.set_rhs(b);

      newton.construct(D_inv);
      ++iter;

      //dit wordt de stap:
      TPM delta;
      delta.convert(newton.gGradient());

      //nog updaten van Z en X
      SUP DZ;
      DZ.fill(delta);

      //DX is B - D^{-1}*DZ*D^{-1}
      SUP DX(B);

      //eerst D^{-1}*DZ*D^{-1} in DX stoppen
      B.L_map(D_inv,DZ);

      DX -= B;

      //voor de zekerheid nog projecteren op juiste subruimte:
      DX.proj_C();

      //welke stapgrootte moet ik nemen?
      if(flag == 0 || flag == 2){//voor centering

         Z += DZ;
         X += DX;

      }
      else{

         //zoek de ideale afstand (geef ook een waarde mee voor de maximale afwijking van het centraal pad):
         a = DZ.line_search(DX,Z,X,2.0);

         Z.daxpy(a,DZ);
         X.daxpy(a,DX);

      }

      //update van enkele belangrijke variabelen
      pd_gap = Z.ddot(X);
      energy = (Z.gI()).ddot(ham);
      center_dev = Z.center_dev(X);

      //keuze voor volgende iteratie:
      if(flag == 0){

         //als hij voldoende gecenterd is, exit.
         if(center_dev < tolerance){

            flag = 1;
            gamma = 0.0;

         }

      }
      else if(flag == 1){

         if(pd_gap < tolerance)//exit when converged
            flag = 3;
         else{//center when not convergence

            flag = 2;
            gamma = 1.0;

         }

      }
      else{//flag == 2: dus na een centering stap

         if(pd_gap < tolerance)//exit when converged
            flag = 3;
         else{//take another step downwards when not converged

            flag = 1;
            gamma = 0;

         }

      }

   }

   cout << endl;
   cout << "FINAL RESULT " << endl;
   cout << endl;
   cout << "E_0 = " << energy << " with accuracy of " << pd_gap << " and a deviation from centrality of " << center_dev << endl;
   cout << endl;
   cout << "Total number of newton steps:\t" << iter << endl;

   Gradient::clear();

   SPSPM::clear();
   TPTPM::clear();

   PPHM::clear();
   DPM::clear();
   PHM::clear();
   TPM::clear();

   Tools::clear();

   return 0;
}

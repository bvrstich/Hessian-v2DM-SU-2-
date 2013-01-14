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

vector< vector<int> > *TPM::t2s;
int ***TPM::s2t;

double **TPM::norm;

/**
 * construct the static lists
 */
void TPM::init(){

   int M = Tools::gM();

   //allocatie van s2t
   s2t = new int ** [2];

   for(int S = 0;S < 2;++S){

      s2t[S] = new int * [M];

      for(int a = 0;a < M;++a)
         s2t[S][a] = new int [M];

   }

   t2s = new vector< vector<int> > [2];

   norm = new double * [M];

   for(int a = 0;a < M;++a)
      norm[a] = new double [M];

   vector<int> v(2);

   //initialisatie van de arrays
   int tp = 0;

   //symmetrical array: S = 0
   for(int a = 0;a < M;++a)
      for(int b = a;b < M;++b){

         v[0] = a;
         v[1] = b;

         t2s[0].push_back(v);

         s2t[0][a][b] = tp;
         s2t[0][b][a] = tp;

         if(a == b)
            norm[a][b] = 1.0/(std::sqrt(2.0));
         else
            norm[a][b] = 1.0;

         norm[b][a] = norm[a][b];

         ++tp;

      }

   //watch it!
   tp = 0;

   //antisymmetrical array: S = 1
   for(int a = 0;a < M;++a)
      for(int b = a + 1;b < M;++b){

         v[0] = a;
         v[1] = b;

         t2s[1].push_back(v);

         s2t[1][a][b] = tp;
         s2t[1][b][a] = tp;

         ++tp;

      }

}

/**
 * deallocate the static lists
 */
void TPM::clear(){

   int M = Tools::gM();

   for(int S = 0;S < 2;++S){

      for(int a = 0;a < M;++a)
         delete [] s2t[S][a];

      delete [] s2t[S];

   }
   
   delete [] s2t;

   delete [] t2s;

   for(int a = 0;a < M;++a)
      delete [] norm[a];

   delete [] norm;

}

/**
 * standard constructor for a spinsymmetrical tp matrix: constructs BlockMatrix object with 2 blocks, for S = 0 or 1,
 * @param M nr of sp orbitals
 * @param N nr of particles
 */
TPM::TPM() : BlockMatrix(2) {

   //set the dimension and degeneracy of the two blocks:
   this->setMatrixDim(0,t2s[0].size(),1);
   this->setMatrixDim(1,t2s[1].size(),3);

}

/**
 * copy constructor: constructs Matrix object of dimension M*(M - 1)/2 and fills it with the content of matrix tpm_c
 * @param tpm_c object that will be copied into this.
 */
TPM::TPM(const TPM &tpm_c) : BlockMatrix(tpm_c){ }

/**
 * destructor: if counter == 1 the memory for the static lists t2s en s2t will be deleted.
 * 
 */
TPM::~TPM(){ }

/**
 * access the elements of the the blocks in sp mode, the symmetry or antisymmetry of the blocks is automatically accounted for:\n\n
 * Antisymmetrical for S = 1, symmetrical in the sp orbitals for S = 0\n\n
 * @param S The spinquantumnumber that identifies the block
 * @param a first sp index that forms the tp row index i of spin S, together with b
 * @param b second sp index that forms the tp row index i of spin S, together with a
 * @param c first sp index that forms the tp column index j of spin S, together with d
 * @param d second sp index that forms the tp column index j of spin S, together with c
 * @return the number on place TPM(B,i,j) with the right phase.
 */
double TPM::operator()(int S,int a,int b,int c,int d) const{

   if(S == 0){

      int i = s2t[0][a][b];
      int j = s2t[0][c][d];

      return (*this)(S,i,j);

   }
   else{

      if( (a == b) || (c == d) )
         return 0;
      else{

         int i = s2t[1][a][b];
         int j = s2t[1][c][d];

         int phase = 1;

         if(a > b)
            phase *= -1;
         if(c > d)
            phase *= -1;

         return phase*(*this)(S,i,j);

      }

   }

}

ostream &operator<<(ostream &output,const TPM &tpm_p){

   for(int S = 0;S < 2;++S){

      output << S << "\t" << tpm_p.gdim(S) << "\t" << tpm_p.gdeg(S) << std::endl;
      output << std::endl;

      for(int i = 0;i < tpm_p.gdim(S);++i)
         for(int j = 0;j < tpm_p.gdim(S);++j){

            output << i << "\t" << j << "\t|\t" << tpm_p.t2s[S][i][0] << "\t" << tpm_p.t2s[S][i][1]

               << "\t" << tpm_p.t2s[S][j][0] << "\t" << tpm_p.t2s[S][j][1] << "\t" << tpm_p(S,i,j) << endl;

         }

      std::cout << std::endl;

   }

   return output;

}

/**
 * construct the spinsymmetrical hubbard hamiltonian with on site repulsion U
 * @param U onsite repulsion term
 */
void TPM::hubbard(double U){

   int N = Tools::gN();
   int M = Tools::gM();

   int a,b,c,d;//sp (lattice sites here) orbitals

   double ward = 1.0/(N - 1.0);

   int sign;

   for(int S = 0;S < 2;++S){

      sign = 1 - 2*S;

      for(int i = 0;i < this->gdim(S);++i){

         a = t2s[S][i][0];
         b = t2s[S][i][1];

         for(int j = i;j < this->gdim(S);++j){

            c = t2s[S][j][0];
            d = t2s[S][j][1];

            (*this)(S,i,j) = 0;

            //eerst hopping
            if( (a == c) && ( ( (b + 1)%M == d ) || ( b == (d + 1)%M ) ) )
               (*this)(S,i,j) -= ward;

            if( (b == c) && ( ( (a + 1)%M == d ) || ( a == (d + 1)%M ) ) )
               (*this)(S,i,j) -= sign*ward;

            if( (a == d) && ( ( (b + 1)%M == c ) || ( b == (c + 1)%M ) ) )
               (*this)(S,i,j) -= sign*ward;

            if( (b == d) && ( ( (a + 1)%M == c ) || ( a == (c + 1)%M ) ) )
               (*this)(S,i,j) -= ward;

            //only on-site interaction for singlet tp states:
            if(S == 0)
               if(i == j && a == b)
                  (*this)(S,i,j) += 2.0*U;

            (*this)(S,i,j) *= norm[a][b] * norm[c][d];

         }
      }

   }

   this->symmetrize();

}

/**
 * The spincoupled Q map
 * @param option = 1, regular Q map , = -1 inverse Q map
 * @param tpm_d the TPM of which the Q map is taken and saved in this.
 */
void TPM::Q(int option,const TPM &tpm_d){

   double a = 1;
   double b = 1.0/(Tools::gN()*(Tools::gN() - 1.0));
   double c = 1.0/(Tools::gN() - 1.0);

   this->Q(option,a,b,c,tpm_d);

}

/**
 * The spincoupled Q-like map: see primal-dual.pdf for more info (form: Q^S(A,B,C)(TPM) )
 * @param option = 1, regular Q-like map , = -1 inverse Q-like map
 * @param A factor in front of the two particle piece of the map
 * @param B factor in front of the no particle piece of the map
 * @param C factor in front of the single particle piece of the map
 * @param tpm_d the TPM of which the Q-like map is taken and saved in this.
 */
void TPM::Q(int option,double A,double B,double C,const TPM &tpm_d){

   int M = Tools::gM();

   //for inverse
   if(option == -1){

      B = (B*A + 2.0*B*C*M - 2.0*C*C)/( A * (2.0*C*(M - 1.0) -  A) * ( A + 2.0*B*M*(2.0*M - 1.0) - 2.0*C*(2.0*M - 1.0) ) );
      C = C/(A*(2.0*C*(M - 1.0) - A));
      A = 1.0/A;

   }

   SPM spm;
   spm.bar(C,tpm_d);

   //de trace*2 omdat mijn definitie van trace in berekeningen over alle (alpha,beta) loopt
   double ward = B*tpm_d.trace()*2.0;

   int sign;

   //loop over the spinblocks
   for(int S = 0;S < 2;++S){

      //symmetry or antisymmetry?
      sign = 1 - 2*S;

      for(int i = 0;i < this->gdim(S);++i){

         int a = t2s[S][i][0];
         int b = t2s[S][i][1];

         for(int j = i;j < this->gdim(S);++j){

            int c = t2s[S][j][0];
            int d = t2s[S][j][1];

            //here starts the Q-map

            //the tp part
            (*this)(S,i,j) = A*tpm_d(S,i,j);

            //the np part
            if(i == j)
               (*this)(S,i,i) += ward;

            //and four sp parts:
            if(a == c)
               (*this)(S,i,j) -= norm[a][b] * norm[c][d] * spm(b,d);

            if(b == c)
               (*this)(S,i,j) -= sign * norm[a][b] * norm[c][d] * spm(a,d);

            if(a == d)
               (*this)(S,i,j) -= sign * norm[a][b] * norm[c][d] * spm(b,c);

            if(b == d)
               (*this)(S,i,j) -= norm[a][b] * norm[c][d] * spm(a,c);

         }
      }

   }

   this->symmetrize();

}

/**
 * initialize this onto the unitmatrix with trace N*(N - 1)/2
 */
void TPM::unit(){

   double ward = Tools::gN()*(Tools::gN() - 1.0)/(2.0*Tools::gM()*(2.0*Tools::gM() - 1.0));

   for(int S = 0;S < 2;++S){

      for(int i = 0;i < this->gdim(S);++i){

         (*this)(S,i,i) = ward;

         for(int j = i + 1;j < this->gdim(S);++j)
            (*this)(S,i,j) = (*this)(S,j,i) = 0.0;

      }
   }

}

/**
 * fill the TPM object with the S^2 matrix
 */
void TPM::set_S_2(){

   *this = 0.0;

   for(int i = 0;i < this->gdim(0);++i)
      (*this)(0,i,i) = -1.5 * (Tools::gN() - 2.0)/(Tools::gN() - 1.0);

   for(int i = 0;i < this->gdim(1);++i)
      (*this)(1,i,i) = -1.5 * (Tools::gN() - 2.0)/(Tools::gN() - 1.0) + 2.0;

}

/**
 * convert a 2DM object from vector to matrix/TPM form
 * @param grad input Gradient object
 */
void TPM::convert(const Gradient &grad){

   int tpmm;

   for(int S = 0;S < 2;++S)
      for(int i = 0;i < gdim(S);++i)
         for(int j = i;j < gdim(S);++j){

            tpmm = TPTPM::gt2tpmm(S,i,j);

            (*this)(S,i,j) = grad[tpmm]/ ( 2.0 * Gradient::gnorm(tpmm) * (2.0 * S + 1.0) );

         }

   this->symmetrize();

}

/**
 * perform a line search what step size in along the Newton direction is ideal.
 * @param t potential scaling factor
 * @param P SUP matrix containing the inverse of the constraints (carrier space matrices)
 * @param ham Hamiltonian of the problem
 * @return the steplength
 */
double TPM::line_search(double t,SUP &P,TPM &ham){

   double tolerance = 1.0e-5*t;

   if(tolerance < 1.0e-12)
      tolerance = 1.0e-12;

   //neem de wortel uit P
   P.sqrt(1);

   //maak eerst een SUP van delta
   SUP S_delta;

   S_delta.fill(*this);

   //hulpje om dingskes in te steken:
   SUP hulp;

   hulp.L_map(P,S_delta);

   EIG eigen(hulp);

   double a = 0;

   double b = -1.0/eigen.min();

   double c(0);

   double ham_delta = ham.ddot(*this);

   while(b - a > tolerance){

      c = (b + a)/2.0;

      if( (ham_delta - t*eigen.lsfunc(c)) < 0.0)
         a = c;
      else
         b = c;

   }

   return c;

}

/**
 * perform a line search what step size in along the Newton direction is ideal, this one is used for extrapolation.
 * @param t potential scaling factor
 * @param rdm TPM containing the current approximation of the rdm.
 * @param ham Hamiltonian of the problem
 * @return the steplength
 */
double TPM::line_search(double t,TPM &rdm,TPM &ham){

   SUP P;
 
   P.fill(rdm);

   P.invert();

   return this->line_search(t,P,ham);

}

/**
 * @return The expectation value of the total spin for the TPM.
 */
double TPM::S_2() const{

   double ward = 0.0;

   for(int i = 0;i < this->gdim(0);++i)
      ward += -1.5 * (Tools::gN() - 2.0)/(Tools::gN() - 1.0) * (*this)(0,i,i);

   for(int i = 0;i < this->gdim(1);++i)
      ward += 3.0 * ( -1.5 * (Tools::gN() - 2.0)/(Tools::gN() - 1.0) + 2.0 ) * (*this)(1,i,i);

   return ward;

}

/**
 * @return the dimension associated with block 'S'
 * @param S the blockindex S
 */
int TPM::gdim(int S){

   return t2s[S].size();

}

/**
 * access to the TPM norms from outside the class
 */
double TPM::gnorm(int a,int b){

   return norm[a][b];

}

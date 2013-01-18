#include <iostream>
#include <fstream>
#include <cmath>

using std::ostream;
using std::endl;

#include "include.h"

int EIG::dim;

/**
 * initialize the statics
 */
void EIG::init(){

   int M = 2*Tools::gM();

   dim = M*(M - 1);

#ifdef __G_CON
   dim += M*M;
#endif

#ifdef __T1_CON
   dim += M*(M - 1)*(M - 2)/6;
#endif

#ifdef __T2_CON
   dim += M*M*(M - 1)/2;
#endif

}

/**
 * standard constructor with initialization on the eigenvalues of a SUP object.
 * @param sup input SUP object that will be destroyed after this function is called. The eigenvectors
 * of the matrix will be stored in the columns of the original SUP matrix.
 */
EIG::EIG(SUP &sup){

   vI = new BlockVector<TPM>(sup.gI());

#ifdef __Q_CON
   vQ = new BlockVector<TPM>(sup.gQ());
#endif

#ifdef __G_CON
   vG = new BlockVector<PHM>(sup.gG());
#endif

#ifdef __T1_CON
   vT1 = new BlockVector<DPM>(sup.gT1());
#endif

#ifdef __T2_CON
   vT2 = new BlockVector<PPHM>(sup.gT2());
#endif

}

/**
 * Copy constructor\n
 * allocates the memory for the eigenvalues of a SUP object and copies the content of eig_c into it.
 * @param eig_c The input EIG that will be copied into this.
 */
EIG::EIG(const EIG &eig_c){

   vI = new BlockVector<TPM>(eig_c.gvI());

#ifdef __Q_CON
   vQ = new BlockVector<TPM>(eig_c.gvQ());
#endif

#ifdef __G_CON
   vG = new BlockVector<PHM>(eig_c.gvG());
#endif

#ifdef __T1_CON
   vT1 = new BlockVector<DPM>(eig_c.gvT1());
#endif

#ifdef __T2_CON
   vT2 = new BlockVector<PPHM>(eig_c.gvT2());
#endif

}

/**
 * overload equality operator
 * @param eig_c object that will be copied into this.
 */
EIG &EIG::operator=(const EIG &eig_c){

   *vI = eig_c.gvI();

#ifdef __Q_CON
   *vQ = eig_c.gvQ();
#endif

#ifdef __G_CON
   *vG = eig_c.gvG();
#endif

#ifdef __T1_CON
   *vT1 = eig_c.gvT1();
#endif

#ifdef __T2_CON
   *vT2 = eig_c.gvT2();
#endif

   return *this;

}

/**
 * Destructor, deallocation of the memory
 */
EIG::~EIG(){

   delete vI;

#ifdef __Q_CON
   delete vQ;
#endif

#ifdef __G_CON
   delete vG;
#endif

#ifdef __T1_CON
   delete vT1;
#endif

#ifdef __T2_CON
   delete vT2;
#endif

}

ostream &operator<<(ostream &output,const EIG &eig_p){

   std::cout << eig_p.gvI() << std::endl;

#ifdef __Q_CON
   std::cout << eig_p.gvQ() << std::endl;
#endif

#ifdef __G_CON
   std::cout << eig_p.gvG() << std::endl;
#endif

#ifdef __T1_CON
   std::cout << eig_p.gvT1() << std::endl;
#endif

#ifdef __T2_CON
   std::cout << eig_p.gvT2() << std::endl;
#endif

   return output;

}

/** 
 * get the BlockVector<TPM> object containing the eigenvalues of the TPM block I
 * @return a BlockVector<TPM> object containing the desired eigenvalues
 */
BlockVector<TPM> &EIG::gvI(){

   return *vI;

}

/** 
 * const version
 * get the BlockVector<TPM> object containing the eigenvalues of the TPM block I
 * @return a BlockVector<TPM> object containing the desired eigenvalues
 */
const BlockVector<TPM> &EIG::gvI() const{

   return *vI;

}

#ifdef __Q_CON

/** 
 * get the BlockVector<TPM> object containing the eigenvalues of the TPM block Q
 * @return a BlockVector<TPM> object containing the desired eigenvalues
 */
BlockVector<TPM> &EIG::gvQ(){

   return *vQ;

}

/** 
 * const version
 * get the BlockVector<TPM> object containing the eigenvalues of the TPM block Q
 * @return a BlockVector<TPM> object containing the desired eigenvalues
 */
const BlockVector<TPM> &EIG::gvQ() const{

   return *vQ;

}

#endif

#ifdef __G_CON

/** 
 * get the BlockVector<PHM> object containing the eigenvalues of the PHM block G
 * @return a BlockVector<PHM> object containing the desired eigenvalues
 */
BlockVector<PHM> &EIG::gvG(){

   return *vG;

}

/** 
 * get the BlockVector<PHM> object containing the eigenvalues of the PHM block G, const version
 * @return a BlockVector<PHM> object containing the desired eigenvalues
 */
const BlockVector<PHM> &EIG::gvG() const{

   return *vG;

}

#endif

#ifdef __T1_CON

/** 
 * get the BlockVector<DPM> object containing the eigenvalues of the DPM block T1 of the SUP matrix
 * @return a BlockVector<DPM> object containing the desired eigenvalues
 */
BlockVector<DPM> &EIG::gvT1(){

   return *vT1;

}

/** 
 * get the BlockVector<DPM> object containing the eigenvalues of the DPM block T1 of the SUP matrix, const version
 * @return a BlockVector<DPM> object containing the desired eigenvalues
 */
const BlockVector<DPM> &EIG::gvT1() const{

   return *vT1;

}

#endif

#ifdef __T2_CON

/** 
 * get the BlockVector<PPHM> object containing the eigenvalues of the PPHM block T2 of the SUP matrix
 * @return a BlockVector<PPHM> object containing the desired eigenvalues
 */
BlockVector<PPHM> &EIG::gvT2(){

   return *vT2;

}

/** 
 * get the BlockVector<PPHM> object containing the eigenvalues of the PPHM block T2 of the SUP matrix, const version
 * @return a BlockVector<PPHM> object containing the desired eigenvalues
 */
const BlockVector<PPHM> &EIG::gvT2() const{

   return *vT2;

}

#endif

/**
 * @return the minimal element present in this EIG object.
 * watch out, only works when EIG is filled with the eigenvalues of a diagonalized SUP matrix
 */
double EIG::min() const{

   //lowest eigenvalue of P block
   double ward = vI->min();

#ifdef __Q_CON
   //lowest eigenvalue of Q block
   if(ward > vQ->min())
      ward = vQ->min();
#endif

#ifdef __G_CON
   //lowest eigenvalue of G block
   if(ward > vG->min())
      ward = vG->min();
#endif

#ifdef __T1_CON
   //lowest eigenvalue of the T1 block
   if(ward > vT1->min())
      ward = vT1->min();
#endif

#ifdef __T2_CON
   //lowest eigenvalue of the T2 block
   if(ward > vT2->min())
      ward = vT2->min();
#endif

   return ward;

}

/**
 * @return the maximum element present in this EIG object.
 * watch out, only works when EIG is filled with the eigenvalues of a diagonalized SUP matrix
 */
double EIG::max() const{

   //highest eigenvalue of P block
   double ward = vI->max();

#ifdef __Q_CON
   //highest eigenvalue of Q block
   if(ward < vQ->max())
      ward = vQ->max();
#endif

#ifdef __G_CON
   //highest eigenvalue of G block
   if(ward < vG->max())
      ward = vG->max();
#endif

#ifdef __T1_CON
   //highest eigenvalue of the T1 block
   if(ward < vT1->max())
      ward = vT1->max();
#endif

#ifdef __T2_CON
   //highest eigenvalue of the T2 block
   if(ward < vT2->max())
      ward = vT2->max();
#endif

   return ward;

}

/**
 * @return The deviation of the central path as calculated with the logarithmic barrierfunction, the EIG object is calculated
 * in SUP::center_dev.
 */
double EIG::center_dev() const{

   double sum = vI->sum() + vQ->sum();

   double log_product = vI->log_product() + vQ->log_product();

#ifdef __G_CON
   sum += vG->sum();

   log_product += vG->log_product();
#endif

#ifdef __T1_CON
   sum += vT1->sum();

   log_product += vT1->log_product();
#endif

#ifdef __T2_CON
   sum += vT2->sum();

   log_product += vT2->log_product();
#endif

   return dim*log(sum/(double)dim) - log_product;

}

/**
 * @return the deviation of the central path measured trough the logarithmic potential barrier (see primal_dual.pdf), when you take a stepsize alpha from
 * the point (S,Z) in the primal dual newton direction (DS,DZ), for which you have calculated the generalized eigenvalues eigen_S and eigen_Z in SUP::line_search.
 * (*this) = eigen_S --> generalized eigenvalues for the DS step
 * @param alpha the stepsize
 * @param eigen_Z --> generalized eigenvalues for the DS step
 * @param c_S = Tr (DS Z)/Tr (SZ): parameter calculated in SUP::line_search
 * @param c_Z = Tr (S DZ)/Tr (SZ): parameter calculated in SUP::line_search
 */
double EIG::centerpot(double alpha,const EIG &eigen_Z,double c_S,double c_Z) const{

   double ward = dim*log(1.0 + alpha*(c_S + c_Z));

   ward -= vI->centerpot(alpha) + eigen_Z.gvI().centerpot(alpha);

   ward -= vQ->centerpot(alpha) + eigen_Z.gvQ().centerpot(alpha);

#ifdef __G_CON
   ward -= vG->centerpot(alpha) + eigen_Z.gvG().centerpot(alpha);
#endif

#ifdef __T1_CON
   ward -= vT1->centerpot(alpha) + eigen_Z.gvT1().centerpot(alpha);
#endif

#ifdef __T2_CON
   ward -= vT2->centerpot(alpha) + eigen_Z.gvT2().centerpot(alpha);
#endif

   return ward;

}

/**
 * @return the total dimension
 */
int EIG::gdim(){

   return dim;

}

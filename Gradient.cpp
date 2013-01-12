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
 * @return the pointer to the vector, for blas and lapack stuff
 */
const double *Gradient::gpointer() const{

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

#ifndef NEWTON_H
#define NEWTON_H

#include <iostream>
#include <cstdlib>
#include <vector>

using std::ostream;
using std::vector;

#include "Vector.h"

/**
 * @author Brecht Verstichel
 * @date 21-11-2012\n\n
 * This is a class written for the construction and solution of the Newton system. 
 */
class Newton{

   public:

      //constructor
      Newton();

      //copy constructor
      Newton(const Newton &);

      //destructor
      virtual ~Newton();

      const Hessian &gH() const;

      void construct(double,const TPM &,const SUP &);

   private:

      //!hessian matrix
      Hessian *H;

      //!input gradient, output delta
      Gradient *gradient;

};

#endif

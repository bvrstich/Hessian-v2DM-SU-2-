#ifndef GRADIENT_H
#define GRADIENT_H

#include <iostream>
#include <fstream>

using std::ostream;

class TPM;
class SUP;

/**
 * @author Brecht Verstichel
 * @date 12-01-2013\n\n
 * This class Gradient is a class written for the construction and the storage of the right-hand side of the Newton equation
 */
class Gradient {

   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param gradient_p the Gradient you want to print
    */
   friend ostream &operator<<(ostream &output,const Gradient &gradient_p);

   public:
      
      //constructor
      Gradient();

      //copy constructor
      Gradient(const Gradient &);

      //destructor
      virtual ~Gradient();

      double &operator[](int);

      const double &operator[](int) const;

      void construct(double,const TPM &,const SUP &);

      const double *gpointer() const;

      double *gpointer();

      void convert(const TPM &);

      static void init();

      static void clear();
      
      static double gnorm(int);

      static double gnorm(int,int);


   private:

      //!vector on 2DM + lagrange multiplier space,
      double *gradient;

      //!norms associated with the vector
      static double *norm;

      //!dimension of the vector
      static int n;

};

#endif

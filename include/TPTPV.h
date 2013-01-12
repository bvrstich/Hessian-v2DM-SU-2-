#ifndef TPTPV_H
#define TPTPV_H

#include <iostream>
#include <fstream>

using std::ostream;

/**
 * @author Brecht Verstichel
 * @date 12-01-2013\n\n
 * This class TPTPV is a class written for the transformation of a 2DM to a vector.
 * It contains the lists relating these two objects and functions that transform them.
 */

class TPTPV {

   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param tpvv_p the TPTPV you want to print
    */
   friend ostream &operator<<(ostream &output,const TPTPV &tpvv_p);

   public:
      
      //constructor
      TPTPV();

      //copy constructor
      TPTPV(const TPTPV &);

      //destructor
      virtual ~TPTPV();

      double &operator[](int);

      const double &operator[](int) const;

      const double *gpointer() const;

      void convert(const TPM &);

      static double gnorm(int);

      static double gnorm(int,int);
 
      static void init();

      static void clear();

   private:

      //!vector on 2DM space, contains the numbers
      double *tpvv;

      //!norms associated with the vector
      static double *norm;

      //!dimension of the vector
      static int n;

};

#endif

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "Colvar.h"
#include "ActionRegister.h"
#include "tools/KernelFunctions.h"
#include "tools/Communicator.h"

#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdio>

using namespace std;

namespace PLMD{
  namespace colvar{
    
//+PLUMEDOC COLVAR FEDERQ_CLUST
/*

*/
//+ENDPLUMEDOC
   
    class Zangy : public Colvar{

    public:

      vector<AtomNumber> center;    // LIST OF MOLECULES NEEDED BY SMACK -> CENTER
      vector<AtomNumber> start;
      vector<AtomNumber> end;
      vector<AtomNumber> all_atoms;
      vector<double> angles;
      vector<double> width;
      unsigned int n_angles;

      unsigned int nmol;   // total number of molecules involved in the calculations
      unsigned int nmol2;  // nmol2 = nmol*2
      unsigned int atoms;  // total number of atoms involved in the calculations

      double zetac, zetal, sigmac, sigmal;  // values for kernel function
      double zetac_c, zetal_c, zbox;
      double lbound, ubound, lbound_c, ubound_c;
      vector<double> kval {0.0, 0.0, 0.0};  // initialize kernel function and its derivative
      vector<double> dkval {0.0, 0.0, 0.0};  // initialize kernel function and its derivative
      vector<double> phi_i {0.0, 0.0, 0.0}; // initialize kernel function and its derivative
      vector<double> dphi_i {0.0, 0.0, 0.0}; // initialize kernel function and its derivative
      vector<double> phi_j {0.0, 0.0, 0.0}; // initialize kernel function and its derivative
      vector<double> dphi_j {0.0, 0.0, 0.0}; // initialize kernel function and its derivative
      void kernel(double);  // kernel function

      double sigmaf, f_ij, df_ij, r_0; // values for switching function f_ij
      void stepfunction(double);  // switching function f_ij



      ofstream fdbg;   // definition of object needed for debugging
      ofstream rdbg;
      
      Zangy(const ActionOptions&);              //    CONSTRUCTOR
      ~Zangy();                                 //    DESTRUCTOR
      // active methods:
      static void registerKeywords( Keywords& keys );  // KEYWORDS
      virtual void calculate();                        // CALCULATE CV 

      double dotprod(Vector, Vector);
      double norm2(Vector);

    };
    
    PLUMED_REGISTER_ACTION(Zangy,"ZANGY")
    
    void Zangy::registerKeywords( Keywords& keys ){

      Colvar::registerKeywords(keys);
      keys.add("atoms","CENTER","the labels of the atoms acting as center of the molecules");
      //keys.add("atoms","CENTER2","the labels of the atoms acting as center of the neighboring molecules");
keys.add("atoms","START","the labels of the atoms acting as start of the intramolecular vector");
      keys.add("atoms","END","the labels of the atoms acting as end  of the intramolecular vector");
      keys.add("compulsory","ANGLES"," Angles that have to be used in the SMAC calculation (need to be associated with width)");
      keys.add("compulsory","N_ANGLES"," Number of angles that have to be used in the SMAC calculation (need to be associated with width)");
      keys.add("compulsory","WIDTH"," Width of the Gaussian used to calculate the angles distribution (need to be associated with angles)");
      keys.add("compulsory","LBOUND","lower bound in z direction below which molecules are not considered");
      keys.add("compulsory","UBOUND","lower bound in z direction above which molecules are not considered");
      keys.add("compulsory","ZETAC","crystal side boundary of the interval switching function (in fractional coordinates)");
      keys.add("compulsory","ZETAL","liquid side boundary of the interval switching function (in fractional coordinates)");
      keys.add("compulsory","SIGMAC","crystal side width of gaussian for the interval switching function (in fractional coordinates)");
      keys.add("compulsory","SIGMAL","liquid side width of gaussian for the interval switching function (in fractional coordinates)");
      keys.add("compulsory","R_0","distance for nearest neighbours");
      keys.add("compulsory","SIGMAF","steepness of switching function for nearest neighbour cutoff function f_ij (in fractional coordinates)");

      keys.remove("NOPBC");

    }
    
    Zangy::Zangy(const ActionOptions&ao):
      PLUMED_COLVAR_INIT(ao)
    {

      parseAtomList("CENTER",center);
      parseAtomList("START",start);
      parseAtomList("END",end);

      parse("N_ANGLES",n_angles);
      parseVector("ANGLES",angles);
      log.printf("ANGLES ARE: %f %f %f %f \n", angles[0], angles[1], angles[2], angles[3]);
      width.resize(n_angles);
      parseVector("WIDTH", width);

      parse("ZETAC", zetac);
      parse("ZETAL", zetal);
      parse("SIGMAC", sigmac);
      parse("SIGMAL", sigmal);
      
      parse("LBOUND", lbound);
      parse("UBOUND", ubound);
      parse("R_0", r_0);
      parse("SIGMAF", sigmaf);

      nmol = center.size();
      nmol2 = nmol*2;

      if (nmol == 0) error("no molecules specified");

      for (unsigned int i = 0; i < nmol; i++) {
        all_atoms.push_back(center[i]);  // vector.push_back(element) adds element to the end of the vector
      }
      for (unsigned int i = 0; i < nmol; i++) {
        all_atoms.push_back(start[i]);
      }
      for (unsigned int i = 0; i < nmol; i++) {
        all_atoms.push_back(end[i]);  // all atoms involved in the calculations are now in a vector of size 3*mols
      }

      atoms = all_atoms.size();

      addValueWithDerivatives();  // informs plumed core that we require space to store the value      

      setNotPeriodic();        // of the CV and that the CV will act on the list of atoms
      requestAtoms(center);    // named center

      checkRead();  // check that everything on the input line has been read properly,
                    // all parse command should follow before checkRead()

    }
    
void Zangy::kernel(double z) {
  // switching function to limit action of Gsmac CV
  double f_lc = 0;
  double f_ll = 0;

  kval[0] = 0;
  kval[1] = 0;
  kval[2] = 0;
  dkval[2] = 0;

  // get value of switching function in dependence of z
  f_lc = 1/(1 + exp(-sigmac*(z-zetac_c)));
  f_ll = 1/(1 + exp(-sigmal*(z-zetal_c)));

  kval[2] = f_lc*(1-f_ll);
  if ( kval[2] < 0.000001 ) { 
    kval[2] = 0.0;
    dkval[2] = 0.0;
  } else if ( kval[2] > 0.999999 ) {
    kval[2] = 1.0;
    dkval[2] = 0.0;
  } else {
    dkval[2] = sigmac*f_lc*(1-f_lc)*(1-f_ll) - sigmal*f_lc*f_ll*(1-f_ll);
  }
  kval[1] = kval[2];
  kval[0] = kval[2];

}


void Zangy::stepfunction(double r_ij) {
  // switching function to determin the number of nearest neighbours

  double f_f = 0;
  f_ij = 0;
  df_ij = 0;

  // get value of switching function in dependence of r_ij
  f_f = 1/(1 + exp(-sigmaf*(r_ij - r_0)));
  f_ij = 1 - f_f;

  df_ij = - sigmaf*f_f*(1-f_f);
  if ( f_ij < 0.000001 ) {
    df_ij = 0.0;
    f_ij = 0.0;
  } else if ( f_ij > 0.999999 ) {
    df_ij = 0.0;
    f_ij = 1.0;
  }

}


double Zangy::norm2(Vector vect){  // CALCULATE THE NORM OF A VECTOR
  return (vect[0]*vect[0]+vect[1]*vect[1]+vect[2]*vect[2]);
}


double Zangy::dotprod(Vector vect1,Vector vect2){  // CALCULATE THE SCALAR PRODUCT OF TWO VECTORS
  return(vect1[0]*vect2[0]+vect1[1]*vect2[1]+vect1[2]*vect2[2]);
}

    
void Zangy::calculate()
{

  zbox = getBox()[2][2];
  zetac_c = zetac*zbox;
  zetal_c = zetal*zbox;

  lbound_c = lbound*zbox;
  ubound_c = ubound*zbox;

  double cv_val;   // CV
  cv_val=0;
  
  Tensor virial;   // VIRIAL
  virial.zero();   // no virial contribution
  vector<Vector> deriv(getNumberOfAtoms());  // DERIVATIVES, vector of customized Plumed vectors
  
  unsigned int stride=comm.Get_size();  // SET THE PARALLELIZATION VARIABLES for the for loops
  unsigned int rank=comm.Get_rank();
  
  stride=comm.Get_size();
  rank=comm.Get_rank();  

  for(unsigned int i=rank;i<nmol;i+=stride) {  // SUM OVER MOLECULES

    Vector pos_i = getPosition(i);  // position of molecule i

    if ( (lbound_c < pos_i[2]) && (pos_i[2] < ubound_c) ) { 
  
      kernel(pos_i[2]);  // calculate value of kernel function kval and its derivative dkval
      phi_i = kval;
      dphi_i = dkval;

      Vector v_i = pbcDistance(getPosition(i+nmol),getPosition(i+nmol2));

      vector<double> f(nmol);      // switching function f
      vector<Vector> df_i(nmol);   // derivative of f with respect to x_i (3 values per entry)
      vector<Vector> df_j(nmol);   // derivative of f with respect to x_j (3 values per entry)

      vector<double> omega(nmol);  // angle gaussians sum
      vector<Vector> domega_i;     // derivative of angle gaussian sum with respect to x_i
      vector<Vector> domega_j;     // derivative of angle gaussian sum with respect to x_j
      vector<Vector> domega_is;    // derivative of angle gaussian sum with respect to x_is
      vector<Vector> domega_ie;    // derivative of angle gaussian sum with respect to x_ie
      vector<Vector> domega_js;    // derivative of angle gaussian sum with respect to x_js
      vector<Vector> domega_je;    // derivative of angle gaussian sum with respect to x_je
      // x_is is vector of atom corresponding to start atom of molecule i
      // x_ie is vector of atom corresponding to end atom of molecule i
      // x_js is vector of atom corresponding to start atom of molecule j
      // x_je is vector of atom corresponding to end atom of molecule j

      Vector dist; // 3D vector
  
      for( int j=0; j < nmol; ++j) {                   // SUM OVER NEIGHBORS

	Vector pos_j = getPosition(j);

     	if ( (j != i) && (lbound_c < pos_j[2]) && (pos_j[2] < ubound_c) ) {
          double modij;
          dist=pbcDistance(pos_i,pos_j);  // distance vector between molecule i and molecule j
          modij=dist.modulo();            // scalar length of the distance vector

          kernel(pos_j[2]);  // calculate value of kernel function kval and its derivative dkval
          phi_j = kval;      // spatial switching function for molecule j
          dphi_j = dkval;    // derivative of spatial switching function j with respect to x_j

          //log << "j: " << j << ", zpos: " << zpos << ", kval[2]: " << kval[2] << "\n";
 
          stepfunction(modij);  // calculate switching function for molecule pair i and j
          f[j] += f_ij*phi_i[2]*phi_j[2];  // coordination number of molecule i
  
          double dfdix = 0;
          for(unsigned int ix=0; ix<3; ix++){
            dfdix = -df_ij*dist[ix]/modij*phi_i[2]*phi_j[2];
            df_i[i][ix] += dfdix + f_ij*dphi_i[ix]*phi_j[2];  // derivative of switching function with respect to x_i
      	    df_j[j][ix] += -dfdix + f_ij*phi_i[2]*dphi_j[ix];  // derivative of switching function with respect to x_j
	  }

	  Vector v_j = pbcDistance(getPosition(j+nmol), getPosition(j+nmol2));
          
	  double alpha;
          alpha = dotprod(v_i, v_j)/sqrt(norm2(v_i)*norm2(v_j));

	  if (alpha > 1.0) {
            alpha = 0.99999;
          } else if (alpha < -1.0) {
            alpha = -0.99999;
          }

          double theta;
	  theta = acos(alpha);
	  
	  omega[j] = 0;
	  double domega = 0;

          for (unsigned int k = 0; k < angles.size(); k++) {
            double gaussians;
            gaussians = exp(-((theta - angles[k])*(theta - angles[k]))/(2*width[k]*width[k]));
            omega[j] += gaussians;
	    domega += - gaussians*(theta - angles[k])/(width[k]*width[k]);
	  }

          double dvac;  // no idea what dvac is
	  dvac = - 1/(sqrt( norm2(v_i)*norm2(v_j))*sqrt(1.0 - alpha*alpha));

          for (unsigned int ix = 0; ix < 3; ix++) {
            domega_i[j][ix] = - domega * dvac * (v_i[ix] - dotprod(v_i, v_j) * v_j[ix]/norm2(v_j));
            domega_j[j][ix] =   domega * dvac * (v_j[ix] - dotprod(v_j, v_j) * v_i[ix]/norm2(v_i));
	  }

          for (unsigned int ix = 0; ix < 3; ix++) {
            
	  }

          cv_val += omega[j]*f[j];

	}
      }
  
//      cv_val += 0;  // SUM THE TOTAL CV
  
    }
  }
 
 
  comm.Sum(deriv);
  comm.Sum(cv_val);
  comm.Sum(virial);

  for(unsigned i=0;i<nmol;i++) {
    setAtomsDerivatives(i,deriv[i]);
  }
  
  setBoxDerivatives(virial);
  setValue(cv_val);

}
  
  
    
    Zangy::~Zangy(){
    }

  }
}

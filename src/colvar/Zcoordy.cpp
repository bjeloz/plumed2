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
   
    class Zcoordy : public Colvar{

    public:

      vector<AtomNumber> center;    // LIST OF MOLECULES NEEDED BY SMACK -> CENTER

      unsigned int nmol;   // total number of molecules involved in the calculations

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

      unsigned int stride;
      unsigned int rank;

      ofstream fdbg;   // definition of object needed for debugging
      ofstream rdbg;
      
      Zcoordy(const ActionOptions&);              //    CONSTRUCTOR
      ~Zcoordy();                                 //    DESTRUCTOR
      // active methods:
      static void registerKeywords( Keywords& keys );  // KEYWORDS
      virtual void calculate();                        // CALCULATE CV 
      
    };
    
    PLUMED_REGISTER_ACTION(Zcoordy,"ZCOORDY")
    
    void Zcoordy::registerKeywords( Keywords& keys ){

      Colvar::registerKeywords(keys);
      keys.add("atoms","CENTER","the labels of the atoms acting as center of the molecules");
      //keys.add("atoms","CENTER2","the labels of the atoms acting as center of the neighboring molecules");
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
    
    Zcoordy::Zcoordy(const ActionOptions&ao):
      PLUMED_COLVAR_INIT(ao)
    {

      parseAtomList("CENTER",center);

      parse("ZETAC", zetac);
      parse("ZETAL", zetal);
      parse("SIGMAC", sigmac);
      parse("SIGMAL", sigmal);
      
      parse("LBOUND", lbound);
      parse("UBOUND", ubound);
      parse("R_0", r_0);
      parse("SIGMAF", sigmaf);

      nmol=center.size();
      
      if(nmol==0) error("no molecules specified");

      addValueWithDerivatives();  // informs plumed core that we require space to store the value      

      setNotPeriodic();        // of the CV and that the CV will act on the list of atoms
      requestAtoms(center);    // named center

      checkRead();  // check that everything on the input line has been read properly,
                    // all parse command should follow before checkRead()

    }
    
void Zcoordy::kernel(double z) {
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


void Zcoordy::stepfunction(double r_ij) {
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


    
void Zcoordy::calculate()
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
  
  //unsigned int stride=comm.Get_size();  // SET THE PARALLELIZATION VARIABLES for the for loops
  //unsigned int rank=comm.Get_rank();
  
  stride=comm.Get_size();
  rank=comm.Get_rank();  

  for(unsigned int i=rank;i<nmol;i+=stride) {  // SUM OVER MOLECULES

    Vector pos_i = getPosition(i);

    if ( (lbound_c < pos_i[2]) && (pos_i[2] < ubound_c) ) { 
  
      kernel(pos_i[2]);  // calculate value of kernel function kval and its derivative dkval
      phi_i = kval;
      dphi_i = dkval;

      double n;
      n=0.;
  
      //vector<double> f;       // SWITCHING FUNCTION
  
      Vector dist; // 3D vector
  
      for(unsigned int j=0; j < nmol; ++j) {                   // SUM OVER NEIGHBORS

	Vector pos_j = getPosition(j);
	
     	if ( (j != i) && (lbound_c < pos_j[2]) && (pos_j[2] < ubound_c) ) {
          double modij;
          dist=pbcDistance(pos_i,pos_j);    // DISTANCE BETWEEN THEM
          modij=dist.modulo();  // scalar length of the distance vector

          kernel(pos_j[2]);  // calculate value of kernel function kval and its derivative dkval
          phi_j = kval;
          dphi_j = dkval;

          //log << "j: " << j << ", zpos: " << zpos << ", kval[2]: " << kval[2] << "\n";
    
          // calculate 
          stepfunction(modij);  // calculate switching function for molecule pair i and j
          n += f_ij*phi_i[2]*phi_j[2];  // coordination number of molecule i
  
          double dfdix = 0;
          for(unsigned int ix=0; ix<3; ix++){
            dfdix = -df_ij*dist[ix]/modij*phi_i[2]*phi_j[2];
            deriv[i][ix] += dfdix + f_ij*dphi_i[ix]*phi_j[2];  // derivative of switching function with respect to x_i
      	    deriv[j][ix] += -dfdix + f_ij*phi_i[2]*dphi_j[ix];  // derivative of switching function with respect to x_j 
            // calculate dn:
          }
        }
      }
    

      Vector position = getPosition(i);
      double zpos = position[2];
      kernel(zpos);  // calculate value of kernel function kval and its derivative dkval
  
      cv_val += n;  // SUM THE TOTAL CV

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
  
  
    
    Zcoordy::~Zcoordy(){
    }

  }
}

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
#include "tools/SwitchingFunction.h"
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
   
    class Zdensy : public Colvar{
      SwitchingFunction distanceswitch;  // definition of switching function for f_ij
      SwitchingFunction coordswitch;     // definition of switching function for rho_i

    public:

      vector<AtomNumber> start;      // LIST OF MOLECULES NEEDED BY SMACK -> START
      vector<AtomNumber> end;        // LIST OF MOLECULES NEEDED BY SMACK -> END
      vector<AtomNumber> center;    // LIST OF MOLECULES NEEDED BY SMACK -> CENTER1
      //vector<AtomNumber> center2;    // LIST OF MOLECULES NEEDED BY SMACK -> CENTER2
      vector<AtomNumber> all_atoms;  // LIST OF MOLECULES NEEDED BY SMACK -> CENTER
      vector<double> angles;
      vector<double> width;

      unsigned int mols;   // total number of molecules involved in the calculations
      unsigned int atoms;  // total number of atoms involved in the calculation

      double zetac, zetal, sigmac, sigmal; // values for kernel function
      vector<double> kval {0.0, 0.0, 0.0};  // initialize kernel function and its derivative
      vector<double> dkval {0.0, 0.0, 0.0};  // initialize kernel function and its derivative
      void kernel(double);  // kernel function

      ofstream fdbg;   // definition of object needed for debugging
      ofstream rdbg;
      
      
      struct Zdensylist_s {   // STRUCTURE NEEDED TO CONSTRUCT THE VERLET LIST
	double rcut;      // CUT OFF OF THE LIST
	double rskin;     // SKIN TO CHECK IF THE LIST HAS TO BE RECALCULATED
	double rskin2;    // SKIN SQUARED
	int step;
	vector<int> nn;       // NUMBERS OF NEIGHBORS
	vector<vector<int> > ni; // List, it's a vector of vectors, for each molecule one vector
	vector<Vector> pos1;  // POSITIONS CENTER MOLECULE
	vector<Vector> pos2;  // POSITIONS NEIGHBORING MOLECULES
      } Zdensylist;  // object of type Zdensylist_s defined
      
      
      Zdensy(const ActionOptions&);              //    CONSTRUCTOR
      ~Zdensy();                                 //    DESTRUCTOR
      // active methods:
      static void registerKeywords( Keywords& keys );  // KEYWORDS
      virtual void calculate();                        // CALCULATE CV 

      void Zdensy_newlist(vector<AtomNumber>&, Zdensylist_s &Zdensylist);   // NEW VERLIST
      void Zdensy_checklist(vector<AtomNumber>&, Zdensylist_s &Zdensylist); // CHECK THE VERLET LIST
      double dotprod(Vector,Vector);
      double norm2(Vector);
      
    };
    
    PLUMED_REGISTER_ACTION(Zdensy,"ZDENSY")
    
    void Zdensy::registerKeywords( Keywords& keys ){

      Colvar::registerKeywords(keys);
      keys.add("atoms","CENTER","the labels of the atoms acting as center of the molecules");
      keys.add("compulsory","ZETAC","crystal side boundary of the interval switching function (in fractional coordinates)");
      keys.add("compulsory","ZETAL","liquid side boundary of the interval switching function (in fractional coordinates)");
      keys.add("compulsory","SIGMAC","crystal side width of gaussian for the interval switching function (in fractional coordinates)");
      keys.add("compulsory","SIGMAL","liquid side width of gaussian for the interval switching function (in fractional coordinates)");

      keys.remove("NOPBC");

    }
    
    Zdensy::Zdensy(const ActionOptions&ao):
      PLUMED_COLVAR_INIT(ao)
    {
      
      parseAtomList("CENTER",center);

      parse("ZETAC",zetac);
      parse("ZETAL",zetal);
      parse("SIGMAC",sigmac);
      parse("SIGMAL",sigmal);

      mols=center.size();
      
      if(mols==0) error("no molecules specified");
      
      for(unsigned int i=0 ; i < mols ; i++){
	all_atoms.push_back(center[i]);  // vector.push_back(element) adds element to the
      }                                   // end of the vector
      
      atoms=all_atoms.size();
      
      addValueWithDerivatives();  // informs plumed core that we require space to store the value      
      setNotPeriodic();           // of the CV and that the CV will act on the list of atoms
      requestAtoms(all_atoms);    // named all_atoms

      checkRead();  // check that everything on the input line has been read properly,

    }
    
void Zdensy::kernel(double z) {
  // switching function to limit action of Gsmac CV

  // transform input values into cartesian coordinates
  double zbox = getBox()[2][2];
  double zetac_c = zetac*zbox;
  double zetal_c = zetal*zbox;

  double f_lc = 0;
  double f_ll = 0;

  // get value of switching function in dependence of z
  f_lc = 1/(1 + exp(-sigmac*(z-zetac_c)));
  f_ll = 1/(1 + exp(-sigmal*(z-zetal_c)));

  kval[2] = f_lc*(1-f_ll);
  if ( kval[2]<0.000001 ) { 
    kval[2]=0.0;
  }
  //kval[1] = kval[2];
  //kval[0] = kval[2];

  dkval[2] = sigmac*f_lc*(1-f_lc)*(1-f_ll) - sigmal*f_lc*f_ll*(1-f_ll);
  if ( kval[2]<0.000001 ) { 
    dkval[2]=0.0;
  }

}

    
void Zdensy::calculate()
{
  
  double cv_val;   // CV
  cv_val=0;
  
  Tensor virial;   // VIRIAL
  virial.zero();   // no virial contribution
  vector<Vector> deriv(getNumberOfAtoms());  // DERIVATIVES, vector of customized Plumed vectors
  
  unsigned int stride=comm.Get_size();  // SET THE PARALLELIZATION VARIABLES for the for loops
  unsigned int rank=comm.Get_rank();
  
  stride=comm.Get_size();
  rank=comm.Get_rank();
  
  for(unsigned int i=rank;i<mols;i+=stride) {                   // SUM OVER MOLECULES

    Vector position = getPosition(i);
    double zpos = position[2];
    kernel(zpos);   // calculate value of kernel function kval and its derivative dkval

    cv_val += kval[2];     // SUM THE TOTAL CV 

    for (unsigned int ix=0; ix<2; ix++) {
      deriv[i][ix]  +=  0;
    }
    
    deriv[i][2]  +=  dkval[2];

  }
 
  comm.Sum(deriv);
  comm.Sum(cv_val);  
  comm.Sum(virial);  

 for(unsigned i=0;i<atoms;i++) {
   setAtomsDerivatives(i,deriv[i]);
 }
  
  setBoxDerivatives(virial);
  setValue(cv_val);
  
}
  
    
Zdensy::~Zdensy(){
}

  }
}

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
   
    class Zdensy2 : public Colvar{

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
      double rho_i, drho_i, n_0; // values for density function rho_i (sigmaf is used here as well)
      void ffunction(double);  // switching function f_ij
      void rhofunction(double);  // switching function f_ij

      vector<double> list;  // inded list of all molecules in defined layer interval
      unsigned int ll;      // length of list
      void layerindices(vector<Vector>&);  // index list of all molecules positioned within the defined layer interval

      unsigned int stride;
      unsigned int rank;

      ofstream fdbg;   // definition of object needed for debugging
      ofstream rdbg;
      
      Zdensy2(const ActionOptions&);              //    CONSTRUCTOR
      ~Zdensy2();                                 //    DESTRUCTOR
      // active methods:
      static void registerKeywords( Keywords& keys );  // KEYWORDS
      virtual void calculate();                        // CALCULATE CV 
      
    };
    
    PLUMED_REGISTER_ACTION(Zdensy2,"ZDENSY2")
    
    void Zdensy2::registerKeywords( Keywords& keys ){

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
      keys.add("compulsory","N_0","number of minimum nearest neighbours");
      keys.add("compulsory","SIGMAF","steepness of switching function for nearest neighbour cutoff function f_ij (in fractional coordinates)");

      keys.remove("NOPBC");

    }
    
    Zdensy2::Zdensy2(const ActionOptions&ao):
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
      parse("N_0", n_0);
      parse("SIGMAF", sigmaf);

      nmol=center.size(); 
      list.resize(nmol);  // allocate enough space for layer molecules list

      if(nmol==0) error("no molecules specified");

      addValueWithDerivatives();  // informs plumed core that we require space to store the value      

      setNotPeriodic();        // of the CV and that the CV will act on the list of atoms
      requestAtoms(center);    // named center

      checkRead();  // check that everything on the input line has been read properly,
                    // all parse command should follow before checkRead()

    }
    
void Zdensy2::kernel(double z) {
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


void Zdensy2::ffunction(double r_ij) {
  // switching function to determine the cutoff range of neighbours

  double f_f;
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


void Zdensy2::rhofunction(double n_i) {
  // switching function to determine the local density of molecule i

  rho_i = 0;
  drho_i = 0;

  // get value of switching function in dependence of n_i
  rho_i = 1/(1 + exp(-sigmaf*(n_i - n_0)));

  drho_i = sigmaf*f_ij*(1-f_ij);
  if ( rho_i < 0.000001 ) {
    drho_i = 0.0;
    rho_i = 0.0;
  } else if ( f_ij > 0.999999 ) {
    drho_i = 0.0;
    rho_i = 1.0;
  }

}


void Zdensy2::layerindices(vector<Vector> &pos) {

  vector<Vector> ** pospointer;  // declare pointer for matix pos
  pospointer = &pos;  // assign address of pos to pospointer
  unsigned int k = 0;

  //unsigned int stride;
  //unsigned int rank;

  //stride = comm.Get_size();
  //rank = comm.Get_rank();

  Vector pos_i;
  //for (unsigned int i=rank; i<nmol; i+=stride) {
  for (unsigned int i=0; i<nmol; i++) {

    pos_i = getPosition(i);

    //cout << "pos_i[0]: " << pos_i[0] << ", pos_i[1]: " << pos_i[1] << ", pos_i[2]: " << pos_i[2] << "\n";

    if ( (lbound_c < pos_i[2]) && (pos_i[2] < ubound_c) ) { 
      list[k] = i;

      for (unsigned int ix=0; i<nmol; i++) {
        **pospointer[k][ix] = pos_i[ix];
      }

      k++;

    }
  }

  ll = k;  // length of neighbour list

}

    
void Zdensy2::calculate()
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

  for (unsigned int i=0; i<nmol; i++) {
    for (unsigned int ix=0; ix<3; ix++) {
      deriv[i][ix] = 0;
    }
  }

  vector<Vector> pos(nmol)
  layerindices(pos);  // construct list of molecules iniside the defined interval

  //for (unsigned int i=0; i<ll; i++) {
  //  cout << "list[i]: " << list[i] << "\n";
  //}
  //cout << "--> ll: " << ll << "\n";

  vector<int> drho(ll);  // DERIVATIVES, vector of customized Plumed vectors
  vector<Vector> df(ll);
  //vector<int> drho(nmol);  // DERIVATIVES, vector of customized Plumed vectors
  //vector<Vector> df(nmol);

  unsigned int stride;  // SET THE PARALLELIZATION VARIABLES for the for loops
  unsigned int rank;

  stride=comm.Get_size();  // SET THE PARALLELIZATION VARIABLES for the for loops
  rank=comm.Get_rank();

  for (unsigned int i=rank; i<ll; i+=stride) {  // SUM OVER MOLECULES
  
    kernel(pos[i]);  // calculate value of kernel function kval and its derivative dkval
    phi_i = kval;
    dphi_i = dkval;

    double n;
    n=0.;
  
    //vector<double> f;       // SWITCHING FUNCTION
  
    Vector dist; // 3D vector
  
    for(unsigned int j=0; j<ll; ++j) {                   // SUM OVER NEIGHBORS

      Vector pos_j = getPosition(j);

      if ( (j != i) && (lbound_c < pos_j[2]) && (pos_j[2] < ubound_c) ) {
        double modij;
        dist=pbcDistance(pos[i],pos[j]);    // DISTANCE BETWEEN THEM
        modij=dist.modulo();  // scalar length of the distance vector

        kernel(pos_j[2]);  // calculate value of kernel function kval and its derivative dkval
        phi_j = kval;
        dphi_j = dkval;

        //log << "j: " << j << ", zpos: " << zpos << ", kval[2]: " << kval[2] << "\n";
    
        // calculate
        ffunction(modij);  // calculate switching function for molecule pair i and j
        n += f_ij*phi_i[2]*phi_j[2];  // coordination number of molecule i
  
        double dfdix = 0;
        for (unsigned int ix=0; ix<3; ix++){
          dfdix = -df_ij*dist[ix]/modij*phi_i[2]*phi_j[2];
          df[i][ix] += dfdix + f_ij*dphi_i[ix]*phi_j[2];  // derivative of switching function with respect to x_i
          df[j][ix] += -dfdix + f_ij*phi_i[2]*dphi_j[ix];  // derivative of switching function with respect to x_j 
        }
      }
    }

    rhofunction(n);
    drho[i] = drho_i;

    // calculate the derivative
    for (unsigned int ix=0; ix<3; ix++){
      deriv[i][ix] += df[i][ix];  // derivative of switching function with respect to x_i
      //deriv[j][ix] += df[j][ix];  // derivative of switching function with respect to x_j 
    }

    // sum the total CV
    cv_val += rho_i;

  }
 
  comm.Sum(deriv);
  comm.Sum(cv_val);
  //comm.Sum(virial);

  for(unsigned i=0;i<nmol;i++) {
    setAtomsDerivatives(i,deriv[i]);
  }
  
  setBoxDerivatives(virial);
  setValue(cv_val);

}
  
  
    
    Zdensy2::~Zdensy2(){
    }

  }
}

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
   
    class GsmacIL : public Colvar{
      SwitchingFunction distanceswitch;  // definition of switching function for f_ij
      SwitchingFunction coordswitch;     // definition of switching function for rho_i

    public:

      vector<AtomNumber> start;      // LIST OF MOLECULES NEEDED BY SMACK -> START
      vector<AtomNumber> end;        // LIST OF MOLECULES NEEDED BY SMACK -> END
      vector<AtomNumber> center;    // LIST OF MOLECULES NEEDED BY SMACK -> CENTER
      vector<AtomNumber> all_atoms;  // LIST OF MOLECULES NEEDED BY SMACK
      vector<double> angles;
      vector<double> width;

      unsigned int nmol;   // total number of molecules involved in the calculations
      unsigned int atoms;  // total number of atoms involved in the calculations

      double zetac, zetal, sigmac, sigmal; // values for kernel function
      double zbox, zetac_c, zetal_c;
      vector<double> kval {0.0, 0.0, 0.0};  // initialize kernel function and its derivative
      vector<double> dkval {0.0, 0.0, 0.0};  // initialize kernel function and its derivative
      void kernel(double);  // kernel function

      double lbound, ubound, lbound_c, ubound_c;
      vector<double> list;  // inded list of all molecules in defined layer interval
      unsigned int ll;      // length of list
      void layerindices(vector<Vector>&, vector<int>&);  // index list of all molecules positioned within the defined layer interval
      vector<Vector> pos;
     
      ofstream fdbg;   // definition of object needed for debugging
      ofstream rdbg;
      
      GsmacIL(const ActionOptions&);              //    CONSTRUCTOR
      ~GsmacIL();                                 //    DESTRUCTOR
      // active methods:
      static void registerKeywords( Keywords& keys );  // KEYWORDS
      virtual void calculate();                        // CALCULATE CV 

      double dotprod(Vector,Vector);
      double norm2(Vector);
      
    };
    
    PLUMED_REGISTER_ACTION(GsmacIL,"GSMACIL")
    
    void GsmacIL::registerKeywords( Keywords& keys ){

      Colvar::registerKeywords(keys);
      keys.add("atoms","CENTER","the labels of the atoms acting as center of the molecules");
      keys.add("atoms","START","the labels of the atoms acting as start of the intramolecular vector");
      keys.add("atoms","END","the labels of the atoms acting as end  of the intramolecular vector");
      keys.add("compulsory","ANGLES"," Angles that have to be used in the SMAC calculation (need to be associated with width)");
      keys.add("compulsory","N_ANGLES"," Number of angles that have to be used in the SMAC calculation (need to be associated with width)");
      keys.add("compulsory","WIDTH"," Width of the Gaussian used to calculate the angles distribution (need to be associated with angles)");
      keys.add("optional","SWITCH","Used if you want to employ an alternative swiching function.");
      keys.add("optional","SWITCH_COORD","Used if you want to employ an alternative swiching function.");
      keys.add("compulsory","R_CUT","Verlet lists, default r_cut = R_max");
      keys.add("compulsory","R_SKIN","Verlet lists, default r_skin = 1.25 R_max");
     
      keys.add("compulsory","ZETAC","crystal side boundary of the interval switching function (in fractional coordinates)");
      keys.add("compulsory","ZETAL","liquid side boundary of the interval switching function (in fractional coordinates)");
      keys.add("compulsory","SIGMAC","crystal side width of gaussian for the interval switching function (in fractional coordinates)");
      keys.add("compulsory","SIGMAL","liquid side width of gaussian for the interval switching function (in fractional coordinates)");
      keys.add("compulsory","LBOUND","lower bound in z direction below which molecules are not considered");
      keys.add("compulsory","UBOUND","lower bound in z direction above which molecules are not considered");

      keys.remove("NOPBC");

    }
    
    GsmacIL::GsmacIL(const ActionOptions&ao):
      PLUMED_COLVAR_INIT(ao)
    {
      
      string sw,errors;        //// SWITCH TO DEFINE WHICH SWF WE ARE GOING TO USE FOR THE PAIR
      parse("SWITCH",sw);
      if(sw.length()>0){
	distanceswitch.set(sw,errors);
	if( errors.length()!=0 ) error("problem reading SWITCH keyword : " + errors );
      } else {                  // IN CASE NO SF ARE PROVIDED, USE THE STANDARD RATIONAL ONE
	int nn=6;
	int mm=12;
	double d0=0.0;
	double r0=0.0;
	parse("R_0",r0);
	if(r0<=0.0) error("R_0 should be explicitly specified and positive");
	parse("D_0",d0);
	parse("NN",nn);
	parse("MM",mm);
	distanceswitch.set(nn,mm,r0,d0);
      }
      
      parse("SWITCH_COORD",sw);   //// SWITCH TO DEFINE WHICH SWF WE ARE GOING TO USE FOR THE DENSITY
      if(sw.length()>0){
	coordswitch.set(sw,errors);
	if( errors.length()!=0 ) error("problem reading SWITCH keyword : " + errors );
      } else {                  // IN CASE NO SF ARE PROVIDED, USE THE STANDARD RATIONAL ONE
	int nn=6;
	int mm=12;
	double d0=0.0;
	double r0=0.0;
	parse("R_0",r0);
	if(r0<=0.0) error("R_0 should be explicitly specified and positive");
	parse("D_0",d0);
	parse("NN",nn);
	parse("MM",mm);
	coordswitch.set(nn,mm,r0,d0);
      }


      unsigned int nAngles;
      parse("N_ANGLES",nAngles);
      angles.resize(nAngles);
      parseVector("ANGLES",angles);
      log.printf("ANGLES ARE: %f %f %f %f \n", angles[0], angles[1], angles[2], angles[3]);
      width.resize(nAngles);
      parseVector("WIDTH",width);

      parseAtomList("CENTER",center);
      parseAtomList("START",start);
      parseAtomList("END",end);

      parse("ZETAC", zetac);
      parse("ZETAL", zetal);
      parse("SIGMAC", sigmac);
      parse("SIGMAL", sigmal);
      parse("LBOUND", lbound);
      parse("UBOUND", ubound);


      nmol=center.size();
      
      
      if(nmol==0) error("no molecules specified");
      
      for(unsigned int i=0 ; i < nmol ; i++){
	all_atoms.push_back(center[i]);  // vector.push_back(element) adds element to the
      }                                   // end of the vector
      for(unsigned int i=0 ; i < nmol ; i++){
	all_atoms.push_back(start[i]);
      }
      for(unsigned int i=0 ; i < nmol ; i++){
	all_atoms.push_back(end[i]);  // all_atoms is now a vector of size 4*nmol
      }
      
      atoms=all_atoms.size();
      
      addValueWithDerivatives();  // informs plumed core that we require space to store the value      
      setNotPeriodic();           // of the CV and that the CV will act on the list of atoms
      requestAtoms(all_atoms);    // named all_atoms

      double r_cut,r_skin;
      parse("R_CUT",r_cut);
      parse("R_SKIN",r_skin);

     
      checkRead();  // check that everything on the input line has been read properly,
                    // all parse command should follow before checkRead()
      
      log<<"  contacts are counted with cutoff "<< distanceswitch.description()<<"\n";
      log<<"  correct density is evaluated above with cutoff "<< coordswitch.description()<<"\n";

    }
    
void GsmacIL::kernel(double z) {
  // switching function to limit action of Gsmac CV

  double f_lc;
  double f_ll;

  // get value of switching function in dependence of z
  f_lc = 1/(1 + exp(-sigmac*(z-zetac_c)));
  f_ll = 1/(1 + exp(-sigmal*(z-zetal_c)));

  kval[2] = f_lc*(1-f_ll);
  if ( kval[2]<0.000001 ) { 
    kval[2]=0.0;
  }
  kval[1] = kval[2];
  kval[0] = kval[2];

  dkval[2] = sigmac*f_lc*(1-f_lc)*(1-f_ll) - sigmal*f_lc*f_ll*(1-f_ll);
  if ( kval[2]<0.000001 ) { 
    dkval[2]=0.0;
  }

}


// parallelization not necessary
void GsmacIL::layerindices(vector<Vector> &pos, vector<int> &list) {

  unsigned int k = 0;

  //unsigned int stride;
  //unsigned int rank;

  //stride = comm.Get_size();
  //rank = comm.Get_rank();

  Vector pos_i;
  for (unsigned int i=0; i<nmol; i++) {
  //for (unsigned int i=rank; i<nmol; i+=stride) {

    pos_i = getPosition(i);

    if ( (lbound_c < pos_i[2]) && (pos_i[2] < ubound_c) ) {
      list[k] = i;  // list is the carrier of the true index i for the molecule. k is the index of the molecules in the lbound_c-ubound_c interval region.

      for (unsigned int ix=0; ix<3; ix++) {
        pos[k][ix] = pos_i[ix];   // pos has meaningful entries only up to ll
      }

      k++;

    }
  }

  //comm.Sum(k);
  ll = k;  // length of neighbour list

}


    
void GsmacIL::calculate()
{

  // transform fractional coordinates into cartesian coordinates for kernel layer function
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
  
  unsigned int stride;  // SET THE PARALLELIZATION VARIABLES for the for loops
  unsigned int rank;
  
  stride=comm.Get_size();
  rank=comm.Get_rank();

  vector<Vector> pos(nmol);
  vector<int> list(nmol);
  layerindices(pos, list);  // construct list of molecules iniside the defined interval, pos, together with the list, list, which carries the true index of the molecule.
  
  //cout << "  ###### ll: " << ll << endl;

  // get the vector of each molecule
  vector<Vector> v(ll); 
  for (unsigned int i = rank; i < ll; i += stride) {  // DEFINE THE INTRAMOLECULAR VECTOR
    unsigned int start = list[i] + nmol;  // true index of vector start atom
    unsigned int end = list[i] + 2*nmol;  // true index of vector end atom
    v[i]=pbcDistance(getPosition(start),getPosition(end));

    //cout << "v[i]: " << v[i][0] << ", " << v[i][1] << ", " << v[i][2] << endl;

  }

  comm.Sum(v);  // MERGING UP

  for (unsigned int i = rank; i < ll; i += stride) {    // SUM OVER MOLECULES
    double n, angtot;
    unsigned int start_i, end_i;
    double S1, S2, S3;
    start_i = start[i].serial();       // serial()??? start_i has the true index of the 
    end_i = end[i].serial();

    vector<double> f(ll);            // SWITCHING FUNCTION
    vector<double> omega(ll);        // SWITCHING FUNCTION
    vector<Vector> domega1(ll);      // ANGULAR PART1
    vector<Vector> domega2(ll);      // ANGULAR PART2
    vector<Vector> df_a(ll);         // ANGULAR PART2
    vector<Vector> df_b(ll);         // ANGULAR PART2
    
    n=0.;
    angtot=0.;
    
    Vector dist;

    for (unsigned int j = 0; j < ll; ++j) {    // SUM OVER NEIGHBORS
      if ( j != i ) {
        double modij;
        dist=pbcDistance(pos[i],pos[j]);    // DISTANCE BETWEEN THEM
        modij=dist.modulo();                                              //
 
        double zpos = pos[j][2];
        kernel(zpos);   // calculate value of kernel function kval and its derivative dkval
     
        double dfunc=0.;         // CALCULATING SWITCHING FUNCTION AND 
        double func=0.;
        func = distanceswitch.calculate(modij,dfunc);
        f[j] = func*kval[2];
        n += f[j];               // CALCULATING THE COORDINATION NUMBER

        for (unsigned int ix=0 ;ix<3 ;ix++) {
          df_a[j][ix] = -dfunc*dist[ix]*kval[2];   // DERIVATIVE OF THE SWITCHING FUNCTION, what the hell?!
          df_b[j][ix] = -dfunc*dist[ix]*kval[2] - func*dkval[ix];   // DERIVATIVE OF THE SWITCHING FUNCTION, the minus is set because of the implementation of the derivatives in the temp variable further below in the code. 
        }

        double alpha;
        alpha = dotprod(v[i],v[j])/sqrt(norm2(v[i])*norm2(v[j]));  //  CHANGE THIS PART IN MOR C++ STYLE

        if ( alpha >= 1.0 ) {  // NOT SURE IF THIS CHECK IS NECESSARY ANYMORE. DOUBLECHECK!
          alpha = 0.99999;
        } else if ( alpha <= -1.0 ){
          alpha = -0.99999;
        }
 
        double phi;
        double domega;
        phi = acos(alpha);
        
        omega[j] = 0.;
        domega = 0.;
        for (unsigned int k=0 ; k < angles.size() ;k++) {  // ANGULAR SERIES! NEED TO BE SUBSTITUTED WITH THE KERNEL
          double e1;	
          e1 = exp(-((phi - angles[k])*(phi - angles[k]))/(2*width[k]*width[k]));  // GAUSSIAN
          omega[j]+=e1;                                            // SUM OVER GAUSSIANS
          domega += - e1*(phi - angles[k])/(width[k]*width[k]);    // DERIVATIVE OF GAUSSIAN
        }
        
        double dvac;  // PART OF THE DERIVATIVES OF THE ANGULAR TERM
        dvac = - 1/(sqrt(norm2(v[i])*norm2(v[j]))*sqrt(1.-alpha*alpha));
        
        for (unsigned int ix=0 ;ix<3 ;ix++) {
          domega1[j][ix] = - domega * dvac * (v[i][ix] - dotprod(v[i],v[j]) * v[j][ix]/norm2(v[j]));  // DERIVATIVE ANGULAR PART 1
          domega2[j][ix] =   domega * dvac * (v[j][ix] - dotprod(v[i],v[j]) * v[i][ix]/norm2(v[i]));  // DERIVATIVE ANGULAR PART 2
        }
 
        angtot += omega[j]*f[j];    // TOTAL OF THE ANGULAR PART

      }
    }

    double drho = 0.0;         
    double rho = 0.0;    // DENSITY SWITCHING FUNCTION

    rho = coordswitch.calculate(n,drho);   // AND DERIVATIVES           

    double zpos = pos[i][2];
    kernel(zpos);   // calculate value of kernel function kval and its derivative dkval

    if (n>0.) {           // THESE PARTS ARE USEFUL IN THE CALCULATION OF THE DERIVATIVES
      // NEW OPERATIONS TO SPEED UP
      S1 = drho * angtot;    
      S2 = -angtot/n;      
      S3 = rho*angtot/n;
      cv_val += S3 * kval[2];     // SUM THE TOTAL CV 
    } else {
      S1 = 0.;
      S2 = 0.;
      S3 = 0.;
    }

    start_i = list[i] + nmol;
    end_i = list[i] + 2*nmol;

    for (unsigned int j=0; j<ll; ++j) {    // SUM OVER NEIGHBORS
      if (j != i) {
        unsigned int start_j, end_j;
  
        start_j = list[j] + nmol;
        end_j   = list[j] + 2*nmol;
  
        unsigned int m=3;
        for (unsigned int ix=0; ix<m; ix++) {    ///  SMAC DERIVATIVES
  	// CENTER
          double temp = 0;
  	if (n > 0.0) {
  	  temp = ( (omega[j]  +  S2) * rho / n + S1) * kval[2]; // * df[j][ix];
        	
          deriv[list[i]][ix] += temp * df_a[j][ix];
          deriv[list[j]][ix] -= temp * df_b[j][ix];

        }
         
          // START, check whether the derivatives are right or not 
          double temp1 = 0;
          double temp2 = 0;
  
          if (n > 0.0) {
            temp1 = rho * f[j] * domega1[j][ix] / n * kval[2];
            temp2 = rho * f[j] * domega2[j][ix] / n * kval[2];
            
            deriv[start_i][ix]   -=   temp2;     // ONLY ANGULAR PART
            deriv[start_j][ix]   +=   temp1;     // ONLY ANGULAR PART
  	  // END
            deriv[end_i][ix]     +=   temp2;     // ONLY ANGULAR PART
            deriv[end_j][ix]     -=   temp1;     // ONLY ANGULAR PART
          }
  
        }
      }
    }

    for (unsigned int ix=0; ix<3; ix++) {
      deriv[list[i]][ix]  +=  S3 * dkval[ix];
    }
    
  }

 
  comm.Sum(deriv);
  comm.Sum(cv_val);  
  comm.Sum(virial);  

  for(unsigned i=0;i<atoms;i++) {
    setAtomsDerivatives(i, deriv[i]);
  }
  
  setBoxDerivatives(virial);
  setValue(cv_val);
  
}


double GsmacIL::norm2(Vector vect){           /// CALCULATE THE NORM OF A VECTOR
  return (vect[0]*vect[0]+vect[1]*vect[1]+vect[2]*vect[2]);
}

double GsmacIL::dotprod(Vector vect1,Vector vect2){           /// CALCULATE THE SCALAR PRODUCT OF TWO VECTORS
  return(vect1[0]*vect2[0]+vect1[1]*vect2[1]+vect1[2]*vect2[2]);
}
    
GsmacIL::~GsmacIL(){
}

  }
}


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

    class GsmacSmooth : public Colvar{
      SwitchingFunction distanceswitch;      // NEED TO DEFINE MULTIPLE SWF
      SwitchingFunction coordswitch;      // NEED TO DEFINE MULTIPLE SWF

    public:

      vector<AtomNumber> start;             // LIST OF MOLECULES NEEDED BY SMACK -> START
      vector<AtomNumber> end;               // LIST OF MOLECULES NEEDED BY SMACK -> END
      vector<AtomNumber> center1;            // LIST OF MOLECULES NEEDED BY SMACK -> CENTER1
      vector<AtomNumber> center2;            // LIST OF MOLECULES NEEDED BY SMACK -> CENTER2
      vector<AtomNumber> all_atoms;            // LIST OF MOLECULES NEEDED BY SMACK -> CENTER
      vector<double> angles;
      vector<double> width;

      unsigned int mols;                           //  total number of molecules
      unsigned int atoms;

      double zetac, zetal, treshold;

      ofstream fdbg;
      ofstream rdbg;



      struct  GsmacSmoothlist_s{       // STRUCTURE NEEDED TO CONSTRUCT THE VERLET LIST

        double rcut;             // CUT OFF OF THE LIST
        double rskin;            // SKIN TO CHECK IF THE LIST HAVE TO BE RECALCULATED
        double rskin2;            // SKIN SQUARED
        int  step;
        vector<int> nn;          // NUMBERS OF NEIGHBORS
        vector<vector<int> > ni; // LIST FOR ATOM I
        vector<Vector> pos1;      // POSITIONS CENTER MOLECULE
        vector<Vector> pos2;      // POSITIONS NEIGHBORING MOLECULES
      } GsmacSmoothlist;


      GsmacSmooth(const ActionOptions&);              //    CONSTRUCTOR
      ~GsmacSmooth();                                 //    DESTRUCTOR
      // active methods:
      static void registerKeywords( Keywords& keys );    // KEYWORDS
      virtual void calculate();                          // CALCULATE CV

      void GsmacSmooth_newlist(vector<AtomNumber>&, GsmacSmoothlist_s &GsmacSmoothlist);     // NEW VERLIST
      void GsmacSmooth_checklist(vector<AtomNumber>&, GsmacSmoothlist_s &GsmacSmoothlist);   // CHECK THE V
      double dotprod(Vector,Vector);
      double norm2(Vector);

    };

    PLUMED_REGISTER_ACTION(GsmacSmooth,"GSMACSMOOTH")

    void GsmacSmooth::registerKeywords( Keywords& keys ){

      //  OrientationSphere::registerKeywords(keys);

      // IN DEFINING THE KEYWORDS I INTENTIONALLY DEFINE EVERYTHING HERE, WITHOUT CALLING ANY OTHER
      // CLASS/FUNCTION/OBJECT SIMPLY BECAUSE i PREFER TO HAVE EVERYTHING UNDER CONTROL.
      Colvar::registerKeywords(keys);
      keys.add("atoms","CENTER1","the labels of the atoms acting as center of the molecules");

      keys.add("atoms","CENTER2","the labels of the atoms acting as center of the neighboring molecules");

      keys.add("atoms","START","the labels of the atoms acting as start of the intramolecular vector");

      keys.add("atoms","END","the labels of the atoms acting as end  of the intramolecular vector");

      keys.add("compulsory","ANGLES"," Angles that have to be used in the SMAC calculation (need to be associated with width)");

      keys.add("compulsory","N_ANGLES"," Number of angles that have to be used in the SMAC calculation (need to be associated with width)");

      keys.add("compulsory","WIDTH"," Width of the Gaussian used to calculate the angles distribution (need to be associated with angles)");

      keys.add("optional","SWITCH","This keyword is used if you want to employ an alternative to the continuous swiching function defined above. "
      "The following provides information on the \\ref switchingfunction that are available. "
      "When this keyword is present you no longer need the NN, MM, D_0 and R_0 keywords.");

      keys.add("optional","SWITCH_COORD","This keyword is used if you want to employ an alternative to the continuous swiching function defined above. "
      "The following provides information on the \\ref switchingfunction that are available. "
      "When this keyword is present you no longer need the NN, MM, D_0 and R_0 keywords.");

      keys.add("compulsory","R_CUT","Verlet lists, default r_cut = R_max");
      keys.add("compulsory","R_SKIN","Verlet lists, default r_skin = 1.25 R_max");

      keys.add("compulsory","ZETAC","crystal side boundary of the interval switching function (in fractional coordinates)");
      keys.add("compulsory","ZETAL","liquid side boundary of the interval switching function (in fractional coordinates)");
      keys.add("compulsory","TRESHOLD","Treshold value above which the molecule is considered crystalline");

      keys.remove("NOPBC");


    }

    GsmacSmooth::GsmacSmooth(const ActionOptions&ao):
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

      parseAtomList("CENTER1",center1);
      parseAtomList("CENTER2",center2);
      parseAtomList("START",start);
      parseAtomList("END",end);

      parse("ZETAC",zetac);
      parse("ZETAL",zetal);
      parse("TRESHOLD",treshold);

      mols=center1.size();
      GsmacSmoothlist.pos1.resize(mols);
      GsmacSmoothlist.pos2.resize(mols);
      GsmacSmoothlist.nn.resize(mols);
      GsmacSmoothlist.ni.resize(mols);
      for(unsigned int i=0; i < mols ; i++){
        GsmacSmoothlist.ni[i].resize(mols);
      }


      if(mols==0) error("no molecules specified");

      for(unsigned int i=0 ; i < mols ; i++){
        all_atoms.push_back(center1[i]);
      }
      for(unsigned int i=0 ; i < mols ; i++){
        all_atoms.push_back(center2[i]);
      }
      for(unsigned int i=0 ; i < mols ; i++){
        all_atoms.push_back(start[i]);
      }
      for(unsigned int i=0 ; i < mols ; i++){
        all_atoms.push_back(end[i]);
      }

      atoms=all_atoms.size();

      addValueWithDerivatives();
      setNotPeriodic();

      // The requested atoms have following shape:
      // atoms = [ center1, center2, start, end ]
      // size(atoms) = 4*mols
      requestAtoms(all_atoms);

      double r_cut,r_skin;
      parse("R_CUT",r_cut);
      parse("R_SKIN",r_skin);

      GsmacSmoothlist.rcut=r_cut;
      GsmacSmoothlist.rskin=r_skin;
      GsmacSmoothlist.rskin2=r_skin*r_skin;

      checkRead();

      log<<"  contacts are counted with cutoff "<< distanceswitch.description()<<"\n";
      log<<"  correct density is evaluated above with cutoff "<< coordswitch.description()<<"\n";

      GsmacSmoothlist.step=0;

    }


    void GsmacSmooth::calculate()
    {
      //Vector c1 = getPosition(1);
      //std::cout << "center1[1]: " << c1[0] << ", " << c1[1] << ", " <<c1[2] << '\n';
      //std::cout << "atoms: " << getNumberOfAtoms() << '\n';

      double cv_val;       // CV
      cv_val=0;
      vector<double> cv_vec(mols);  // (CV[1], ..., CV[mols])

      Tensor virial;              // VIRIAL
      virial.zero();  // no virial contribution

      // no derivatives
      vector<Vector> deriv(getNumberOfAtoms());
      Vector ze;
      ze.zero();
      fill(deriv.begin(), deriv.end(), ze);            

      if(GsmacSmoothlist.step==0) {   // CHECK IF NEIGH LIST HAVE TO BE CONSTRUCTED
        GsmacSmooth_newlist(center1,GsmacSmoothlist);
      }
      GsmacSmoothlist.step=1; GsmacSmooth_checklist(center1,GsmacSmoothlist);   // CALL NEIGH LIST IN CASE

      unsigned int stride=comm.Get_size();                        // SET THE PARALLELIZATION VARIABLES
      unsigned int rank=comm.Get_rank();
      
      stride=comm.Get_size();
      rank=comm.Get_rank();

      // Get molecule vectors for angle calculation in v:
      vector<Vector> v(getNumberOfAtoms());                          //
      for(unsigned int i = rank; i < mols; i += stride) {                   // DEFINE THE INTRAMOLECULAR VECTOR
        int start=i + (int) mols + (int) mols;
        int end=i + (int) mols + (int) mols + (int) mols;
        v[i]=pbcDistance(getPosition(start),getPosition(end));        // THIS IS DONE ONCE PER STEPS FOR ALL MOLECULES
      }
      comm.Sum(v);  // MERGING UP

      vector<vector<double> > f(mols); // SWF of all neighbors j for each mol i
      vector<double> n(mols, 0.);      // Coordination number of each mol i

      // Sum over molecules
      for(unsigned int i=rank;i<mols;i+=stride) {
        double angtot;

        f[i].resize(GsmacSmoothlist.nn[i]);             // SWFs to neighbors of molecule i
        vector<double> omega(GsmacSmoothlist.nn[i]);                            // SWITCHING FUNCTION
        angtot=0.;

        Vector dist;

        for( int j=0; j < GsmacSmoothlist.nn[i]; ++j) { // SUM OVER NEIGHBORS
          int index_j;
          index_j=GsmacSmoothlist.ni[i][j]-mols;        // TRUE INDEX OF THE J NEIGHBOR
          //      if(i==index_j){continue;}
          Vector comp;
          double modij;
          dist=pbcDistance(getPosition(i),getPosition(index_j+mols));    // DISTANCE BETWEEN THEM
          modij=dist.modulo();                                              //

          double dfunc=0.;   // CALCULATING SWITCHING FUNCTION AND
          f[i][j] = distanceswitch.calculate(modij,dfunc);
          n[i] += f[i][j];                // CALCULATING THE COORDINATION NUMBER

          // ANGLE BETWEEN THE DIRECTIVE VECTOR        THE INDEX ARE WRONG IN THE SCALAR PRODUCT! CHANGE THEM!
          double alpha;
          alpha = dotprod(v[i],v[index_j])/sqrt(norm2(v[i])*norm2(v[index_j]));    ///  CHANGE THIS PART IN MOR C++ STYLE

          if(alpha>=1.0){        // NOT SURE IF THIS CHECK IS NECESSARY ANYMORE. DOUBLECHECK!
            alpha=0.99999;
          }else if(alpha<=-1.0){
            alpha=-0.99999;
          }

          double phi;
          phi=acos(alpha);

          omega[j]=0.;
          for(unsigned int k=0 ; k < angles.size() ;k++){   // ANGULAR SERIES! NEED TO BE SUBSTITUTED WITH THE KERNEL
            double e1;
            e1 = exp(-((phi - angles[k])*(phi - angles[k]))/(2*width[k]*width[k]));    // GAUSSIAN
            omega[j]+=e1;         // SUM OVER GAUSSIANS
          }



          angtot += omega[j]*f[i][j];  // TOTAL OF THE ANGULAR PART

        }

        double drho=0.;
        double rho;       // DENSITY SWITCHING FUNCTION
        rho = coordswitch.calculate(n[i],drho);    // AND DERIVATIVES

        if(n[i]>0.){
          cv_vec[i] = rho*angtot/(n[i]);  // CV VALUE OF EACH MOLECULE
        }

      }

      // Smooth CV:
      comm.Sum(cv_vec);
      vector<double> cv_vec_smooth(mols, 0);
      for(unsigned int i=rank; i<mols; i+=stride) {  // SMOOTH CV OVER NEIGHBORS
          Vector position = getPosition(i);
          double zpos = position[2];
          double zbox = getBox()[2][2];
          double zeta = zpos/zbox;
          if ( (zetac < zeta) && (zeta < zetal) ) {
              cv_vec_smooth[i] = cv_vec[i] / n[i];
              for( int j=0; j<GsmacSmoothlist.nn[i]; ++j) {  // SUM OVER NEIGHBORS
                  int index_j;
                  index_j=GsmacSmoothlist.ni[i][j]-mols;
                  cv_vec_smooth[i] += cv_vec[index_j] * f[i][j] / n[i];
              }
          } else {
              cv_vec_smooth[i] = 0;
          }
      }

      comm.Sum(cv_vec_smooth);


      cv_val = 0;
      for(unsigned int i=rank; i<mols; i+=stride) {
          if (cv_vec_smooth[i] > treshold) {
              cv_vec_smooth[i] = 1;
          } else {
              cv_vec_smooth[i] = 0;
          }
        cv_val += cv_vec_smooth[i];
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



    void GsmacSmooth::GsmacSmooth_newlist(vector<AtomNumber> &list1, GsmacSmoothlist_s &GsmacSmoothlist)
    {
      Vector rij, test;
      double mod_rij;


      for(unsigned int j=0; j<list1.size(); ++j){
        //test=getPosition(j);
        GsmacSmoothlist.pos1[j] = getPosition(j);
        GsmacSmoothlist.pos2[j] = getPosition(j+list1.size());
      }

      unsigned int stride;
      unsigned int rank;

      stride=comm.Get_size();  //Number of processes
      rank=comm.Get_rank(); //Rank of pr

      for(unsigned int j=rank; j<list1.size(); j+=stride) {     // sum over grid
        GsmacSmoothlist.nn[j]=0; //reset nlistsize
        for(unsigned int i=0; i<list1.size(); ++i) {                                           // sum over atoms
          if(i!=j){
            rij = pbcDistance(GsmacSmoothlist.pos1[j],GsmacSmoothlist.pos2[i]);
            mod_rij=rij.modulo2();

            if (mod_rij < GsmacSmoothlist.rskin2){ //if distance < rskin
              GsmacSmoothlist.ni[j][GsmacSmoothlist.nn[j]]= i + (unsigned int) list1.size(); //index in neighlist
              GsmacSmoothlist.nn[j]++; //increment nn
            }
          }
        }
      }

      //    comm.Sum(GsmacSmoothlist.nn);
      //    comm.Sum(GsmacSmoothlist.ni[0]);
    }

    void GsmacSmooth::GsmacSmooth_checklist(vector<AtomNumber> &list1, GsmacSmoothlist_s &GsmacSmoothlist)
    {
      unsigned int j;
      double dr=(GsmacSmoothlist.rskin-GsmacSmoothlist.rcut)*0.5;
      Vector rij;

      for (j=0; j<list1.size(); ++j) {
        //check position variations of center molecules
        rij = pbcDistance(getPosition(j),GsmacSmoothlist.pos1[j]);
        if( fabs(rij[0])>dr || fabs(rij[1])>dr || fabs(rij[2])>dr ) {
          GsmacSmooth_newlist(list1,GsmacSmoothlist);
          break;
        }
        //check position variations of neighbor molecules
        rij = pbcDistance(getPosition(j+list1.size()),GsmacSmoothlist.pos2[j]);
        if( fabs(rij[0])>dr || fabs(rij[1])>dr || fabs(rij[2])>dr ) {
          GsmacSmooth_newlist(list1,GsmacSmoothlist);
          break;
        }
      }
    }

    double GsmacSmooth::norm2(Vector vect){           /// CALCULATE THE NORM OF A VECTOR
      return (vect[0]*vect[0]+vect[1]*vect[1]+vect[2]*vect[2]);
    }

    double GsmacSmooth::dotprod(Vector vect1,Vector vect2){           /// CALCULATE THE SCALAR PRODUCT OF TWO VECTORS
      return(vect1[0]*vect2[0]+vect1[1]*vect2[1]+vect1[2]*vect2[2]);
    }

    GsmacSmooth::~GsmacSmooth(){
    }

  }
}

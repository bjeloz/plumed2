/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2017 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

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
#include "MultiColvarBase.h"
#include "AtomValuePack.h"
#include "tools/NeighborList.h"
#include "core/ActionRegister.h"
#include "tools/SwitchingFunction.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD {
namespace multicolvar {

//+PLUMEDOC MCOLVAR SMACS
/*
*/
//+ENDPLUMEDOC


class Smacs : public MultiColvarBase {
private:
  double rcut2;
  SwitchingFunction switchingFunction;
public:
  static void registerKeywords( Keywords& keys );
  explicit Smacs(const ActionOptions&);
// active methods:
  virtual double compute( const unsigned& tindex, AtomValuePack& myatoms ) const ;
/// Returns the number of coordinates of the field
  bool isPeriodic() { return false; }
};

PLUMED_REGISTER_ACTION(Smacs,"SMACS")

void Smacs::registerKeywords( Keywords& keys ) {
  MultiColvarBase::registerKeywords( keys );
  //keys.use("SPECIES"); keys.use("SPECIESA"); keys.use("SPECIESB");
  keys.add("atoms","CENTER","the labels of the atoms acting as center of the molecules");
  keys.add("atoms","START","the labels of the atoms acting as start of the molecules");
  keys.add("atoms","END","the labels of the atoms acting as end of the molecules");
  keys.add("compulsory","NN","6","The n parameter of the switching function ");
  keys.add("compulsory","MM","0","The m parameter of the switching function; 0 implies 2*NN");
  keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
  keys.add("compulsory","R_0","The r_0 parameter of the switching function");
  //keys.add("optional","R_POWER","Multiply the coordination number function by a power of r, "
  //         "as done in White and Voth (see note above, default: no)");
  keys.add("optional","SWITCH","This keyword is used if you want to employ an alternative to the continuous swiching function defined above. "
           "The following provides information on the \\ref switchingfunction that are available. "
           "When this keyword is present you no longer need the NN, MM, D_0 and R_0 keywords.");
  keys.add("optional","SWITCH_COORD","This keyword is used if you want to employ an alternative to the continuous swiching function defined above. "
           "The following provides information on the \\ref switchingfunction that are available. "
           "When this keyword is present you no longer need the NN, MM, D_0 and R_0 keywords.");
  // Use actionWithDistributionKeywords
  keys.use("MEAN"); keys.use("MORE_THAN"); keys.use("LESS_THAN"); keys.use("MAX");
  keys.use("MIN"); keys.use("BETWEEN"); keys.use("HISTOGRAM"); keys.use("MOMENTS");
  keys.use("ALT_MIN"); keys.use("LOWEST"); keys.use("HIGHEST");
}

Smacs::Smacs(const ActionOptions&ao):
  Action(ao),
  MultiColvarBase(ao)
{
  std::vector<AtomNumber> all_atoms;
  std::vector<AtomNumber> center_atoms;
  std::vector<AtomNumber> start_atoms;
  std::vector<AtomNumber> end_atoms;
  ActionAtomistic::parseAtomList("CENTER", center_atoms );
  ActionAtomistic::parseAtomList("START", start_atoms );
  ActionAtomistic::parseAtomList("END", end_atoms );
  all_atoms.reserve ( center_atoms.size() + start_atoms.size() + end_atoms.size() );
  all_atoms.insert ( all_atoms.end(), center_atoms.begin(), center_atoms.end() );
  all_atoms.insert ( all_atoms.end(), start_atoms.begin(), start_atoms.end() );
  all_atoms.insert ( all_atoms.end(), end_atoms.begin(), end_atoms.end() );
  setupMultiColvarBase( all_atoms );
  std::vector<bool> catom_ind(all_atoms.size(),false);
  for(unsigned i=0; i<center_atoms.size(); ++i) {
    catom_ind[i] = true;
  }
  setAtomsForCentralAtom(catom_ind);

  // Read in the switching function
  std::string sw, errors; parse("SWITCH",sw);
  if(sw.length()>0) {
    switchingFunction.set(sw,errors);
    if( errors.length()!=0 ) error("problem reading SWITCH keyword : " + errors );
  } else {
    double r_0=-1.0, d_0; int nn, mm;
    parse("NN",nn); parse("MM",mm);
    parse("R_0",r_0); parse("D_0",d_0);
    if( r_0<0.0 ) error("you must set a value for R_0");
    switchingFunction.set(nn,mm,r_0,d_0);

  }
  log.printf("  coordination of central atom and those within %s\n",( switchingFunction.description() ).c_str() );

  //get cutoff of switching function
  double rcut = switchingFunction.get_dmax();

  // Set the link cell cutoff
  setLinkCellCutoff( rcut );
  rcut2 = rcut * rcut;

  // And setup the ActionWithVessel
  //std::vector<AtomNumber> all_atoms; setupMultiColvarBase( all_atoms );
  checkRead();
}

double Smacs::compute( const unsigned& tindex, AtomValuePack& myatoms ) const {
  // Calculate the coordination number
  double dfunc, d2, sw, d, raised;
  for(unsigned i=1; i<myatoms.getNumberOfAtoms(); ++i) {
    Vector& distance=myatoms.getPosition(i);
    if ( (d2=distance[0]*distance[0])<rcut2 &&
         (d2+=distance[1]*distance[1])<rcut2 &&
         (d2+=distance[2]*distance[2])<rcut2 &&
         d2>epsilon ) {

      //sw = switchingFunction.calculateSqr( d2, dfunc );
      log.printf("distance %d %f",i,d2);
      }
    }
  return d2;
}

}
}

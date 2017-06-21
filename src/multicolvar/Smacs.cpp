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
  int r_power;
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
  keys.use("SPECIES"); keys.use("SPECIESA"); keys.use("SPECIESB");
  keys.add("compulsory","NN","6","The n parameter of the switching function ");
  keys.add("compulsory","MM","0","The m parameter of the switching function; 0 implies 2*NN");
  keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
  keys.add("compulsory","R_0","The r_0 parameter of the switching function");
  keys.add("optional","R_POWER","Multiply the coordination number function by a power of r, "
           "as done in White and Voth (see note above, default: no)");
  keys.add("optional","SWITCH","This keyword is used if you want to employ an alternative to the continuous swiching function defined above. "
           "The following provides information on the \\ref switchingfunction that are available. "
           "When this keyword is present you no longer need the NN, MM, D_0 and R_0 keywords.");
  // Use actionWithDistributionKeywords
  keys.use("MEAN"); keys.use("MORE_THAN"); keys.use("LESS_THAN"); keys.use("MAX");
  keys.use("MIN"); keys.use("BETWEEN"); keys.use("HISTOGRAM"); keys.use("MOMENTS");
  keys.use("ALT_MIN"); keys.use("LOWEST"); keys.use("HIGHEST");
}

Smacs::Smacs(const ActionOptions&ao):
  Action(ao),
  MultiColvarBase(ao),
  r_power(0)
{

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

  //parse power
  parse("R_POWER", r_power);
  if(r_power > 0) {
    log.printf("  Multiplying switching function by r^%d\n", r_power);
    double offset = switchingFunction.calculate(rcut*0.9999, rcut2) * pow(rcut*0.9999, r_power);
    log.printf("  You will have a discontinuous jump of %f to 0 near the cutoff of your switching function. "
               "Consider setting D_MAX or reducing R_POWER if this is large\n", offset);
  }

  // Set the link cell cutoff
  setLinkCellCutoff( rcut );
  rcut2 = rcut * rcut;

  // And setup the ActionWithVessel
  std::vector<AtomNumber> all_atoms; setupMultiColvarBase( all_atoms ); checkRead();
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

      sw = switchingFunction.calculateSqr( d2, dfunc );
      if(r_power > 0) {
        d = sqrt(d2); raised = pow( d, r_power - 1 );
        accumulateSymmetryFunction( 1, i, sw * raised * d,
                                    (dfunc * d * raised + sw * r_power) * distance,
                                    (-dfunc * d * raised - sw * r_power) * Tensor(distance, distance),
                                    myatoms );
      } else {
        accumulateSymmetryFunction( 1, i, sw, (dfunc)*distance, (-dfunc)*Tensor(distance,distance), myatoms );
      }
    }
  }

  return myatoms.getValue(1);
}

}
}

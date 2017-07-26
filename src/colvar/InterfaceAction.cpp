/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2017 The plumed team
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
#include "Colvar.h"
#include "ActionRegister.h"

#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

using namespace std;

namespace PLMD {
  namespace colvar {

    //+PLUMEDOC COLVAR INTERFACEACTION
    /*
    Calculate the solid liquid interface in the z direction.

    Only the solvent atoms are needed to find the interface, therefore for this CV only the solvent atoms are considered.

    */
    //+ENDPLUMEDOC

    class InterfaceAction : public Colvar {
        bool isnotscaled;
        int N_ato_mol, N_ato_tot, N_mol_tot, binsz, threshold;
        double zetac, zetal;

      public:
        InterfaceAction(const ActionOptions&);
        virtual void calculate();
        static void registerKeywords( Keywords& keys );
    };

    PLUMED_REGISTER_ACTION(InterfaceAction,"INTERFACEACTION")

    void InterfaceAction::registerKeywords( Keywords& keys ) {
      Colvar::registerKeywords(keys);
      keys.add("atoms","SOLVENTGROUP","group of solvent atoms involved in the calculation");
      keys.add("compulsory","SOLVENTATOMOL","number of atoms in one solvent molecule");
      keys.add("compulsory","BINSZ","number of equidistant bins in the z direction for the calculation of the z dependent solvent density");
      keys.add("compulsory","THRESHOLD","threshold of number of solvent molecules per bin; value should corresponds to the one where the solid liquid interface is located");
      keys.add("compulsory","ZETAC","action distance from interface ZI into the crystal");
      keys.add("compulsory","ZETAL","action distance from interface ZI into the liquid");
    }

    InterfaceAction::InterfaceAction(const ActionOptions&ao):
      PLUMED_COLVAR_INIT(ao),
      isnotscaled(false)
    {

      /////////////////////////////////////////////
      // ~~ read and calculate all input data ~~ //
      /////////////////////////////////////////////

      // atom group
      vector<AtomNumber> atom_list;
      parseAtomList("SOLVENTGROUP", atom_list);
      
      // number of atoms per solvent molecule
      parse("SOLVENTATOMOL", N_ato_mol);

      // total number of solvent atoms
      N_ato_tot = atom_list.size();

      // total number of solvent molecules
      N_mol_tot = N_ato_tot/N_ato_mol;

      // write to log
      log.printf("Total number of solvent atoms:\t %d\n", N_ato_tot);
      log.printf("Total number of solvent molecules:\t %d\n", N_mol_tot);


      // number of bins for the construction of the solvent density histogram in the z direction
      parse("BINSZ", binsz);

      // threshold number of molecules per bin where the interface should be located
      parse("THRESHOLD", threshold);

      // write to log
      log.printf("Bin size in z direction:\t %d\n", binsz);
      log.printf("Number of solvent molecules per bin threshold:\t %d\n", threshold);


      // interval limit in crystal
      parse("ZETAC", zetac);
      
      // interval limit in liquid
      parse("ZETAL", zetal);

      // write to log
      log.printf("Interval limit in crystal:\t %d\n", zetac);
      log.printf("Interval limit in liquid:\t %d\n", zetal);
      log.printf("  \n");
      log.flush();

    }

    //////////////////////
    // ~~ calculator ~~ //
    //////////////////////
    void InterfaceAction::calculate() {
      
      Tensor virial;
      virial.zero();  // no virial contribution

    }

  }
}

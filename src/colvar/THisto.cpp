/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   See http://www.plumed-code.org for more information.

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

   CmuMD method file, 
   see (and cite) Perego, Salvalaglio, Parrinello J. Chem. Phys. 142 144113 (2015) 
   http://scitation.aip.org/content/aip/journal/jcp/142/14/10.1063/1.4917200
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

namespace PLMD{
  namespace colvar{

    //+PLUMEDOC COLVAR THISTO 
    /*
      Calculates the solute and solvent concentration in a planar shell of the box
      WARNING!!! The derivative contains only the outer boundary terms, the resulting bias force is non-conservative!!! NO METADYNAMICS, use only for CmuMD!!!
    */
    //+ENDPLUMEDOC
   
    class THisto : public Colvar {
      //      SwitchingFunction switchingFunction; //instance sw.f
      bool isnotscaled;
      bool isFirstStep;
      int  N_su, N_sv, Na_sv_permol, Na_su_permol, Na_sv, Na_su, N_mol, com_sv, com_su, nbin;
      double  nint;
      
    public:
      THisto(const ActionOptions&);
      virtual void calculate();
      static void registerKeywords(Keywords& keys);
      ofstream histo;      // write histogram
      ofstream interface;  // write interface position 
      ofstream fdbg;
    };

    PLUMED_REGISTER_ACTION(THisto,"THISTO")

// calculate solvent density, its histogram, solute density, and print histogram (4 blocks)

    void THisto::registerKeywords(Keywords& keys){
      Colvar::registerKeywords(keys);
      keys.add("atoms","GROUP","the group of atoms involved in the calculation"); //*
      keys.add("compulsory","NSV","Solvent atoms"); //*
      keys.add("optional","SOLUTE","Solute tot atoms"); //*
      keys.add("optional","NST","Solute atoms"); //*
      keys.add("optional","NZ","Interface localization: zbin");
      keys.add("optional","NINT","Interface localization: int density");
      keys.add("optional","COMST","solute COM");
      keys.add("optional","COMSV","solvent COM");
      keys.addFlag("NOSCALE",false,"use absolute length units");
      keys.remove("NOPBC");
    }

    THisto::THisto(const ActionOptions&ao):
      PLUMED_COLVAR_INIT(ao),
      //init bool parameters
      isnotscaled(false)
    {

      
      //Read atom group
      vector<AtomNumber> at_list;
      
      parseAtomList("GROUP",at_list);
      
      Na_sv_permol=1; //default
      parse("NSV",Na_sv_permol); //get number of atoms per molecule
      
      N_su=0; //default
      Na_su=0;
      Na_su_permol=1;
      
      parse("SOLUTE",Na_su); //get number of solute atoms
      parse("NST",Na_su_permol);
      
      //Solution numbers
      N_su=(int)(Na_su/Na_su_permol); //Number of solute atoms
      Na_sv=at_list.size()-Na_su; //Number of solvent atoms
      N_sv=(int)(Na_sv/Na_sv_permol); //Number of solvent molecules
      N_mol=N_sv+N_su; //Number of total molecules
      
      log.printf("Number of atoms:\tw %d\tu %d\n",Na_sv,Na_su);
      log.printf("Number of molecules:\ttot %d\t w %d\tu %d\n",N_mol,N_sv,N_su);
      
      
   
      //COM flags
      com_sv=-1;
      com_su=-1;
      parse("COMSV",com_sv); 
      parse("COMST",com_su); 
      
      nint=0.0; //default values
      nbin=100;
      parse("NINT",nint); //interface boundary concentration
      parse("NZ",nbin); //z histogram bins
      log.printf("Histogram:\tnint %lf\tnbin %d\n",nint,nbin);
      log.flush(); //DBG
      //other bool parameters 
      parseFlag("NOSCALE",isnotscaled);
      checkRead();
      addValueWithDerivatives(); 
      setNotPeriodic();
      
      //log atom lists
      log.printf("  of %d atoms\n",at_list.size());
      for(unsigned int i=0;i<at_list.size();++i){
	log.printf("  %d", at_list[i].serial());
      }
      log.printf("  \n");
      if(N_su>0){ 
	log.printf("of which the first %d are solute atoms\n",N_su);
      }
      requestAtoms(at_list);
      log.printf("  \n");
      log.flush();       
      isFirstStep=true;

      // open the histo and interface files
      histo.open("histo", std::ofstream::out | std::ofstream::trunc); // trunc deletes the old histo file if there was one
      histo << "#! bin_number bin_position solvent_conc solute_conc\n";
      histo.close();

      // open the histo and interface files
      interface.open("interface", std::ofstream::out | std::ofstream::trunc);
      interface << "#! time left_interface_pos left_interface_fract_pos right_interface_pos rigth_interface_fract_pos\n";
      interface.close();

    }

 
    // calculator
    void THisto::calculate()    
    {

      // open the histo and interface files
      histo.open("histo", std::ofstream::out | std::ofstream::app); // app appends data to existing histo file
      interface.open("interface", std::ofstream::out | std::ofstream::app);

      // center of mass vectors for the solvent and solute
      vector<Vector> com_solv(N_sv);
      vector<Vector> com_solu(N_su);


      //Parallel parameters
      
      unsigned int stride;
      unsigned int rank; 
      
      stride=comm.Get_size();  //Number of processes
      rank=comm.Get_rank(); //Rank of present process
     

      //Box size
      
      double LBC[3];
      
      for(int i=0;i<3;++i) LBC[i]=getBox()[i][i]; 
    
      
      //Histogram settings (for interface localization)
    
      //histz-array allocation
      vector<int> histz_sv(nbin,0.0);  // solvent
      vector<int> histz_su(nbin,0.0);  // solute
      int nz=0;
      //bins width (time dependent)
      double dz=LBC[2]/nbin;
      double Vbin=LBC[0]*LBC[1]*dz; //Bin volume [nm^3], needed to find the interface

      // write historgram of solute
      for(int i=rank; i<N_su; i+=stride){
        com_solu[i]=getPosition(i);
        nz=(int)(com_solu[i][2]/dz); //fill histogram
        histz_su[nz]+=1;
      }
   
      // write historgram of solvent
      for(int i=rank; i<N_sv; i+=stride){
        com_solv[i]=getPosition(Na_su+i*Na_sv_permol+com_sv);
        nz=(int)(com_solv[i][2]/dz); //fill histogram
        histz_sv[nz]+=1;
      }
 
      // Smooth histogram
      vector<double> smooth_histz_sv(nbin,0.0);
      vector<double> smooth_histz_su(nbin,0.0);
      for(int i=0; i<nbin; i+=1){
        if (i==0) {
          smooth_histz_sv[i]=(histz_sv[i]+histz_sv[i+1])/2.;
          smooth_histz_su[i]=(histz_su[i]+histz_su[i+1])/2.;
        } else if (i==(nbin-1)) {
          smooth_histz_sv[i]=(histz_sv[i]+histz_sv[i-1])/2.;
          smooth_histz_su[i]=(histz_su[i]+histz_su[i-1])/2.;
        } else {
          smooth_histz_sv[i]=(histz_sv[i]+histz_sv[i-1]+histz_sv[i+1])/3.;
          smooth_histz_su[i]=(histz_su[i]+histz_su[i-1]+histz_su[i+1])/3.;
        }
        histo << i << " " << i*dz+dz/2 << " " << smooth_histz_sv[i] << " " << smooth_histz_su[i] << "\n"; //useful to plot the solute distribution across the box
      }

      for(int i=0; i<nbin; i+=1){
        histz_sv[i] = smooth_histz_sv[i];
      }


      //Get the liquid-crystal interfaces
      double halfbin, ileft, iright, zleft, zright, zleft_frac, zright_frac, time;

      //interface finder
      isFirstStep=false;
      halfbin=(int)(LBC[2]/(2*dz));  // get the middle bin, LBC[2] corresponds to the z length of the box
      int p=0;
      int pmone=0;

      //find the crystal if it's not at the half, it finds the crystal before halfbin exceeds the limits 
      //3 adjacent bins with water concentration < than nint/3
      while((histz_sv[halfbin]+histz_sv[halfbin+1]+histz_sv[halfbin-1]) > nint*Vbin){
        p++;
        pmone=2*(p%2)-1;
        halfbin=halfbin+p*pmone; //Move through the bins
      }

      // put halfbin inside the crystal volume (3 bins, WARNING!!! parameter dependent)
      // the search algorithm starts inside the crystal
      ileft=halfbin;
      while(histz_sv[ileft] < nint*Vbin){
        ileft=ileft-1;
        if(ileft<0) ileft=ileft+nbin; //pbc on left
      }

      iright=halfbin; //WARNING!!! parameter dependent
      if(iright>=nbin) iright=iright-nbin; //pbc on right
      while(histz_sv[iright]< nint*Vbin){
        iright=iright+1;
        if(iright>=nbin) iright=iright-nbin; //pbc on right
      }

      zleft=dz*(ileft+1); //left interface coordinate
      zright=dz*(iright); //right interface coordinate

      zleft_frac = (ileft+1)/nbin;
      zright_frac = (iright)/nbin;

      time = getTime();

      interface << time << " " << zleft << " " << zleft_frac << " " << zright << " " << zright_frac << "\n";

    // close open files
    histo.close();
    interface.close();

    }

  }
}

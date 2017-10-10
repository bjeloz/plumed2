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

    //+PLUMEDOC COLVAR TNSHELL 
    /*
      Calculates the solute and solvent concentration in a planar shell of the box
      WARNING!!! The derivative contains only the outer boundary terms, the resulting bias force is non-conservative!!! NO METADYNAMICS, use only for CmuMD!!!
    */
    //+ENDPLUMEDOC
   
    class TNshell : public Colvar {
      //      SwitchingFunction switchingFunction; //instance sw.f
      bool issolute,isdelta;
      bool isnotscaled;
      bool isFirstStep;
      int storeHalfBin;
      int  N_st, N_sv, Na_sv_permol, Na_st_permol, Na_sv, Na_st, N_mol, com_sv, com_st, nbin, asymm;
      double  iD_CR, iD_F, iCR_Size, iw_force, iw_in, iw_out, co_out, co_in, co_f, nint, fixi;
      
    public:
      TNshell(const ActionOptions&);
      virtual void calculate();
      static void registerKeywords(Keywords& keys);
      double sigmon(double z, double Coff);
      double sigmoff(double z, double Coff);
      double dsig(double z, double Coff);
      ofstream fdbg;
    };

    PLUMED_REGISTER_ACTION(TNshell,"TNSHELL")

    void TNshell::registerKeywords(Keywords& keys){
      Colvar::registerKeywords(keys);
      keys.add("atoms","GROUP","the group of atoms involved in the calculation");
      keys.add("compulsory","NSV","Solvent atoms");
      keys.add("optional","SOLUTE","Solute tot atoms");
      keys.add("optional","NST","Solute atoms");
      keys.add("compulsory","DCR","CR distance");
      keys.add("compulsory","CRSIZE","CR size");
      keys.add("optional","DF","Force distance");
      keys.add("compulsory","WF","force sigma length");
      keys.add("optional","COF","force sigma cutoff");
      keys.add("optional","WIN","in sigma length");
      keys.add("optional","COIN","in sigma cutoff");
      keys.add("optional","WOUT","out sigma length");
      keys.add("optional","COOUT","out sigma cutoff");
      keys.add("optional","FIXED","fixed interface");      
      keys.add("optional","NZ","Interface localization: zbin");
      keys.add("optional","NINT","Interface localization: int density");
      keys.add("optional","COMST","solute COM");
      keys.add("optional","COMSV","solvent COM");
      keys.addFlag("NOSCALE",false,"use absolute length units");
      keys.add("optional","ASYMM","only left(smaller z) or right (larger z) considered");
      keys.addFlag("DELTA",false,"concentration gradient");
      keys.remove("NOPBC");
    }

    TNshell::TNshell(const ActionOptions&ao):
      PLUMED_COLVAR_INIT(ao),
      //init bool parameters
      isnotscaled(false)
    {
      
      //Read atom group
      vector<AtomNumber> at_list;
      
      parseAtomList("GROUP",at_list);
      
      Na_sv_permol=1; //default
      parse("NSV",Na_sv_permol); //get number of atoms per molecule
      
      N_st=0; //default
      Na_st=0;
      Na_st_permol=1;
      
      parse("SOLUTE",Na_st); //get number of solute atoms
      parse("NST",Na_st_permol);
      
      //Solution numbers
      N_st=(int)(Na_st/Na_st_permol); //Number of solute atoms
      Na_sv=at_list.size()-Na_st; //Number of solvent atoms
      N_sv=(int)(Na_sv/Na_sv_permol); //Number of solvent molecules
      N_mol=N_sv+N_st; //Number of total molecules
      
      log.printf("Number of atoms:\tw %d\tu %d\n",Na_sv,Na_st);
      log.printf("Number of molecules:\ttot %d\t w %d\tu %d\n",N_mol,N_sv,N_st);
      
      //Parameters (force position and switching function temperature)
      
      parse("DCR",iD_CR); //CR distance from interface 
      parse("CRSIZE",iCR_Size); //CR Size
      iD_F=iD_CR+iCR_Size; //initialize D_F: force distance from interface
      parse("DF",iD_F); 
      if(iD_F<iD_CR+iCR_Size){ //re-initialize D_F if inside CR
	iD_F=iD_CR+iCR_Size; 
	log.printf("D_F inside CR region, reset at the boundary");
      }
      parse("WF",iw_force); //Fermi Fun T at DF      
      co_f=20.0; //initialize cut-off in
      parse("COF",co_f); //cut-off for Fermi f
      iw_in=iw_force; //initialize w_in
      parse("WIN",iw_in); //Fermi Fun T at CRin
      co_in=co_f; //initialize cut-off in
      parse("COIN",co_in); //cut-off for Fermi f
      iw_out=iw_force; //initialize w_out
      parse("WOUT",iw_out); //Fermi Fun T at CRout
      co_out=co_f; //initialize cut-off in
      parse("COOUT",co_out); //cut-off for Fermi f
      
      log.printf("Geometry:\tD_CR %lf\tCR_size %lf\tD_F %lf\n",iD_CR,iCR_Size,iD_F);
      log.flush();
      
      fixi=-1.; //default fixed inactive
      parse("FIXED",fixi); //fixed interface coordinate (always scaled)
      if(fixi>=0){
	log.printf("Fixed interface at:\t %lf\n L_box",fixi);
      }
      
      //Asymmetry
      
      asymm=0; //default no asymmetry
      parse("ASYMM",asymm); //cut-off for Fermi f
      if(asymm<0){
	log.printf("Only left CR considered");
      }else if(asymm>0){
	log.printf("Only right CR considered");
      }
      
      parseFlag("DELTA",isdelta);
      if(isdelta){
	log.printf("Difference between right and left CR calculated");
      }
   
      //COM flags
      com_sv=-1;
      com_st=-1;
      parse("COMSV",com_sv); 
      parse("COMST",com_st); 
      
      nint=0.0; //default values
      nbin=100;
      parse("NINT",nint); //interface boundary concentration
      parse("NZ",nbin); //z histogram bins
      if(fixi<0){
	log.printf("Histogram:\tnint %lf\tnbin %d\n",nint,nbin);
      }
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
      if(N_st>0){ 
	log.printf("of which the first %d are solute atoms\n",N_st);
      }
      requestAtoms(at_list);
      log.printf("  \n");
      log.flush();       
      isFirstStep=true;
    }

    double TNshell::sigmon(double z, double Coff){
      double sig;
      if( z < -Coff){
	sig=0.0;
      }else if(z > Coff){
	sig=1.0;
      }else{
	sig=1.0/(exp(-z)+1.0);
      }
      return(sig);
    }
  

    double TNshell::sigmoff(double z, double Coff){
      double sig;
      if( z < -Coff){
	sig=1.0;
      }else if(z > Coff){
	sig=0.0;
      }else{
	sig=1.0/(exp(z)+1.0);
      }
      return(sig);
    }
    
    double TNshell::dsig(double z, double Coff){
      double dsig;
      if(fabs(z) > Coff){
	dsig=0.0;
      }else{
	dsig=0.5/(1.0+cosh(z));
      }
      return(dsig);
    }


    
    // calculator
    void TNshell::calculate()    
    {
      
      double n_CR,n_CRr,n_CRl;
      Tensor virial;
      
      virial.zero();  //no virial contribution
      
      //Vector deriv;
      
      vector<Vector> deriv(getNumberOfAtoms());
      vector<Vector> com_solv(N_sv);
      Vector diff;

      Vector ze;
      ze.zero();
      //init derivatives
      fill(deriv.begin(), deriv.end(), ze);

      //Parallel parameters
      
      unsigned int stride;
      unsigned int rank; 
      
      stride=comm.Get_size();  //Number of processes
      rank=comm.Get_rank(); //Rank of present process
      
      //Solvent position matrix allocation 
      vector<Vector> solve_x(Na_sv_permol);
      //Solvent mass array allocation 
      vector<double> solve_m(Na_sv_permol);
      
      //Solvent masses and total mass
      double M_sv=0.0;
      for(int i=0;i<Na_sv_permol;++i){
	solve_m[i]=getMass(Na_st+i); //the first Na_st are skipped
	M_sv += solve_m[i];
      }

      //Box size
      
      double LBC[3];
      
      for(int i=0;i<3;++i) LBC[i]=getBox()[i][i]; 
    
      double D_CR,CR_Size,D_F,w_force,w_in,w_out;
      
      double fix_int=0;
      if(fixi>=0.0) fix_int=LBC[2]*fixi; //fixed interface

      if(!isnotscaled){ //rescale input distances
	
	D_CR=LBC[2]*iD_CR;
	CR_Size=LBC[2]*iCR_Size;
	D_F=LBC[2]*iD_F;
	w_force=LBC[2]*iw_force;
	w_in=LBC[2]*iw_in;
	w_out=LBC[2]*iw_out;
	
      }
      
      //rescale the cut-offs
      
      double VCR;
      if(asymm==0){
	VCR=2*LBC[0]*LBC[1]*CR_Size; //CR volume 
      }else{
	VCR=LBC[0]*LBC[1]*CR_Size; //CR volume 
      }
  
     
      
      //Histogram settings (for interface localization)
    
      //histz-array allocation
      vector<int> histz(nbin,0.0);
      int nz=0;
      //bins width (time dependent)
      double dz=LBC[2]/nbin;
      double Vbin=LBC[0]*LBC[1]*dz; //Bin volume [nm^3]  
    
      //center of mass vector
      for(int i=rank; i<N_sv; i+=stride){
	com_solv[i].zero();
	//center of mass
	if (com_sv<0) {
	  solve_x[0] = getPosition(Na_st+i*Na_sv_permol);
	  for(int j=1;j<Na_sv_permol;++j){
	    solve_x[j] = getPosition(Na_st+i*Na_sv_permol+j);
	    diff = pbcDistance(solve_x[0],solve_x[j]);
	    com_solv[i] += solve_m[j]*diff;
	  }      
	  com_solv[i] = com_solv[i] / M_sv + solve_x[0]; 
	  //impose PBC on com (useless on x and y for now!!!)
	  if(com_solv[i][2]<0) com_solv[i][2]=com_solv[i][2]+LBC[2];
	  if(com_solv[i][2]>=LBC[2]) com_solv[i][2]=com_solv[i][2]-LBC[2];
	}else{
	  //no com
	  com_solv[i]=getPosition(Na_st+i*Na_sv_permol+com_sv);
	}
	
	if(fixi<0){
	  nz=(int)(com_solv[i][2]/dz); //fill histogram
	  histz[nz]+=1;
	}
       // log.printf("nz=%d\n",nz);
      }
 
      // Smooth histogram
      vector<double> smooth_histz(nbin,0.0);
       for(int i=0; i<nbin; i+=1){
         if (i==0) {
           smooth_histz[i]=(histz[i]+histz[i+1])/2.;
         } else if (i==(nbin-1)) {
           smooth_histz[i]=(histz[i]+histz[i-1])/2.;
         } else {
           smooth_histz[i]=(histz[i]+histz[i-1]+histz[i+1])/3.;
         }
         //log.printf("Histogram= %d %f\n",i,smooth_histz[i]); //useful to plot the solute distribution across the box
       }

// include code block here between 314-346 to calcutate solute distribution in z direction, average every 5 ps; on the fly calculation

       for(int i=0; i<nbin; i+=1){
            histz[i] = smooth_histz[i];
       }

      //communicate
      comm.Sum(histz);
      comm.Sum(com_solv);

      //Get the liquid-crystal interfaces
      double halfbin, ileft, iright, zleft, zright;
      
       //interface finder
       if(fixi<0){
         if (isFirstStep) {
            isFirstStep=false;
            halfbin=(int)(LBC[2]/(2*dz));
            int p=0;
            int pmone=0;
      
	    //find the crystal if it's not at the half, it finds the crystal before halfbin exceeds the limits 
     	    //3 adjacent bins with water concentration < than nint/3
   	    while((histz[halfbin]+histz[halfbin+1]+histz[halfbin-1]) > nint*Vbin){
	         p++;
	         pmone=2*(p%2)-1;
	         halfbin=halfbin+p*pmone; //Move through the bins
            }
         } else {
           halfbin=storeHalfBin;
         }

	
	//put halfbin inside the crystal volume (3 bins, WARNING!!! parameter dependent)
	
	ileft=halfbin;
	while(histz[ileft] < nint*Vbin){
	  ileft=ileft-1;
	  if(ileft<0) ileft=ileft+nbin; //pbc on left
	}
	
        iright=halfbin; //WARNING!!! parameter dependent
	if(iright>=nbin) iright=iright-nbin; //pbc on right
	while(histz[iright]< nint*Vbin){
	  iright=iright+1;
	  if(iright>=nbin) iright=iright-nbin; //pbc on right
	}
      
        storeHalfBin=(ileft+iright)/2;
	zleft=dz*(ileft+1); //left interface coordinate
	zright=dz*(iright); //right interface coordinate
      }else{
	zleft=fix_int;
	zright=fix_int;
      }
     
      //Fermi function parameters
      double ZCRrin, ZCRrout, ZCRlin, ZCRlout, ZFright, ZFleft;
      ZCRlin=zleft-D_CR;
      ZCRlout=zleft-D_CR-CR_Size;
      ZFleft=zleft-D_F;
      
      ZCRrin=zright+D_CR;
      ZCRrout=zright+D_CR+CR_Size;
      ZFright=zright+D_F;
      
      //Evaluate concentration and derivatives
      //if isolute is true C counts the solute molecules, else the solvent ones
      
      n_CR=0.0;
      n_CRr=0.0;
      n_CRl=0.0;
      
      double zin,zout,n_lx,n_rx,n_x,zl,zr,dfunc,dl,dr;
      int k;
      if(N_st == 0){ //if solvent specie is restrained
	for(int i=rank; i<N_sv; i+=stride){
	  //Fermi-like weighting
	  dfunc=0;
	  dl=0;
	  dr=0;
	  n_lx=0;
	  n_rx=0;
	  //left-side sigma
	  if(asymm<=0){
	    zin=(com_solv[i][2]-ZCRlin)/w_in;
	    zout=(com_solv[i][2]-ZCRlout)/w_out;
	    //with periodic image, sigma on at zout, off at zin
	    n_lx=sigmon(zout,co_out)*sigmoff(zin,co_in)+sigmon(zout-LBC[2]/w_out,co_out)*sigmoff(zin-LBC[2]/w_in,co_in);
	    
	    //Derivatives (only outer boundary derivatives!!!)
	    zl=(com_solv[i][2]-ZFleft)/w_force;
	    dl=(dsig(zl,co_f)+dsig(zl-LBC[2]/w_force,co_f))/w_force;	
	  }
	  //right-side sigma
	  if(asymm>=0){
	    zin=(com_solv[i][2]-ZCRrin)/w_in;
	    zout=(com_solv[i][2]-ZCRrout)/w_out;
	    //with periodic image, sigma on at zin, off at zout
	    n_rx=sigmon(zin,co_in)*sigmoff(zout,co_out)+sigmon(zin+LBC[2]/w_in,co_in)*sigmoff(zout+LBC[2]/w_out,co_out);
	    
	    zr=(com_solv[i][2]-ZFright)/w_force;
	    dr=(-dsig(zr,co_f)-dsig(zr+LBC[2]/w_force,co_f))/w_force;
	  }

	  //update CV (for now this is the number of molcules)
	  n_CRr+=n_rx;
	  n_CRl+=n_lx;
	  //n_CR+=n_x;
	  
	  if(isdelta){
	    dfunc=dr-dl;
	  }else{
	    dfunc=dr+dl;
	  }
	  
	  
	  if(com_sv<0){ //com coordinates
	    for(int l=0; l<Na_sv_permol; ++l){
	      k=Na_st+i*Na_sv_permol+l; //atom counter
	      deriv[k][2]=getMass(k)/M_sv*(dfunc/VCR); //com affects the derivatives
	    }
	  }else{//single atom coordinates
	    k=Na_st+i*Na_sv_permol+com_sv ; //atom counter (just the derivatives with respect to "com" atom coordinates)	   
	    deriv[k][2]=dfunc/VCR;
	  }
	}
	vector<Vector>().swap(com_solv);
      }else{ //if solute specie is restrained
	
	vector<Vector> com_solut(N_st);
	//Solute position matrix allocation 
	vector<Vector> solut_x(Na_st_permol);
	//Solute mass array allocation 
	vector<double> solut_m(Na_st_permol);
	//Solute masses and total mass
	double M_st=0.0;
	for(int i=0;i<Na_st_permol;++i){
	  solut_m[i]=getMass(i);
	  M_st += solut_m[i];
	}
	
	for(int i=rank; i<N_st; i+=stride){

	  dfunc=0;
	  dl=0;
	  dr=0;
	  n_lx=0;
	  n_rx=0;
	  com_solut[i].zero();
	  //center of mass
	  if (com_st<0) {
	    solut_x[0] = getPosition(i*Na_st_permol);
	    for(int j=1; j<Na_st_permol; ++j){
	      solut_x[j] = getPosition(i*Na_st_permol+j);
	      diff = pbcDistance(solut_x[0],solut_x[j]);
	      com_solut[i] += solut_m[j]*diff;
	    }      
	    com_solut[i] = com_solut[i] / M_st + solut_x[0]; 
	    //PBC (only orthorhombic!!!) 
	    //Only on z
	    if(com_solut[i][2]<0) com_solut[i][2]=com_solut[i][2]+LBC[2];
	    if(com_solut[i][2]>=LBC[2]) com_solut[i][2]=com_solut[i][2]-LBC[2];
	  }else{
	    com_solut[i]=getPosition(i*Na_st_permol+com_st);
	  }
	  
	  //Fermi-like weighting
	  if(asymm<=0){
	    //left-side sigma
	    zin=(com_solut[i][2]-ZCRlin)/w_in;
	    zout=(com_solut[i][2]-ZCRlout)/w_out;
	    //with periodic image
	    n_lx=sigmon(zout,co_out)*sigmoff(zin,co_in)+sigmon(zout-LBC[2]/w_out,co_out)*sigmoff(zin-LBC[2]/w_in,co_in);
	    
	    zl=(com_solut[i][2]-ZFleft)/w_force;
	    dl=(dsig(zl,co_f)+dsig(zl-LBC[2]/w_force,co_f))/w_force;	
	  }
	  
	  if(asymm>=0){
	    //right-side sigma
	    zin=(com_solut[i][2]-ZCRrin)/w_in;
	    zout=(com_solut[i][2]-ZCRrout)/w_out;
	    
	    //with periodic image
	    n_rx=sigmon(zin,co_in)*sigmoff(zout,co_out)+sigmon(zin+LBC[2]/w_in,co_in)*sigmoff(zout+LBC[2]/w_out,co_out);
	  
	    zr=(com_solut[i][2]-ZFright)/w_force;
	    dr=(-dsig(zr,co_f)-dsig(zr+LBC[2]/w_force,co_f))/w_force;
	  }

	  //update CV (for now this is the number of molcules)
	  n_CRr+=n_rx;
	  n_CRl+=n_lx;
	  //n_CR+=n_x;
	  
	  if(isdelta){
	    dfunc=dr-dl;
	  }else{
	    dfunc=dr+dl;
	  }
	  
	  
	  if(com_st<0){ //com coordinates
	    for(int l=0; l<Na_st_permol; ++l){
	      k=i*Na_st_permol+l; //atom counter
	      deriv[k][2] = getMass(k)/M_st*(dfunc/VCR);
	    }
	  }else{//single atom coordinates
	    k=i*Na_st_permol+com_st ; //atom counter (just the derivatives with respect to "com" atom coordinates)
	    deriv[k][2] = dfunc/VCR;
	  }
	  
	}
	vector<Vector>().swap(com_solut);
      }
      
      comm.Sum(deriv);
      comm.Sum(n_CRr);
      comm.Sum(n_CRl);
      comm.Sum(virial);
      int Natot=Na_st+Na_sv;
      for(int i=0; i< Natot; ++i){
	setAtomsDerivatives(i, deriv[i]);
      }
      
      vector<Vector>().swap(deriv);
      if(isdelta){
	n_CR=n_CRr-n_CRl;
      }else{
	n_CR=n_CRr+n_CRl;
      }
      setValue(n_CR/VCR);
      setBoxDerivatives(virial);
      //setBoxDerivativesNoPbc();
    }
  }  
}    

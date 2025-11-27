// Need to write a PDF about this biz.
// All quantities derived from Fetter and Walecka Ch.1 Example
#include "electronGasSPBasis.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>

// r_0 is the interparticle distance radius V/N = (4/3)pi r_0^3
// r_s is the Wigner-Seitz radius  r_s = r_0/r_bohr
electronGasSPBasis::electronGasSPBasis(double r_sIn, int shellMaxIn, int NparticlesIn){
  this->r_s = r_sIn;  
  this->nMax = ceil(shellMaxIn/3.);
  this->shellMax = shellMaxIn;
  this->Nparticles = NparticlesIn;

  int E;
  int count = 0;
  int iEnergy = 0;
  int nonZeroShell;
  int shell = 0;
  int nMaxActual = 0;
  // This is related to a cool problem in number theory,
  // I do not want to increment shell if it cannot be represented
  // by the sum of three squares. MOST numbers can be, except 7, 15, 28... 
  // numbers of the form (4^m)*(8k+7). However,
  // ALL numbers can be represented by the sum of four squares!
  // I think this while loop is a good enough solution, since not 
  // many iterations will be wasted. Think harder about an algorithm later.

  while( shell <= this->shellMax ){
    nonZeroShell = 0;
    for( int nx = -this->nMax; nx <= this->nMax; nx++){    
      for( int ny = -this->nMax; ny <= this->nMax; ny++){	
	for( int nz = -this->nMax; nz <= this->nMax; nz++){	  
	  E = nx*nx + ny*ny + nz*nz;
	  if( E == iEnergy ){
	    nonZeroShell = 1;
	    if( nx > nMaxActual ){
	      nMaxActual = nx;
	    }
	    for( int sz = -1; sz <= 1; sz = sz+2){
	      // for( int tz = -1; tz <= tzMaxTemp; tz = tz+2){	
	      count++;
		// } // end tz loop
	    } // end sz loop    
	  } // end if
	} // end nz loop 
      } // end ny loop
    } // end nx loop
    iEnergy++;
    if( nonZeroShell == 1 ){
      shell++;
    }
  } // end shell loop  
  
  this->nMax = nMaxActual;
  this->Nspstates = count;
}


void electronGasSPBasis::generateIndexMap(){

  this->indexMap = new int* [Nspstates];
  for(int i = 0; i < this->Nspstates; i++){
    this->indexMap[i] = new int[4];
  } // maybe delete this at some point.
  
  int E;
  int count = 0;
  int iEnergy = 0;;
  int nonZeroShell;
  int shell = 0;

  while( shell <= this->shellMax ){
    nonZeroShell = 0;
    for( int nx = -this->nMax; nx <= this->nMax; nx++){    
      for( int ny = -this->nMax; ny <= this->nMax; ny++){	
	for( int nz = -this->nMax; nz <= this->nMax; nz++){	  
	  E = nx*nx + ny*ny + nz*nz;
	  if( E == iEnergy ){
	    nonZeroShell = 1;
	    for( int sz = -1; sz <= 1; sz = sz+2){	      
	      this->indexMap[count][0] = nx;
	      this->indexMap[count][1] = ny;
	      this->indexMap[count][2] = nz;
	      this->indexMap[count][3] = sz;
	      
	      count++;
	      
	    } // end sz loop    
	  } // end if
	} // end nz loop 
      } // end ny loop
    } // end nx loop
    iEnergy++;
    if( nonZeroShell == 1 ){
      shell++;
    }
  } // end shell loop  

} // end generateIndexMap
 
   
// requires generate index map first
void electronGasSPBasis::generateBasis(){
  this->spEnergy = new double[this->Nspstates];
         
  // Want things in Hartrees and Angstroms.
  // 27.21138505 eV/Hartree, 10,000 Angstroms/micron
  double hbarc = 0.1973269788; // eV micron from NIST site
  double eVs_in_a_Hartree = 27.21138505; // eV Wikipedia
  hbarc = hbarc*10000/eVs_in_a_Hartree; // Hartree Angstrom
  double massc2 = 0.5109989461; // MeV wikipedia electron mass 
  massc2 = massc2*1000000/eVs_in_a_Hartree; // Hartree
  double fine_struc = 1.0/137.035999139; // Wikipedia


  double r_bohr = hbarc/(massc2*fine_struc); // Angstrom
  double r_0 = this->r_s*r_bohr;
  double prefactor = hbarc*hbarc/(2.*massc2);
  this->L = pow( (4./3.)*M_PI*this->Nparticles, 1./3.)*r_0; // Angstrom
  this->e_sq = fine_struc*hbarc; // Hartree Angstrom

 
  std::cout << "mass*c^2: " << massc2 << " Hartree" << std::endl;
  std::cout << "hbar^2/(2m) = " << prefactor << " Hartree*angstrom^2" << std::endl;
  std::cout << "r_0 = " << r_0 << " angstrom, V/N= " << (4./3.)*M_PI*r_0*r_0*r_0 << " angstrom^3" << std::endl;
  std::cout << "k_fermi " << std::setprecision(16) << pow(9.*M_PI/4.,1./3.)/r_0 << " angstrom^-1" << std::endl;

  double E;
  
  // k = 2*pi*n/L
  // L = V^(1/3)
 
  for(int i = 0; i < this->Nspstates; i++){
    E = ( prefactor*4*M_PI*M_PI/(this->L*this->L) ) * 
      (this->indexMap[i][0]*this->indexMap[i][0] + 
       this->indexMap[i][1]*this->indexMap[i][1] + 
       this->indexMap[i][2]*this->indexMap[i][2]);       
    this->spEnergy[i] = E;
  }  
} // end generateBasis


int electronGasSPBasis::checkSympqrs(int p, int q, int r, int s){
  int result = 0;
  //int ** indexMap = SPbasis.indexMap;
  // check spin
  if( (this->indexMap[p][3] + this->indexMap[q][3] == 
       this->indexMap[r][3] + this->indexMap[s][3]) ){
    // check momentum conservation
    if( (this->indexMap[p][0] + this->indexMap[q][0] == 
	 this->indexMap[r][0] + this->indexMap[s][0]) && 
	(this->indexMap[p][1] + this->indexMap[q][1] == 
	 this->indexMap[r][1] + this->indexMap[s][1]) && 
	(this->indexMap[p][2] + this->indexMap[q][2] == 
	 this->indexMap[r][2] + this->indexMap[s][2] ) ) {
      result = 1;
    }    
  }
  return result;
} // end checkSympqrs

int electronGasSPBasis::checkModSympqrs(int p, int q, int r, int s){
  int result = 0;
  // check spin and isospin
  if( (this->indexMap[p][3] - this->indexMap[q][3] == 
       this->indexMap[r][3] - this->indexMap[s][3]) ){
    // check momentum conservation
    if( (this->indexMap[p][0] - this->indexMap[q][0] == 
	 this->indexMap[r][0] - this->indexMap[s][0]) && 
	(this->indexMap[p][1] - this->indexMap[q][1] == 
	 this->indexMap[r][1] - this->indexMap[s][1]) && 
	(this->indexMap[p][2] - this->indexMap[q][2] == 
	 this->indexMap[r][2] - this->indexMap[s][2] ) ) {
      result = 1;
    }    
  }
  return result;
} // end checkSympqrs

int electronGasSPBasis::checkChanSym(int p, int q, int ichan){
  int result = 0;

  if( (this->indexMap[p][0] + this->indexMap[q][0] == 
       this->chanValue[ichan].chanNx) &&
      (this->indexMap[p][1] + this->indexMap[q][1] == 
       this->chanValue[ichan].chanNy) &&
      (this->indexMap[p][2] + this->indexMap[q][2] == 
       this->chanValue[ichan].chanNz) &&
      (this->indexMap[p][3] + this->indexMap[q][3] == 
       this->chanValue[ichan].chanSz) ){
    result = 1;
  }
  return result;
} // end checkChanSym

int electronGasSPBasis::checkChanSym2(int p, int q, int chanNx, int chanNy, int chanNz, int chanSz){
  int result = 0;
 
  if( (this->indexMap[p][0] + this->indexMap[q][0] == chanNx) &&
      (this->indexMap[p][1] + this->indexMap[q][1] == chanNy) &&
      (this->indexMap[p][2] + this->indexMap[q][2] == chanNz) &&
      (this->indexMap[p][3] + this->indexMap[q][3] == chanSz) ){
   
    result = 1;
  }
  return result;
} // end checkChanSym

int electronGasSPBasis::checkChanModSym2(int p, int q, int chanNx, int chanNy, int chanNz, int chanSz){
  int result = 0;
   
  if( (this->indexMap[p][0] - this->indexMap[q][0] == chanNx) &&
      (this->indexMap[p][1] - this->indexMap[q][1] == chanNy) &&
      (this->indexMap[p][2] - this->indexMap[q][2] == chanNz) &&
      (this->indexMap[p][3] - this->indexMap[q][3] == chanSz) ){
   
    result = 1;
  }
  return result;
} // end checkChanSym


int electronGasSPBasis::checkChanModSym(int p, int q, int ichan){
  int result = 0;

  if( (this->indexMap[p][0] - this->indexMap[q][0] == 
       this->chanModValue[ichan].chanNx) &&
      (this->indexMap[p][1] - this->indexMap[q][1] == 
       this->chanModValue[ichan].chanNy) &&
      (this->indexMap[p][2] - this->indexMap[q][2] == 
       this->chanModValue[ichan].chanNz) &&
      (this->indexMap[p][3] - this->indexMap[q][3] == 
       this->chanModValue[ichan].chanSz) ){
    result = 1;
  }
  return result;
} // end checkChanModSym



void electronGasSPBasis::setUpTwoStateChannels(){
  int channelNmax = 2*this->nMax;  
  int channelSzMax = 2;
    
  // if(this->tzMax == 1){
  //   channelTzMax = -2; 
  //   modChannelTzMax = 0;
  // } else if(this->tzMax == 2){
  //   channelTzMax = 2;    
  //   modChannelTzMax = 2;
  // } else {
  //   channelTzMax = 8000;   
  //   modChannelTzMax = 8000;
  //   std::cout << "wrong tzMax" << std::endl;
  // }
  
  int channelCount = 0;
  for(int ichanNx = -channelNmax; ichanNx <= channelNmax; ichanNx++){
    for(int ichanNy = -channelNmax; ichanNy <= channelNmax; ichanNy++){
      for(int ichanNz = -channelNmax; ichanNz <= channelNmax; ichanNz++){
	for(int ichanSz = -channelSzMax; ichanSz <= channelSzMax; ichanSz = ichanSz +2){
	  // for(int ichanTz = -2; ichanTz <= channelTzMax; ichanTz = ichanTz +2){   
	    for(int p = 0; p < this->Nspstates; p++){	      
	      for(int q = 0; q < this->Nspstates; q++){
		if( this->checkChanSym2(p,q,ichanNx,ichanNy,ichanNz,ichanSz) == 1){
		  channelCount++;
		  goto loopbreak1;
		} // end if
	      } // end q
	    } // end p
	  loopbreak1:;
	   
	    // } // end ichantz
	} // end ichansz
      } // end ichannz
    } // end ichanny
  } // end ichannx

  //int nonZeroChan;   
  this->Nchannels = channelCount;
  this->chanValue = new channelBundle[channelCount];
  this->chanModValue = new channelBundle[channelCount];

  channelCount = 0;
  for(int ichanNx = -channelNmax; ichanNx <= channelNmax; ichanNx++){
    for(int ichanNy = -channelNmax; ichanNy <= channelNmax; ichanNy++){
      for(int ichanNz = -channelNmax; ichanNz <= channelNmax; ichanNz++){
	for(int ichanSz = -channelSzMax; ichanSz <= channelSzMax; ichanSz = ichanSz +2){
	  //for(int ichanTz = -2; ichanTz <= channelTzMax; ichanTz = ichanTz +2){
 	     
	    for(int p = 0; p < this->Nspstates; p++){
	      for(int q = 0; q < this->Nspstates; q++){
		if( this->checkChanSym2(p,q,ichanNx,ichanNy,ichanNz,ichanSz) == 1){
		  this->chanValue[channelCount].chanNx = ichanNx;
		  this->chanValue[channelCount].chanNy = ichanNy;
		  this->chanValue[channelCount].chanNz = ichanNz;
		  this->chanValue[channelCount].chanSz = ichanSz;
		  // this->chanValue[channelCount].chanTz = ichanTz;
		  channelCount++; 
		  goto loopbreak2;
		} // end if
	      } // end q
	    } // end p
	  loopbreak2:;  
	    
	    //} // end ichantz
	} // end ichansz
      } // end ichannz
    } // end ichanny
  } // end ichannx
   
      
  channelCount = 0;
  for(int ichanNx = -channelNmax; ichanNx <= channelNmax; ichanNx++){
    for(int ichanNy = -channelNmax; ichanNy <= channelNmax; ichanNy++){
      for(int ichanNz = -channelNmax; ichanNz <= channelNmax; ichanNz++){
	for(int ichanSz = -channelSzMax; ichanSz <= channelSzMax; ichanSz = ichanSz +2){
	  //for(int imodChanTz = -modChannelTzMax; imodChanTz <= modChannelTzMax; imodChanTz = imodChanTz +2){	    
	     
	    for(int p = 0; p < this->Nspstates; p++){
	      for(int q = 0; q < this->Nspstates; q++){
		if( this->checkChanModSym2(p,q,ichanNx,ichanNy,ichanNz,ichanSz) == 1){
		  this->chanModValue[channelCount].chanNx = ichanNx;
		  this->chanModValue[channelCount].chanNy = ichanNy;
		  this->chanModValue[channelCount].chanNz = ichanNz;
		  this->chanModValue[channelCount].chanSz = ichanSz;
		  // this->chanModValue[channelCount].chanTz = imodChanTz;
		  channelCount++;
		  goto loopbreak3;
		} // end if
	      } // end q
	    } // end p
	  loopbreak3:;
	    
	    // } // end ichantz
	} // end ichansz
      } // end ichannz
    } // end ichanny
  } // end ichannx
} // end calcNumChannels

void electronGasSPBasis::printBasis(){
  
  std::cout << "SPBasis:" << std::endl;
  std::cout << "i nx ny nz sz E" << std::endl;  
  for( int p = 0; p < this->Nspstates; p++ ) {    
    std::cout << p << " " << this->indexMap[p][0] << " " << this->indexMap[p][1] << " " << this->indexMap[p][2] << " " << this->indexMap[p][3] << " " << this->spEnergy[p] << std::endl;    
  }
} // end printBasis

void electronGasSPBasis::deallocate(){

  for(int i = 0; i < this->Nspstates; i++){
    delete [] this->indexMap[i];
  }
  delete [] this->indexMap;
  delete [] this->spEnergy;
  delete [] this->chanValue;   
  delete [] this->chanModValue;
} // end deallocate

double electronGasSPBasis::calc_TBME(int p, int q, int r, int s){
  int *qi = this->indexMap[p];  
  int *qj = this->indexMap[q];  
  int *qk = this->indexMap[r];
  int *ql = this->indexMap[s];
  double vout = 0.;

  // spin deltas  
  if( qi[3] == qk[3] && qj[3] == ql[3] ){
  
  int initialMomx, initialMomy, initialMomz;
  int finalMomx, finalMomy, finalMomz;
 
  initialMomx = qi[0] + qj[0];
  initialMomy = qi[1] + qj[1];
  initialMomz = qi[2] + qj[2];
  finalMomx = qk[0] + ql[0];
  finalMomy = qk[1] + ql[1];
  finalMomz = qk[2] + ql[2];
  
   // Momentum Conservation Checks.
  if( initialMomx == finalMomx && initialMomy == finalMomy && initialMomz == finalMomz ){
   
    double prefactor = this->e_sq/(this->L*this->L*this->L); // hard code this later?

    double relMom = 0.0; 
    double term;
    double mu = 0.0;
    for( int i = 0; i < 3; i++){
      //ki[i] = 2*M_PI*qi[i]/L;
      //kj[i] = 2*M_PI*qj[i]/L;
      // kk[i] = 2*M_PI*qk[i]/L;
      // kl[i] = 2*M_PI*ql[i]/L;
      term = 2*M_PI*(qi[i] - qk[i])/this->L;
      relMom += term*term;
    }
    
    // get rid of 0 momentum transfer modes
    if(relMom < 1.e-5){
      vout = 0.0;
    } else {
      vout = prefactor*4*M_PI/(relMom + mu*mu);
    }

  } // end Momentum Conversation Check
  } // end spin check
  return vout;
}

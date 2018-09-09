// Need to write a PDF about this biz.

#include "infMatterSPBasis.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>

infMatterSPBasis::infMatterSPBasis(double densityIn, int tzMaxIn, int shellMaxIn, int NparticlesIn){
  this->density = densityIn;
  this->tzMax = tzMaxIn;
  this->nMax = ceil(shellMaxIn/3.);
  this->shellMax = shellMaxIn;
  this->Nparticles = NparticlesIn;

  int tzMaxTemp = this->tzMax - 1; 
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
	      for( int tz = -1; tz <= tzMaxTemp; tz = tz+2){	
		count++;
	      } // end tz loop
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


void infMatterSPBasis::generateIndexMap(){

  this->indexMap = new int* [Nspstates];
  for(int i = 0; i < this->Nspstates; i++){
    this->indexMap[i] = new int[5];
  } // maybe delete this at some point.

  int tzMaxTemp = this->tzMax - 1;
  //int shellMax = 3*this->nMax*this->nMax;
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
	      for( int tz = -1; tz <= tzMaxTemp; tz = tz+2){	
		this->indexMap[count][0] = nx;
		this->indexMap[count][1] = ny;
		this->indexMap[count][2] = nz;
		this->indexMap[count][3] = sz;
		this->indexMap[count][4] = tz;
		count++;
	      } // end tz loop
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
void infMatterSPBasis::generateBasis(){
  this->spEnergy = new double[this->Nspstates];
  double hbarc = 197.3269788; // MeVfm
  //double m_neutronc2 = 939.565378; // MeV // accurate value
  double m_neutronc2 = 939.565; // MeV Gaute's value
  //double m_protonc2 = 938.272046; // MeV
  //double m_symmc2 = (m_neutronc2 + m_protonc2)*0.5;
  double massc2;
  if( tzMax == 1){
    massc2 = m_neutronc2;
  } else if( tzMax == 2){
    //massc2 = m_symmc2;
    massc2 = m_neutronc2;
  } else {
    std::cout << "tzMax = 1 or 2 ONLY" << std::endl;
  }
 
  double prefactor = hbarc*hbarc/(2.*massc2);
  double L = pow(Nparticles/density, 1./3.);
  std::cout << "mass*c^2: " << massc2 << "MeV" << std::endl;
  std::cout << "hbar^2/(2m) = " << prefactor << "MeV*fm^2" << std::endl;
  std::cout << "L = " << L << "fm V= " << L*L*L << "fm^3" << std::endl;
  std::cout << "k_fermi " << std::setprecision(16) << pow((6.*M_PI*M_PI*density/(2.*tzMax)),1./3.) << "fm^-1" << std::endl;

  double E;
  // vector<state> psi; // vector of sp states
  for(int i = 0; i < this->Nspstates; i++){
    E = ( prefactor*4*M_PI*M_PI/(L*L) ) * 
      (this->indexMap[i][0]*this->indexMap[i][0] + 
       this->indexMap[i][1]*this->indexMap[i][1] + 
       this->indexMap[i][2]*this->indexMap[i][2]);       
    this->spEnergy[i] = E;
  }  
} // end generateBasis


int infMatterSPBasis::checkSympqrs(int p, int q, int r, int s){
  int result = 0;
  //int ** indexMap = SPbasis.indexMap;
  // check spin and isospin
  if( (this->indexMap[p][3] + this->indexMap[q][3] == 
       this->indexMap[r][3] + this->indexMap[s][3]) && 
      (this->indexMap[p][4] + this->indexMap[q][4] == 
       this->indexMap[r][4] + this->indexMap[s][4]) ){
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

int infMatterSPBasis::checkModSympqrs(int p, int q, int r, int s){
  int result = 0;
  // check spin and isospin
  if( (this->indexMap[p][3] - this->indexMap[q][3] == 
       this->indexMap[r][3] - this->indexMap[s][3]) && 
      (this->indexMap[p][4] - this->indexMap[q][4] == 
       this->indexMap[r][4] - this->indexMap[s][4]) ){
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

int infMatterSPBasis::checkChanSym(int p, int q, int ichan){
  int result = 0;

  if( (this->indexMap[p][0] + this->indexMap[q][0] == 
       this->chanValue[ichan].chanNx) &&
      (this->indexMap[p][1] + this->indexMap[q][1] == 
       this->chanValue[ichan].chanNy) &&
      (this->indexMap[p][2] + this->indexMap[q][2] == 
       this->chanValue[ichan].chanNz) &&
      (this->indexMap[p][3] + this->indexMap[q][3] == 
       this->chanValue[ichan].chanSz) &&
      (this->indexMap[p][4] + this->indexMap[q][4] == 
       this->chanValue[ichan].chanTz) ){
    result = 1;
  }
  return result;
} // end checkChanSym

int infMatterSPBasis::checkChanSym2(int p, int q, int chanNx, int chanNy, int chanNz, int chanSz, int chanTz){
  int result = 0;
 
  if( (this->indexMap[p][0] + this->indexMap[q][0] == chanNx) &&
      (this->indexMap[p][1] + this->indexMap[q][1] == chanNy) &&
      (this->indexMap[p][2] + this->indexMap[q][2] == chanNz) &&
      (this->indexMap[p][3] + this->indexMap[q][3] == chanSz) &&
      (this->indexMap[p][4] + this->indexMap[q][4] == chanTz) ){
   
    result = 1;
  }
  return result;
} // end checkChanSym

int infMatterSPBasis::checkChanModSym2(int p, int q, int chanNx, int chanNy, int chanNz, int chanSz, int chanTz){
  int result = 0;
   
  if( (this->indexMap[p][0] - this->indexMap[q][0] == chanNx) &&
      (this->indexMap[p][1] - this->indexMap[q][1] == chanNy) &&
      (this->indexMap[p][2] - this->indexMap[q][2] == chanNz) &&
      (this->indexMap[p][3] - this->indexMap[q][3] == chanSz) &&
      (this->indexMap[p][4] - this->indexMap[q][4] == chanTz) ){
   
    result = 1;
  }
  return result;
} // end checkChanSym


int infMatterSPBasis::checkChanModSym(int p, int q, int ichan){
  int result = 0;

  if( (this->indexMap[p][0] - this->indexMap[q][0] == 
       this->chanModValue[ichan].chanNx) &&
      (this->indexMap[p][1] - this->indexMap[q][1] == 
       this->chanModValue[ichan].chanNy) &&
      (this->indexMap[p][2] - this->indexMap[q][2] == 
       this->chanModValue[ichan].chanNz) &&
      (this->indexMap[p][3] - this->indexMap[q][3] == 
       this->chanModValue[ichan].chanSz) &&
      (this->indexMap[p][4] - this->indexMap[q][4] == 
       this->chanModValue[ichan].chanTz) ){
    result = 1;
  }
  return result;
} // end checkChanModSym



void infMatterSPBasis::setUpTwoStateChannels(){
  int channelNmax = 2*this->nMax;  
  int channelTzMax;
  int modChannelTzMax;
  int channelSzMax = 2;
    
  if(this->tzMax == 1){
    channelTzMax = -2; 
    modChannelTzMax = 0;
  } else if(this->tzMax == 2){
    channelTzMax = 2;    
    modChannelTzMax = 2;
  } else {
    channelTzMax = 8000;   
    modChannelTzMax = 8000;
    std::cout << "wrong tzMax" << std::endl;
  }
  
  int channelCount = 0;
  for(int ichanNx = -channelNmax; ichanNx <= channelNmax; ichanNx++){
    for(int ichanNy = -channelNmax; ichanNy <= channelNmax; ichanNy++){
      for(int ichanNz = -channelNmax; ichanNz <= channelNmax; ichanNz++){
	for(int ichanSz = -channelSzMax; ichanSz <= channelSzMax; ichanSz = ichanSz +2){
	  for(int ichanTz = -2; ichanTz <= channelTzMax; ichanTz = ichanTz +2){
	     
   
	    for(int p = 0; p < this->Nspstates; p++){	      
	      for(int q = 0; q < this->Nspstates; q++){
		if( this->checkChanSym2(p,q,ichanNx,ichanNy,ichanNz,ichanSz,ichanTz) == 1){
		  channelCount++;
		  goto loopbreak1;
		} // end if
	      } // end q
	    } // end p
	  loopbreak1:;
	   
	  } // end ichantz
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
	  for(int ichanTz = -2; ichanTz <= channelTzMax; ichanTz = ichanTz +2){
 	     
	    for(int p = 0; p < this->Nspstates; p++){
	      for(int q = 0; q < this->Nspstates; q++){
		if( this->checkChanSym2(p,q,ichanNx,ichanNy,ichanNz,ichanSz,ichanTz) == 1){
		  this->chanValue[channelCount].chanNx = ichanNx;
		  this->chanValue[channelCount].chanNy = ichanNy;
		  this->chanValue[channelCount].chanNz = ichanNz;
		  this->chanValue[channelCount].chanSz = ichanSz;
		  this->chanValue[channelCount].chanTz = ichanTz;
		  channelCount++; 
		  goto loopbreak2;
		} // end if
	      } // end q
	    } // end p
	  loopbreak2:;  
	    
	  } // end ichantz
	} // end ichansz
      } // end ichannz
    } // end ichanny
  } // end ichannx
   
      
  channelCount = 0;
  for(int ichanNx = -channelNmax; ichanNx <= channelNmax; ichanNx++){
    for(int ichanNy = -channelNmax; ichanNy <= channelNmax; ichanNy++){
      for(int ichanNz = -channelNmax; ichanNz <= channelNmax; ichanNz++){
	for(int ichanSz = -channelSzMax; ichanSz <= channelSzMax; ichanSz = ichanSz +2){
	  for(int imodChanTz = -modChannelTzMax; imodChanTz <= modChannelTzMax; imodChanTz = imodChanTz +2){	    
	     
	    for(int p = 0; p < this->Nspstates; p++){
	      for(int q = 0; q < this->Nspstates; q++){
		if( this->checkChanModSym2(p,q,ichanNx,ichanNy,ichanNz,ichanSz,imodChanTz) == 1){
		  this->chanModValue[channelCount].chanNx = ichanNx;
		  this->chanModValue[channelCount].chanNy = ichanNy;
		  this->chanModValue[channelCount].chanNz = ichanNz;
		  this->chanModValue[channelCount].chanSz = ichanSz;
		  this->chanModValue[channelCount].chanTz = imodChanTz;
		  channelCount++;
		  goto loopbreak3;
		} // end if
	      } // end q
	    } // end p
	  loopbreak3:;
	    
	  } // end ichantz
	} // end ichansz
      } // end ichannz
    } // end ichanny
  } // end ichannx
} // end calcNumChannels

void infMatterSPBasis::printBasis(){
  
  std::cout << "SPBasis:" << std::endl;
  std::cout << "i nx ny nz sz tz E" << std::endl;  
  for( int p = 0; p < this->Nspstates; p++ ) {    
    std::cout << p << " " << this->indexMap[p][0] << " " << this->indexMap[p][1] << " " << this->indexMap[p][2] << " " << this->indexMap[p][3] << " " << this->indexMap[p][4] << " " << this->spEnergy[p] << std::endl;    
  }
} // end printBasis

void infMatterSPBasis::deallocate(){

  for(int i = 0; i < this->Nspstates; i++){
    delete [] this->indexMap[i];
  }
  delete [] this->indexMap;
  delete [] this->spEnergy;
  delete [] this->chanValue;   
  delete [] this->chanModValue;
} // end deallocate



double infMatterSPBasis::calc_TBME(int p, int q, int r, int s){
   double vout = 0.;
   int *qi = this->indexMap[p];  
   int *qj = this->indexMap[q];  
   int *qk = this->indexMap[r];
   int *ql = this->indexMap[s];

   int initialMomx, initialMomy, initialMomz;
   int finalMomx, finalMomy, finalMomz;
   // Momentum Conservation Checks.
  initialMomx = qi[0] + qj[0];
  initialMomy = qi[1] + qj[1];
  initialMomz = qi[2] + qj[2];
  finalMomx = qk[0] + ql[0];
  finalMomy = qk[1] + ql[1];
  finalMomz = qk[2] + ql[2];
  
  // maybe add a conservation of spin if here.
  if( initialMomx == finalMomx && initialMomy == finalMomy && initialMomz == finalMomz ){
    double L = pow(this->Nparticles/this->density,1./3.);
    double V_R, V_T, V_S;
    double V_0R, V_0T, V_0S;
    double kappa_R, kappa_T, kappa_S;
    double *ki = new double[3];  
    double *kj = new double[3];  
    double *kk = new double[3];
    double *kl = new double[3];
    double *relMomBra = new double[3];
    double *relMomKet = new double[3];
    double *relMomTransf = new double[3];
    double qSquared, spinEx, isoSpinEx;
    double IsIt, PsIt, PsPt, IsPt;
    V_0R = 200; //MeV
    V_0T = 178; //MeV
    V_0S = 91.85; //MeV
    kappa_R = 1.487; //fm^-2
    kappa_T = 0.639; //fm^-2
    kappa_S = 0.465; //fm^-2
    
    qSquared = 0.;
    for( int i = 0; i < 3; i++){
      ki[i] = 2*M_PI*qi[i]/L;
      kj[i] = 2*M_PI*qj[i]/L;
      kk[i] = 2*M_PI*qk[i]/L;
      kl[i] = 2*M_PI*ql[i]/L;
      relMomBra[i] = 0.5*(ki[i] - kj[i]);
      relMomKet[i] = 0.5*(kk[i] - kl[i]);
      relMomTransf[i] = relMomBra[i] - relMomKet[i];
      qSquared += relMomTransf[i]*relMomTransf[i];
    }

    V_R = V_0R/(L*L*L)*pow(M_PI/kappa_R,1.5)*exp(-qSquared/(4*kappa_R));
    V_T = -V_0T/(L*L*L)*pow(M_PI/kappa_T,1.5)*exp(-qSquared/(4*kappa_T));
    V_S = -V_0S/(L*L*L)*pow(M_PI/kappa_S,1.5)*exp(-qSquared/(4*kappa_S));
    
    spinEx = spinExchangeMtxEle(qi[3],qj[3],qk[3],ql[3]);
    isoSpinEx = spinExchangeMtxEle(qi[4],qj[4],qk[4],ql[4]);
    
    // 4 terms, IsIt, PsIt, PsPt, IsPt
    // identity spin, identity isospin
    IsIt = kron_del(qi[3],qk[3])*kron_del(qj[3],ql[3])*kron_del(qi[4],qk[4])*kron_del(qj[4],ql[4]);
    // Exchange spin, identity isospin
    PsIt = spinEx*kron_del(qi[4],qk[4])*kron_del(qj[4],ql[4]);
    // Exchange spin, Exchange isospin
    PsPt = spinEx*isoSpinEx;
    // identity spin, Exchange isospin
    IsPt = kron_del(qi[3],qk[3])*kron_del(qj[3],ql[3])*isoSpinEx;
    
    vout = 0.5*(V_R + 0.5*V_T + 0.5*V_S)*IsIt
      + 0.25*(V_T - V_S)*PsIt
      - 0.5*(V_R + 0.5*V_T + 0.5*V_S)*PsPt
      - 0.25*(V_T - V_S)*IsPt;
    
    delete[] ki;  
    delete[] kj;
    delete[] kk;
    delete[] kl;
    delete[] relMomBra;
    delete[] relMomKet;
    delete[] relMomTransf;
  } // end Momentum Conversation Check

  return vout;
}


int infMatterSPBasis::spinExchangeMtxEle(int i, int j, int k, int l){
  if( i == l && j == k ){
    return 1;
  } else {
    return 0;
  }
} // end spinEx



int infMatterSPBasis::kron_del(int i, int j){
  if(i != j){
    return 0;
  }
  return 1;
} // end kron_del



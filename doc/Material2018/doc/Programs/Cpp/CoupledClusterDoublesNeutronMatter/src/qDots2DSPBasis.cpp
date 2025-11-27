// Need to write a PDF about this biz.

#include "qDots2DSPBasis.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>

qDots2DSPBasis::qDots2DSPBasis(double hbarOmega_In, int shellMax_In, int Nparticles_In){
  this->hbarOmega = hbarOmega_In;
  this->shellMax = shellMax_In;
  this->Nspstates = (shellMax_In+1)*(shellMax_In+2);
  this->Nparticles = Nparticles_In;
  // three 2-body spin states, 2*mlMax projection states
  this->Nchannels = 3*(4*shellMax_In+1);
}


void qDots2DSPBasis::generateIndexMap(){

  this->indexMap = new int* [Nspstates];
  for(int i = 0; i < this->Nspstates; i++){
    this->indexMap[i] = new int[3];
  } // maybe delete this at some point.
 
  
  // for(int i=0; i<numSpstates; i++){
  //   for(int j=0; j<3; j++){
  //     indexMap[i][j] = 1;
  //   } // end j loop
  // } // end i loop

  int n,nMax,count,nDecreaseFlag;  
  // mlMax = this->shellMax;  
  nDecreaseFlag = 0;
  count = 0; // count starts at zero for array purposes
  for(int shell = 0; shell<=this->shellMax; shell++){    
    nMax = shell/2;
    n=0;
    nDecreaseFlag = 0;
    for(int ml=-shell; ml<=shell; ml+=2){       
      indexMap[count][0] = n;
      indexMap[count][1] = ml;
      indexMap[count][2] = +1;      
      count++;      
      indexMap[count][0] = n;
      indexMap[count][1] = ml;
      indexMap[count][2] = -1;      
      count++;

      // increase n until nMax, then decrease n. Sometimes theres only pair of states an n=nMax, sometimes there's two.
      if( nDecreaseFlag == 0 && n < nMax ){
	n++;
      }
      else if( nDecreaseFlag == 0 && n == nMax ){
	nDecreaseFlag++;
	if( shell % 2 == 0 ) {
	  n--;
	}
      }
      else if( nDecreaseFlag == 1 ){
	n--;
      }
      else {
	std::cout << "nDecreaseFlag broken" << std::endl;
      } // end if chain      
    } // end ml loop
  } // end shell loop  
} // end generateIndexMap
 
   
// requires generate index map first
void qDots2DSPBasis::generateBasis(){
  this->spEnergy = new double[this->Nspstates];  
  std::cout << "hbar*Omega: " << this->hbarOmega << "MeV?" << std::endl; 

  double E;  
  // vector<state> psi; // vector of sp states
  for(int i = 0; i < this->Nspstates; i++){
    E = (2*this->indexMap[i][0] + abs(this->indexMap[i][1]) + 1)*this->hbarOmega;    
    this->spEnergy[i] = E;
  }  
} // end generateBasis

// Double check this
// Is adding spins too large of a channel?
// Are each particle's spin conserved individually?
int qDots2DSPBasis::checkSympqrs(int p, int q, int r, int s){
  int result = 0;
  //int ** indexMap = SPbasis.indexMap;
  // check spin 
  if( this->indexMap[p][2] + this->indexMap[q][2] == 
       this->indexMap[r][2] + this->indexMap[s][2] ){
    // check ml conservation
    if( this->indexMap[p][1] + this->indexMap[q][1] == 
	 this->indexMap[r][1] + this->indexMap[s][1] ) {
      result = 1;
    }    
  }
  return result;
} // end checkSympqrs

int qDots2DSPBasis::checkModSympqrs(int p, int q, int r, int s){
  int result = 0;
   // check spin 
  if( this->indexMap[p][2] - this->indexMap[q][2] == 
       this->indexMap[r][2] - this->indexMap[s][2] ){
    // check ml conservation
    if( this->indexMap[p][1] - this->indexMap[q][1] == 
	 this->indexMap[r][1] - this->indexMap[s][1] ) {
      result = 1;
    }    
  }
  return result;
} // end checkSympqrs

int qDots2DSPBasis::checkChanSym(int p, int q, int ichan){
  int result = 0;

  if( (this->indexMap[p][1] + this->indexMap[q][1] == 
       this->chanValue[ichan].chanMl) &&
      (this->indexMap[p][2] + this->indexMap[q][2] == 
       this->chanValue[ichan].chanSz) ){
    result = 1;
  }
  return result;
} // end checkChanSym

int qDots2DSPBasis::checkChanSym2(int p, int q, int chanMl, int chanSz){
  int result = 0;
 
  if( (this->indexMap[p][1] + this->indexMap[q][1] == chanMl) &&
      (this->indexMap[p][2] + this->indexMap[q][2] == chanSz) ){
   
    result = 1;
  }
  return result;
} // end checkChanSym

int qDots2DSPBasis::checkChanModSym2(int p, int q, int chanMl, int chanSz){
  int result = 0;
   
  if( (this->indexMap[p][1] - this->indexMap[q][1] == chanMl) &&
      (this->indexMap[p][2] - this->indexMap[q][2] == chanSz) ){
   
    result = 1;
  }
  return result;
} // end checkChanSym


int qDots2DSPBasis::checkChanModSym(int p, int q, int ichan){
  int result = 0;

  if( (this->indexMap[p][1] - this->indexMap[q][1] == 
       this->chanModValue[ichan].chanMl) &&
      (this->indexMap[p][2] - this->indexMap[q][2] == 
       this->chanModValue[ichan].chanSz) ){
    result = 1;
  }
  return result;
} // end checkChanModSym



void qDots2DSPBasis::setUpTwoStateChannels(){
  int channelMlMax = 2*this->shellMax;  
  int channelSzMax = 2;    
  
  this->chanValue = new channelBundle[this->Nchannels];
  this->chanModValue = new channelBundle[this->Nchannels];

  int channelCount = 0;
  for(int ichanMl = -channelMlMax; ichanMl <= channelMlMax; ichanMl++){   
    for(int ichanSz = -channelSzMax; ichanSz <= channelSzMax; ichanSz = ichanSz +2){  	   
      //std::cout << "channelCount = " << channelCount << "out of: " << this->Nchannels << std::endl;
      this->chanValue[channelCount].chanMl = ichanMl;
      this->chanValue[channelCount].chanSz = ichanSz;
      this->chanModValue[channelCount].chanMl = ichanMl;
      this->chanModValue[channelCount].chanSz = ichanSz;
      channelCount++;    
      
    } // end ichansz     
  } // end ichannml
  
} // end setUpTwoStateChannels

void qDots2DSPBasis::printBasis(){
  
  std::cout << "SPBasis:" << std::endl;
  std::cout << "i n ml sz E" << std::endl;  
  for( int p = 0; p < this->Nspstates; p++ ) {    
    std::cout << p << " " << this->indexMap[p][0] << " " << this->indexMap[p][1] << " " << this->indexMap[p][2] << " " << this->spEnergy[p] << std::endl;    
  }
} // end printBasis

void qDots2DSPBasis::deallocate(){

  for(int i = 0; i < this->Nspstates; i++){
    delete [] this->indexMap[i];
  }
  delete [] this->indexMap;
  delete [] this->spEnergy;
  delete [] this->chanValue;   
  delete [] this->chanModValue;
} // end deallocate



double qDots2DSPBasis::calc_TBME(int p, int q, int r, int s){

  // Code was given to me by Nathan Parzuchowski, I converted from .f90

  // Anisimovas, Matulis. J. Pys.: Condens. Matter 10, 601 (1998)
 
  // For some reason the k and l are switched. I don't know why.  
  // so v_int( i , j , l, k  ) gives V_{ijkl} (Not ANTI-SYMMETRIZED) 
  
  //  qi is a three integer array:  ( n_i , ml_i , ms_i ) 
  
   int *qi = this->indexMap[p];  
   int *qj = this->indexMap[q];
   // k and l flipped 
   int *qk = this->indexMap[s];
   int *ql = this->indexMap[r];



  int gs [4];
  //int qi[3],qj[3],qk[3],ql[3]
  double G,LAM, sm, sm_int, vout;
     
  vout = 0.;
  sm = 0.;
  
  // check symmetries
  // kron_dels can be removed if accessed properly
  if( kron_del(qi[1]+qj[1],qk[1]+ql[1]) * kron_del(qi[2], ql[2]) * kron_del(qk[2], qj[2]) == 1){  
    for(int j1 = 0; j1 <= qi[0]; j1++){
      for(int j2 = 0; j2 <= qj[0]; j2++){
	for(int j3 = 0; j3 <= qk[0]; j3++){
	  for(int j4 = 0; j4 <= ql[0]; j4++){

	    gs[0]=j1+j4+(abs(qi[1])+qi[1])/2. + (abs(ql[1])-ql[1])/2.;
            gs[3]=j1+j4+(abs(qi[1])-qi[1])/2. + (abs(ql[1])+ql[1])/2.;
            gs[1]=j2+j3+(abs(qj[1])+qj[1])/2. + (abs(qk[1])-qk[1])/2.;
            gs[2]=j2+j3+(abs(qj[1])-qj[1])/2. + (abs(qk[1])+qk[1])/2.;
	    
            G = gs[0] + gs[1] + gs[2] + gs[3];
	    sm_int = 0.;
	    for(int l1 = 0; l1 <= gs[0]; l1++){
	      for(int l2 = 0; l2 <= gs[1]; l2++){
		for(int l3 = 0; l3 <= gs[2]; l3++){
		  for(int l4 = 0; l4 <= gs[3]; l4++){
		    LAM = l1+l2+l3+l4;
		    
		    sm_int = sm_int + kron_del(l1+l2,l3+l4) * pow(-1, gs[1]+gs[2]-l2-l3) * bin_coef(gs[0],l1) * bin_coef(gs[1],l2) * bin_coef(gs[2],l3) * bin_coef(gs[3],l4) * specialGamma(1.+ LAM * 0.5) * specialGamma((G-LAM+1.) * 0.5);
		  } // end l4 loop
		} // end l3 loop
	      } // end l2 loop
	    } // end l1 loop
	    
	    
	    sm = sm + pow(-1, j1+j2+j3+j4)/ (factorial(j1)*factorial(j2)*factorial(j3)*factorial(j4)) * bin_coef(qi[0]+abs(qi[1]),qi[0]-j1) *  bin_coef(qj[0]+abs(qj[1]),qj[0]-j2) * bin_coef(qk[0]+abs(qk[1]),qk[0]-j3)  * bin_coef(ql[0]+abs(ql[1]),ql[0]-j4) * pow(0.5,(G+1.) * 0.5) * sm_int;
	    
	  } // end j4 loop
	} // end j3 loop
      } // end j2 loop
    } // end j1 loop
    
    vout = sqrt(this->hbarOmega) * sm * sqrt( fac_over_fac( qi[0], qi[0]+abs(qi[1]) ) * fac_over_fac( qj[0], qj[0]+abs(qj[1]) ) *  fac_over_fac( qk[0], qk[0]+abs(qk[1]) ) * fac_over_fac( ql[0], ql[0]+abs(ql[1]) ) );
    
    
  } // end if

  return vout;
  
} // end v_int


int qDots2DSPBasis::kron_del(int i, int j){
  if(i != j){
    return 0;
  }
  return 1;
} // end kron_del

double qDots2DSPBasis::specialGamma(double x){
  // only works for half or whole ints
  double sqpi, result; 
  int x_c;
  sqpi = sqrt(M_PI);
  x_c = floor(x);

  if( abs(float(x_c) - x) < 1e-3 ){
    result = factorial(x_c - 1);
    return result;
  } else {
    x_c = floor(x-0.5);
    result = fac_over_fac(2*x_c, x_c)/pow(4.0,x_c)*sqpi;
    return result;
  } // end if
} // end specialGamma

double qDots2DSPBasis::bin_coef(int n, int k){
  double result;
  result = fac_over_fac(n,k)/factorial(n-k);
  return result;
} // end bin_coef

double qDots2DSPBasis::fac_over_fac(int a1, int a2){
  double result = 1.;
  
  if(a1 > a2){
    for(int i = a2 + 1; i <= a1; i++){
      result = result*i;
    } // end i loop
    return result;
  } else {
    for(int i = a1 + 1; i <= a2; i++){
      result = result*i;
    } // end i loop
    result = 1./result;
    return result;    
  } // end if
} // end fac_over_fac

int qDots2DSPBasis::factorial(int n){
  // maybe change this to a long?
  int product = 1;
  
  for(int i=1; i<=n; i++){
    product = product*i;
  }
  return product;
} // end factorial


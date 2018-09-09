// Need to write a PDF about this biz.

#include "abstractSPbasis.hpp"
#include "pairingSPBasis.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

pairingSPBasis::pairingSPBasis(double xiIn, double gIn, int numPairStatesIn, int NparticlesIn){
  this->xi = xiIn;
  this->g = gIn;  
  this->Nspstates = 2*numPairStatesIn;
  this->Nparticles = NparticlesIn;  
  this->Nchannels = 3;
}


void pairingSPBasis::generateIndexMap(){

  this->indexMap = new int* [Nspstates];
  for(int i = 0; i < this->Nspstates; i++){
    this->indexMap[i] = new int[2];
  } 

  int count = 0;
  for( int i=0; i<Nspstates/2; i++){
    this->indexMap[count][0] = i+1;
    this->indexMap[count][1] = +1;
    count++;
    this->indexMap[count][0] = i+1;
    this->indexMap[count][1] = -1;
    count++;    
  }   

} // end generateIndexMap
 
   
// requires generate index map first
void pairingSPBasis::generateBasis(){
  
  this->spEnergy = new double[this->Nspstates];

  double E;
  // vector<state> psi; // vector of sp states
  for(int i = 0; i < this->Nspstates; i++){
    E = (indexMap[i][0] - 1.0)*xi;     
    this->spEnergy[i] = E;
  }  
} // end generateBasis


int pairingSPBasis::checkSympqrs(int p, int q, int r, int s){
  int result = 0;
  if( (this->indexMap[p][1] + this->indexMap[q][1] == this->indexMap[r][1] + this->indexMap[s][1]) ){
    result = 1;
  }  
  return result;
} // end checkSympqrs

int pairingSPBasis::checkModSympqrs(int p, int q, int r, int s){
  int result = 0;
  if( (this->indexMap[p][1] - this->indexMap[q][1] == this->indexMap[r][1] - this->indexMap[s][1]) ){
    result = 1;
  }  
  return result;
} // end checkSympqrs

int pairingSPBasis::checkChanSym(int p, int q, int ichan){
  int result = 0;
  if( indexMap[p][1] + indexMap[q][1] == this->chanValue[ichan].chanSz ){
    result = 1;
  }
  return result;
} // end checkChanSym

int pairingSPBasis::checkChanSym2(int p, int q, int chanSz){
  int result = 0;
 
  if( indexMap[p][1] + indexMap[q][1] == chanSz ){   
    result = 1;
  }
  return result;
} // end checkChanSym

int pairingSPBasis::checkChanModSym2(int p, int q, int chanSz){
  int result = 0;
   
 if( indexMap[p][1] - indexMap[q][1] == chanSz ){   
    result = 1;
  }
  return result;
} // end checkChanSym


int pairingSPBasis::checkChanModSym(int p, int q, int ichan){
  int result = 0;
  if( indexMap[p][1] - indexMap[q][1] == this->chanValue[ichan].chanSz ){
    result = 1;
  }
  return result;
} // end checkChanModSym



void pairingSPBasis::setUpTwoStateChannels(){
  
  this->chanValue = new channelBundle[this->Nchannels];
  this->chanModValue = new channelBundle[this->Nchannels];

  int channelSzMax = 2;
  int channelCount = 0;
  for(int ichanSz = -channelSzMax; ichanSz <= channelSzMax; ichanSz = ichanSz +2){
    this->chanValue[channelCount].chanSz = ichanSz;
    this->chanModValue[channelCount].chanSz = ichanSz;
    channelCount++;
  }  
} // end calcNumChannels

void pairingSPBasis::printBasis(){
  std::cout << "# p sz E" << std::endl;
  for(int i = 0; i < this->Nspstates; i++){    
    std::cout << i << " " << this->indexMap[i][0] << " " << this->indexMap[i][1] << " " << this->spEnergy[i] << std::endl;
  }
}

void pairingSPBasis::deallocate(){

  for(int i = 0; i < this->Nspstates; i++){
    delete [] this->indexMap[i];
  }
  delete [] this->indexMap;
  delete [] this->spEnergy;
  delete [] this->chanValue;   
  delete [] this->chanModValue;
} // end deallocate

// for antisymmetrized, do -0.5*g
double pairingSPBasis::calc_TBME(int p, int q, int r, int s){
 double vout = 0.;

  if(kronecker_del(this->indexMap[p][0],this->indexMap[q][0]) * kronecker_del(this->indexMap[r][0],this->indexMap[s][0]) * kronecker_del(this->indexMap[p][1], -this->indexMap[q][1]) * kronecker_del(-this->indexMap[r][1], this->indexMap[s][1]) == 1){
      vout = -0.25*this->g;
      
      if( this->indexMap[p][1] == -this->indexMap[r][1] ){
	vout = -vout;
      }
  }
  return vout;
}

int pairingSPBasis::kronecker_del(int i, int j){
  if(i != j){
    return 0;
  }
  return 1;
} // end kronecker_del


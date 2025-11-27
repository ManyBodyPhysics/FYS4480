#ifndef PAIRING_HPP
#define PAIRING_HPP

#include "abstractSPbasis.hpp"

// struct channelBundle{
//   int chanNx;
//   int chanNy;
//   int chanNz;
//   int chanSz;
//   int chanTz; 
// };

class pairingSPBasis: public abstractSPbasis{
public:
  double xi;
  double g;
  // int Nspstates;
  // int Nparticles;
  // int Nchannels;
  // int ** indexMap;
  // double * spEnergy;
  // channelBundle * chanValue;
  // channelBundle * chanModValue;

  pairingSPBasis(double xiIn, double gIn, int numPairStatesIn, int NparticlesIn);
  void generateIndexMap();
  void generateBasis();
  int checkSympqrs(int p, int q, int r, int s);
  int checkModSympqrs(int p, int q, int r, int s);
  int checkChanSym(int p, int q, int ichan);
  int checkChanSym2(int p, int q, int chanSz);
int checkChanModSym2(int p, int q, int chanSz);
  int checkChanModSym(int p, int q, int ichan);  
  void setUpTwoStateChannels();
  void printBasis();
  void deallocate();

  double calc_TBME(int p, int q, int r, int s);
  int kronecker_del(int i, int j);
};

#endif /* PAIRING_HPP */

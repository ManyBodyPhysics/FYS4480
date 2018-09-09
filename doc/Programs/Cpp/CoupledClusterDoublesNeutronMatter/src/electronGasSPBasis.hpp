#ifndef ELECTRONGAS_HPP
#define ELECTRONGAS_HPP

#include "abstractSPbasis.hpp"

class electronGasSPBasis: public abstractSPbasis{
public:
  

  int nMax;
  int shellMax;
  double r_s;
  double L;
  double e_sq;
  // int Nspstates;
  // int Nparticles;
  // int Nchannels;
  // int ** indexMap;
  // double * spEnergy;
  // channelBundle * chanValue;
  // channelBundle * chanModValue;
  

  electronGasSPBasis(double r_sIn, int shellMaxIn, int NparticlesIn);
  void generateIndexMap();
  void generateBasis();
  int checkSympqrs(int p, int q, int r, int s);
  int checkModSympqrs(int p, int q, int r, int s);
  int checkChanSym(int p, int q, int ichan);
  int checkChanSym2(int p, int q, int chanNx, int chanNy, int chanNz, int chanSz);
  int checkChanModSym2(int p, int q, int chanNx, int chanNy, int chanNz, int chanSz);
  int checkChanModSym(int p, int q, int ichan);  
  void setUpTwoStateChannels();
  void printBasis();
  void deallocate();

  double calc_TBME(int p, int q, int r, int s);
};

#endif /* ELECTRONGAS_HPP */

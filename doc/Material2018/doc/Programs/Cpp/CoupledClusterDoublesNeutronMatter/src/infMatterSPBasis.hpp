#ifndef INFMATTER_HPP
#define INFMATTER_HPP

#include "abstractSPbasis.hpp"

// struct channelBundle{
//   int chanNx;
//   int chanNy;
//   int chanNz;
//   int chanSz;
//   int chanTz; 
// };

class infMatterSPBasis: public abstractSPbasis{
public:
   double density;
   int tzMax;
   int nMax;
   int shellMax;
  // int Nspstates;
  // int Nparticles;
  // int Nchannels;
  // int ** indexMap;
  // double * spEnergy;
  // channelBundle * chanValue;
  // channelBundle * chanModValue;

  infMatterSPBasis(double densityIn, int tzMaxIn, int shellMaxIn, int NparticlesIn);// : abstractSPbasis(densityIn, tzMaxIn, shellMaxIn, NparticlesIn) {}
  void generateIndexMap();
  void generateBasis();
  int checkSympqrs(int p, int q, int r, int s);
  int checkModSympqrs(int p, int q, int r, int s);
  int checkChanSym(int p, int q, int ichan);
  int checkChanSym2(int p, int q, int chanNx, int chanNy, int chanNz, int chanSz, int chanTz);
int checkChanModSym2(int p, int q, int chanNx, int chanNy, int chanNz, int chanSz, int chanTz);
  int checkChanModSym(int p, int q, int ichan);  
  void setUpTwoStateChannels();
  void printBasis();
  void deallocate();

  double calc_TBME(int p, int q, int r, int s);
  int spinExchangeMtxEle(int i, int j, int k, int l);
  int kron_del(int i, int j);

};

#endif /* INFMATTER_HPP */

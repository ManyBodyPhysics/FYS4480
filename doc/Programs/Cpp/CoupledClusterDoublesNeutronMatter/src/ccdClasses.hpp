#include "abstractSPbasis.hpp"
#include "SymBlock.hpp"

class IndexBundle{
public:
  int channel;
  int row;
  int col;
  IndexBundle ();
  IndexBundle (int,int,int);
};


class Amplitudes{
public:
  int numChannels;

  SymmetryBlock * pphh;
  SymmetryBlock * phhp;
  SymmetryBlock * p_phh;
  SymmetryBlock * hpp_h;

  IndexBundle **** MegaMap;

  void updateMods();
  void rotateSPChannels();

  void generateMegaMap(int numChannels, abstractSPbasis * SPbasis);
  void deallocateMegaMap(int numChannels);
  IndexBundle get_Mod_Bundle(int p, int q, int r, int s, int numChannels, SymmetryBlock* O_pqrs, abstractSPbasis * SPbasis);
  
};


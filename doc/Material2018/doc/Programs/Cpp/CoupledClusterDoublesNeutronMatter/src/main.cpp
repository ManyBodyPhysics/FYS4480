// CCD code for arbitrary basis
// Preliminary code version, Morten Hjorth-Jensen and Justin Lietz, 2016


#include <iostream>
#include <cmath>
#include <iomanip>
#include <cblas.h>
#include <omp.h>
#include "SymBlock.hpp"
#include "abstractSPbasis.hpp"
#include "infMatterSPBasis.hpp"
#include "minnesotaInfMatter.hpp"
#include "electronGasSPBasis.hpp"
#include "electronGasInteraction.hpp"
#include "pairingSPBasis.hpp"
#include "pairingInteraction.hpp"
#include "qDots2DSPBasis.hpp"
#include "ccdClasses.hpp"

using namespace std;


// Two Body operator using channels
struct TB_OpChannels{
  SymmetryBlock * hhhh;
  SymmetryBlock * hhhp;
  SymmetryBlock * hhpp;  
  SymmetryBlock * hpph;  
  SymmetryBlock * hppp;
  SymmetryBlock * pppp;  
  SymmetryBlock * phh_p;
  SymmetryBlock * h_hpp; 
  SymmetryBlock * hpph_hphp_mod;
  SymmetryBlock * hhpp_hpph_mod; 

  //For Xnn
  double ** pp;
  double ** hh;
};


// Pass vnn by ref to not create an image of the whole thing.
void loadChannelBlocks1(TB_OpChannels &vnn, abstractSPbasis * SPbasis, double memRequested);
void loadChannelBlocks2(TB_OpChannels &vnn, abstractSPbasis * SPbasis, double memRequested);
void setUpAmplitudes(Amplitudes &t2, Amplitudes &t2_temp, abstractSPbasis * SPbasis);
double MBPT2corr(TB_OpChannels &vnn, int numChannels, abstractSPbasis * SPbasis);
void CCDstep(TB_OpChannels &vnn, TB_OpChannels &Xnn, Amplitudes &t2, Amplitudes &t2_temp, int numChannels, abstractSPbasis * SPbasis);
double CCDcorr(TB_OpChannels &vnn, Amplitudes &t2, int numChannels);
double referenceEnergy(abstractSPbasis * SPbasis);
double calcVpqrs(int p, int q, int r, int s, abstractSPbasis * SPbasis);


///////////////////////////////////////////////////
// begin of main program
///////////////////////////////////////////////////

int cheatingBasisParameter = -1;

int main(int argc, char * argv[])
{  
  
  double hw,xi,g,memRequested,density,r_s;
  int tzMax,shellMax,Nparticles, numPairStates;
  
  abstractSPbasis * SPbasis;
  
  if( argc != 7){
    cout << "Please select : 0 xi g numPairStates Nparticles memRequested (Gb): values at command line" << endl;
    cout << "or : 1 density tzMax shellMax Nparticles memRequested (Gb): values at command line" << endl;
    cout << "or : 2 0 r_s shellMax Nparticles memRequested (Gb): values at command line" << endl;
     cout << "or : 3 0 hbarOmega shellMax Nparticles memRequested (Gb): values at command line" << endl;
    return 0;    
  } else if ( atoi(argv[1]) == 0 ){
    cheatingBasisParameter = 0;
    xi = atof(argv[2]);
    g = atof(argv[3]);
    numPairStates = atoi(argv[4]);
    Nparticles = atoi(argv[5]);
    memRequested = atof(argv[6]);    
    SPbasis = new pairingSPBasis (xi,g,numPairStates,Nparticles);
  } else if ( atoi(argv[1]) == 1 ){
    cheatingBasisParameter = 1;
    density = atof(argv[2]);
    tzMax = atoi(argv[3]);
    shellMax = atoi(argv[4]);
    Nparticles = atoi(argv[5]);
    memRequested = atof(argv[6]);
    SPbasis = new infMatterSPBasis (density,tzMax,shellMax,Nparticles);
  } else if ( atoi(argv[1]) == 2 ){
    cheatingBasisParameter = 2;
    r_s = atof(argv[3]);    
    shellMax = atoi(argv[4]);
    Nparticles = atoi(argv[5]);
    memRequested = atof(argv[6]);
    SPbasis = new electronGasSPBasis (r_s,shellMax,Nparticles);
   } else if ( atoi(argv[1]) == 3 ){
    cheatingBasisParameter = 3;
    hw = atof(argv[3]);    
    shellMax = atoi(argv[4]);
    Nparticles = atoi(argv[5]);
    memRequested = atof(argv[6]);
    SPbasis = new qDots2DSPBasis (hw,shellMax,Nparticles);
  } else {

    cout << "bases: 0 for pairing, 1 for infMatter, 2 for electronGas, 3 for 2d qdots." << endl;
    return 0;
  }
   
  
  
  //omp_set_num_threads(1);  
  double wallTimeStart, wallTimeEnd;
  wallTimeStart = omp_get_wtime();


  int fermiLevel;
  fermiLevel = Nparticles;

  
  //abstractSPbasis * SPbasis = new infMatterSPBasis (density,tzMax,shellMax,Nparticles);

  
  SPbasis->generateIndexMap();
  SPbasis->generateBasis();
  
  
  int Nspstates = SPbasis->Nspstates;
 
  //Check build
  SPbasis->printBasis();

  //return 0;
 
  double Eref;
  Eref = referenceEnergy(SPbasis);
  

  // rotate sp energies to HF energies
  double ei,ea;
  for(int i=0; i<fermiLevel; i++){
    ei = SPbasis->spEnergy[i];
    for(int j=0; j<fermiLevel; j++){
      ei += calcVpqrs(i,j,i,j,SPbasis);
    }
    SPbasis->spEnergy[i] = ei;
  } 
  for(int a=fermiLevel; a<Nspstates; a++){
    ea = SPbasis->spEnergy[a];
    for(int j=0; j<fermiLevel; j++){
      ea += calcVpqrs(a,j,a,j,SPbasis);
    }
    SPbasis->spEnergy[a] = ea;
  }

  //Check rotation 
  SPbasis->printBasis();


  cout << endl << "CHANNELS" << endl;
  SPbasis->setUpTwoStateChannels();  
  int numChannels = SPbasis->Nchannels;
  cout << "numChannels: " << numChannels << endl;
  //int channelNmax = 2*SPbasis->nMax;
  // int MaxNumChannels = (2*channelNmax+1)*(2*channelNmax+1)*(2*channelNmax+1)*3*(2*tzMax-1);
  // cout << "MaxNumChannels: " << MaxNumChannels << endl;

  wallTimeEnd = omp_get_wtime();
  cout << endl << "Setup time: " << wallTimeEnd - wallTimeStart << endl << endl;

  //return 0;
   

  // Move this initializing to a function probably
  Amplitudes t2,t2_temp;
  t2.numChannels = numChannels;
  t2_temp.numChannels = numChannels;

  TB_OpChannels vnn,Xnn;
  vnn.hhhh = new SymmetryBlock[numChannels];
  vnn.hhhp = new SymmetryBlock[numChannels];
  vnn.hhpp = new SymmetryBlock[numChannels];
  vnn.hhpp_hpph_mod = new SymmetryBlock[numChannels];
  vnn.phh_p = new SymmetryBlock[Nspstates-Nparticles];
  vnn.h_hpp = new SymmetryBlock[Nparticles];
  vnn.hpph = new SymmetryBlock[numChannels]; 
  vnn.hpph_hphp_mod = new SymmetryBlock[numChannels]; 
  vnn.hppp = new SymmetryBlock[numChannels];
  vnn.pppp = new SymmetryBlock[numChannels];
  t2.pphh = new SymmetryBlock[numChannels];
  t2_temp.pphh = new SymmetryBlock[numChannels];
  t2.phhp = new SymmetryBlock[numChannels]; 
  t2_temp.phhp = new SymmetryBlock[numChannels]; 
  t2.p_phh = new SymmetryBlock[Nspstates-Nparticles]; 
  t2_temp.p_phh = new SymmetryBlock[Nspstates-Nparticles];
  t2.hpp_h = new SymmetryBlock[Nparticles];
  t2_temp.hpp_h = new SymmetryBlock[Nparticles];
  

  loadChannelBlocks2(vnn,SPbasis,memRequested);
  setUpAmplitudes(t2,t2_temp,SPbasis);

  
 
  t2.generateMegaMap(numChannels, SPbasis);

      
  int hhDim,hpDim;
  int numParticleStates = Nspstates - Nparticles;  
  
  Xnn.pp = new double*[numParticleStates];
  Xnn.hh = new double*[Nparticles];
  Xnn.hhhh = new SymmetryBlock[numChannels];
  Xnn.hpph = new SymmetryBlock[numChannels];
  Xnn.hpph_hphp_mod = new SymmetryBlock[numChannels];

  for(int a = 0; a < numParticleStates; a++){
    Xnn.pp[a] = new double[numParticleStates];
  }
  for(int i = 0; i < Nparticles; i++){
    Xnn.hh[i] = new double[Nparticles];
  }
  for(int ichan = 0; ichan < numChannels; ichan++){
    hhDim = vnn.hhhh[ichan].getRowNum();    
    hpDim = vnn.hpph[ichan].getRowNum();

    Xnn.hhhh[ichan].allocate(hhDim,hhDim);
    Xnn.hhhh[ichan].zeros();
    Xnn.hpph[ichan].allocate(hpDim,hpDim);
    Xnn.hpph[ichan].zeros();  
    Xnn.hpph_hphp_mod[ichan].allocate(hpDim,hpDim);
    Xnn.hpph_hphp_mod[ichan].zeros();    
  }
 
  
  double corrCCD,corrMBPT2;
  corrCCD = 0.;
  corrMBPT2 = MBPT2corr(vnn, numChannels, SPbasis);

  cout << "MBPT2corr = " << setprecision(16) << corrMBPT2 << endl;
  cout << "Eref: " << Eref << endl;
  cout << "Empbt2 = " << (Eref + corrMBPT2) << endl;

  cout << "MBPT2corr/A = " << setprecision(16) << corrMBPT2/Nparticles << endl;
  cout << "Eref/A: " << Eref/Nparticles << endl;
  cout << "Empbt2/A = " << (Eref + corrMBPT2)/Nparticles << endl;

  double tolerance = 1.e-6;
  double energyDiff = 1.0;
  double corrEnergyPrev;
  double value;  
   
  corrEnergyPrev = corrCCD;
   
  int II=0;
  corrCCD = 100.;
  corrEnergyPrev = corrCCD;
  energyDiff = 100.0;

  wallTimeEnd = omp_get_wtime();
  cout << endl << "Pre-iteration time: " << wallTimeEnd - wallTimeStart << endl << endl;

  //return 0;

  while(abs(energyDiff) > tolerance){   
    CCDstep(vnn,Xnn,t2,t2_temp,numChannels,SPbasis);
    corrCCD = CCDcorr(vnn, t2, numChannels);

    energyDiff = corrCCD - corrEnergyPrev;
    corrEnergyPrev = corrCCD;
    for(int ichan = 0; ichan < numChannels; ichan++){
      for(int pp = 0; pp < t2.pphh[ichan].getRowNum(); pp++){
	for(int hh = 0; hh < t2.pphh[ichan].getColNum(); hh++){
	  value = t2.pphh[ichan].getElement(pp,hh);
	  t2_temp.pphh[ichan].setElement(pp,hh,value);
	}
      }
    }
    t2.updateMods();
    cout << "i=" << II << ", corrCCD = " << corrCCD << ", energy diff = " << energyDiff << endl;
    cout << "i=" << II << ", corrCCD/A = " << corrCCD/Nparticles << ", energy diff/A = " << energyDiff/Nparticles << endl;
     II++;       
  }

  cout << "Using " << omp_get_max_threads() << " threads maximum." << endl;
  cout << setprecision(16) << "E_CCD= " << (Eref + corrCCD) << endl << endl;  
  cout << setprecision(16) << "E_CCD/A= " << (Eref + corrCCD)/Nparticles << endl << endl;  
  wallTimeEnd = omp_get_wtime();
  cout << "Final time: " << wallTimeEnd - wallTimeStart << endl;
   
  // cout << "Gaute Benchmarks for SNM, density = 0.16, nmax=2, A=28" << endl;
  // cout << "E0/A = -23.33771, MBPT2/A = -26.34531, CCD/A = -28.26330" << endl;
  // cout << "Gaute Benchmarks for PNM, density = 0.08, nmax=2, A=14" << endl;
  // cout << "E0/A = 10.33372, MBPT2/A = 9.70395, CCD/A = 9.69991" << endl;
  // cout << "My Benchmarks for SNM, density = 0.16, nmax=1, A=28" << endl;
  // cout << "E0/A = -23.33771, MBPT2/A = -26.34531, CCD/A = -26.644310" << endl;
  // cout << "My Benchmarks for PNM, density = 0.08, nmax=1, A=14" << endl;
  // cout << "E0/A = 10.33372, MBPT2/A = 10.15807, CCD/A = 10.155078" << endl;


  t2.deallocateMegaMap(numChannels);
  
   // Free all that memory
   for(int ichan = 0; ichan < numChannels; ichan++){
     vnn.hhhh[ichan].deallocate();
     vnn.hhhp[ichan].deallocate();
     vnn.hhpp[ichan].deallocate();   
     vnn.hhpp_hpph_mod[ichan].deallocate(); 
     vnn.hpph[ichan].deallocate();  
     vnn.hpph_hphp_mod[ichan].deallocate();   
     vnn.hppp[ichan].deallocate();
     vnn.pppp[ichan].deallocate();
     t2.pphh[ichan].deallocate();
     t2_temp.pphh[ichan].deallocate();
     t2.phhp[ichan].deallocate();
     t2_temp.phhp[ichan].deallocate();
   }

    
   for(int p = 0; p < Nspstates-Nparticles; p++){    
     vnn.phh_p[p].deallocate();  
     t2.p_phh[p].deallocate();    
     t2_temp.p_phh[p].deallocate();
   }

  
   for(int i = 0; i < Nparticles; i++){
     vnn.h_hpp[i].deallocate();
     t2.hpp_h[i].deallocate();
     t2_temp.hpp_h[i].deallocate();
   }
   
   
   // Free all that memory
   for(int ichan = 0; ichan < numChannels; ichan++){       
     Xnn.hhhh[ichan].deallocate();
     Xnn.hpph[ichan].deallocate();    
     Xnn.hpph_hphp_mod[ichan].deallocate();     
   }
  
   for(int a = 0; a < Nspstates-Nparticles; a++){
     delete [] Xnn.pp[a];
   }
   
   for(int i = 0; i < Nparticles; i++){
     delete [] Xnn.hh[i];
   }
   
   delete [] vnn.hhhh;
   delete [] vnn.hhhp;
   delete [] vnn.hhpp;
   delete [] vnn.hhpp_hpph_mod;
   delete [] vnn.phh_p;
   delete [] vnn.hpph; 
   delete [] vnn.hpph_hphp_mod;  
   delete [] vnn.hppp;
   delete [] vnn.pppp;
   delete [] vnn.h_hpp;
   delete [] t2.pphh;
   delete [] t2_temp.pphh;
   delete [] t2.phhp;
   delete [] t2_temp.phhp;  
   delete [] t2.p_phh;  
   delete [] t2_temp.p_phh;
   delete [] t2.hpp_h;
   delete [] t2_temp.hpp_h;
   
   delete [] Xnn.pp;
   delete [] Xnn.hh;
   delete [] Xnn.hhhh;
   delete [] Xnn.hpph;
   delete [] Xnn.hpph_hphp_mod;
   //} // end run loop
  
   SPbasis->deallocate();   

   return 0;  
 }

 double MBPT2corr(TB_OpChannels &vnn, int numChannels, abstractSPbasis * SPbasis){
   double corr = 0.;  
   double energyDenom;
   int h1,h2,p1,p2;
   for(int ichan = 0; ichan < numChannels; ichan++){    
     for(int hhIndex = 0; hhIndex < vnn.hhpp[ichan].getRowNum(); hhIndex++){
       for(int ppIndex = 0; ppIndex < vnn.hhpp[ichan].getColNum(); ppIndex++){
	 h1 = vnn.hhpp[ichan].rowMap[hhIndex][0];
	 h2 = vnn.hhpp[ichan].rowMap[hhIndex][1];
	 p1 = vnn.hhpp[ichan].colMap[ppIndex][0];
	 p2 = vnn.hhpp[ichan].colMap[ppIndex][1];
	 energyDenom = SPbasis->spEnergy[h1] + SPbasis->spEnergy[h2] - SPbasis->spEnergy[p1] - SPbasis->spEnergy[p2];
	 corr += vnn.hhpp[ichan].getElement(hhIndex,ppIndex)*vnn.hhpp[ichan].getElement(hhIndex,ppIndex)/energyDenom;	
       }
     }
   }
   corr = 0.25*corr;
   return corr;
 } // end MBPT2corr


 // Going with the form of CCD equations in Gaute's Lecture Notes
void CCDstep(TB_OpChannels &vnn, TB_OpChannels &Xnn, Amplitudes &t2, Amplitudes &t2_temp, int numChannels, abstractSPbasis * SPbasis){  
   int Nparticles = SPbasis->Nparticles;
   int Nspstates = SPbasis->Nspstates;

   int i,j,a,b; 
   double vijab, value, energyDenom, mixing;  
   double ALPHA, BETA;
   int M,K,N,LDA,LDB,LDC;  
   
   int chanTemp0,rowTemp0,colTemp0;
   int chanTemp1,rowTemp1,colTemp1;
   int chanTemp2,rowTemp2,colTemp2;
   int chanTemp3,rowTemp3,colTemp3;
   int chanTemp4,rowTemp4,colTemp4;
   int chanTemp5,rowTemp5,colTemp5;
   int chanTemp6,rowTemp6,colTemp6;
   int chanTemp7,rowTemp7,colTemp7;     
   double startSplitTime = omp_get_wtime();

   double endSplitTime = omp_get_wtime();
   double splitTicksTaken = endSplitTime - startSplitTime;
   double splitDuration = splitTicksTaken / (double) CLOCKS_PER_SEC;
         

  endSplitTime = omp_get_wtime();
  splitTicksTaken = endSplitTime - startSplitTime;
  splitDuration = splitTicksTaken / (double) CLOCKS_PER_SEC;
  cout << "t2_hpp_h, t2.p_phh, t2.phhp calc: " << splitDuration << endl;

  startSplitTime = omp_get_wtime();    
  // Xnn.pp
  for(int b = 0; b < Nspstates-Nparticles; b++){
    for(int c = 0; c < Nspstates-Nparticles; c++){       
      value = 0.;      
      if( b == c){
	for(int DKL = 0; DKL < vnn.phh_p[c].getRowNum(); DKL++){
	  value += t2.p_phh[b].getElement(0,DKL)*vnn.phh_p[c].getElement(DKL,0);
	}
      }      
      Xnn.pp[b][c] = -0.5*value;      
    } // end c loop
  } // end b loop
 
  endSplitTime = omp_get_wtime();
  splitTicksTaken = endSplitTime - startSplitTime;
  splitDuration = splitTicksTaken / (double) CLOCKS_PER_SEC;
  cout << "Xnn.pp calc: " << splitDuration << endl;

  startSplitTime = omp_get_wtime();
  // load Xnn.hh
  for(int k = 0; k < Nparticles; k++){
    for(int j = 0; j < Nparticles; j++){
      value = 0.;
      if( k == j){
	for(int iBra = 0; iBra < vnn.h_hpp[k].getColNum(); iBra++){
	  value += vnn.h_hpp[k].getElement(0,iBra)*t2.hpp_h[j].getElement(iBra,0);
	}
      }      
      Xnn.hh[k][j] = 0.5*value;      
    } // end j loop
  } // end k loop
  
  endSplitTime = omp_get_wtime();
  splitTicksTaken = endSplitTime - startSplitTime;
  splitDuration = splitTicksTaken / (double) CLOCKS_PER_SEC;
  cout << "Xnn.hh calc: " << splitDuration << endl;

  startSplitTime = omp_get_wtime();   
  value = 0.;  

  // load Xnn.hhhh
  for(int ichan = 0; ichan < numChannels; ichan++){      
    ALPHA = 0.5;
    BETA = 1.0;
    M = vnn.hhpp[ichan].getRowNum();
    K = vnn.hhpp[ichan].getColNum();
    N = t2_temp.pphh[ichan].getColNum();
    LDA = K;
    LDB = N;
    LDC = N;  
   
    for(int IJ = 0; IJ < vnn.hhhh[ichan].getRowNum(); IJ++){
      for(int KL = 0; KL < vnn.hhhh[ichan].getColNum(); KL++){
	value = vnn.hhhh[ichan].getElement(IJ,KL);// + Xnn.hhhh[ichan].getElement(IJ,KL);
	Xnn.hhhh[ichan].setElement(IJ,KL,value);
      }
    }
    if( M != 0 && K != 0 && N !=0 ){  
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, ALPHA, vnn.hhpp[ichan].pMatElements, LDA, t2_temp.pphh[ichan].pMatElements, LDB, BETA, Xnn.hhhh[ichan].pMatElements, LDC);
    }    
  }  
 
  endSplitTime = omp_get_wtime();
  splitTicksTaken = endSplitTime - startSplitTime;
  splitDuration = splitTicksTaken / (double) CLOCKS_PER_SEC;
  cout << "Xnn.hhhh calc: " << splitDuration << endl;

  startSplitTime = omp_get_wtime();
  // load Xnn.hpph_hphp_mod
  for(int ichan = 0; ichan < numChannels; ichan++){      
    ALPHA = 0.5;
    BETA = 1.0;
    M = vnn.hhpp_hpph_mod[ichan].getRowNum();
    K = vnn.hhpp_hpph_mod[ichan].getColNum();
    N = t2.phhp[ichan].getColNum();
    LDA = K;
    LDB = N;
    LDC = N;  
   
    for(int KC = 0; KC < vnn.hpph_hphp_mod[ichan].getRowNum(); KC++){
      for(int JB = 0; JB < vnn.hpph_hphp_mod[ichan].getColNum(); JB++){
  	value = vnn.hpph_hphp_mod[ichan].getElement(KC,JB);// + Xnn.hhhh[ichan].getElement(IJ,KL);
  	Xnn.hpph_hphp_mod[ichan].setElement(KC,JB,value);
      }
    }
    if( M != 0 && K != 0 && N !=0 ){  
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, ALPHA, vnn.hhpp_hpph_mod[ichan].pMatElements, 
		  LDA, t2.phhp[ichan].pMatElements, LDB, BETA, Xnn.hpph_hphp_mod[ichan].pMatElements, LDC);
    }    
  }    

  endSplitTime = omp_get_wtime();
  splitTicksTaken = endSplitTime - startSplitTime;
  splitDuration = splitTicksTaken / (double) CLOCKS_PER_SEC;
  cout << "Xnn.hpph_hphp_mod calc: " << splitDuration << endl;

  startSplitTime = omp_get_wtime();

  //rotate t2.p_phh & t2.hpp_h
  t2.rotateSPChannels();

 

  endSplitTime = omp_get_wtime();
  splitTicksTaken = endSplitTime - startSplitTime;
  splitDuration = splitTicksTaken / (double) CLOCKS_PER_SEC;
  cout << "t2.p_phh & t2.hpp_h rotate: " << splitDuration << endl;

  startSplitTime = omp_get_wtime();  
  // KLsum
  for(int ichan = 0; ichan < numChannels; ichan++){      
    ALPHA = 0.5;
    BETA = 0.0;
    M = t2_temp.pphh[ichan].getRowNum();
    K = t2_temp.pphh[ichan].getColNum();
    N = Xnn.hhhh[ichan].getColNum();
    LDA = K;
    LDB = N;
    LDC = N;  
   
    if( M != 0 && K != 0 && N !=0 ){  
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, ALPHA, t2_temp.pphh[ichan].pMatElements, LDA, Xnn.hhhh[ichan].pMatElements, LDB, BETA, t2.pphh[ichan].pMatElements, LDC);
    }
  }

  endSplitTime = omp_get_wtime();
  splitTicksTaken = endSplitTime - startSplitTime;
  splitDuration = splitTicksTaken / (double) CLOCKS_PER_SEC;
  cout << "KLsum calc: " << splitDuration << endl;

  startSplitTime = omp_get_wtime();

  // KCsum
  for(int ichan = 0; ichan < numChannels; ichan++){      
    ALPHA = 1.0;
    BETA = 0.0;
    M = t2.phhp[ichan].getRowNum();
    K = t2.phhp[ichan].getColNum();
    N = Xnn.hpph_hphp_mod[ichan].getColNum();
    LDA = K;
    LDB = N;
    LDC = N;      

    if( M != 0 && K != 0 && N !=0 ){  
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, ALPHA, t2.phhp[ichan].pMatElements, LDA, Xnn.hpph_hphp_mod[ichan].pMatElements, LDB, BETA, t2_temp.phhp[ichan].pMatElements, LDC);
    }
  }  

  endSplitTime = omp_get_wtime();
  splitTicksTaken = endSplitTime - startSplitTime;
  splitDuration = splitTicksTaken / (double) CLOCKS_PER_SEC;
  cout << "KCsum calc: " << splitDuration << endl;

  startSplitTime = omp_get_wtime();

  // CD sum
  for(int ichan = 0; ichan < numChannels; ichan++){      
    ALPHA = 0.5;
    BETA = 1.0;
    M = vnn.pppp[ichan].getRowNum();
    K = vnn.pppp[ichan].getColNum();
    N = t2_temp.pphh[ichan].getColNum();
    LDA = K;
    LDB = N;
    LDC = N;      

    if( M != 0 && K != 0 && N !=0 ){  
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, ALPHA, vnn.pppp[ichan].pMatElements, LDA, t2_temp.pphh[ichan].pMatElements, LDB, BETA, t2.pphh[ichan].pMatElements, LDC);
    }
  }

  endSplitTime = omp_get_wtime();
  splitTicksTaken = endSplitTime - startSplitTime;
  splitDuration = splitTicksTaken / (double) CLOCKS_PER_SEC;
  cout << "CDsum calc: " << splitDuration << endl;

  startSplitTime = omp_get_wtime();
    
  // Csum     
  for(int B = 0; B < Nspstates-Nparticles; B++){    
    for(int AIJ = 0; AIJ < t2_temp.p_phh[B].getColNum(); AIJ++){
      value = Xnn.pp[B][B]*t2.p_phh[B].getElement(0,AIJ);      
      t2_temp.p_phh[B].setElement(0,AIJ,value);      
    }
  }

  endSplitTime = omp_get_wtime();
  splitTicksTaken = endSplitTime - startSplitTime;
  splitDuration = splitTicksTaken / (double) CLOCKS_PER_SEC;
  cout << "Csum calc: " << splitDuration << endl;

  startSplitTime = omp_get_wtime(); 
  
  // Ksum
  for(int J = 0; J < Nparticles; J++){    
    for(int IAB = 0; IAB < t2.hpp_h[J].getRowNum(); IAB++){      
      value = t2.hpp_h[J].getElement(IAB,0)*Xnn.hh[J][J];
      t2_temp.hpp_h[J].setElement(IAB,0,value);
    }
  }
  
  endSplitTime = omp_get_wtime();
  splitTicksTaken = endSplitTime - startSplitTime;
  splitDuration = splitTicksTaken / (double) CLOCKS_PER_SEC;
  cout << "Ksum calc: " << splitDuration << endl;

  startSplitTime = omp_get_wtime();  

  for(int ichan = 0; ichan < numChannels; ichan++){
    for(int AB = 0; AB < vnn.hhpp[ichan].getColNum(); AB++){
      for(int IJ = 0; IJ < vnn.hhpp[ichan].getRowNum(); IJ++){
      	
	a = vnn.hhpp[ichan].colMap[AB][0];
  	b = vnn.hhpp[ichan].colMap[AB][1];
	i = vnn.hhpp[ichan].rowMap[IJ][0];
  	j = vnn.hhpp[ichan].rowMap[IJ][1];
  	
  	energyDenom = SPbasis->spEnergy[i] + SPbasis->spEnergy[j] - SPbasis->spEnergy[a] - SPbasis->spEnergy[b];	
  	vijab = vnn.hhpp[ichan].getElement(IJ,AB);

	chanTemp0 = t2.MegaMap[0][ichan][AB][IJ].channel;
	rowTemp0 = t2.MegaMap[0][ichan][AB][IJ].row;
	colTemp0 = t2.MegaMap[0][ichan][AB][IJ].col;

	chanTemp1 = t2.MegaMap[1][ichan][AB][IJ].channel;
	rowTemp1 = t2.MegaMap[1][ichan][AB][IJ].row;
	colTemp1 = t2.MegaMap[1][ichan][AB][IJ].col;

	chanTemp2 = t2.MegaMap[2][ichan][AB][IJ].channel;
	rowTemp2 = t2.MegaMap[2][ichan][AB][IJ].row;
	colTemp2 = t2.MegaMap[2][ichan][AB][IJ].col;

	chanTemp3 = t2.MegaMap[3][ichan][AB][IJ].channel;
	rowTemp3 = t2.MegaMap[3][ichan][AB][IJ].row;
	colTemp3 = t2.MegaMap[3][ichan][AB][IJ].col;

	chanTemp4 = t2.MegaMap[4][ichan][AB][IJ].channel;
	rowTemp4 = t2.MegaMap[4][ichan][AB][IJ].row;
	colTemp4 = t2.MegaMap[4][ichan][AB][IJ].col;

	chanTemp5 = t2.MegaMap[5][ichan][AB][IJ].channel;
	rowTemp5 = t2.MegaMap[5][ichan][AB][IJ].row;
	colTemp5 = t2.MegaMap[5][ichan][AB][IJ].col;

	chanTemp6 = t2.MegaMap[6][ichan][AB][IJ].channel;
	rowTemp6 = t2.MegaMap[6][ichan][AB][IJ].row;
	colTemp6 = t2.MegaMap[6][ichan][AB][IJ].col;

	chanTemp7 = t2.MegaMap[7][ichan][AB][IJ].channel;
	rowTemp7 = t2.MegaMap[7][ichan][AB][IJ].row;
	colTemp7 = t2.MegaMap[7][ichan][AB][IJ].col;

  	value = vijab + t2.pphh[ichan].getElement(AB,IJ) + t2_temp.phhp[chanTemp4].getElement(rowTemp4,colTemp4) 
	  - t2_temp.phhp[chanTemp5].getElement(rowTemp5,colTemp5) - t2_temp.phhp[chanTemp6].getElement(rowTemp6,colTemp6)
	  + t2_temp.phhp[chanTemp7].getElement(rowTemp7,colTemp7)
	  + t2_temp.p_phh[chanTemp3].getElement(rowTemp3,colTemp3) - t2_temp.p_phh[chanTemp2].getElement(rowTemp2,colTemp2)
	  + t2_temp.hpp_h[chanTemp1].getElement(rowTemp1,colTemp1) - t2_temp.hpp_h[chanTemp0].getElement(rowTemp0,colTemp0);

		
  	value = value/energyDenom;
	mixing = 1.0;
	value = mixing*value + (1.-mixing)*t2_temp.pphh[ichan].getElement(AB,IJ);
  	t2.pphh[ichan].setElement(AB,IJ,value);	
      } // end AB loop
    } // end IJ loop
  } // end ichan loop
  
  endSplitTime = omp_get_wtime();
  splitTicksTaken = endSplitTime - startSplitTime;
  splitDuration = splitTicksTaken / (double) CLOCKS_PER_SEC;
  cout << "t_pphh calc: " << splitDuration << endl;

} // end CCDstep


double CCDcorr(TB_OpChannels &vnn, Amplitudes &t2, int numChannels){
  double corr = 0.;   
  for(int ichan = 0; ichan < numChannels; ichan++){    
    for(int IJ = 0; IJ < vnn.hhpp[ichan].getRowNum(); IJ++){
      for(int AB = 0; AB < vnn.hhpp[ichan].getColNum(); AB++){
	corr += vnn.hhpp[ichan].getElement(IJ,AB)*t2.pphh[ichan].getElement(AB,IJ);
      }
    }
  }
  corr = 0.25*corr;
  return corr;
} // end CCDcorr


double referenceEnergy(abstractSPbasis * SPbasis){  
  int Nparticles = SPbasis->Nparticles;
  int fermiLevel = Nparticles;  
  double vijij = 0.;
  double Eref = 0.;
  
  for(int i = 0; i < fermiLevel-1; i++){
    for(int j = i+1; j < fermiLevel; j++){
      if( SPbasis->checkSympqrs(i,j,i,j) == 1 ){
  	vijij = calcVpqrs(i,j,i,j,SPbasis);
	Eref += vijij;
      } // end if  
    } // end j    
  } // end i
  
  for(int i = 0; i < fermiLevel; i++){
    Eref += SPbasis->spEnergy[i];
  }
  
  return Eref;
}


double calcVpqrs(int p, int q, int r, int s, abstractSPbasis * SPbasis){
  double vpqrs = 0.;
  vpqrs = SPbasis->calc_TBME(p,q,r,s) - SPbasis->calc_TBME(p,q,s,r);
  // if( SPbasis->checkSympqrs(p,q,r,s) == 1){
  //   if( cheatingBasisParameter == 0 ){
  //     //double g = SPbasis->g;
  //     //vpqrs = vint_pairing(SPbasis->indexMap[p],SPbasis->indexMap[q],SPbasis->indexMap[r],SPbasis->indexMap[s],g);
  //     vpqrs = SPbasis->calc_TBME(p,q,r,s);
  //   } else if( cheatingBasisParameter == 1 ){
  //     double L = pow(SPbasis->Nparticles/SPbasis->density,1./3.);  
  //     vpqrs = vint_Minnesota_Momentum(SPbasis->indexMap[p],SPbasis->indexMap[q],SPbasis->indexMap[r],SPbasis->indexMap[s],L) - vint_Minnesota_Momentum(SPbasis->indexMap[p],SPbasis->indexMap[q],SPbasis->indexMap[s],SPbasis->indexMap[r],L);
  //   } else if( cheatingBasisParameter == 2 ){
  //     double hbarc = 0.1973269788; // eV micron
  //     hbarc = hbarc*10000/27.21138505; // Hartree Angstrom
  //     double fine_struc = 1.0/137.035999139;
  //     double massc2 = 0.5109989461; // MeV wikipedia electron mass 
  //     massc2 = massc2*1000000/27.21138505; // Hartree

  //     double r_bohr = hbarc/(massc2*fine_struc);
  //     double L = pow( (4./3.)*M_PI*SPbasis->Nparticles, 1./3.)*SPbasis->r_s*r_bohr;
  //     vpqrs = vint_electronGas_Momentum(SPbasis->indexMap[p],SPbasis->indexMap[q],SPbasis->indexMap[r],SPbasis->indexMap[s],L) - vint_electronGas_Momentum(SPbasis->indexMap[p],SPbasis->indexMap[q],SPbasis->indexMap[s],SPbasis->indexMap[r],L);
  //   } else {
  //     cout << "global cheat parameter broken." << endl;
  //   }    
  // }
  
  return vpqrs;
} // end calcVpqrs

void setUpAmplitudes(Amplitudes &t2, Amplitudes &t2_temp, abstractSPbasis * SPbasis){
  
  int Nparticles = SPbasis->Nparticles;
  int Nspstates = SPbasis->Nspstates;    
  int fermiLevel = Nparticles;   
  
  
  ////////////////////////////////////////////////////////////////////////
  cout << "Setting up t2.pphh,t2_temp.pphh" << endl;
  ////////////////////////////////////////////////////////////////////////
  
  int memoryUsed = 0;     
#pragma omp parallel for reduction(+:memoryUsed) 
  for(int ichan = 0; ichan < SPbasis->Nchannels; ichan++){
    int hhDim,ppDim,hhIndex,ppIndex;
              
     hhDim = 0;	   	     	    
     for(int i = 0; i < fermiLevel; i++){
       for(int j = 0; j < fermiLevel; j++){	
	 if( SPbasis->checkChanSym(i,j,ichan) == 1 ){		 
	   hhDim++;		  
	 } // end if
       } // end j loop
     } // end i loop

     ppDim = 0;
     for(int a = fermiLevel; a < Nspstates; a++){
       for(int b = fermiLevel; b < Nspstates; b++){
	 if( SPbasis->checkChanSym(a,b,ichan) == 1 ){			 
	   ppDim++;
	 } // end if
       } // end b loop
     } // end a loop

     memoryUsed+= hhDim*ppDim;
     t2.pphh[ichan].allocate(ppDim,hhDim);
     t2.pphh[ichan].zeros();
     t2_temp.pphh[ichan].allocate(ppDim,hhDim);
     t2_temp.pphh[ichan].zeros();
	     	    
     ppIndex = 0;
     for(int a = fermiLevel; a < Nspstates; a++){
       for(int b = fermiLevel; b < Nspstates; b++){
	 if( SPbasis->checkChanSym(a,b,ichan) == 1 ){
	   hhIndex = 0;
	   for(int i = 0; i < fermiLevel; i++){
	     for(int j = 0; j < fermiLevel; j++){
	       if( SPbasis->checkChanSym(i,j,ichan) == 1 ){
		 t2.pphh[ichan].rowMap[ppIndex][0] = a;
		 t2.pphh[ichan].rowMap[ppIndex][1] = b;
		 t2.pphh[ichan].colMap[hhIndex][0] = i;
		 t2.pphh[ichan].colMap[hhIndex][1] = j;

		 t2_temp.pphh[ichan].rowMap[ppIndex][0] = a;
		 t2_temp.pphh[ichan].rowMap[ppIndex][1] = b;
		 t2_temp.pphh[ichan].colMap[hhIndex][0] = i;
		 t2_temp.pphh[ichan].colMap[hhIndex][1] = j;
		 hhIndex++;
	       } // end if		      
	     } // end b loop
	   } // end a loop
	   ppIndex++;
	 } // end if
       } // end a loop
     } // end i loop

   } // end ichan
     
   cout << "t2.pphh requires: " << memoryUsed*sizeof(double)/1.e9 << " Gb." << endl;
   cout << "t2_temp.pphh requires: " << memoryUsed*sizeof(double)/1.e9 << " Gb." << endl;
 
 ////////////////////////////////////////////////////////////////////////
   cout << "Setting up t2.p_phh, t2_temp.p_phh" << endl;
   ////////////////////////////////////////////////////////////////////////

   memoryUsed = 0;     
#pragma omp parallel for reduction(+:memoryUsed)   
   for(int a = fermiLevel; a < Nspstates; a++){
     int pDim = 1;
     int phhDim = 0;
     for(int b = fermiLevel; b < Nspstates; b++){
       for(int i = 0; i < fermiLevel; i++){
	 for(int j = 0; j < fermiLevel; j++){
	   if( SPbasis->checkSympqrs(a,b,i,j) == 1 ){
	     phhDim++;
	   }
	 } // end j loop
       } // end i loop
     } // end b loop
     memoryUsed += pDim*phhDim;
     t2.p_phh[a-Nparticles].allocateWithLabels(pDim,phhDim,1,3);
     t2.p_phh[a-Nparticles].zeros();     
     t2_temp.p_phh[a-Nparticles].allocateWithLabels(pDim,phhDim,1,3);
     t2_temp.p_phh[a-Nparticles].zeros();
     //cout << "a: " << a << " phhDim: " << phhDim << endl;
   } // end a loop
     
#pragma omp parallel for
   for(int a = fermiLevel; a < Nspstates; a++){
     int phhIndex = 0;
     for(int b = fermiLevel; b < Nspstates; b++){
       for(int i = 0; i < fermiLevel; i++){
	 for(int j = 0; j < fermiLevel; j++){
	   if( SPbasis->checkSympqrs(a,b,i,j) == 1 ){	   	   
	     t2.p_phh[a-Nparticles].rowMap[0][0] = a;
	     t2.p_phh[a-Nparticles].colMap[phhIndex][0] = b;
	     t2.p_phh[a-Nparticles].colMap[phhIndex][1] = i;
	     t2.p_phh[a-Nparticles].colMap[phhIndex][2] = j;

	     t2_temp.p_phh[a-Nparticles].rowMap[0][0] = a;
	     t2_temp.p_phh[a-Nparticles].colMap[phhIndex][0] = b;
	     t2_temp.p_phh[a-Nparticles].colMap[phhIndex][1] = i;
	     t2_temp.p_phh[a-Nparticles].colMap[phhIndex][2] = j;
	     phhIndex++;
	   }
	 } // end j loop
       } // end i loop
     } // end b loop   
   } // end a loop

   cout << "t2.p_phh requires: " << memoryUsed*sizeof(double)/1.e9 << " Gb." << endl;
   
   ////////////////////////////////////////////////////////////////////////
   cout << "Setting up t2.hpp_h, t2_temp.hpp_h" << endl;
   ////////////////////////////////////////////////////////////////////////

  memoryUsed = 0;     
#pragma omp parallel for reduction(+:memoryUsed)     
   for(int i = 0; i < fermiLevel; i++){
     int hDim = 1;
     int hppDim = 0;
     for(int j = 0; j < fermiLevel; j++){
       for(int a = fermiLevel; a < Nspstates; a++){
	 for(int b = fermiLevel; b < Nspstates; b++){
	   if( SPbasis->checkSympqrs(i,j,a,b) == 1 ){
	     hppDim++;
	   }
	 } // end b loop
       } // end a loop
     } // end j loop
     memoryUsed += hDim*hppDim;
     t2.hpp_h[i].allocateWithLabels(hppDim,hDim,3,1);
     t2.hpp_h[i].zeros();
     t2_temp.hpp_h[i].allocateWithLabels(hppDim,hDim,3,1);
     t2_temp.hpp_h[i].zeros();     
   } // end i loop
   
#pragma omp parallel for
   for(int i = 0; i < fermiLevel; i++){
     int hppIndex = 0;
     for(int j = 0; j < fermiLevel; j++){
       for(int a = fermiLevel; a < Nspstates; a++){
	 for(int b = fermiLevel; b < Nspstates; b++){
	   if( SPbasis->checkSympqrs(i,j,a,b) == 1 ){	   	   
	     t2.hpp_h[i].rowMap[hppIndex][0] = j;
	     t2.hpp_h[i].rowMap[hppIndex][1] = a;
	     t2.hpp_h[i].rowMap[hppIndex][2] = b;
	     t2.hpp_h[i].colMap[0][0] = i;

	     t2_temp.hpp_h[i].rowMap[hppIndex][0] = j;
	     t2_temp.hpp_h[i].rowMap[hppIndex][1] = a;
	     t2_temp.hpp_h[i].rowMap[hppIndex][2] = b;
	     t2_temp.hpp_h[i].colMap[0][0] = i;
	     hppIndex++;
	   }
	 } // end b loop
       } // end a loop
     } // end j loop   
   } // end i loop

   cout << "t2.hpp_h requires: " << memoryUsed*sizeof(double)/1.e9 << " Gb." << endl;

   ///////////////////////////////////////////////////////////////////////////////////////////////

 
    
   ////////////////////////////////////////////////////////////////////////
   cout << "Setting up t2.phhp,t2_temp.phhp" << endl;
   ////////////////////////////////////////////////////////////////////////
   
   memoryUsed = 0;   
#pragma omp parallel for reduction(+:memoryUsed) 
   for(int ichan = 0; ichan < SPbasis->Nchannels; ichan++){
     int hpModDim,phModDim,phIndex,hpIndex;
	 
    
     hpModDim = 0;
     for(int i = 0; i < fermiLevel; i++){
       for(int a = fermiLevel; a < Nspstates; a++){
	 if( SPbasis->checkChanModSym(i,a,ichan) == 1 ){	  
	   hpModDim++;		  
	 } // end if		
       } // end a loop
     } // end i loop

     phModDim = 0;
     for(int a = fermiLevel; a < Nspstates; a++){
       for(int i = 0; i < fermiLevel; i++){
	 if( SPbasis->checkChanModSym(a,i,ichan) == 1 ){	  
	   phModDim++;		  
	 } // end if		
       } // end a loop
     } // end i loop


     memoryUsed+= phModDim*hpModDim;
     t2.phhp[ichan].allocate(phModDim,hpModDim);
     t2.phhp[ichan].zeros();
     t2_temp.phhp[ichan].allocate(phModDim,hpModDim);
     t2_temp.phhp[ichan].zeros();
	     	    
     phIndex = 0;
     for(int a = fermiLevel; a < Nspstates; a++){
       for(int i = 0; i < fermiLevel; i++){
	 if( SPbasis->checkChanModSym(a,i,ichan) == 1 ){
	   hpIndex = 0;
	   for(int j = 0; j < fermiLevel; j++){
	     for(int b = fermiLevel; b < Nspstates; b++){
	       if( SPbasis->checkChanModSym(j,b,ichan) == 1 ){		
		 t2.phhp[ichan].rowMap[phIndex][0] = a;		
		 t2.phhp[ichan].rowMap[phIndex][1] = i;
		 t2.phhp[ichan].colMap[hpIndex][0] = j;
		 t2.phhp[ichan].colMap[hpIndex][1] = b;

		 t2_temp.phhp[ichan].rowMap[phIndex][0] = a;		
		 t2_temp.phhp[ichan].rowMap[phIndex][1] = i;
		 t2_temp.phhp[ichan].colMap[hpIndex][0] = j;
		 t2_temp.phhp[ichan].colMap[hpIndex][1] = b;
		 hpIndex++;			
	       } // end if		      
	     } // end a loop
	   } // end j loop		  	  
	   phIndex++;
	 } // end if		
       } // end b loop
     } // end i loop
	  
   } // end ichan
   
   cout << "t2.phhp requires: " << memoryUsed*sizeof(double)/1.e9 << " Gb." << endl;
   cout << "t2_temp.phhp requires: " << memoryUsed*sizeof(double)/1.e9 << " Gb." << endl;      
}

void loadChannelBlocks1(TB_OpChannels &vnn, abstractSPbasis * SPbasis, double memRequested){
           
   int Nparticles = SPbasis->Nparticles;
   int Nspstates = SPbasis->Nspstates;
      
   int memoryUsed,totalMemoryUsed,ppppMemory,ppppIndexMemory;
   int fermiLevel = Nparticles;   
   int numChannels = SPbasis->Nchannels;
      
   totalMemoryUsed = 0;
   ppppIndexMemory = 0;
   ////////////////////////////////////////////////////////////////////////
   cout << "Setting up vnn.pppp" << endl;
   ////////////////////////////////////////////////////////////////////////
   
   
   memoryUsed = 0;
   totalMemoryUsed = 0;
#pragma omp parallel for reduction(+:memoryUsed,ppppIndexMemory) 
   for(int ichan = 0; ichan < numChannels; ichan++){
     double vpqrs;
     int ppDim,ppIndex,pp2Index;


     ppDim = 0;
     for(int a = fermiLevel; a < Nspstates; a++){
       for(int b = fermiLevel; b < Nspstates; b++){
	 if( SPbasis->checkChanSym(a,b,ichan) == 1 ){			 
	   ppDim++;
	 } // end if
		 
       } // end b loop
     } // end a loop

     memoryUsed += ppDim*ppDim;
     ppppIndexMemory += 4*ppDim;
     vnn.pppp[ichan].allocate(ppDim,ppDim);
     vnn.pppp[ichan].zeros();
	     
	     
     ppIndex = 0;
     for(int a = fermiLevel; a < Nspstates; a++){
       for(int b = fermiLevel; b < Nspstates; b++){
	 if( SPbasis->checkChanSym(a,b,ichan) == 1 ){
	   pp2Index = 0;
	   for(int c = fermiLevel; c < Nspstates; c++){
	     for(int d = fermiLevel; d < Nspstates; d++){
	       if( SPbasis->checkChanSym(c,d,ichan) == 1 ){
		 vpqrs = calcVpqrs(a,b,c,d,SPbasis);
		 vnn.pppp[ichan].setElement(ppIndex,pp2Index,vpqrs);
		 vnn.pppp[ichan].rowMap[ppIndex][0] = a;
		 vnn.pppp[ichan].rowMap[ppIndex][1] = b;
		 vnn.pppp[ichan].colMap[pp2Index][0] = c;
		 vnn.pppp[ichan].colMap[pp2Index][1] = d;			
		 pp2Index++;
	       } // end if 
	     } // end d loop
	   } // end c loop
	   ppIndex++;
	 } // end if		
       } // end b loop
     } // end a loop
   } // end ichan
      
   ppppMemory = memoryUsed + ppppIndexMemory;
   totalMemoryUsed += memoryUsed + ppppIndexMemory;   
   cout << "vnn.pppp requires: " << memoryUsed*sizeof(double)/1.e9 << " Gb." << endl;
   cout << "vnn.pppp requires: " << ppppIndexMemory*sizeof(int)/1.e9 << " Gb in index memory." << endl;
   
   if( totalMemoryUsed/1.e9 > memRequested ){
     cout << "Not requesting enough memory!" << endl;
   }
   
   ////////////////////////////////////////////////////////////////////////
   cout << "Setting up vnn.hhhh" << endl;
   ////////////////////////////////////////////////////////////////////////

   memoryUsed = 0;   
#pragma omp parallel for reduction(+:memoryUsed) 
   for(int ichan = 0; ichan < numChannels; ichan++){
     double vpqrs;
     int hhDim, hhIndex, hh2Index;

     hhDim = 0;	   
     for(int i = 0; i < fermiLevel; i++){
       for(int j = 0; j < fermiLevel; j++){	
	 if( SPbasis->checkChanSym(i,j,ichan) == 1 ){		 
	   hhDim++;		  
	 } // end if

       } // end j loop
     } // end i loop

     memoryUsed+= hhDim*hhDim;
     vnn.hhhh[ichan].allocate(hhDim,hhDim);
     vnn.hhhh[ichan].zeros();
	          
     hhIndex = 0;	    
     for(int i = 0; i < fermiLevel; i++){
       for(int j = 0; j < fermiLevel; j++){
	 if( SPbasis->checkChanSym(i,j,ichan) == 1 ){			   
	   hh2Index = 0;
	   for(int k = 0; k < fermiLevel; k++){
	     for(int l = 0; l < fermiLevel; l++){
	       if( SPbasis->checkChanSym(k,l,ichan) == 1 ){
		 vpqrs = calcVpqrs(i,j,k,l,SPbasis);
		 vnn.hhhh[ichan].setElement(hhIndex,hh2Index,vpqrs);
		 vnn.hhhh[ichan].rowMap[hhIndex][0] = i;
		 vnn.hhhh[ichan].rowMap[hhIndex][1] = j;
		 vnn.hhhh[ichan].colMap[hh2Index][0] = k;
		 vnn.hhhh[ichan].colMap[hh2Index][1] = l;			
		 hh2Index++;
	       } // end if 
	     } // end l loop
	   } // end k loop
	   hhIndex++;
	 } // end if		
       } // end j loop
     } // end i loop

   } // end ichan
   
   totalMemoryUsed += memoryUsed;
   cout << "vnn.hhhh requires: " << memoryUsed*sizeof(double)/1.e9 << " Gb." << endl;
   
   ////////////////////////////////////////////////////////////////////////
   cout << "Setting up vnn.hhhp" << endl;
   ////////////////////////////////////////////////////////////////////////
   
   memoryUsed = 0;     
#pragma omp parallel for reduction(+:memoryUsed) 
   for(int ichan = 0; ichan < numChannels; ichan++){
     int hhDim,hpDim,hhIndex,hpIndex;
     double vpqrs;
  
     hhDim = 0;	   
     // find sizes of the blocks, allocate memory, zero out.	    
     for(int i = 0; i < fermiLevel; i++){
       for(int j = 0; j < fermiLevel; j++){	
	 if( SPbasis->checkChanSym(i,j,ichan) == 1 ){		 
	   hhDim++;		  
	 } // end if

       } // end j loop
     } // end i loop
     hpDim = 0;
     for(int i = 0; i < fermiLevel; i++){
       for(int a = fermiLevel; a < Nspstates; a++){
	 if( SPbasis->checkChanSym(i,a,ichan) == 1 ){		  
	   hpDim++;		  
	 } // end if		 
       } // end a loop
     } // end i loop

     memoryUsed+= hhDim*hpDim;
     vnn.hhhp[ichan].allocate(hhDim,hpDim);
     vnn.hhhp[ichan].zeros();
	     
	     	    
     hhIndex = 0;	    
     for(int i = 0; i < fermiLevel; i++){
       for(int j = 0; j < fermiLevel; j++){
	 if( SPbasis->checkChanSym(i,j,ichan) == 1 ){	
		   
	   hpIndex = 0;
	   for(int k = 0; k < fermiLevel; k++){
	     for(int a = fermiLevel; a < Nspstates; a++){
	       if( SPbasis->checkChanSym(k,a,ichan) == 1 ){
		 vpqrs = calcVpqrs(i,j,k,a,SPbasis);
		 vnn.hhhp[ichan].setElement(hhIndex,hpIndex,vpqrs);
		 vnn.hhhp[ichan].rowMap[hhIndex][0] = i;
		 vnn.hhhp[ichan].rowMap[hhIndex][1] = j;
		 vnn.hhhp[ichan].colMap[hpIndex][0] = k;
		 vnn.hhhp[ichan].colMap[hpIndex][1] = a;			
		 hpIndex++;
	       } // end if 
	     } // end a loop
	   } // end k loop

	   hhIndex++;
	 } // end if		
       } // end j loop
     } // end i loop
   } // end ichan
   
   totalMemoryUsed += memoryUsed;
   cout << "vnn.hhhp requires: " << memoryUsed*sizeof(double)/1.e9 << " Gb." << endl;

   ////////////////////////////////////////////////////////////////////////
   cout << "Setting up vnn.hhpp" << endl;
   ////////////////////////////////////////////////////////////////////////
    
   memoryUsed = 0;     
#pragma omp parallel for reduction(+:memoryUsed) 
   for(int ichan = 0; ichan < numChannels; ichan++){
     int hhDim,ppDim,hhIndex,ppIndex;
     double vpqrs;
  	  
     hhDim = 0;	   
     // find sizes of the blocks, allocate memory, zero out.	    
     for(int i = 0; i < fermiLevel; i++){
       for(int j = 0; j < fermiLevel; j++){	
	 if( SPbasis->checkChanSym(i,j,ichan) == 1 ){		 
	   hhDim++;		  
	 } // end if
       } // end j loop
     } // end i loop

     ppDim = 0;
     for(int a = fermiLevel; a < Nspstates; a++){
       for(int b = fermiLevel; b < Nspstates; b++){
	 if( SPbasis->checkChanSym(a,b,ichan) == 1 ){			 
	   ppDim++;
	 } // end if
       } // end b loop
     } // end a loop

     memoryUsed+= hhDim*ppDim;
     vnn.hhpp[ichan].allocate(hhDim,ppDim);
     vnn.hhpp[ichan].zeros();
	     
	     	    
     hhIndex = 0;	    
     for(int i = 0; i < fermiLevel; i++){
       for(int j = 0; j < fermiLevel; j++){
	 if( SPbasis->checkChanSym(i,j,ichan) == 1 ){		   
	   ppIndex = 0;
	   for(int a = fermiLevel; a < Nspstates; a++){
	     for(int b = fermiLevel; b < Nspstates; b++){
	       if( SPbasis->checkChanSym(a,b,ichan) == 1 ){
		 vpqrs = calcVpqrs(i,j,a,b,SPbasis);
		 vnn.hhpp[ichan].setElement(hhIndex,ppIndex,vpqrs);
		 vnn.hhpp[ichan].rowMap[hhIndex][0] = i;
		 vnn.hhpp[ichan].rowMap[hhIndex][1] = j;
		 vnn.hhpp[ichan].colMap[ppIndex][0] = a;
		 vnn.hhpp[ichan].colMap[ppIndex][1] = b;
		 ppIndex++;
	       } // end if		      
	     } // end b loop
	   } // end a loop
	   hhIndex++;
	 } // end if		
       } // end j loop
     } // end i loop
   } // end ichan
   
   totalMemoryUsed += memoryUsed;
   cout << "vnn.hhpp requires: " << memoryUsed*sizeof(double)/1.e9 << " Gb." << endl;
   
   ////////////////////////////////////////////////////////////////////////
   cout << "Setting up vnn.hpph" << endl;
   ////////////////////////////////////////////////////////////////////////
      
   memoryUsed = 0;     
#pragma omp parallel for reduction(+:memoryUsed) 
   for(int ichan = 0; ichan < numChannels; ichan++){
     int hpDim,phDim,hpIndex,phIndex;
     double vpqrs;

     hpDim = 0;
     for(int i = 0; i < fermiLevel; i++){
       for(int a = fermiLevel; a < Nspstates; a++){
	 if( SPbasis->checkChanSym(i,a,ichan) == 1 ){	  
	   hpDim++;		  
	 } // end if		
       } // end a loop
     } // end i loop

     phDim = 0;
     for(int a = fermiLevel; a < Nspstates; a++){
       for(int i = 0; i < fermiLevel; i++){
	 if( SPbasis->checkChanSym(a,i,ichan) == 1 ){	  
	   phDim++;		  
	 } // end if		
       } // end a loop
     } // end i loop

     memoryUsed+= hpDim*phDim;
     vnn.hpph[ichan].allocate(hpDim,phDim);
     vnn.hpph[ichan].zeros();
	     	    
     hpIndex = 0;
     for(int i = 0; i < fermiLevel; i++){
       for(int a = fermiLevel; a < Nspstates; a++){
	 if( SPbasis->checkChanSym(i,a,ichan) == 1 ){
	   phIndex = 0;
	   for(int b = fermiLevel; b < Nspstates; b++){
	     for(int j = 0; j < fermiLevel; j++){	    	    
	       if( SPbasis->checkChanSym(b,j,ichan) == 1 ){
		 vpqrs = calcVpqrs(i,a,b,j,SPbasis);
		 vnn.hpph[ichan].setElement(hpIndex,phIndex,vpqrs);
		 vnn.hpph[ichan].rowMap[hpIndex][0] = i;		
		 vnn.hpph[ichan].rowMap[hpIndex][1] = a;
		 vnn.hpph[ichan].colMap[phIndex][0] = b;
		 vnn.hpph[ichan].colMap[phIndex][1] = j;
		 phIndex++;			
	       } // end if		      
	     } // end j loop
	   } // end b loop
	   hpIndex++;
	 } // end if
       } // end a loop
     } // end i loop
  
   } // end ichan
  
   totalMemoryUsed += memoryUsed;
   cout << "vnn.hpph requires: " << memoryUsed*sizeof(double)/1.e9 << " Gb." << endl;
   ///////////////////////////////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////////
   cout << "Setting up vnn.hppp" << endl;
   ////////////////////////////////////////////////////////////////////////
   
   memoryUsed = 0;     
#pragma omp parallel for reduction(+:memoryUsed) 
   for(int ichan = 0; ichan < numChannels; ichan++){
     int hpDim,ppDim,hpIndex,ppIndex;
     double vpqrs;

     hpDim = 0;
     for(int i = 0; i < fermiLevel; i++){
       for(int a = fermiLevel; a < Nspstates; a++){
	 if( SPbasis->checkChanSym(i,a,ichan) == 1 ){	  
	   hpDim++;		  
	 } // end if		
       } // end a loop
     } // end i loop

     ppDim = 0;
     for(int a = fermiLevel; a < Nspstates; a++){
       for(int b = fermiLevel; b < Nspstates; b++){
	 if( SPbasis->checkChanSym(a,b,ichan) == 1 ){			 
	   ppDim++;
	 } // end if
       } // end b loop
     } // end a loop

     memoryUsed+= hpDim*ppDim;
     vnn.hppp[ichan].allocate(hpDim,ppDim);
     vnn.hppp[ichan].zeros();
	     	    
     hpIndex = 0;
     for(int i = 0; i < fermiLevel; i++){
       for(int a = fermiLevel; a < Nspstates; a++){
	 if( SPbasis->checkChanSym(i,a,ichan) == 1 ){
	   ppIndex = 0;
	   for(int b = fermiLevel; b < Nspstates; b++){
	     for(int c = fermiLevel; c < Nspstates; c++){
	       if( SPbasis->checkChanSym(b,c,ichan) == 1 ){
		 vpqrs = calcVpqrs(i,a,b,c,SPbasis);
		 vnn.hppp[ichan].setElement(hpIndex,ppIndex,vpqrs);
		 vnn.hppp[ichan].rowMap[hpIndex][0] = i;
		 vnn.hppp[ichan].rowMap[hpIndex][1] = a;
		 vnn.hppp[ichan].colMap[ppIndex][0] = b;
		 vnn.hppp[ichan].colMap[ppIndex][1] = c;
		 ppIndex++;			
	       } // end if		      
	     } // end c loop
	   } // end b loop	   
	   hpIndex++;
	 } // end if
       } // end a loop
     } // end i loop	     
	    
   } // end ichan
   
   totalMemoryUsed += memoryUsed;
   cout << "vnn.hppp requires: " << memoryUsed*sizeof(double)/1.e9 << " Gb." << endl;

  
   
  
   ////////////////////////////////////////////////////////////
   cout << "Setting up vnn.phh_p" << endl;
   /////////////////////////////////////////////////////////////
   memoryUsed = 0;     
#pragma omp parallel for reduction(+:memoryUsed)   
   for(int a = fermiLevel; a < Nspstates; a++){
     int pDim = 1;     
     int phhDim = 0;
     

     for(int b = fermiLevel; b < Nspstates; b++){
       for(int i = 0; i < fermiLevel; i++){
   	 for(int j = 0; j < fermiLevel; j++){
   	   if( SPbasis->checkSympqrs(a,b,i,j) == 1 ){
   	     phhDim++;
   	   }
   	 } // end j loop
       } // end i loop
     } // end b loop
     vnn.phh_p[a-Nparticles].allocateWithLabels(phhDim,pDim,3,1);
     vnn.phh_p[a-Nparticles].zeros();  
     memoryUsed += phhDim*pDim;
   } // end a loop
   
#pragma omp parallel for 
   for(int a = fermiLevel; a < Nspstates; a++){
     int phhIndex = 0;
     double vpqrs;
     for(int b = fermiLevel; b < Nspstates; b++){
       for(int i = 0; i < fermiLevel; i++){
   	 for(int j = 0; j < fermiLevel; j++){
   	   if( SPbasis->checkSympqrs(a,b,i,j) == 1 ){	   	   
   	     vpqrs = calcVpqrs(i,j,a,b,SPbasis);
   	     vnn.phh_p[a-Nparticles].setElement(phhIndex,0,vpqrs);
   	     //cout << "a: " << a << " b: " << b << " i: " << i << " j: " << j << " vpqrs: " << vpqrs << endl;
   	     vnn.phh_p[a-Nparticles].rowMap[phhIndex][0] = b;
   	     vnn.phh_p[a-Nparticles].rowMap[phhIndex][1] = i;
   	     vnn.phh_p[a-Nparticles].rowMap[phhIndex][2] = j;
   	     vnn.phh_p[a-Nparticles].colMap[0][0] = a;	    
   	     phhIndex++;
   	   }
   	 } // end j loop
       } // end i loop
     } // end b loop   
   } // end a loop
 
   totalMemoryUsed += memoryUsed;
   cout << "vnn.phh_p requires: " << memoryUsed*sizeof(double)/1.e9 << " Gb." << endl;

   ////////////////////////////////////////////////////////////
   cout << "Setting up vnn.h_hpp" << endl;
   /////////////////////////////////////////////////////////////

   memoryUsed = 0;     
#pragma omp parallel for reduction(+:memoryUsed)  
   for(int i = 0; i < fermiLevel; i++){
     int hDim = 1;
     int hppDim = 0;
     for(int j = 0; j < fermiLevel; j++){
       for(int a = fermiLevel; a < Nspstates; a++){
	 for(int b = fermiLevel; b < Nspstates; b++){
	   if( SPbasis->checkSympqrs(i,j,a,b) == 1 ){
	     hppDim++;
	   }
	 } // end b loop
       } // end a loop
     } // end j loop
     vnn.h_hpp[i].allocateWithLabels(hDim,hppDim,1,3);
     vnn.h_hpp[i].zeros();  
     memoryUsed += hDim*hppDim;
   } // end i loop

#pragma omp parallel for
   for(int i = 0; i < fermiLevel; i++){
     int hppIndex = 0;
     double vpqrs;
     for(int j = 0; j < fermiLevel; j++){
       for(int a = fermiLevel; a < Nspstates; a++){
	 for(int b = fermiLevel; b < Nspstates; b++){
       	   if( SPbasis->checkSympqrs(i,j,a,b) == 1 ){	   	   
	     vpqrs = calcVpqrs(i,j,a,b,SPbasis);
	     vnn.h_hpp[i].setElement(0,hppIndex,vpqrs);
	     //cout << "a: " << a << " b: " << b << " i: " << i << " j: " << j << " vpqrs: " << vpqrs << endl;
	     vnn.h_hpp[i].rowMap[0][0] = i;
	     vnn.h_hpp[i].colMap[hppIndex][0] = j;
	     vnn.h_hpp[i].colMap[hppIndex][1] = a;
	     vnn.h_hpp[i].colMap[hppIndex][2] = b;
	    	    
	     hppIndex++;
	   }
	 } // end b loop
       } // end a loop
     } // end j loop   
   } // end i loop

   totalMemoryUsed += memoryUsed;
   cout << "vnn.h_hpp requires: " << memoryUsed*sizeof(double)/1.e9 << " Gb." << endl;

  
   
   ////////////////////////////////////////////////////////////////////////
   cout << "Setting up vnn.hhpp_hpph_mod" << endl;
   ////////////////////////////////////////////////////////////////////////
  

   memoryUsed = 0;   
#pragma omp parallel for reduction(+:memoryUsed) 
   for(int ichan = 0; ichan < numChannels; ichan++){
     int hpModDim,phModDim,phIndex,hpIndex;
     double vpqrs;

     hpModDim = 0;
     for(int i = 0; i < fermiLevel; i++){
       for(int a = fermiLevel; a < Nspstates; a++){
	 if( SPbasis->checkChanModSym(i,a,ichan) == 1 ){	  
	   hpModDim++;		  
	 } // end if		
       } // end a loop
     } // end i loop

     phModDim = 0;
     for(int a = fermiLevel; a < Nspstates; a++){
       for(int i = 0; i < fermiLevel; i++){
	 if( SPbasis->checkChanModSym(a,i,ichan) == 1 ){	  
	   phModDim++;		  
	 } // end if		
       } // end a loop
     } // end i loop

     memoryUsed+= hpModDim*phModDim;
     vnn.hhpp_hpph_mod[ichan].allocate(hpModDim,phModDim);
     vnn.hhpp_hpph_mod[ichan].zeros();
	     
	     	    
     hpIndex = 0;
     for(int i = 0; i < fermiLevel; i++){
       for(int a = fermiLevel; a < Nspstates; a++){
	 if( SPbasis->checkChanModSym(i,a,ichan) == 1 ){
	   phIndex = 0;
	   for(int b = fermiLevel; b < Nspstates; b++){
	     for(int j = 0; j < fermiLevel; j++){		      
	       if( SPbasis->checkChanModSym(b,j,ichan) == 1 ){		
		 vpqrs = calcVpqrs(i,j,a,b,SPbasis); 
		 vnn.hhpp_hpph_mod[ichan].setElement(hpIndex,phIndex,vpqrs);
		 vnn.hhpp_hpph_mod[ichan].rowMap[hpIndex][0] = i;		
		 vnn.hhpp_hpph_mod[ichan].rowMap[hpIndex][1] = a;
		 vnn.hhpp_hpph_mod[ichan].colMap[phIndex][0] = b;
		 vnn.hhpp_hpph_mod[ichan].colMap[phIndex][1] = j;
		 phIndex++;
	       } // end if		      
	     } // end j loop
	   } // end b loop			  
	   hpIndex++;
	 } // end if		
       } // end b loop
     } // end i loop	     

   } // end ichan

   
   totalMemoryUsed += memoryUsed;
   cout << "vnn.hhpp_hpph_mod requires: " << memoryUsed*sizeof(double)/1.e9 << " Gb." << endl;

   ///////////////////////////////////////////////////////////////////////////////////////////////


   ////////////////////////////////////////////////////////////////////////
   cout << "Setting up vnn.hpph_hphp_mod" << endl;
   ////////////////////////////////////////////////////////////////////////
   
   memoryUsed = 0;   
#pragma omp parallel for reduction(+:memoryUsed) 
   for(int ichan = 0; ichan < numChannels; ichan++){
     int hpModDim,hpIndex,hp2Index;
     double vpqrs;

	  
     hpModDim = 0;
     for(int i = 0; i < fermiLevel; i++){
       for(int a = fermiLevel; a < Nspstates; a++){
	 if( SPbasis->checkChanModSym(i,a,ichan) == 1 ){	  
	   hpModDim++;		  
	 } // end if		
       } // end a loop
     } // end i loop

	    
     memoryUsed+= hpModDim*hpModDim;
     vnn.hpph_hphp_mod[ichan].allocate(hpModDim,hpModDim);
     vnn.hpph_hphp_mod[ichan].zeros();	
	     	     	    
     hpIndex = 0;
     for(int i = 0; i < fermiLevel; i++){
       for(int b = fermiLevel; b < Nspstates; b++){
	 if( SPbasis->checkChanModSym(i,b,ichan) == 1 ){
	   hp2Index = 0;		   
	   for(int j = 0; j < fermiLevel; j++){
	     for(int a = fermiLevel; a < Nspstates; a++){		      
	       if( SPbasis->checkChanModSym(j,a,ichan) == 1 ){		
		 vpqrs = calcVpqrs(i,a,b,j,SPbasis); 
		 vnn.hpph_hphp_mod[ichan].setElement(hpIndex,hp2Index,vpqrs);
		 vnn.hpph_hphp_mod[ichan].rowMap[hpIndex][0] = i;		
		 vnn.hpph_hphp_mod[ichan].rowMap[hpIndex][1] = b;
		 vnn.hpph_hphp_mod[ichan].colMap[hp2Index][0] = j;
		 vnn.hpph_hphp_mod[ichan].colMap[hp2Index][1] = a;
		 
		 hp2Index++;
			 
	       } // end if	
		       	      
	     } // end b loop
	   } // end j loop
	   hpIndex++;
	 } // end if
       } // end a loop
     }	// end i loop
	     
   } // end ichan
   
   totalMemoryUsed += memoryUsed;
   cout << "vnn.hpph_hphp_mod requires: " << memoryUsed*sizeof(double)/1.e9 << " Gb." << endl;
   
   
   
   cout << endl << "total Memory required: " << totalMemoryUsed *sizeof(double)/1.e9 << " Gb." << endl;
   cout << "Vnn.pppp takes up " << 100.*ppppMemory/totalMemoryUsed << "% of that." << endl;
 } // end parallel channel blocks1

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////


void loadChannelBlocks2(TB_OpChannels &vnn, abstractSPbasis * SPbasis, double memRequested){
           
   int Nparticles = SPbasis->Nparticles;
   int Nspstates = SPbasis->Nspstates;
      
   int memoryUsed,totalMemoryUsed,ppppMemory,ppppIndexMemory;
   int fermiLevel = Nparticles;   
   int numChannels = SPbasis->Nchannels;
      
   totalMemoryUsed = 0;
   ppppIndexMemory = 0;
   ////////////////////////////////////////////////////////////////////////
   cout << "Setting up vnn.pppp" << endl;
   ////////////////////////////////////////////////////////////////////////
   
   
   memoryUsed = 0;
   totalMemoryUsed = 0;
#pragma omp parallel for reduction(+:memoryUsed,ppppIndexMemory) 
   for(int ichan = 0; ichan < numChannels; ichan++){
     double vpqrs;
     int ppDim,ppIndex,pp2Index;


     ppDim = 0;
     for(int a = fermiLevel; a < Nspstates; a++){
       for(int b = fermiLevel; b < Nspstates; b++){
	 if( SPbasis->checkChanSym(a,b,ichan) == 1 ){			 
	   ppDim++;
	 } // end if		 
       } // end b loop
     } // end a loop

     memoryUsed += ppDim*ppDim;
     ppppIndexMemory += 4*ppDim;
     vnn.pppp[ichan].allocate(ppDim,ppDim);
     vnn.pppp[ichan].zeros();
	     
	     
     ppIndex = 0;
     for(int a = fermiLevel; a < Nspstates; a++){
       for(int b = fermiLevel; b < Nspstates; b++){
	 if( SPbasis->checkChanSym(a,b,ichan) == 1 ){
	   pp2Index = 0;
	   for(int c = fermiLevel; c < Nspstates; c++){
	     for(int d = fermiLevel; d < Nspstates; d++){
	       if( SPbasis->checkChanSym(c,d,ichan) == 1 ){
		 vpqrs = calcVpqrs(a,b,c,d,SPbasis);
		 vnn.pppp[ichan].setElement(ppIndex,pp2Index,vpqrs);
		 vnn.pppp[ichan].rowMap[ppIndex][0] = a;
		 vnn.pppp[ichan].rowMap[ppIndex][1] = b;
		 vnn.pppp[ichan].colMap[pp2Index][0] = c;
		 vnn.pppp[ichan].colMap[pp2Index][1] = d;			
		 pp2Index++;
	       } // end if 
	     } // end d loop
	   } // end c loop
	   ppIndex++;
	 } // end if		
       } // end b loop
     } // end a loop
   } // end ichan
      
   ppppMemory = memoryUsed + ppppIndexMemory;
   totalMemoryUsed += memoryUsed + ppppIndexMemory;   
   cout << "vnn.pppp requires: " << memoryUsed*sizeof(double)/1.e9 << " Gb." << endl;
   cout << "vnn.pppp requires: " << ppppIndexMemory*sizeof(int)/1.e9 << " Gb in index memory." << endl;
   
   if( totalMemoryUsed/1.e9 > memRequested ){
     cout << "Not requesting enough memory!" << endl;
   }
   
   ////////////////////////////////////////////////////////////////////////
   cout << "Setting up vnn.hhhh" << endl;
   ////////////////////////////////////////////////////////////////////////

   memoryUsed = 0;   
#pragma omp parallel for reduction(+:memoryUsed) 
   for(int ichan = 0; ichan < numChannels; ichan++){
     double vpqrs;
     int hhDim, hhIndex, hh2Index;

     hhDim = 0;	   
     for(int i = 0; i < fermiLevel; i++){
       for(int j = 0; j < fermiLevel; j++){	
	 if( SPbasis->checkChanSym(i,j,ichan) == 1 ){		 
	   hhDim++;		  
	 } // end if

       } // end j loop
     } // end i loop

     memoryUsed+= hhDim*hhDim;
     vnn.hhhh[ichan].allocate(hhDim,hhDim);
     vnn.hhhh[ichan].zeros();
	          
     hhIndex = 0;	    
     for(int i = 0; i < fermiLevel; i++){
       for(int j = 0; j < fermiLevel; j++){
	 if( SPbasis->checkChanSym(i,j,ichan) == 1 ){			   
	   hh2Index = 0;
	   for(int k = 0; k < fermiLevel; k++){
	     for(int l = 0; l < fermiLevel; l++){
	       if( SPbasis->checkChanSym(k,l,ichan) == 1 ){
		 vpqrs = calcVpqrs(i,j,k,l,SPbasis);
		 vnn.hhhh[ichan].setElement(hhIndex,hh2Index,vpqrs);
		 vnn.hhhh[ichan].rowMap[hhIndex][0] = i;
		 vnn.hhhh[ichan].rowMap[hhIndex][1] = j;
		 vnn.hhhh[ichan].colMap[hh2Index][0] = k;
		 vnn.hhhh[ichan].colMap[hh2Index][1] = l;			
		 hh2Index++;
	       } // end if 
	     } // end l loop
	   } // end k loop
	   hhIndex++;
	 } // end if		
       } // end j loop
     } // end i loop

   } // end ichan
   
   totalMemoryUsed += memoryUsed;
   cout << "vnn.hhhh requires: " << memoryUsed*sizeof(double)/1.e9 << " Gb." << endl;
   
   ////////////////////////////////////////////////////////////////////////
   cout << "Setting up vnn.hhhp" << endl;
   ////////////////////////////////////////////////////////////////////////
   
   memoryUsed = 0;     
#pragma omp parallel for reduction(+:memoryUsed) 
   for(int ichan = 0; ichan < numChannels; ichan++){
     int hhDim,hpDim,hhIndex,hpIndex;
     double vpqrs;
  
     hhDim = 0;	   
     // find sizes of the blocks, allocate memory, zero out.	    
     for(int i = 0; i < fermiLevel; i++){
       for(int j = 0; j < fermiLevel; j++){	
	 if( SPbasis->checkChanSym(i,j,ichan) == 1 ){		 
	   hhDim++;		  
	 } // end if

       } // end j loop
     } // end i loop
     hpDim = 0;
     for(int i = 0; i < fermiLevel; i++){
       for(int a = fermiLevel; a < Nspstates; a++){
	 if( SPbasis->checkChanSym(i,a,ichan) == 1 ){		  
	   hpDim++;		  
	 } // end if		 
       } // end a loop
     } // end i loop

     memoryUsed+= hhDim*hpDim;
     vnn.hhhp[ichan].allocate(hhDim,hpDim);
     vnn.hhhp[ichan].zeros();
	     
	     	    
     hhIndex = 0;	    
     for(int i = 0; i < fermiLevel; i++){
       for(int j = 0; j < fermiLevel; j++){
	 if( SPbasis->checkChanSym(i,j,ichan) == 1 ){	
		   
	   hpIndex = 0;
	   for(int k = 0; k < fermiLevel; k++){
	     for(int a = fermiLevel; a < Nspstates; a++){
	       if( SPbasis->checkChanSym(k,a,ichan) == 1 ){
		 vpqrs = calcVpqrs(i,j,k,a,SPbasis);
		 vnn.hhhp[ichan].setElement(hhIndex,hpIndex,vpqrs);
		 vnn.hhhp[ichan].rowMap[hhIndex][0] = i;
		 vnn.hhhp[ichan].rowMap[hhIndex][1] = j;
		 vnn.hhhp[ichan].colMap[hpIndex][0] = k;
		 vnn.hhhp[ichan].colMap[hpIndex][1] = a;			
		 hpIndex++;
	       } // end if 
	     } // end a loop
	   } // end k loop

	   hhIndex++;
	 } // end if		
       } // end j loop
     } // end i loop
   } // end ichan
   
   totalMemoryUsed += memoryUsed;
   cout << "vnn.hhhp requires: " << memoryUsed*sizeof(double)/1.e9 << " Gb." << endl;

   ////////////////////////////////////////////////////////////////////////
   cout << "Setting up vnn.hhpp" << endl;
   ////////////////////////////////////////////////////////////////////////
    
   memoryUsed = 0;     
#pragma omp parallel for reduction(+:memoryUsed) 
   for(int ichan = 0; ichan < numChannels; ichan++){
     int hhDim,ppDim,hhIndex,ppIndex;
     double vpqrs;
  	  
     hhDim = 0;	   
     // find sizes of the blocks, allocate memory, zero out.	    
     for(int i = 0; i < fermiLevel; i++){
       for(int j = 0; j < fermiLevel; j++){	
	 if( SPbasis->checkChanSym(i,j,ichan) == 1 ){		 
	   hhDim++;		  
	 } // end if
       } // end j loop
     } // end i loop

     ppDim = 0;
     for(int a = fermiLevel; a < Nspstates; a++){
       for(int b = fermiLevel; b < Nspstates; b++){
	 if( SPbasis->checkChanSym(a,b,ichan) == 1 ){			 
	   ppDim++;
	 } // end if
       } // end b loop
     } // end a loop

     memoryUsed+= hhDim*ppDim;
     vnn.hhpp[ichan].allocate(hhDim,ppDim);
     vnn.hhpp[ichan].zeros();
	     
	     	    
     hhIndex = 0;	    
     for(int i = 0; i < fermiLevel; i++){
       for(int j = 0; j < fermiLevel; j++){
	 if( SPbasis->checkChanSym(i,j,ichan) == 1 ){		   
	   ppIndex = 0;
	   for(int a = fermiLevel; a < Nspstates; a++){
	     for(int b = fermiLevel; b < Nspstates; b++){
	       if( SPbasis->checkChanSym(a,b,ichan) == 1 ){
		 vpqrs = calcVpqrs(i,j,a,b,SPbasis);
		 vnn.hhpp[ichan].setElement(hhIndex,ppIndex,vpqrs);
		 vnn.hhpp[ichan].rowMap[hhIndex][0] = i;
		 vnn.hhpp[ichan].rowMap[hhIndex][1] = j;
		 vnn.hhpp[ichan].colMap[ppIndex][0] = a;
		 vnn.hhpp[ichan].colMap[ppIndex][1] = b;
		 ppIndex++;
	       } // end if		      
	     } // end b loop
	   } // end a loop
	   hhIndex++;
	 } // end if		
       } // end j loop
     } // end i loop
   } // end ichan
   
   totalMemoryUsed += memoryUsed;
   cout << "vnn.hhpp requires: " << memoryUsed*sizeof(double)/1.e9 << " Gb." << endl;
   
   ////////////////////////////////////////////////////////////////////////
   cout << "Setting up vnn.hpph" << endl;
   ////////////////////////////////////////////////////////////////////////
      
   memoryUsed = 0;     
#pragma omp parallel for reduction(+:memoryUsed) 
   for(int ichan = 0; ichan < numChannels; ichan++){
     int hpDim,phDim,hpIndex,phIndex;
     double vpqrs;

     hpDim = 0;
     for(int i = 0; i < fermiLevel; i++){
       for(int a = fermiLevel; a < Nspstates; a++){
	 if( SPbasis->checkChanSym(i,a,ichan) == 1 ){	  
	   hpDim++;		  
	 } // end if		
       } // end a loop
     } // end i loop

     phDim = 0;
     for(int a = fermiLevel; a < Nspstates; a++){
       for(int i = 0; i < fermiLevel; i++){
	 if( SPbasis->checkChanSym(a,i,ichan) == 1 ){	  
	   phDim++;		  
	 } // end if		
       } // end a loop
     } // end i loop

     memoryUsed+= hpDim*phDim;
     vnn.hpph[ichan].allocate(hpDim,phDim);
     vnn.hpph[ichan].zeros();
	     	    
     hpIndex = 0;
     for(int i = 0; i < fermiLevel; i++){
       for(int a = fermiLevel; a < Nspstates; a++){
	 if( SPbasis->checkChanSym(i,a,ichan) == 1 ){
	   phIndex = 0;
	   for(int b = fermiLevel; b < Nspstates; b++){
	     for(int j = 0; j < fermiLevel; j++){	    	    
	       if( SPbasis->checkChanSym(b,j,ichan) == 1 ){
		 vpqrs = calcVpqrs(i,a,b,j,SPbasis);
		 vnn.hpph[ichan].setElement(hpIndex,phIndex,vpqrs);
		 vnn.hpph[ichan].rowMap[hpIndex][0] = i;		
		 vnn.hpph[ichan].rowMap[hpIndex][1] = a;
		 vnn.hpph[ichan].colMap[phIndex][0] = b;
		 vnn.hpph[ichan].colMap[phIndex][1] = j;
		 phIndex++;			
	       } // end if		      
	     } // end j loop
	   } // end b loop
	   hpIndex++;
	 } // end if
       } // end a loop
     } // end i loop
  
   } // end ichan
  
   totalMemoryUsed += memoryUsed;
   cout << "vnn.hpph requires: " << memoryUsed*sizeof(double)/1.e9 << " Gb." << endl;
   ///////////////////////////////////////////////////////////////////////////////////////////////

   ////////////////////////////////////////////////////////////////////////
   cout << "Setting up vnn.hppp" << endl;
   ////////////////////////////////////////////////////////////////////////
   
   memoryUsed = 0;     
#pragma omp parallel for reduction(+:memoryUsed) 
   for(int ichan = 0; ichan < numChannels; ichan++){
     int hpDim,ppDim,hpIndex,ppIndex;
     double vpqrs;

     hpDim = 0;
     for(int i = 0; i < fermiLevel; i++){
       for(int a = fermiLevel; a < Nspstates; a++){
	 if( SPbasis->checkChanSym(i,a,ichan) == 1 ){	  
	   hpDim++;		  
	 } // end if		
       } // end a loop
     } // end i loop

     ppDim = 0;
     for(int a = fermiLevel; a < Nspstates; a++){
       for(int b = fermiLevel; b < Nspstates; b++){
	 if( SPbasis->checkChanSym(a,b,ichan) == 1 ){			 
	   ppDim++;
	 } // end if
       } // end b loop
     } // end a loop

     memoryUsed+= hpDim*ppDim;
     vnn.hppp[ichan].allocate(hpDim,ppDim);
     vnn.hppp[ichan].zeros();
	     	    
     hpIndex = 0;
     for(int i = 0; i < fermiLevel; i++){
       for(int a = fermiLevel; a < Nspstates; a++){
	 if( SPbasis->checkChanSym(i,a,ichan) == 1 ){
	   ppIndex = 0;
	   for(int b = fermiLevel; b < Nspstates; b++){
	     for(int c = fermiLevel; c < Nspstates; c++){
	       if( SPbasis->checkChanSym(b,c,ichan) == 1 ){
		 vpqrs = calcVpqrs(i,a,b,c,SPbasis);
		 vnn.hppp[ichan].setElement(hpIndex,ppIndex,vpqrs);
		 vnn.hppp[ichan].rowMap[hpIndex][0] = i;
		 vnn.hppp[ichan].rowMap[hpIndex][1] = a;
		 vnn.hppp[ichan].colMap[ppIndex][0] = b;
		 vnn.hppp[ichan].colMap[ppIndex][1] = c;
		 ppIndex++;			
	       } // end if		      
	     } // end c loop
	   } // end b loop	   
	   hpIndex++;
	 } // end if
       } // end a loop
     } // end i loop	     
	    
   } // end ichan
   
   totalMemoryUsed += memoryUsed;
   cout << "vnn.hppp requires: " << memoryUsed*sizeof(double)/1.e9 << " Gb." << endl;

  
   
  
   ////////////////////////////////////////////////////////////
   cout << "Setting up vnn.phh_p" << endl;
   /////////////////////////////////////////////////////////////
   memoryUsed = 0;     
#pragma omp parallel for reduction(+:memoryUsed)   
   for(int a = fermiLevel; a < Nspstates; a++){
     int pDim = 1;     
     int phhDim = 0;
     

     for(int b = fermiLevel; b < Nspstates; b++){
       for(int i = 0; i < fermiLevel; i++){
   	 for(int j = 0; j < fermiLevel; j++){
   	   if( SPbasis->checkSympqrs(a,b,i,j) == 1 ){
   	     phhDim++;
   	   }
   	 } // end j loop
       } // end i loop
     } // end b loop
     vnn.phh_p[a-Nparticles].allocateWithLabels(phhDim,pDim,3,1);
     vnn.phh_p[a-Nparticles].zeros();  
     memoryUsed += phhDim*pDim;
   } // end a loop
   
#pragma omp parallel for 
   for(int a = fermiLevel; a < Nspstates; a++){
     int phhIndex = 0;
     double vpqrs;
     for(int b = fermiLevel; b < Nspstates; b++){
       for(int i = 0; i < fermiLevel; i++){
   	 for(int j = 0; j < fermiLevel; j++){
   	   if( SPbasis->checkSympqrs(a,b,i,j) == 1 ){	   	   
   	     vpqrs = calcVpqrs(i,j,a,b,SPbasis);
   	     vnn.phh_p[a-Nparticles].setElement(phhIndex,0,vpqrs);
   	     //cout << "a: " << a << " b: " << b << " i: " << i << " j: " << j << " vpqrs: " << vpqrs << endl;
   	     vnn.phh_p[a-Nparticles].rowMap[phhIndex][0] = b;
   	     vnn.phh_p[a-Nparticles].rowMap[phhIndex][1] = i;
   	     vnn.phh_p[a-Nparticles].rowMap[phhIndex][2] = j;
   	     vnn.phh_p[a-Nparticles].colMap[0][0] = a;	    
   	     phhIndex++;
   	   }
   	 } // end j loop
       } // end i loop
     } // end b loop   
   } // end a loop
 
   totalMemoryUsed += memoryUsed;
   cout << "vnn.phh_p requires: " << memoryUsed*sizeof(double)/1.e9 << " Gb." << endl;

   ////////////////////////////////////////////////////////////
   cout << "Setting up vnn.h_hpp" << endl;
   /////////////////////////////////////////////////////////////

   memoryUsed = 0;     
#pragma omp parallel for reduction(+:memoryUsed)  
   for(int i = 0; i < fermiLevel; i++){
     int hDim = 1;
     int hppDim = 0;
     for(int j = 0; j < fermiLevel; j++){
       for(int a = fermiLevel; a < Nspstates; a++){
	 for(int b = fermiLevel; b < Nspstates; b++){
	   if( SPbasis->checkSympqrs(i,j,a,b) == 1 ){
	     hppDim++;
	   }
	 } // end b loop
       } // end a loop
     } // end j loop
     vnn.h_hpp[i].allocateWithLabels(hDim,hppDim,1,3);
     vnn.h_hpp[i].zeros();  
     memoryUsed += hDim*hppDim;
   } // end i loop

#pragma omp parallel for
   for(int i = 0; i < fermiLevel; i++){
     int hppIndex = 0;
     double vpqrs;
     for(int j = 0; j < fermiLevel; j++){
       for(int a = fermiLevel; a < Nspstates; a++){
	 for(int b = fermiLevel; b < Nspstates; b++){
       	   if( SPbasis->checkSympqrs(i,j,a,b) == 1 ){	   	   
	     vpqrs = calcVpqrs(i,j,a,b,SPbasis);
	     vnn.h_hpp[i].setElement(0,hppIndex,vpqrs);
	     //cout << "a: " << a << " b: " << b << " i: " << i << " j: " << j << " vpqrs: " << vpqrs << endl;
	     vnn.h_hpp[i].rowMap[0][0] = i;
	     vnn.h_hpp[i].colMap[hppIndex][0] = j;
	     vnn.h_hpp[i].colMap[hppIndex][1] = a;
	     vnn.h_hpp[i].colMap[hppIndex][2] = b;
	    	    
	     hppIndex++;
	   }
	 } // end b loop
       } // end a loop
     } // end j loop   
   } // end i loop

   totalMemoryUsed += memoryUsed;
   cout << "vnn.h_hpp requires: " << memoryUsed*sizeof(double)/1.e9 << " Gb." << endl;

  
   
   ////////////////////////////////////////////////////////////////////////
   cout << "Setting up vnn.hhpp_hpph_mod" << endl;
   ////////////////////////////////////////////////////////////////////////
  

   memoryUsed = 0;   
#pragma omp parallel for reduction(+:memoryUsed) 
   for(int ichan = 0; ichan < numChannels; ichan++){
     int hpModDim,phModDim,phIndex,hpIndex;
     double vpqrs;

     hpModDim = 0;
     for(int i = 0; i < fermiLevel; i++){
       for(int a = fermiLevel; a < Nspstates; a++){
	 if( SPbasis->checkChanModSym(i,a,ichan) == 1 ){	  
	   hpModDim++;		  
	 } // end if		
       } // end a loop
     } // end i loop

     phModDim = 0;
     for(int a = fermiLevel; a < Nspstates; a++){
       for(int i = 0; i < fermiLevel; i++){
	 if( SPbasis->checkChanModSym(a,i,ichan) == 1 ){	  
	   phModDim++;		  
	 } // end if		
       } // end a loop
     } // end i loop

     memoryUsed+= hpModDim*phModDim;
     vnn.hhpp_hpph_mod[ichan].allocate(hpModDim,phModDim);
     vnn.hhpp_hpph_mod[ichan].zeros();
	     
	     	    
     hpIndex = 0;
     for(int i = 0; i < fermiLevel; i++){
       for(int a = fermiLevel; a < Nspstates; a++){
	 if( SPbasis->checkChanModSym(i,a,ichan) == 1 ){
	   phIndex = 0;
	   for(int b = fermiLevel; b < Nspstates; b++){
	     for(int j = 0; j < fermiLevel; j++){		      
	       if( SPbasis->checkChanModSym(b,j,ichan) == 1 ){		
		 vpqrs = calcVpqrs(i,j,a,b,SPbasis); 
		 vnn.hhpp_hpph_mod[ichan].setElement(hpIndex,phIndex,vpqrs);
		 vnn.hhpp_hpph_mod[ichan].rowMap[hpIndex][0] = i;		
		 vnn.hhpp_hpph_mod[ichan].rowMap[hpIndex][1] = a;
		 vnn.hhpp_hpph_mod[ichan].colMap[phIndex][0] = b;
		 vnn.hhpp_hpph_mod[ichan].colMap[phIndex][1] = j;
		 phIndex++;
	       } // end if		      
	     } // end j loop
	   } // end b loop			  
	   hpIndex++;
	 } // end if		
       } // end b loop
     } // end i loop	     

   } // end ichan

   
   totalMemoryUsed += memoryUsed;
   cout << "vnn.hhpp_hpph_mod requires: " << memoryUsed*sizeof(double)/1.e9 << " Gb." << endl;

   ///////////////////////////////////////////////////////////////////////////////////////////////


   ////////////////////////////////////////////////////////////////////////
   cout << "Setting up vnn.hpph_hphp_mod" << endl;
   ////////////////////////////////////////////////////////////////////////
   
   memoryUsed = 0;   
#pragma omp parallel for reduction(+:memoryUsed) 
   for(int ichan = 0; ichan < numChannels; ichan++){
     int hpModDim,hpIndex,hp2Index;
     double vpqrs;

	  
     hpModDim = 0;
     for(int i = 0; i < fermiLevel; i++){
       for(int a = fermiLevel; a < Nspstates; a++){
	 if( SPbasis->checkChanModSym(i,a,ichan) == 1 ){	  
	   hpModDim++;		  
	 } // end if		
       } // end a loop
     } // end i loop

	    
     memoryUsed+= hpModDim*hpModDim;
     vnn.hpph_hphp_mod[ichan].allocate(hpModDim,hpModDim);
     vnn.hpph_hphp_mod[ichan].zeros();	
	     	     	    
     hpIndex = 0;
     for(int i = 0; i < fermiLevel; i++){
       for(int b = fermiLevel; b < Nspstates; b++){
	 if( SPbasis->checkChanModSym(i,b,ichan) == 1 ){
	   hp2Index = 0;		   
	   for(int j = 0; j < fermiLevel; j++){
	     for(int a = fermiLevel; a < Nspstates; a++){		      
	       if( SPbasis->checkChanModSym(j,a,ichan) == 1 ){		
		 vpqrs = calcVpqrs(i,a,b,j,SPbasis); 
		 vnn.hpph_hphp_mod[ichan].setElement(hpIndex,hp2Index,vpqrs);
		 vnn.hpph_hphp_mod[ichan].rowMap[hpIndex][0] = i;		
		 vnn.hpph_hphp_mod[ichan].rowMap[hpIndex][1] = b;
		 vnn.hpph_hphp_mod[ichan].colMap[hp2Index][0] = j;
		 vnn.hpph_hphp_mod[ichan].colMap[hp2Index][1] = a;
		 
		 hp2Index++;
			 
	       } // end if	
		       	      
	     } // end b loop
	   } // end j loop
	   hpIndex++;
	 } // end if
       } // end a loop
     }	// end i loop
	     
   } // end ichan
   
   totalMemoryUsed += memoryUsed;
   cout << "vnn.hpph_hphp_mod requires: " << memoryUsed*sizeof(double)/1.e9 << " Gb." << endl;
   
   
   
   cout << endl << "total Memory required: " << totalMemoryUsed *sizeof(double)/1.e9 << " Gb." << endl;
   cout << "Vnn.pppp takes up " << 100.*ppppMemory/totalMemoryUsed << "% of that." << endl;
 } // end parallel channel blocks2













































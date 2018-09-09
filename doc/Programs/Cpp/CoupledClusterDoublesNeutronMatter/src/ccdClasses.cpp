#include "ccdClasses.hpp"
#include <iostream>

//////////////////////////////////////////////
// IndexBundle
//////////////////////////////////////////////

IndexBundle::IndexBundle (){
  this->channel = 0;
  this->row = 0;
  this->col = 0;
}

IndexBundle::IndexBundle(int channelIn, int rowIn, int colIn){
  this->channel = channelIn;
  this->row = rowIn;
  this->col = colIn;
}

//////////////////////////////////////////////
// Amplitudes
//////////////////////////////////////////////


void Amplitudes::updateMods(){
  double value;
  int chanTemp1,rowTemp1,colTemp1;
  int chanTemp2,rowTemp2,colTemp2;
  int chanTemp4,rowTemp4,colTemp4;

  // set up t2.hpp_h, t2.p_phh, and t2.phhp
  for(int ichan = 0; ichan < this->numChannels; ichan++){
    for(int AB = 0; AB < this->pphh[ichan].getRowNum(); AB++){
      for(int IJ = 0; IJ < this->pphh[ichan].getColNum(); IJ++){
	value = this->pphh[ichan].getElement(AB,IJ);
	chanTemp2 = this->MegaMap[2][ichan][AB][IJ].channel;
	rowTemp2 = this->MegaMap[2][ichan][AB][IJ].row;
	colTemp2 = this->MegaMap[2][ichan][AB][IJ].col;
	this->p_phh[chanTemp2].setElement(rowTemp2,colTemp2,value);

	chanTemp1 = this->MegaMap[1][ichan][AB][IJ].channel;
	rowTemp1 = this->MegaMap[1][ichan][AB][IJ].row;
	colTemp1 = this->MegaMap[1][ichan][AB][IJ].col;
	this->hpp_h[chanTemp1].setElement(rowTemp1,colTemp1,value);

	chanTemp4 = this->MegaMap[4][ichan][AB][IJ].channel;
	rowTemp4 = this->MegaMap[4][ichan][AB][IJ].row;
	colTemp4 = this->MegaMap[4][ichan][AB][IJ].col;
	this->phhp[chanTemp4].setElement(rowTemp4,colTemp4,value);
      }
    }
  }
} // end updateMods

void Amplitudes::rotateSPChannels(){
  double value;
  int chanTemp0,rowTemp0,colTemp0;
  int chanTemp3,rowTemp3,colTemp3;
 
  // rotate t2.p_phh & t2.hpp_h
  for(int ichan = 0; ichan < this->numChannels; ichan++){
    for(int AB = 0; AB < this->pphh[ichan].getRowNum(); AB++){
      for(int IJ = 0; IJ < this->pphh[ichan].getColNum(); IJ++){
	value = this->pphh[ichan].getElement(AB,IJ);
	chanTemp3 = this->MegaMap[3][ichan][AB][IJ].channel;
	rowTemp3 = this->MegaMap[3][ichan][AB][IJ].row;
	colTemp3 = this->MegaMap[3][ichan][AB][IJ].col;
	this->p_phh[chanTemp3].setElement(rowTemp3,colTemp3,value);
	
	chanTemp0 = this->MegaMap[0][ichan][AB][IJ].channel;
	rowTemp0 = this->MegaMap[0][ichan][AB][IJ].row;
	colTemp0 = this->MegaMap[0][ichan][AB][IJ].col;
	this->hpp_h[chanTemp0].setElement(rowTemp0,colTemp0,value);	
      }
    }
  }
} // end rotateSPChannels

//////////////////////////////////////////////
// MegaMap
//////////////////////////////////////////////


void Amplitudes::generateMegaMap(int numChannels, abstractSPbasis * SPbasis){
  int Nparticles = SPbasis->Nparticles;  
  
  this->MegaMap = new IndexBundle***[8];
  for(int i = 0; i < 8; i++){
    this->MegaMap[i] = new IndexBundle**[numChannels];
  }

  IndexBundle temp;

  int a,b,i,j;
  int aa,bb,ii,jj;
  int rowNum,colNum;
  int chanTemp0,rowTemp0,colTemp0;
  int chanTemp1,rowTemp1,colTemp1;
  int chanTemp2,rowTemp2,colTemp2;
  int chanTemp3,rowTemp3,colTemp3;
  int chanTemp4,rowTemp4,colTemp4;
  int chanTemp5,rowTemp5,colTemp5;
  int chanTemp6,rowTemp6,colTemp6;
  int chanTemp7,rowTemp7,colTemp7;
  int memory = 0;

  // First Index: 8 partitions
  // 0: hpp_h: iab_j
  // 1: P(ij) hpp_h: jab_i
  // 2: p_phh: a_bij
  // 3: P(ab) p_phh: b_aij
  // 4: phhp: ai_jb
  // 5: P(ij) phhp: aj_ib
  // 6: P(ab) phhp: bi_ja
  // 7: P(ij)P(ab) phhp: bj_ia


  for(int ichan = 0; ichan < numChannels; ichan++){

    rowNum = this->pphh[ichan].getRowNum();
    colNum = this->pphh[ichan].getColNum();
    memory += rowNum*colNum;

    for(int i = 0; i < 8; i++){
      this->MegaMap[i][ichan] = new IndexBundle*[rowNum];
    }
    
    for(int AB = 0; AB < rowNum; AB++){      

      
      for(int i = 0; i < 8; i++){
	this->MegaMap[i][ichan][AB] = new IndexBundle[colNum];
      }
            
      for(int IJ = 0; IJ < colNum; IJ++){
	a = this->pphh[ichan].rowMap[AB][0];
	b = this->pphh[ichan].rowMap[AB][1];
	i = this->pphh[ichan].colMap[IJ][0];
	j = this->pphh[ichan].colMap[IJ][1];

	/////////////////////////////////
	// 0: hpp_h: iab_j
	/////////////////////////////////
	chanTemp0 = -1;
	rowTemp0 = -1;
	colTemp0 = -1;

	//if( checkSympqrs(p,q,r,s,indexMap) == 1){
	for(int IAB = 0; IAB < this->hpp_h[j].getRowNum(); IAB++){
	  ii = this->hpp_h[j].rowMap[IAB][0];
	  aa = this->hpp_h[j].rowMap[IAB][1];
	  bb = this->hpp_h[j].rowMap[IAB][2];
	  jj = this->hpp_h[j].colMap[0][0];
            
	  if( (i == ii) && (a == aa) && (b == bb) && (j == jj) ){
	    chanTemp0 = j;
	    rowTemp0 = IAB;
	    colTemp0 = 0;
	  } // end if	  
	} // end IAB

	this->MegaMap[0][ichan][AB][IJ].channel = chanTemp0;
	this->MegaMap[0][ichan][AB][IJ].row = rowTemp0;
	this->MegaMap[0][ichan][AB][IJ].col = colTemp0;

	/////////////////////////////////
	// 1: P(ij) hpp_h: jab_i
	/////////////////////////////////
	chanTemp1 = -1;
	rowTemp1 = -1;
	colTemp1 = -1;	

	for(int JAB = 0; JAB < this->hpp_h[i].getRowNum(); JAB++){
	  jj = this->hpp_h[i].rowMap[JAB][0];
	  aa = this->hpp_h[i].rowMap[JAB][1];
	  bb = this->hpp_h[i].rowMap[JAB][2];
	  ii = this->hpp_h[i].colMap[0][0];
      
	  if( (i == ii) && (a == aa) && (b == bb) && (j == jj) ){
	    chanTemp1 = i;
	    rowTemp1 = JAB;
	    colTemp1 = 0;
	  } // end if	  
	} // end IAB
	

	this->MegaMap[1][ichan][AB][IJ].channel = chanTemp1;
	this->MegaMap[1][ichan][AB][IJ].row = rowTemp1;
	this->MegaMap[1][ichan][AB][IJ].col = colTemp1;
	
	/////////////////////////////////
	// 2: p_phh: a_bij
	/////////////////////////////////
	chanTemp2 = -1;
	rowTemp2 = -1;
	colTemp2 = -1;

	for(int BIJ = 0; BIJ < this->p_phh[a-Nparticles].getColNum(); BIJ++){
	  aa = this->p_phh[a-Nparticles].rowMap[0][0];
	  bb = this->p_phh[a-Nparticles].colMap[BIJ][0];
	  ii = this->p_phh[a-Nparticles].colMap[BIJ][1];
	  jj = this->p_phh[a-Nparticles].colMap[BIJ][2];
          
      
	  if( (i == ii) && (a == aa) && (b == bb) && (j == jj) ){
	    chanTemp2 = a-Nparticles;
	    rowTemp2 = 0;
	    colTemp2 = BIJ;
	  } // end if	  
	} // end IAB
	

	this->MegaMap[2][ichan][AB][IJ].channel = chanTemp2;
	this->MegaMap[2][ichan][AB][IJ].row = rowTemp2;
	this->MegaMap[2][ichan][AB][IJ].col = colTemp2;

	/////////////////////////////////
	// 3: P(ab) p_phh: b_aij
	/////////////////////////////////
	chanTemp3 = -1;
	rowTemp3 = -1;
	colTemp3 = -1;

	for(int AIJ = 0; AIJ < this->p_phh[b-Nparticles].getColNum(); AIJ++){
	  bb = this->p_phh[b-Nparticles].rowMap[0][0];
	  aa = this->p_phh[b-Nparticles].colMap[AIJ][0];
	  ii = this->p_phh[b-Nparticles].colMap[AIJ][1];
	  jj = this->p_phh[b-Nparticles].colMap[AIJ][2];
      
	        
	  if( (i == ii) && (a == aa) && (b == bb) && (j == jj) ){
	    chanTemp3 = b-Nparticles;
	    rowTemp3 = 0;
	    colTemp3 = AIJ;
	  } // end if	  
	} // end IAB
	

	this->MegaMap[3][ichan][AB][IJ].channel = chanTemp3;
	this->MegaMap[3][ichan][AB][IJ].row = rowTemp3;
	this->MegaMap[3][ichan][AB][IJ].col = colTemp3;

	/////////////////////////////////
	// 4: phhp: ai_jb
	/////////////////////////////////
	chanTemp4 = -1;
	rowTemp4 = -1;
	colTemp4 = -1;

	temp = get_Mod_Bundle(a,i,j,b,numChannels,this->phhp,SPbasis);
	chanTemp4 = temp.channel;
	rowTemp4 = temp.row;
	colTemp4 = temp.col;
	
	this->MegaMap[4][ichan][AB][IJ].channel = chanTemp4;
	this->MegaMap[4][ichan][AB][IJ].row = rowTemp4;
	this->MegaMap[4][ichan][AB][IJ].col = colTemp4;

	/////////////////////////////////
	// 5: P(ij) phhp: aj_ib
	/////////////////////////////////
	chanTemp5 = -1;
	rowTemp5 = -1;
	colTemp5 = -1;

	temp = get_Mod_Bundle(a,j,i,b,numChannels,this->phhp,SPbasis);
	chanTemp5 = temp.channel;
	rowTemp5 = temp.row;
	colTemp5 = temp.col;
	
	this->MegaMap[5][ichan][AB][IJ].channel = chanTemp5;
	this->MegaMap[5][ichan][AB][IJ].row = rowTemp5;
	this->MegaMap[5][ichan][AB][IJ].col = colTemp5;

	/////////////////////////////////
	// 6: P(ab) phhp: bi_ja
	/////////////////////////////////
	chanTemp6 = -1;
	rowTemp6 = -1;
	colTemp6 = -1;

	temp = get_Mod_Bundle(b,i,j,a,numChannels,this->phhp,SPbasis);
	chanTemp6 = temp.channel;
	rowTemp6 = temp.row;
	colTemp6 = temp.col;
	
	this->MegaMap[6][ichan][AB][IJ].channel = chanTemp6;
	this->MegaMap[6][ichan][AB][IJ].row = rowTemp6;
	this->MegaMap[6][ichan][AB][IJ].col = colTemp6;

	/////////////////////////////////
	// 7: P(ij)P(ab) phhp: bj_ia
	/////////////////////////////////
	chanTemp7 = -1;
	rowTemp7 = -1;
	colTemp7 = -1;

	temp = get_Mod_Bundle(b,j,i,a,numChannels,this->phhp,SPbasis);
	chanTemp7 = temp.channel;
	rowTemp7 = temp.row;
	colTemp7 = temp.col;
	
	this->MegaMap[7][ichan][AB][IJ].channel = chanTemp7;
	this->MegaMap[7][ichan][AB][IJ].row = rowTemp7;
	this->MegaMap[7][ichan][AB][IJ].col = colTemp7;

	

	if( chanTemp7 == -1 || rowTemp7 == -1 || colTemp7 == -1){
	  std::cout << "from generate:" << std::endl;
	  std::cout << "temp.channel: " << temp.channel << " temp.row: " << temp.row << " temp.col: " << temp.col << std::endl;
	  
	  std::cout << "ichan: " << ichan << " AB: " << AB << " IJ: " << IJ << std::endl;
	  std::cout << "a: " << a << " b: " << b << " i: " << i << " j: " << j << std::endl;
	  std::cout << "chan: " << this->MegaMap[4][ichan][AB][IJ].channel 
	       << " row: " << this->MegaMap[4][ichan][AB][IJ].row 
	       << " col: " << this->MegaMap[4][ichan][AB][IJ].col << std::endl;
	}
	
      } // end IJ
    } // end AB
  } // end ichan

  memory = 8*sizeof(IndexBundle)*memory;
  std::cout << "T2Map requires " << memory/1.e9 << " Gb of memory." << std::endl;
} // end generate

// maybe broken
void Amplitudes::deallocateMegaMap(int numChannels){
  
  int rowNum;  

  for(int i = 0; i < 8; i++){
    for(int ichan = 0; ichan < numChannels; ichan++){
      rowNum = this->pphh[ichan].getRowNum();      
      for(int irow = 0; irow < rowNum; irow++){	
	  delete [] this->MegaMap[i][ichan][irow];	  
      } // end rowNum
      delete [] this->MegaMap[i][ichan];
    } // end ichan
    delete [] this->MegaMap[i];
  } // end i
  delete [] this->MegaMap;
}



IndexBundle Amplitudes::get_Mod_Bundle(int p, int q, int r, int s, int numChannels, SymmetryBlock* O_pqrs, abstractSPbasis * SPbasis){
  
  IndexBundle temp;
  temp.channel = -1;
  temp.row = -1;
  temp.col = -1;
  
  int p1,q1,r1,s1;
  
  // find the channel
  for(int ichan = 0; ichan < numChannels; ichan++){
    int chanDim = O_pqrs[ichan].getRowNum();
    if( chanDim == 0){
      continue;
    } 
    chanDim = O_pqrs[ichan].getColNum();
    if( chanDim == 0){
      continue;
    } 
    
    if( SPbasis->checkChanModSym(p,q,ichan) == 1 ){
       
      //now find the element
      for(int irow = 0; irow < O_pqrs[ichan].getRowNum(); irow++){
	for(int icol = 0; icol < O_pqrs[ichan].getColNum(); icol++){
	  p1 = O_pqrs[ichan].rowMap[irow][0];
	  q1 = O_pqrs[ichan].rowMap[irow][1];
	  r1 = O_pqrs[ichan].colMap[icol][0];
	  s1 = O_pqrs[ichan].colMap[icol][1];
	  if( (p == p1) && (q == q1) && (r == r1) && (s == s1) ){
	   
	    temp.channel = ichan;	 
	    temp.row = irow;	       
	    temp.col = icol;

	    return temp;
	  } // end if
	} // end col loop
      } // end row loop
    } // end if
  } // end ichan loop  
  
  return temp;
} // get_mod_Bundle


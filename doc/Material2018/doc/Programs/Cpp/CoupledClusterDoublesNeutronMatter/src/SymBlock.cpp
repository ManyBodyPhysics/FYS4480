#include "SymBlock.hpp"
#include <iostream>
#include <iomanip>

SymmetryBlock::SymmetryBlock () {
  this->Nrows = 0;
  this->Ncols = 0;
  this->rowLabelNum = 0;
  this->colLabelNum = 0;
} // end default constructor

SymmetryBlock::SymmetryBlock (int NrowsInput, int NcolsInput){
  this->allocate(NrowsInput, NcolsInput);
} // end functional constructor

SymmetryBlock::~SymmetryBlock () {
  if( this->Nrows != 0 || this->Ncols !=0 ){
    std::cout << "Destructor called before deallocate!" << std::endl;
  }     
} // end destructor

void SymmetryBlock::allocate(int NrowsInput, int NcolsInput){
  this->Nrows = NrowsInput;
  this->Ncols = NcolsInput;
  this->rowLabelNum = 2;
  this->colLabelNum = 2;
  this->pMatElements = new double[NrowsInput*NcolsInput];
  // for(int i=0; i<NrowsInput; i++) {
  //   this->pMatElements[i] = new double[NcolsInput];
  // }

  this->rowMap = new int*[NrowsInput];
  for(int i=0; i<NrowsInput; i++){
    this->rowMap[i] = new int[2];
  }
  
  this->colMap = new int*[NcolsInput];
  for(int j=0; j<NcolsInput; j++){
    this->colMap[j] = new int[2];
  }
} // end allocate

void SymmetryBlock::allocateWithLabels(int NrowsInput, int NcolsInput, int rowLabelNumInput, int colLabelNumInput){
  this->Nrows = NrowsInput;
  this->Ncols = NcolsInput;
  this->rowLabelNum = rowLabelNumInput;
  this->colLabelNum = colLabelNumInput;
  this->pMatElements = new double[NrowsInput*NcolsInput];
  // for(int i=0; i<NrowsInput; i++) {
  //   this->pMatElements[i] = new double[NcolsInput];
  // }

  this->rowMap = new int*[this->Nrows];
  for(int i=0; i<NrowsInput; i++){
    this->rowMap[i] = new int[this->rowLabelNum];
  }
  
  this->colMap = new int*[this->Ncols];
  for(int j=0; j<NcolsInput; j++){
    this->colMap[j] = new int[this->colLabelNum];
  }
} // end allocateWithLabels

void SymmetryBlock::deallocate(){
  
  // for(int i=0; i < this->Nrows; i++){
  //    delete [] this->pMatElements[i];
  // }
  
  for(int i=0; i < this->Nrows; i++){
    delete [] this->rowMap[i];
  }

  if( this->Nrows != 0 ){
    delete [] this->rowMap;
  }

  for(int j=0; j < this->Ncols; j++){
    delete [] this->colMap[j];
  }
  
  if( this->Ncols != 0 ){
    delete [] this->colMap;
  }

  if( this->Nrows != 0 && this->Ncols != 0){
    delete [] this->pMatElements;
  }

  this->Nrows = 0;
  this->Ncols = 0;
  this->rowLabelNum = 0;
  this->colLabelNum = 0;

} // end deallocate

int SymmetryBlock::getRowNum(){
  return this->Nrows;
} // end getRowNum

int SymmetryBlock::getColNum(){
  return this->Ncols;
} // end getColNum

double SymmetryBlock::getElement(int row, int col){
  return this->pMatElements[row*this->Ncols + col];
} // end getElement

void SymmetryBlock::setElement(int row, int col, double value){
  this->pMatElements[row*this->Ncols + col] = value;
} // end setElement

void SymmetryBlock::zeros(){
  for(int row=0; row < this->Nrows; row++){
    for(int col=0; col < this->Ncols; col++){
      this->pMatElements[row*this->Ncols + col] = 0.0;
    }
  }
  for(int row=0; row < this->Nrows; row++){
    for(int iRowLabel = 0; iRowLabel < this->rowLabelNum; iRowLabel++){
      this->rowMap[row][iRowLabel] = 0;
    }
    //this->rowMap[row][0] = 0;
    //this->rowMap[row][1] = 0;
  }
  for(int col=0; col < this->Ncols; col++){
    for(int iColLabel = 0; iColLabel < this->colLabelNum; iColLabel++){
      this->colMap[col][iColLabel] = 0;
    }
    //this->colMap[col][0] = 0;
    //this->colMap[col][1] = 0;
  }
} // end zeros

void SymmetryBlock::print(){
  std::cout << std::setprecision(4); 
  for(int irow=0; irow < this->Nrows; irow++){
    for(int icol=0; icol < this-> Ncols; icol++){
      std::cout << std::setw(8) << this->getElement(irow,icol) << " ";
    }
    std::cout << std::endl;
  }
} // end print

void SymmetryBlock::printRowMap(){  
  
  for(int irow=0; irow < this->Nrows; irow++){  
    std::cout << irow << ": <";  
    for(int iRowLabel = 0; iRowLabel < this->rowLabelNum; iRowLabel++){
	std::cout << this->rowMap[irow][iRowLabel] << ",";
      }       
    std::cout << "|" << std::endl;
  }
} // end printRowMap

void SymmetryBlock::printColMap(){  
  for(int icol=0; icol < this->Ncols; icol++){   
    std::cout << icol << ": |"; 
    for(int iColLabel = 0; iColLabel < this->colLabelNum; iColLabel++){
	std::cout << this->colMap[icol][iColLabel] << ",";
      }       
    std::cout << ">" << std::endl;     
  }
} // end printColMap

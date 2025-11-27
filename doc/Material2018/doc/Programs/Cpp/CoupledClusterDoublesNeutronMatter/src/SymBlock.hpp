#ifndef SYMBLOCK_HPP
#define SYMBLOCK_HPP

class SymmetryBlock {
  // private by default  
  int Nrows;
  int Ncols;  
  int rowLabelNum;
  int colLabelNum;
public:
  double * pMatElements;
  int ** rowMap;
  int ** colMap;
  SymmetryBlock ();
  SymmetryBlock (int,int);
  ~SymmetryBlock ();

  void allocate(int,int);
  //void allocate(int,int,int,int);
  void allocateWithLabels(int,int,int,int);
  void deallocate();
  int getRowNum();
  int getColNum();
  double getElement(int,int);
  void setElement(int,int,double);
  void zeros();
  void print();
  void printRowMap();
  void printColMap();
};

#endif /* SYMBLOCK_HPP */

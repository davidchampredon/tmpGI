//
//  dcMatrix.h
//  Epidemic_Models
//
//  Created by David Champredon on 12-05-27.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef dcMatrix_h
#define dcMatrix_h

#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>
#include <vector>
#include <assert.h>
#include <cstdlib>


using namespace std;

class Matrix
{
public:
	unsigned long nbRows;	//Dimensions
	unsigned long nbCols;
    
	vector<double> val;	//Value for each element
    
	double & operator()(unsigned long,unsigned long);
	double & operator[](unsigned long);
    
    unsigned long getNbRows()     {return nbRows;}
    unsigned long getNbCols()     {return nbCols;}
    
    
    // Constructors
	
	//Matrix(){nbRows=0;nbCols=0;}
	Matrix(){}
	
	Matrix(unsigned long h,unsigned long l) 
	{
		nbRows=h; nbCols=l; 
		vector<double> tmp(l*h,0.0);
		val=tmp;
	}
	
	Matrix(unsigned long n) 
	{
		nbRows=n; nbCols=n;
		vector<double> tmp(n*n, 0.0);
		val=tmp;
	}
	
    Matrix(string pathFile);
	
	Matrix(vector<double> v);	// Matrix (n,1) from a single vector(n)
	
	
    
    void    resize(unsigned long n)  
	{
		nbRows=n;nbCols=n; 
		val.resize(n*n);
	}
    
	void    resize(unsigned long nRows, unsigned long nCols)  
	{
		nbRows=nRows; nbCols=nCols; 
		this->val.resize(nRows*nCols);
	}
	
	void	clear() {val.clear();nbCols=0; nbRows=0;}
    
	void    display();
        
	void    RandomInit();
    
    // Files operations
	void    FromFile(const char*);
    void    FromFile(string);
	void	FromFile_Rows(string fileName, unsigned long nrow);
	
    void    WriteToFile(string);
	void    WriteToFileCSV(string);
	void    WriteToFileCSV(string fileName, vector<string> header);
	
	
    // Operations on embedded vectors
    vector<double>  extractColumn(unsigned long j_col);
	vector<double>	extractRow(unsigned long i_row);
	
	void            addRowVector(vector<double> v);
	void            addRowVector(vector<unsigned long> v);
	
    void            addColVector(vector<double> v);
    
	void			removeRow(unsigned long i_row);	// removes row 'i_row' and resize Matrix
	void			removeCol(unsigned long j_col);	// removes column 'j_col' and resize Matrix

	
	
	// Extract the row #i of the matrix
	// "i" is calculated such that it is
	// the smallest element of column "j_col"
  	vector<double>	extractRow_cond_minElement(unsigned long j_col);
	
    // Operations on elements
    
    void    setAllValues(double value);
	void	setValueFromMatrix(Matrix M);
	
	double  sumAllElements();
	double  sumLine(unsigned long i);		// sum all elements of line #i
	double  sumColumn(unsigned long j);	// sum all elements of line #i
	
	// conditional sum 
	double	sumColumnIf(unsigned long colToSum, unsigned long colToTest,
						double lowerBound, double upperBound);

	// counts nb of elements which are lower<element<upper  
	unsigned long		countColumnIf(unsigned long colToTest,
						  double lowerBound, double upperBound);
    
    void    setColumnValues(unsigned long colNb, vector<double> v);
	void	setRowValues(unsigned long rowNb_start0, vector<double> v);
	void	setRowValues(unsigned long rowNb_start0, vector<unsigned long> v);
	
	unsigned long countNonZeroElements();
	unsigned long countNonZeroElements_line(unsigned long i);
	
    Matrix  transpose();
    
    bool    isSymetric();

    double  determinant();
	
	Matrix	getMinor(unsigned long row, unsigned long col);
	
	Matrix	inverse();
    
    Matrix  Cholesky();
    
    double  getMinimumValue();
    double  getMaximumValue();
	
};

Matrix operator + (Matrix &A,Matrix &B);
Matrix operator - (Matrix &A,Matrix &B);
Matrix operator * (Matrix &A,Matrix &B);
Matrix operator * (double a,Matrix &A);

Matrix Id(unsigned long n);

Matrix power(Matrix A,int n);

Matrix cholesky(Matrix A);	//renvoi la Matrix triangul L tq://si A symetrique,carre //L*transpo(L)=A

double distance_Matrix(Matrix A, Matrix B, double power);	// Euclidian distance b/w two matrices

Matrix rowBind(Matrix A, Matrix B);

Matrix	transpo(Matrix A);        // FIX ME : to delete if not used elsewhere
double	Det(Matrix A);           // FIX ME : to delete if not used elsewhere
unsigned long		test_sym(Matrix A);       // FIX ME : to delete if not used elsewhere



#endif

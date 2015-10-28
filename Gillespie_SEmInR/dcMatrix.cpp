//
//  dcMatrix.cpp
//  Epidemic_Models
//
//  Created by David Champredon on 12-05-27.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include "dcMatrix.h"

 

// ////////////////////////////////////////////////////////////////////
//              CONSTRUCTORS
// ////////////////////////////////////////////////////////////////////


Matrix::Matrix(string fileName)
{
    ifstream thefile (fileName.c_str());
	
	for(unsigned long i=0;i<nbRows;i++)		
		for(unsigned long j=0;j<nbCols;j++)
		{
			thefile >> val[nbCols*i+j];
		}
	
}



Matrix::Matrix(vector<double> v)
{
	unsigned long nCol = v.size();
	nbRows = nCol;
	nbCols = 1;
	
	val.resize(nCol);
	
	for (unsigned long i=0; i<nCol; i++) 
	{
		val[i] = v[i];
	}
}


// ////////////////////////////////////////////////////////////////////
//              OPERATORS
// ////////////////////////////////////////////////////////////////////

double & Matrix::operator () (unsigned long i,unsigned long j)
{
	/// Retrieve matrix element
	/// top left element is M(0,0)
	/// bottom right element is M(n-1,m-1)
	
    assert(i>=0);assert(i<nbRows);
    assert(j>=0);assert(j<nbCols);
    
    return val[nbCols*i+j];
}

double & Matrix::operator [] (unsigned long i)
{
    assert(i>=0);assert(i<nbRows);assert(nbCols==1);
    return val[i];
}


// ////////////////////////////////////////////////////////////////////
//              Set Functions
// ////////////////////////////////////////////////////////////////////



void Matrix::RandomInit()
{
	srand ( time(NULL) );
	
	for(unsigned long i=0;i<nbRows;i++)		
		for(unsigned long j=0;j<nbCols;j++)
		{
			val[nbCols*i+j]=rand() % 11 - 10;
		}
    
}

void Matrix::setAllValues(double value)
{
    for(unsigned long i=0;i<nbRows;i++)		
		for(unsigned long j=0;j<nbCols;j++)
		{
			val[nbCols*i+j] = value;
		}
    
}

void Matrix::setValueFromMatrix(Matrix M)
{
	for(unsigned long i=0;i<nbRows;i++)		
		for(unsigned long j=0;j<nbCols;j++)
		{
			val[nbCols*i+j] = M(i,j);
		}
}


// ////////////////////////////////////////////////////////////////////
//              File Functions
// ////////////////////////////////////////////////////////////////////


void Matrix::FromFile(const char* fileName)
{
	ifstream thefile (fileName);
	
	for(unsigned long i=0;i<nbRows;i++)		
		for(unsigned long j=0;j<nbCols;j++)
		{
			thefile >> val[nbCols*i+j];
		}
	
}

void Matrix::FromFile(string fileName)
{
	// READ A FILE AND FILL THE 
	// ELEMENT VALUES TO THE MATRIX
	// *** WARNING *** MUST BE SAME SIZE!!!
	
	ifstream thefile (fileName.c_str());
	
	for(unsigned long i=0;i<nbRows;i++)		
		for(unsigned long j=0;j<nbCols;j++)
		{
			thefile >> val[nbCols*i+j];
		}	
}

void Matrix::FromFile_Rows(string fileName, unsigned long nrow)
{
	// READ A FILE AND FILL THE 
	// ELEMENT VALUES TO THE MATRIX

	// NUMBER OF ROWS PRE-SPECIFIED
	
	
	ifstream thefile (fileName.c_str());
	
	if (!thefile)
	{
		cout<<endl<<" ERROR [FromFile_Rows]: This file is not found: "<<fileName<<endl;
		exit(1);
	}
	
	vector<double> x;	
	double file_value;
	
	x.clear();
	
	while (!thefile.eof())   
	{
		thefile >> file_value;
		if (thefile.eof()) break;
		x.push_back (file_value);	
	}
	
	thefile.close();
	
	unsigned long n = x.size();
	
	if (n%nrow != 0 )
	{
		cout << endl << "ERROR [FromFile_Rows]: number of rows ("<< nrow<<") does not divide number of data("<<n<<")!"<<endl;
		exit(1);
	}

	this->resize(nrow, n/nrow);
	
	unsigned long cnt=0;
	
	for(unsigned long i=0;i<nbRows;i++)		
		for(unsigned long j=0;j<nbCols;j++)
		{
			val[nbCols*i+j] = x[cnt];
			cnt++;
		}	
	
}



void Matrix::WriteToFile(string fileName)
{
	ofstream thefile (fileName.c_str());
	
	for(unsigned long i=0;i<nbRows;i++)
    {
		for(unsigned long j=0;j<nbCols;j++)
		{
			thefile << val[nbCols*i+j] << "\t";
		}
        thefile << endl;
    }
}

void Matrix::WriteToFileCSV(string fileName)
{
	ofstream thefile (fileName.c_str());
	
	for(unsigned long i=0;i<nbRows;i++)
    {
		for(unsigned long j=0;j<nbCols;j++)
		{
			thefile << val[nbCols*i+j] ;
			if (j<nbCols-1) thefile<< ",";
		}
        thefile << endl;
    }
}

void Matrix::WriteToFileCSV(string fileName, vector<string> headers)
{
	ofstream thefile (fileName.c_str());
	
	if (nbCols!=headers.size())
	{
		cout << "ERROR -- WriteToFileCSV: Headers size ("<<headers.size()<<") does not match Matrix columns size ("<<nbCols<<")!"<<endl;
		cout << "Cannot write matrix to this file: "<<fileName<<endl;
		exit(1);
	}
	
	// Write headers on the first line
	for (unsigned long j=0; j<nbCols; j++) 
	{
		thefile << headers[j];
		if (j<nbCols-1) thefile<< ",";
	}
	thefile << endl;
	
	// Then write data
	for(unsigned long i=0;i<nbRows;i++)
    {
		for(unsigned long j=0;j<nbCols;j++)
		{
			thefile << val[nbCols*i+j] ;
			if (j<nbCols-1) thefile<< ",";
		}
        thefile << endl;
    }
}



// ////////////////////////////////////////////////////////////////////
//              Vector Functions
// ////////////////////////////////////////////////////////////////////



void Matrix::addRowVector(vector<double> v)
{    
	Matrix res;
	
	if (nbRows==0 && nbCols==0)
	{
		// If Matrix is empty, then this vector initialize the Matrix
		// as a 1 row matrix
		
		Matrix tmp(1,v.size());
		
		for (unsigned long j=0; j<v.size(); j++) 
		{
			tmp(0,j) = v[j];
		}
		res = tmp;
	}
	
	else 
	{
		// If Matrix is NOT empty
		
		if(nbCols != v.size())
		{
			cout<<"CAN'T ADD ROW VECTOR TO MATRIX, SIZES DO NOT MATCH: Matrix cols = "<< nbCols
			<<" vs Vector size = "<< v.size() << endl;
			exit(1);
		}
		
		Matrix tmp(nbRows+1,nbCols);
		
		for(unsigned long i=0; i<nbRows; i++)
		{
			for(unsigned long j=0; j<nbCols;j++)
			{
				tmp(i,j) = val[i*nbCols+j];//this->val[i*ncol+j];
			}
		}
		
		for(unsigned long j=0; j<nbCols;j++)
		{
			tmp(nbRows,j) = v[j];
		}
		
		res=tmp;
	}


    
	*this = res;


		
	/* WHY DOESN'T THIS WORK????
	
	
	vector<double> tmp = val;
	val.resize(nbCols*(nbRows+1));
	
	for (unsigned long j=0; j<nbCols; j++) 
	{
		val[nbRows*nbCols+j]=v[j];
	}
	*/
	
	/*
	 WHY DOESN'T THIS WORK????
	 
	 resize(nbRows+1, nbCols);
	 this->setRowValues(nbRows, v);
	
	 */
}


void Matrix::addRowVector(vector<unsigned long> v)
{
	vector<double> tmp;
	
	for (unsigned long i=0; i<v.size(); i++)
		tmp.push_back((double)(v[i]));
	
	addRowVector(tmp);
}


void Matrix::addColVector(vector<double> v)
{
	/// Add a column vector to a Matrix.
	/// If the matrix is empty, create a nx1 matrix
	
    unsigned long nrow = this->nbRows;
    unsigned long ncol = this->nbCols;
	
	
    if(nrow != v.size() && nrow>0)
	{
		cout<<"CAN'T ADD Col VECTOR TO MATRIX, SIZES DO NOT MATCH: Matrix rows = "<<nrow
		<<" vs Vector size = "<<v.size()<<endl;
		exit(1);
	}
    
    Matrix tmp(v.size(),ncol+1);
    
    for(unsigned long i=0; i<v.size(); i++)
    {
        for(unsigned long j=0; j<ncol;j++)
        {
            tmp(i,j) = this->val[i*ncol+j];
        }
    }
    
    for(unsigned long i=0; i<v.size(); i++)
    {
        tmp(i,ncol) = v[i];
    }
    
    *this = tmp;
}

void Matrix::removeRow(unsigned long i_row)
{
	for (unsigned long j=nbCols-1; j>=0; j--) 
	{
		val.erase(val.begin() + i_row*nbCols+j);
	}
	nbRows--;
}


void Matrix::removeCol(unsigned long j_col)
{
	for (unsigned long i=nbRows-1; i>=0; i--)
	{
		val.erase(val.begin() + i*nbCols+j_col);
	}
	nbCols--;
}


vector<double> Matrix::extractColumn(unsigned long j_col)
{
	if (j_col >= nbCols) {
		cout << endl << " ERROR [Matrix::extractColumn]:cannot extract col("<< j_col;
		cout << ") greater than size (0--"<< nbCols-1<< ")!"<<endl;
		exit(1);
	}
	
	vector<double> v(nbRows);
	
    for(unsigned long i=0; i<nbRows; i++) v[i]= val[i*nbCols+ j_col];
    
	return v;
}

vector<double> Matrix::extractRow(unsigned long i_row)
{
	if (i_row >= nbRows) {
		cout << endl << " ERROR [Matrix::extractRow]:cannot extract row("<< i_row;
		cout << ") greater than size (0--"<< nbRows-1<< ")!"<<endl;
		exit(1);
	}
	
	vector<double> v(nbCols);
	
    for(unsigned long j=0; j<nbCols; j++) 
		v[j]= val[i_row*nbCols+ j];
    
	return v;
}


vector<double>	Matrix::extractRow_cond_minElement(unsigned long j_col)
{
	if (j_col >= nbCols) {
		cout << endl << " ERROR [Matrix::extractRow_cond_minElement]:cannot extract row("<< j_col;
		cout << ") greater than size (0--"<< nbCols-1<< ")!"<<endl;
		exit(1);
	}
	
	vector<double> col_select = extractColumn(j_col);
	
	//argminElementVector(col_select); DOESN'T WORK WHEN #include "dcTools.h" ???
	
	double min = 999999999;
	unsigned long argmin = 0;
    
    for (unsigned long i=0; i<col_select.size(); i++) 
	{
        if (min > col_select[i]) 
		{
            min = col_select[i];
			argmin = i;
        }
    }
	
	return extractRow(argmin);
}




void Matrix::setColumnValues(unsigned long colNb_start0, vector<double> v)
{
    assert(v.size() == nbRows);
    
    for(unsigned long i=0; i<nbRows; i++) 
        val[i*nbCols+ colNb_start0] = v[i] ;
}


void Matrix::setRowValues(unsigned long rowNb_start0, vector<double> v)
{
	/// Set the values of a given row
	
    assert(v.size() == nbCols);
    
    for(unsigned long j=0; j<nbCols; j++) 
        val[rowNb_start0*nbCols+ j] = v[j] ;
}


void Matrix::setRowValues(unsigned long rowNb_start0, vector<unsigned long> v)
{
	/// Set the values of a given row
	
	assert(v.size() == nbCols);
	
	for(unsigned long j=0; j<nbCols; j++)
		val[rowNb_start0*nbCols+ j] = (double)(v[j]) ;
}



unsigned long Matrix::countNonZeroElements()
{
	/// Counts the number of non-zeros elements
	
	unsigned long cnt = 0;
	
	for (unsigned long i=0; i<nbRows; i++) {
		for (unsigned long j=0; j<nbCols; j++) {
			
			double x = 	val[nbCols*i+j];
			if(x>0 || x<0) cnt++;
		}
	}
	return cnt;
}


unsigned long Matrix::countNonZeroElements_line(unsigned long i)
{
	/// Counts the number of non-zeros elements for a given line only
	
	unsigned long cnt = 0;
	
	for (unsigned long j=0; j<nbCols; j++) {
		
		double x = 	val[nbCols*i+j];
		if(x>0 || x<0) cnt++;
	}
	return cnt;
}


// ////////////////////////////////////////////////////////////////////
//              Usual Matrix Functions
// ////////////////////////////////////////////////////////////////////



void Matrix::display()
{
    cout<<endl;
	cout << "Matrix dimension: "<<nbRows<<"x"<<nbCols<<endl;
    for(unsigned long i=0;i<nbRows;i++)
    {
        cout<<"[ ";
        for(unsigned long j=0;j<nbCols;j++)
        {
            cout<<val[i*nbCols+j];
			if (j<nbCols-1) cout<<"\t";
        }
        cout<<"]"<<endl;
    }
}

bool Matrix::isSymetric()
{
	if (nbCols!=nbRows) return false;
	unsigned long nCol = nbRows;
    unsigned long i=1,j=0;
	
	if (nCol==1) return 1;
	
    
	for (i=1;i<nCol;i++)
        for (j=0;j<i;j++)
            if ( val[i*nbCols+j]!=val[j*nbCols+i] ) {break;}
	if ( (i==nCol) && (j==nCol-1) ) return true;
	else return false;
}


Matrix Matrix::transpose()
{
    Matrix B(nbCols,nbRows);
    
    for (unsigned long i=0;i<nbCols;i++)
        for (unsigned long j=0;j<nbRows;j++)
        {
            B(i,j)=val[j*nbCols+i];
        }
    return B;
}


double Matrix::determinant()
{
    assert (nbCols==nbRows);

    double s = 0;
    unsigned long nCol = nbCols;

    if (nCol==2) return val[0*nbCols+0]*val[1*nbCols+1]
                        - val[0*nbCols+1]*val[1*nbCols+0]; //A(0,0)*A(1,1)-A(0,1)*A(1,0);
    
    else
    {
        Matrix B(nCol-1);
        for(unsigned long k=0;k<nCol;k++)
        {
            for(unsigned long i=0;i<nCol;i++)
                for(unsigned long j=1;j<nCol;j++) 
                { 
                    if (i<k) B(i,j-1)=val[i*nbCols+j];
                    if (i>k) B(i-1,j-1)=val[i*nbCols+j];
                }
            s += val[k*nbCols+0]*B.determinant()*pow (-1,k);
        } 
        return s;
    }
}

Matrix Matrix::getMinor(unsigned long row, unsigned long col)
{
	if (nbCols!=nbRows)
	{
		cout << endl<< "ERROR [Matrix::getMinor] : non square matrix" <<endl;
		exit(1);
	}
	
	unsigned long n = nbCols;
	
	Matrix B(n-1);
	// calculate the cofactor of element (row,col)
	
	// indicate which col and row is being copied to dest
	unsigned long colCount=0,rowCount=0;
	
	for(unsigned long i = 0; i < n; i++ )
	{
		if( i != row )
		{
			colCount = 0;
			for(unsigned long j = 0; j < n; j++ )
			{
				// when j is not the element
				if( j != col )
				{
					B(rowCount,colCount) = val[i*n+j];
					colCount++;
				}
			}
			rowCount++;
		}
	}

	return B;
}


Matrix Matrix::inverse()
{	
	if (nbRows!=nbCols)
	{
		cout << "ERROR [Matrix::inverse] : cannot inverse a non-square matrix!"<<endl;
		exit(1);
	}
	
	unsigned long n = nbCols;
	
	// Copy this Matrix unsigned longo A
	Matrix A(n);
	for (unsigned long i=0; i<n; i++) {
		for (unsigned long j=0; j<n; j++) {
			A(i,j) = val[n*i+j];
		}
	}
	
	// get the determinant of A
	double det = A.determinant();
	
	if (det==0)
	{
		cout << "ERROR [Matrix::inverse] : cannot inverse matrix (determinant	=0)!"<<endl;
		exit(1);
	}
	
    double oneOverDet = 1.0/det;
	
	Matrix AInverse(n);
	
	if (n==2)
	{
		AInverse(0,0) = A(1,1)*oneOverDet;
		AInverse(0,1) = -A(0,1)*oneOverDet;
		AInverse(1,0) = -A(1,0)*oneOverDet;
		AInverse(1,1) = A(0,0)*oneOverDet;
	}
	
	if (n>2)
	{
		Matrix Aminor(n-1);
		
		for(unsigned long j=0;j<n;j++)
		{
			for(unsigned long i=0;i<n;i++)
			{
				// get the co-factor (matrix) of A(j,i)
				Aminor = A.getMinor(j,i);
				
				AInverse(i,j) = oneOverDet*Aminor.determinant(); 
				
				if( (i+j)%2 == 1)
					AInverse(i,j) = -AInverse(i,j);
			}
		}
	}
	
	return AInverse;
}



Matrix Matrix::Cholesky()
{				
    if(nbCols!=nbRows) 		
    { cout<<"Cholesky: non square Matrix !"<<endl;
        exit(1);
    }
    
    if( ! this->isSymetric() ) 
    {  cout<<"Cholesky: non symetric Matrix !"<<endl;
        exit(1);
    }
	
    unsigned long nCol = nbCols;
    Matrix L(nCol);			
    float s;
    
    for(unsigned long i=0;i<nCol;i++)
    {
        s=0;
        for(unsigned long k=0;k<i;k++) s+=L(i,k)*L(i,k);
        L(i,i) = sqrt(val[i*nbCols+i]-s);
        
        for(unsigned long j=i+1;j<nCol;j++) 
        {
            s=0;
            for(unsigned long k=0;k<i;k++) s+=L(j,k)*L(i,k);
            L(j,i)=(val[i*nbCols+j]-s)/L(i,i); 
        }
    }
    return L;
} 


double Matrix::getMinimumValue()
{
    double min = 1e36;
    
    for (unsigned long i=0;i<nbCols;i++)
        for (unsigned long j=0;j<nbRows;j++)
        {
            double tmp = val[i*nbCols+j];
            if (tmp<min) min=tmp;
            
        }
    return min;
}


double Matrix::getMaximumValue()
{
    double maxi = -1e36;
    
    for (unsigned long i=0;i<nbCols;i++)
        for (unsigned long j=0;j<nbRows;j++)
        {
            double tmp = val[i*nbCols+j];
            if (maxi<tmp) maxi=tmp;
            
        }
    return maxi;
}

double Matrix::sumAllElements()
{
	double s=0;
	
	for(unsigned long i=0;i<nbRows;i++)
	{
		for(unsigned long j=0;j<nbCols;j++)
		{
			s=s+val[i*nbCols+j];
		}
        
	}
	
	return s;
}

double Matrix::sumLine(unsigned long i)
{
	double s=0.0;
	
	for(unsigned long k=0;k<nbCols;k++)
	{
        s += val[i*nbCols+k];
	}
	
	return s;
}

double Matrix::sumColumn(unsigned long j)
{
	double s=0.0;
	
	for(unsigned long k=0;k<nbRows ;k++)
	{
        s += val[k*nbCols+j];
	}
	
	return s;
}


double Matrix::sumColumnIf(unsigned long colToSum, unsigned long colToTest,
						   double lowerBound, double upperBound)
{
	double s=0.0;
	
	for(unsigned long k=0;k<nbRows ;k++)
	{
		if (val[k*nbCols+colToTest]>lowerBound &&
			val[k*nbCols+colToTest]<upperBound)
			
			s += val[k*nbCols+colToSum];
	}
	
	return s;
}

unsigned long Matrix::countColumnIf(unsigned long colToTest, double lowerBound, double upperBound)
{
	unsigned long c = 0;
	
	for(unsigned long k=0;k<nbRows ;k++)
	{
		if (val[k*nbCols+colToTest]>lowerBound &&
			val[k*nbCols+colToTest]<upperBound)
			
			c ++;
	}
	
	return c;
}
	

/* ************************************************************ */


Matrix operator + (Matrix &A,Matrix &B)
{
    assert(A.nbCols==B.nbCols && A.nbRows==B.nbRows);
    Matrix C(A.nbRows,A.nbCols);
    for (unsigned long i=0;i<A.nbRows;i++)
        for (unsigned long j=0;j<A.nbCols;j++)
        {
            C(i,j)=A(i,j)+B(i,j);
        }
    return C;   
}

Matrix operator - (Matrix &A,Matrix &B)
{
    assert(A.nbCols==B.nbCols && A.nbRows==B.nbRows);
    Matrix C(A.nbRows,A.nbCols);
    for (unsigned long i=0;i<A.nbRows;i++)
        for (unsigned long j=0;j<A.nbCols;j++)
        {
            C(i,j)=A(i,j)-B(i,j);
        }
    return C;   
}


Matrix operator * (Matrix &A,Matrix &B)
{
    assert(A.nbCols==B.nbRows);
    Matrix C(A.nbRows,B.nbCols);
    double s=0;
    for (unsigned long i=0;i<A.nbRows;i++)
        for (unsigned long j=0;j<B.nbCols;j++)
        {
            s=0;
            for(unsigned long k=0;k<A.nbCols;k++) s+=A(i,k)*B(k,j);
            C(i,j)=s;
        }
    return C;   
}

Matrix operator * (double a,Matrix &A)
{
    Matrix C(A.nbRows,A.nbCols);
    for (unsigned long i=0;i<A.nbRows;i++)
        for (unsigned long j=0;j<A.nbCols;j++)
        {
            C(i,j)=a*A(i,j);
        }
    return C;   
}





Matrix transpo(Matrix A)
{
    Matrix B(A.nbCols,A.nbRows);
    for (unsigned long i=0;i<A.nbCols;i++)
        for (unsigned long j=0;j<A.nbRows;j++)
        {
            B(i,j)=A(j,i);
        }
    return B;
}

unsigned long test_sym(Matrix A)
{
	if (A.nbCols!=A.nbRows) return 0;
	unsigned long nCol=A.nbRows,i=0,j=0;
	
	if (nCol==1) return 1;
	
    
	for (i=1;i<nCol;i++)
        for (j=0;j<i;j++)
            if ( A(i,j)!=A(j,i) ) {break;}
	if ( (i==nCol) && (j==nCol-1) ) return 1;
	else return 0;
}

Matrix Id(unsigned long nCol)
{
    Matrix I(nCol);
    for (unsigned long i=0;i<nCol;i++) 
        for (unsigned long j=0;j<=i;j++)
        {
            if(i==j) I(i,j)=1;
            else { I(i,j)=0; I(j,i)=0; }
        } 
    return I;
}



Matrix power(Matrix A,unsigned long nCol)
{
    if (A.nbCols!=A.nbRows)
    {cout<<"Power over non square matrix impossible!"<<endl;
        exit(1);}
    else
    {
        Matrix B(A.nbCols);
        B=Id(A.nbCols);
        for(unsigned long i=0;i<nCol;i++) {B=B*A;}
        return B;
    }
}

Matrix cholesky(Matrix A)	//renvoi la Matrix triangul L tq:
{				//si A symetrique,carre
    if(A.nbCols!=A.nbRows) 		//L*transpo(L)=A
    { cout<<"Cholesky(non memb): non square Matrix !"<<endl;
        exit(1);}
    if(!test_sym(A)) 
    {  cout<<"Cholesky(non memb): non symetric Matrix !"<<endl;
        exit(1);
    }
	
    unsigned long nCol=A.nbCols;
    Matrix L(nCol);			
    float s;
    for(unsigned long i=0;i<nCol;i++)
    {
        s=0;
        for(unsigned long k=0;k<i;k++) s+=L(i,k)*L(i,k);
        L(i,i)=sqrt(A(i,i)-s);
        
        for(unsigned long j=i+1;j<nCol;j++) 
        {
            s=0;
            for(unsigned long k=0;k<i;k++) s+=L(j,k)*L(i,k);
            L(j,i)=(A(i,j)-s)/L(i,i); 
        }
    }
    return L;
} 


bool same_size(Matrix A, Matrix B)
{
	bool res = true;
	unsigned long ra = A.getNbRows();
	unsigned long rb = B.getNbRows();
	unsigned long ca = A.getNbCols();
	unsigned long cb = B.getNbCols();
	
	if (ra!=rb || ca!=cb) res = false;
	
	return res;
}

double distance_Matrix(Matrix A, Matrix B, double power)
{
	if (!same_size(A,B))
	{
		cout << "ERROR [distance_Matrix]: matrix not the same size!"<<endl;
		exit(1);
	}
	
	double dist = 0.0;
	
	for (unsigned long i=0; i<A.getNbRows(); i++)
	{
		for (unsigned long j=0; j<A.getNbCols(); j++) {
			dist += pow(A(i,j)-B(i,j),power);
		}
	}
	
	return pow(dist,1.0/power);
}


Matrix rowBind(Matrix A, Matrix B)
{
	if (A.getNbCols()!=B.getNbCols())
	{
		cout << "ERROR [rowBind]: matrix not the same column size!"<<endl;
		exit(1);
	}
	
	Matrix M(A.getNbRows()+B.getNbRows(),A.getNbCols());
	unsigned long i=0;
	
	for (i=0; i<A.getNbRows(); i++) {
		for (unsigned long j=0; j<A.getNbCols(); j++) {
			M(i,j) = A(i,j);
		}
	}
	for (unsigned long ii=0; ii<B.getNbRows(); ii++) {
		for (unsigned long j=0; j<A.getNbCols(); j++) {
			M(i+ii,j) = B(ii,j);
		}
	}
	
	return M;	
}



double Det(Matrix A)
{
    assert (A.nbCols==A.nbRows);
    double s=0;
    unsigned long nCol=A.nbCols;
    if (nCol==2) return A(0,0)*A(1,1)-A(0,1)*A(1,0);
    else
    {
        Matrix B(nCol-1);
        for(unsigned long k=0;k<nCol;k++)
        {
            for(unsigned long i=0;i<nCol;i++)
                for(unsigned long j=1;j<nCol;j++) 
                { 
                    if (i<k) B(i,j-1)=A(i,j);
                    if (i>k) B(i-1,j-1)=A(i,j);
                }
            s+=A(k,0)*Det(B)*pow (-1,k);
        } 
        return s;
    }
}


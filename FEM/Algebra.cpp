#include "stdafx.h"
#include "algebra.h"
#include "symbol.h"
#include "math.h"

/////////////////////////////////////////////////////////////////////////////
//#ifdef _DEBUG 
//#define new DEBUG_NEW
//#undef THIS_FILE
//static char THIS_FILE[] = __FILE__;
//#endif
/////////////////////////////////////////////////////////////////////////////

void initVector(long n, double *v, int incx, double alfa)
{
	double *p=v;
	for(long i=0; i<n; i++)
	{
		*p=alfa; p+=incx;
	}
}
void initVector(long n, int *v, int incx, int alfa)
{
	int *p=v;
	for(long i=0; i<n; i++)
	{
		*p=alfa; p+=incx;
	}
}
void initVector0(long n, double *v, int incx)
{
	double *p=v;
	for(long i=0; i<n; i++)
	{
		*p=0.0; p+=incx;
	}
}
void initVector0(long n, double *v)
{
	double *p=v;
	for(long i=0; i<n; i++)
	{
		*p=0.0; p++;
	}
}
void initVector0(long n, int *v)
{
	int *p=v;
	for(long i=0; i<n; i++)
	{
		*p=0; p++;
	}
}
bool LUDcmpAndSolve( int n, double**A, double* b )
{
	bool rtn=false;
	rtn=LUDecomposition(n, A);
	LUSolve(n, A, b);

	return rtn;
}
bool LUDecomposition(int n, double**A)
{
	int i, j ,k;
	double big, sum, temp;
	double eps=1.0e-12;

	for(i=0; i<n; i++)
	{
		big=0.0;
		for(j=0; j<n; j++)
			if((temp=fabs(A[j][j])) > big) big=temp;
		if(big <= eps) return FALSE;
	}
	for(j=0; j<n; j++)
	{
		for(i=0; i<=j; i++)
		{
			sum=A[i][j];
			for(k=0; k<i; k++) sum-=A[i][k]*A[k][j];
			A[i][j]=sum;
		}
		for(i=j+1; i<n; i++)
		{
			sum=A[i][j];
			for(k=0; k<j; k++)
				sum-=A[i][k]*A[k][j];
			A[i][j]=sum/A[j][j];
		}
	}

	return TRUE;	
}
void LUSolve(int n, double**A, double* b )
{
	LowerSolve(n, A, b);
	UpperSolve(n, A, b);
}
void LowerSolve(int n, double**A, double* b )
{
	int i, j;
	double sum;
	
	for(i=0; i<n; i++)
	{
		sum=b[i];
		for(j=0; j<=i-1; j++) sum-=A[i][j]*b[j];
		b[i]=sum;
	}
}
void UpperSolve(int n, double**A, double* b )
{
	int i, j;
	double sum;

	for(i=n-1; i>=0; i--)
	{
		sum=b[i];
		for(j=i+1; j<n; j++) sum-=A[i][j]*b[j];
		b[i]=sum/A[i][i];
	}
}
double det(double **D,int n) //输入代表矩阵的二维数组、矩阵阶数，返回矩阵的行列式
{
	double d=0;

	// 一阶二阶直接计算
	if(n==1)d=D[0][0];
	if(n==2)d=D[0][0]*D[1][1]-D[0][1]*D[1][0];
	else{
		for(int k=0;k<n;k++){
			// 为代数余子式申请内存
			double **M;
			M=(double**)malloc((n-1)*sizeof(double*));
			for(int i=0;i<n-1;i++)
				M[i]=(double*)malloc((n-1)*sizeof(double));

			// 为代数余子式赋值
			for(int i=0;i<n-1;i++)
				for(int j=0;j<n-1;j++)
					M[i][j]=D[i+1][j<k?j:j+1];

			// 按第一行展开，递归计算行列式，注意元素0则不展开可以加快计算速度
			if(D[0][k])
				d+=D[0][k]*det(M,n-1)*(((2+k)%2)?-1:1);

			// 释放内存
			for(int i=0;i<n-1;i++)free(M[i]);
			free(M);
		}
	}
	return d;                   
}












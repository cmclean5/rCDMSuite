//define guards, so headers are declare only once.
#ifndef SPECTRALMODULARITY_H
#define SPECTRALMODULARITY_H

/*
  Use of the stack for memory allocation. This 
 should be faster for large networks, but will need to reset NSIZE 
 large enough for your network size, and then re-Make. 

 */

#include "Headers.h"
#include "network.h"
#include "nr3.h"

class SpectralModularity {

 public:
   SpectralModularity();
   SpectralModularity(network *, edgelist *, double *, int, int);
  ~SpectralModularity();
  int calculateSpectralModularity();
  void setMinCn( int );
  void settol  ( double );
  void setPrint( bool );
  
 private:
  void calculateB( double *, int );
  int  delta( int, int );
  void split( double *, int, int *, const char* );
  void deltaModularityMax( int, double & );
  void maxModularity( double & );
  void modifySplit( int );
  void deltaModularity( double & );
  void maximiseIndexVectors();
  void calculateEigenVectors();
  void setupMatrices();
  void assignSpace();
  void freeSpace();

  
  //global variable values  
  //static const bool PRINT = true;
  bool PRINT;
  static const int dummy  = -1000;
  static const int NSIZE  = 5000;
  static const int MSIZE  = NSIZE * NSIZE;

  double tol;//the tolerance value, 10^-5; eigenvalues below this threshold are not used
  int MINCn;//The minimum cluster size

  network *gg;
  double *A;
  double *Bgi;  //The Modularity matrix, Bgi
  int    NR_Bgi;//Number of rows of Bgi
  int    NC_Bgi;//Number of cols of Bgi
  int    M;//number of edges
  bool   usedBgi;
  
  double specQ;
  double NORM;
  double betai;

  int    MAXK;  //Counter storing the maximum community number so far
  
  double u       [NSIZE];

  int    SI      [NSIZE];
  int    si      [NSIZE];
  int    visited [NSIZE];
  int    keys_p  [NSIZE];
  int    keys_n  [NSIZE];  

  double Bgi_temp[MSIZE];
  
};

struct Symmeig {
	Int n;
	MatDoub z;
	VecDoub d,e;
	Bool yesvecs;

	Symmeig(MatDoub_I &a, Bool yesvec=true) : n(a.nrows()), z(a), d(n),
		e(n), yesvecs(yesvec)
	{
		tred2();
		tqli();
		sort();
	}
	Symmeig(VecDoub_I &dd, VecDoub_I &ee, Bool yesvec=true) :
		n(dd.size()), d(dd), e(ee), z(n,n,0.0), yesvecs(yesvec)
	{
		for (Int i=0;i<n;i++) z[i][i]=1.0;
		tqli();
		sort();
	}
	void sort() {
		if (yesvecs)
			eigsrt(d,&z);
		else
			eigsrt(d);
	}

        void eigsrt(VecDoub_IO &d, MatDoub_IO *v=NULL){
	Int k;
	Int n=d.size();
	for (Int i=0;i<n-1;i++) {
	  Doub p=d[k=i];
	  for (Int j=i;j<n;j++)
	    if (d[j] >= p) p=d[k=j];
	  if (k != i) {
	    d[k]=d[i];
	    d[i]=p;
	    if (v != NULL)
	      for (Int j=0;j<n;j++) {
		p=(*v)[j][i];
		(*v)[j][i]=(*v)[j][k];
		(*v)[j][k]=p;
	      }
	  }
	}
	}
  
        void tred2(){

	  Int l,k,j,i;
	  Doub scale,hh,h,g,f;
	  for (i=n-1;i>0;i--) {
	    l=i-1;
	    h=scale=0.0;
	    if (l > 0) {
	      for (k=0;k<i;k++)
		scale += abs(z[i][k]);
	      if (scale == 0.0)
		e[i]=z[i][l];
	      else {
		for (k=0;k<i;k++) {
		  z[i][k] /= scale;
		  h += z[i][k]*z[i][k];
		}
		f=z[i][l];
		g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
		e[i]=scale*g;
		h -= f*g;
		z[i][l]=f-g;
		f=0.0;
		for (j=0;j<i;j++) {
		  if (yesvecs)
		    z[j][i]=z[i][j]/h;
		  g=0.0;
		  for (k=0;k<j+1;k++)
		    g += z[j][k]*z[i][k];
		  for (k=j+1;k<i;k++)
		    g += z[k][j]*z[i][k];
		  e[j]=g/h;
		  f += e[j]*z[i][j];
		}
		hh=f/(h+h);
		for (j=0;j<i;j++) {
		  f=z[i][j];
		  e[j]=g=e[j]-hh*f;
		  for (k=0;k<j+1;k++)
		    z[j][k] -= (f*e[k]+g*z[i][k]);
		}
	      }
	    } else
	      e[i]=z[i][l];
	    d[i]=h;
	  }
	  if (yesvecs) d[0]=0.0;
	  e[0]=0.0;
	  for (i=0;i<n;i++) {
	    if (yesvecs) {
	      if (d[i] != 0.0) {
		for (j=0;j<i;j++) {
		  g=0.0;
		  for (k=0;k<i;k++)
		    g += z[i][k]*z[k][j];
		  for (k=0;k<i;k++)
		    z[k][j] -= g*z[k][i];
		}
	      }
	      d[i]=z[i][i];
	      z[i][i]=1.0;
	      for (j=0;j<i;j++) z[j][i]=z[i][j]=0.0;
	    } else {
	      d[i]=z[i][i];
	    }
	  }
	}


  void tqli(){

    Int m,l,iter,i,k;
    Doub s,r,p,g,f,dd,c,b;
    const Doub EPS=numeric_limits<Doub>::epsilon();
    for (i=1;i<n;i++) e[i-1]=e[i];
    e[n-1]=0.0;
    for (l=0;l<n;l++) {
	  iter=0;
	  do {
	    for (m=l;m<n-1;m++) {
	      dd=abs(d[m])+abs(d[m+1]);
	      if (abs(e[m]) <= EPS*dd) break;
	    }
	    if (m != l) {
	      if (iter++ == 30) ;//cout << "Warning: Too many iterations in tqli" << endl;
	      g=(d[l+1]-d[l])/(2.0*e[l]);
	      r=pythag(g,1.0);
	      g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
	      s=c=1.0;
	      p=0.0;
	      for (i=m-1;i>=l;i--) {
		f=s*e[i];
		b=c*e[i];
		e[i+1]=(r=pythag(f,g));
		if (r == 0.0) {
		  d[i+1] -= p;
		  e[m]=0.0;
		  break;
		}
		s=f/r;
		c=g/r;
		g=d[i+1]-p;
		r=(d[i]-g)*s+2.0*c*b;
		d[i+1]=g+(p=s*r);
		g=c*r-b;
		if (yesvecs) {
		  for (k=0;k<n;k++) {
		    f=z[k][i+1];
		    z[k][i+1]=s*z[k][i]+c*f;
		    z[k][i]=c*z[k][i]-s*f;
		  }
		}
	      }
	      if (r == 0.0 && i >= l) continue;
	      d[l] -= p;
	      e[l]=g;
	      e[m]=0.0;
	    }
	  } while (m != l);
    }
  }

  
  Doub pythag(const Doub a, const Doub b){
    Doub absa=abs(a), absb=abs(b);
    return (absa > absb ? absa*sqrt(1.0+SQR(absb/absa)) :
	    (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb))));
  }

};

#endif

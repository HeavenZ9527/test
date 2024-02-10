/*
 * svdSolveEquation.h
 *   利用SVD分解解方程Ax=b,  其解为x=V*diag(1/wi)*U'*b
 *   只能解非齐次方程组
*/

#include <iostream.h>
#include<math.h>

#define NRANSI
#include"nrutil.h"

void svdcmp(float **a, int m, int n, float w[], float **v)
{
	float pythag(float a, float b);
	int flag,i,its,j,jj,k,l,nm;
	float anorm,c,f,g,h,s,scale,x,y,z,*rv1;

	rv1=vector(1,n);
	g=scale=anorm=0.0;
	for (i=1;i<=n;i++) {
		l=i+1;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i <= m) {
			for (k=i;k<=m;k++) scale += (float)fabs(a[k][i]);
			if (scale) {
				for (k=i;k<=m;k++) {
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}
				f=a[i][i];
				g =  (float)-SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][i]=f-g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
					f=s/h;
					for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
				}
				for (k=i;k<=m;k++) a[k][i] *= scale;
			}
		}
		w[i]=scale *g;
		g=s=scale=0.0;
		if (i <= m && i != n) {
			for (k=l;k<=n;k++) scale += (float) fabs(a[i][k]);
			if (scale) {
				for (k=l;k<=n;k++) {
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g =  (float)-SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][l]=f-g;
				for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
				for (j=l;j<=m;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
					for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
				}
				for (k=l;k<=n;k++) a[i][k] *= scale;
			}
		}
		anorm=FMAX(anorm, (float)(fabs(w[i])+fabs(rv1[i])));
	}
	for (i=n;i>=1;i--) {
		if (i < n) {
			if (g) {
				for (j=l;j<=n;j++)
					v[j][i]=(a[i][j]/a[i][l])/g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
					for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
				}
			}
			for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=IMIN(m,n);i>=1;i--) {
		l=i+1;
		g=w[i];
		for (j=l;j<=n;j++) a[i][j]=0.0;
		if (g) {
			g= (float)1.0/g;
			for (j=l;j<=n;j++) {
				for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
				f=(s/a[i][i])*g;
				for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
			}
			for (j=i;j<=m;j++) a[j][i] *= g;
		} else for (j=i;j<=m;j++) a[j][i]=0.0;
		++a[i][i];
	}
	for (k=n;k>=1;k--) {
		for (its=1;its<=30;its++) {
			flag=1;
			for (l=k;l>=1;l--) {
				nm=l-1;
				if ((float)(fabs(rv1[l])+anorm) == anorm) {
					flag=0;
					break;
				}
				if ((float)(fabs(w[nm])+anorm) == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if ((float)(fabs(f)+anorm) == anorm) break;
					g=w[i];
					h=pythag(f,g);
					w[i]=h;
					h= (float)1.0/h;
					c=g*h;
					s = -f*h;
					for (j=1;j<=m;j++) {
						y=a[j][nm];
						z=a[j][i];
						a[j][nm]=y*c+z*s;
						a[j][i]=z*c-y*s;
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j=1;j<=n;j++) v[j][k] = -v[j][k];
				}
				break;
			}
			if (its == 30) nrerror("no convergence in 30 svdcmp iterations");
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f= (float)(((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y));
			g=pythag(f,1.0);
			f= (float)((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=pythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=1;jj<=n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=pythag(f,h);
				w[j]=z;
				if (z) {
					z= (float)1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=1;jj<=m;jj++) {
					y=a[jj][j];
					z=a[jj][i];
					a[jj][j]=y*c+z*s;
					a[jj][i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
	free_vector(rv1,1,n);
}


float pythag(float a, float b)
{
	float absa,absb;
	absa= (float)fabs(a);
	absb= (float)fabs(b);
	if (absa > absb) return  (float)(absa*sqrt(1.0+SQR(absb/absa)));
	else return  (float)((absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb))));
}

#undef NRANSI

/*
void main(void)
{
//	extern void svdcmp(float **a, int m, int n, float w[], float **v);

	int i,j,m=3,n=2;/ *确定为m*n的矩阵a[1..m,1..n]* /
	float **a,*w,**v,x[2],xT[2];
	float aT[3][2]={ 
		1,1,
		2,1,
		2,-1};
	float b[3]= {
		2,
		3,
		1.1f};
	a=matrix(1,m,1,n);
	v=matrix(1,n,1,n);
	w=vector(1,n);

	cout<<"the equation is:"<<endl;
	for(i=1;i<=m;i++)
	{
		for(j=1;j<=n;j++)
		{
			a[i][j]=aT[i-1][j-1];
			cout<<a[i][j]<<" \t  ";
		}
		cout<<" = "<<b[i-1]<<endl;
	}
	for(i=1;i<=n;i++)
	{
		for(j=1;j<=n;j++)
			v[i][j]=1.0;
		w[i]=0.0;
		xT[i-1]=0.0;
		x[i-1]=0.0;
	}

   svdcmp(a,  m, n, w, v);
  /////////////////////////////////
   cout<<"\nU is :"<<endl;
   for(j=1;j<=n;j++)
   {
	   for(i=1;i<=m;i++)
	   {
		   cout<<a[j][i]<<"   ";
	   }
	   cout<<endl;

   }

   cout<<"\n\nS:"<<endl;
   for(i=1;i<=n;i++)cout<<w[i]<<"   ";

   cout<<"\n\nV is: "<<endl;

   for(j=1;j<=n;j++)
   {
	   for(i=1;i<=n;i++)
	   {
		   cout<<v[j][i]<<"  ";
	   }
	   cout<<endl;
   }
////////////////////////////////////////

    //x=(v*diag(1/wi))*(u'*b)

   for(j=1;j<=n;j++)
   {
	   for(i=1;i<=n;i++)
	   {
		   if(fabs(w[j])<1e-6)
			   v[i][j]=v[i][j]*0;
		   else
			   v[i][j]=v[i][j]/w[j];
		}
   }
/ *   ///////////////////////////////////V*inv(s)
   cout<<"\n\nV is: "<<endl;

   for(j=1;j<=n;j++)
   {
	   for(i=1;i<=n;i++)
	   {
		   cout<<v[j][i]<<"  ";
	   }
	   cout<<endl;
   }
   /////////////////////////////////
* /
 
   
   for(j=1;j<=n;j++)
   {
	   for(i=1;i<=m;i++)
	   {
		   xT[j-1]=xT[j-1]+a[i][j]*b[i-1];
	   }

   }
/ *
   cout<<"\n\nUT*b"<<endl;
   for(i=1;i<=n;i++)
	   cout<<xT[i-1]<<endl;
* /
   
   cout<<"\n\nthe solution is:"<<endl;
   for(i=1;i<=n;i++)
   {
	   for(j=1;j<=n;j++)
	   {
		   x[i-1]=x[i-1]+v[i][j]*xT[j-1];
	   }
	 cout<<"x["<<i<<"]="<<x[i-1]<<endl;

   }

  cout<<"\n  press enter for exit!"<<endl;
  getchar();

}
*/
//g++ ZDS_R.cpp -O2 -lqd -o ZDS_R
#include<iostream>
#include <iomanip>
#include"qdMatrix.h"
using std::endl;
using std::cout;

#define ERR_TOL 1.0e-30
//#define real float
//#define real double
//#define real dd_real
#define real qd_real

real exact(real t,real lambda)
{
real Z;
Z=exp(lambda*t); //linear case
//Z=1./(1+exp(-lambda*t)); //logistic case
return Z;
}

real funD(real Z,real t,real lambda)
{
real D;
D=lambda*Z;// linear
//D=lambda*(1-Z)*Z;//logistic
return D;
}

real funS(real Z,real D,real t,real lambda)
{
real S;
S=lambda*D;// linear
//S=lambda*D*(1-2*Z);//logistic
return S;
}



real fixe_point(real Z0,real t, real Dt,real lambda,int R,qdMat<real> &coef)
{
// R structural equation and 2R physical equation
// coef(:,r) is the r+1 strutural equation, r=0:R-1
// given r, coef(i,r) are
//    i=0;R, coef for Z(0)...Z(R), i=R+1:2R+1, coef for D(0) to D(R), i=2R+2:3R+2
qdVec<real> Z(R+1),D(R+1),S(R+1),auxZ(R+1);
qdVec<real> aux(R),x(R),b(R);
qdMat<real> A(R,R),Ainv(R,R);
Z(0)=Z0;D(0)=funD(Z(0),t,lambda);S(0)=funS(Z(0),D(0),t,lambda);
for(unsigned int r=1;r<=R;r++)
    {
    Z(r)=Z(r-1)+Dt*D(r-1)+0.5*Dt*Dt*S(r-1);
    D(r)=funD(Z(r),t+r*Dt,lambda);S(r)=funS(Z(r),D(r),t+r*Dt,lambda);
    }
for(unsigned int r=0;r<R;r++)
    {
    aux(r)=(coef(0,r)*Z(0)+coef(R+1,r)*D(0)+coef(2*R+2,r)*S(0));
    for(int s=1;s<=R;s++)
        A(r,s-1)=coef(s,r);
    }

//A.disp();
Ainv=inverse(A);
real err=1;

while (err>ERR_TOL)
    {
    auxZ=Z;
    for(int r=0;r<R;r++) //use SEr to isolate Z1...ZR
        {
        b(r)=-aux(r);// start at zero
        for(int s=1;s<=R;s++)
            {
            b(r)-=(coef(R+1+s,r)*D(s)+coef(2*R+2+s,r)*S(s));
            }
        }
    x=Ainv*b;
    err=0.;
    for(int r=1;r<=R;r++) //update the values
        {
        Z(r)=x(r-1);D(r)=funD(Z(r),t+r*Dt,lambda);S(r)=funS(Z(r),D(r),t+r*Dt,lambda);
        err+=Abs((auxZ(r)-Z(r)));
        }
    //cout<<"error="<< std::scientific<<err<<endl;
    }
return Z(R);
}

real mypow(unsigned int a,unsigned int b)
{
real aux=1.;
real aa=a*1.0;
if(b==0) return aux;
for(unsigned i=0;i<b;i++)
   aux*=aa;
return aux;
}

real mypow(real a,unsigned int b)
{
real aux=1.;
if(b==0) return aux;
for(unsigned i=0;i<b;i++)
   aux*=a;
return aux;
}

real computeA(int k,int m)
    {
    if (m<k) return (real)0.;
    real aux=1.;
    for(int u=0;u<k;u++)
        aux*=(real)(m-u);
    //printf("k=%d, m=%d, A=%f\n",k,m,aux[0]);
    return aux;
    }

real computeQ(unsigned int k,unsigned int m,unsigned int r,real Dt)
    {
    real aux=0;
    if(m<k) return (real)0.;
    if (m==k) return (real)1./mypow(Dt,k);
    aux=mypow(r,m-k);
    aux/=mypow(Dt,k);
    return aux;
    }

qdMat<real>  computeCoef(int K,int R,int M,int S,real Dt)
    {
    assert(S>0);
    qdMat<real> MC(M+1,M+1);
    for (int m=0;m<=M;m++)
        for(int k=0;k<=K;k++)
            for(int r=0;r<=R;r++)
                MC(m,k*(R+1)+r)=computeA(k,m)*computeQ(k,m,r,Dt);
    //MC.disp();
    qdMat<real> MN= MC.extractRow(0,M-S);
    //MN.disp();
    qdMat<real> coef=kern(MN);
    //(MN*coef).disp();
    return coef;
    }



int main()
    {
    unsigned int oldcw;
    fpu_fix_start(&oldcw);
    // [[Z0 Z1,...ZR] [D0 D1,... DR] [S0 S1 SR], ....]
    int K=2; //mandatory for this code
    int R=2;// number of stage
    int M=(K+1)*(R+1)-1; // number of structural equation
    real T = 1.0,N = R*5.0,Dt=T/N;
    // build the structure equations
    qdMat<real> coef=computeCoef(K,R,M,R,Dt);
    //coef.disp();
    real lambda = -1.0;
    real x = 1.0,t=0,err;
    cout.precision(10);
    for (unsigned int n=0;n<N;n=n+R)
        {
        t=n*Dt;
        x=fixe_point(x,t,Dt,lambda,R,coef);
        //cout<<"x="<<std::scientific<<x<<endl;
        }
    err=Abs(  (x-exact(T,lambda)) )/exact(T,lambda);
    cout<<"error="<< std::scientific<<err<<endl;
    //cout<<"delta<<"<<std::scientific<<1.-exact(T,lambda)<<endl;
    fpu_fix_end(&oldcw);
    return 0;
    }


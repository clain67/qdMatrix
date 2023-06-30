//g++ ZDS_sys.cpp -O2 -lqd -o ZDS_sys
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


//system of P equation, v(p) is the vector of the variable p, for the different r and s
// v(p) can be "quase" treated as a scalar cases. the only par in the system is when we call funD and funS.
// coefficient matrices together with inverse matrice are the same (the structural equations are the same whatever the system)
typedef qdVec<qdVec<real>> qdvvect;

real mypow(unsigned int a,unsigned int b)
    {
    real aux=1.;
    real aa=a*1.0;
    if(b==0) return aux;
    for(unsigned int i=0;i<b;i++)
        aux*=aa;
    return aux;
    }

real mypow(real a,unsigned int b)
    {
    real aux=1.;
    if(b==0) return aux;
    for(unsigned int i=0;i<b;i++)
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

/*
real computeQ(int k,int m,int r)
    {
    real aux=0;
    if(m<k) return (real)0.;
    if (m==k) return (real)1.;
    aux=mypow(r,m-k);
    return aux;
    }*/

real computeQ(unsigned int k,unsigned int m,unsigned int r,real Dt)
    {
    real aux=0;
    if(m<k) return (real)0.;
    if (m==k) return (real)1./mypow(Dt,k);
    aux=mypow(r,m-k);
    aux/=mypow(Dt,k);
    return aux;
    }

    qdMat<real> computeCoef(int K,int R,int M,int S,real Dt)
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


/* ============================ EDO data =======================*/
qdVec<real>  exact(unsigned int Dim, real t,qdMat<real> &Lambda)
    {
    qdVec<real> Z(Dim);
    for(unsigned int i=0;i<Dim;i++)
        Z(i)=exp(Lambda(i,i)*t); //linear case
    //Z=1./(1+exp(-lambda*t)); //logistic case
    return Z;
    }

qdVec<real>  funD(unsigned int Dim, qdVec<real>  Z,real t,qdMat<real> &Lambda)
    {
    qdVec<real> D(Dim);
    D=(Lambda*Z);// linear
    //D=lambda*(1-Z)*Z;//logistic
    return D;
    }

qdVec<real>  funS(unsigned int Dim, qdVec<real>  Z,qdVec<real>  D,real t,qdMat<real> &Lambda)
    {
    qdVec<real> S(Dim);
    S=Lambda*D;// linear
    //S=lambda*D*(1-2*Z);//logistic
    return S;
    }

/* ======================== fixe point ==================================*/
qdVec<real> fixe_point(qdVec<real> Z0,real t, real Dt,qdMat<real> &Lambda,int R,qdMat<real> &coef,qdMat<real> &Ainv)
    {
    // R structural equation and 2R physical equation
    // coef(:,r) is the r+1 strutural equation, r=0:R-1
    // given r, coef(i,r) are
    // Z(i),D(i), S(i) a vector of size Dim
    //    i=0;R, coef for Z(0)...Z(R), i=R+1:2R+1, coef for D(0) to D(R), i=2R+2:3R+2
// -------------------- define and size the vectors -----------------
    unsigned int Dim=Z0.size();
    qdvvect Z(R+1),D(R+1),S(R+1),aux(R+1),auxZ(R+1);
    for(unsigned int r=0;r<=R;r++)  //[R+1][Dim]
        {
        Z(r).resize(Dim);
        aux(r).resize(Dim);
        auxZ(r).resize(Dim);
        D(r).resize(Dim);
        S(r).resize(Dim);
        }
    qdvvect x(Dim),b(Dim);
    for(unsigned int d=0;d<Dim;d++)// [Dim][R]
        {
        x(d).resize(R);
        b(d).resize(R);
        }
// ------------------ initialization -----------------------
    Z(0)=Z0;D(0)=funD(Dim,Z(0),t,Lambda);S(0)=funS(Dim,Z(0),D(0),t,Lambda);// [Dim]
    for(unsigned int r=1;r<=R;r++)
        {
        Z(r)=Z(r-1)+Dt*D(r-1)+0.5*Dt*Dt*S(r-1); //[Dim]
        D(r)=funD(Dim,Z(r),t+r*Dt,Lambda);S(r)=funS(Dim,Z(r),D(r),t+r*Dt,Lambda); //[R+1][Dim]
        }
    for(unsigned int r=0;r<R;r++)
        {
        aux(r)=(coef(0,r)*Z(0)+coef(R+1,r)*D(0)+coef(2*R+2,r)*S(0)); //[0:R-1][Dim]
        }
// ---------------------- the fix point loop ------------------------
    real err=1;
    while (err>ERR_TOL)
        {
        for(unsigned int r=1;r<=R;r++) // save the last step
            {
            auxZ(r)=Z(r);
            }
        for(unsigned int d=0;d<Dim;d++)    // ALERT index permutation!
            for(int r=0;r<R;r++) //use SEr to isolate Z1...ZR
                {
                b(d)(r)=-aux(r)(d);// start at zero //b[0:R-1][Dim]
                for(int s=1;s<=R;s++)
                    b(d)(r)-=(coef(R+1+s,r)*D(s)(d)+coef(2*R+2+s,r)*S(s)(d)); //b[0:R-1][Dim]
                }
        // we have to solve x[:,d]=A*b[:,d],
        for(unsigned int d=0;d<Dim;d++)
            x(d)=Ainv*b(d);
        err=0.;
        for(unsigned int d=0;d<Dim;d++)
            for(int r=1;r<=R;r++) //update the values
                Z(r)(d)=x(d)(r-1);
        for(int r=1;r<=R;r++)
            {
            D(r)=funD(Dim,Z(r),t+r*Dt,Lambda);S(r)=funS(Dim,Z(r),D(r),t+r*Dt,Lambda);
            err+=norm1(auxZ(r)-Z(r));
            }
        //cout<<"error="<< std::scientific<<err<<endl;
        }
    //cout<<"----------"<<endl;
    return Z(R);
    }





/* ======================= MAIN ============================*/
int main()
    {
    unsigned int oldcw;
    fpu_fix_start(&oldcw);
    // [[Z0 Z1,...ZR] [D0 D1,... DR] [S0 S1 SR], ....]
    unsigned int Dim=2; // size of the system
    unsigned int K=2; //mandatory for this code
    unsigned int R=2;// number of stage
    unsigned int M=(K+1)*(R+1)-1; // number of structural equation

// ---- build the structure equations and associated stuff ----------------
    real T = 1.0,N = R*10.0,Dt=T/N;
    qdMat<real> coef=computeCoef(K,R,M,R,Dt);
    //coef.disp();
    qdMat<real> A(R,R),Ainv(R,R); //compute the inverse matrix
    for(unsigned int r=0;r<R;r++)
        {
        for(unsigned int s=1;s<=R;s++)
            A(r,s-1)=coef(s,r);
        }
    //A.disp();
    Ainv=inverse(A);

//------------------------ parametrize the EDO problem and discretization -------------------

    real lambda = -1.0;
    qdMat<real> Lambda;Lambda.eyes(Dim,Dim);Lambda*=lambda;Lambda(1,1)=-2;
    real t=0,err;
    qdVec<real> x(Dim); x = (real)1.0;// qdVec<real> x(Dim)=(real) 1.0; does not work
    cout.precision(10);
//--------------------- main loop on time -------------------------
    for (unsigned int n=0;n<N;n=n+R)
        {
        t=n*Dt;
        x=fixe_point(x,t,Dt,Lambda,R,coef,Ainv);
        }
// --------------------- compute the error ----------------------------
    err=0.;
    x.disp();
    qdVec<real> sol=exact(Dim,T,Lambda);
    for(unsigned int i=0;i<Dim;i++)
        err=norminf(x-sol)/norminf(sol);
    cout<<"error="<< std::scientific<<err<<endl;
    fpu_fix_end(&oldcw);
    return 0;
    }

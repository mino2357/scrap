#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <quadmath.h>

using T = __float128;
using state_type = std::array<T,12>;

state_type operator+(const state_type& a,const state_type& b){state_type r;for(size_t i=0;i<a.size();++i)r[i]=a[i]+b[i];return r;}
state_type operator-(const state_type& a,const state_type& b){state_type r;for(size_t i=0;i<a.size();++i)r[i]=a[i]-b[i];return r;}
state_type operator*(const T& c,const state_type& a){state_type r;for(size_t i=0;i<a.size();++i)r[i]=c*a[i];return r;}

struct three_body_rhs{
    void operator()(const state_type& y,state_type& dydt,T)const{
        const T m1=3,m2=4,m3=5;
        T x1=y[0],y1=y[1],vx1=y[2],vy1=y[3];
        T x2=y[4],y2=y[5],vx2=y[6],vy2=y[7];
        T x3=y[8],y3=y[9],vx3=y[10],vy3=y[11];
        dydt[0]=vx1; dydt[1]=vy1; dydt[4]=vx2; dydt[5]=vy2; dydt[8]=vx3; dydt[9]=vy3;
        auto acc=[](T xi,T yi,T xj,T yj,T mj,T& axi,T& ayi){
            T dx=xj-xi,dy=yj-yi; T r2=dx*dx+dy*dy; T r3=r2*sqrtq(r2);
            axi+=mj*dx/r3; ayi+=mj*dy/r3;};
        T ax1=0,ay1=0,ax2=0,ay2=0,ax3=0,ay3=0;
        acc(x1,y1,x2,y2,m2,ax1,ay1); acc(x1,y1,x3,y3,m3,ax1,ay1);
        acc(x2,y2,x1,y1,m1,ax2,ay2); acc(x2,y2,x3,y3,m3,ax2,ay2);
        acc(x3,y3,x1,y1,m1,ax3,ay3); acc(x3,y3,x2,y2,m2,ax3,ay3);
        dydt[2]=ax1; dydt[3]=ay1; dydt[6]=ax2; dydt[7]=ay2; dydt[10]=ax3; dydt[11]=ay3;
    }
};

static std::string to_stringq(T x){char buf[128];quadmath_snprintf(buf,sizeof(buf),"%.40Qe",x);return buf;}

struct recorder{
    std::ofstream& ofs;
    void operator()(const state_type& y,T t)const{
        ofs<<to_stringq(t)<<' '<<to_stringq(y[0])<<' '<<to_stringq(y[1])
           <<' '<<to_stringq(y[4])<<' '<<to_stringq(y[5])
           <<' '<<to_stringq(y[8])<<' '<<to_stringq(y[9])<<'\n';}
};

void modified_midpoint(const three_body_rhs& rhs,const state_type& y,const state_type& dydx,T t,T htot,int nstep,state_type& yout){
    state_type ym=y,yn=y;T h=htot/nstep;for(size_t i=0;i<y.size();++i)yn[i]=y[i]+h*dydx[i];
    state_type dydxt;T x=t+h;for(int n=1;n<nstep;n++){rhs(yn,dydxt,x);for(size_t i=0;i<y.size();++i){T temp=ym[i]+2*h*dydxt[i];ym[i]=yn[i];yn[i]=temp;}x+=h;}
    rhs(yn,dydxt,t+htot);for(size_t i=0;i<y.size();++i)yout[i]=0.5*(ym[i]+yn[i]+h*dydxt[i]);}

bool bulirsch_stoer_step(const three_body_rhs& rhs,state_type& y,T& t,T& h,T tol){
    const int KMAX=8; const int nseq[KMAX]={2,4,6,8,10,12,14,16}; std::array<state_type,KMAX> table; T x[KMAX];
    state_type dydx;rhs(y,dydx,t);state_type yscale;for(size_t i=0;i<y.size();++i)yscale[i]=fabsq(y[i])+fabsq(h*dydx[i])+static_cast<T>(1e-30);
    T htemp=h;const T safety=0.9;while(true){for(int k=0;k<KMAX;k++){int n=nseq[k];state_type ytemp;modified_midpoint(rhs,y,dydx,t,htemp,n,ytemp);T xx=htemp/n;x[k]=xx*xx;table[k]=ytemp;if(k==0)continue;for(int j=k-1;j>=0;--j){T factor=x[k]/x[j]-1;for(size_t i=0;i<y.size();++i)table[j][i]=table[j+1][i]+(table[j+1][i]-table[j][i])/factor;}state_type yerr;for(size_t i=0;i<y.size();++i)yerr[i]=table[0][i]-table[1][i];T errmax=0;for(size_t i=0;i<y.size();++i){T err=fabsq(yerr[i]/yscale[i]);if(err>errmax)errmax=err;}if(errmax<tol){y=table[0];t+=htemp;double fac=pow((double)(tol/errmax),1.0/(2*k+1));if(fac>5)fac=5;h=htemp*T(safety*fac);return true;}}htemp*=0.5;if(htemp<static_cast<T>(1e-12)){y=table[0];t+=htemp;h=htemp;return true;}}
}

void integrate(const three_body_rhs& rhs,state_type& y,T t0,T t1,T h,T tol,recorder& obs){T t=t0;obs(y,t);while(t<t1){if(t+h>t1)h=t1-t;bulirsch_stoer_step(rhs,y,t,h,tol);obs(y,t);}}

int main(){
    state_type y={T(1),T(3),T(0),T(0),T(-2),T(-1),T(0),T(0),T(1),T(-1),T(0),T(0)};
    T t_end=15;std::ofstream ofs("case05.dat");recorder obs{ofs};three_body_rhs rhs;integrate(rhs,y,0,t_end,static_cast<T>(1e-3),static_cast<T>(1e-12),obs);
    char buf[256];quadmath_snprintf(buf,sizeof(buf),"%.15Qe %.15Qe %.15Qe %.15Qe %.15Qe %.15Qe\n",y[0],y[1],y[4],y[5],y[8],y[9]);std::cout<<buf;return 0;}

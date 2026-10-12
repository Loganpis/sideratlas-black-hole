// Original Kerr null-ray integrator, in M=G=c=1 units.
// Potentials: Gralla & Lupsasca 2020, PRD101 044032, equations (7),(8),(18).
// The independent coordinate is s=B*tau, B=sqrt(alpha^2+beta^2+1),
// where tau increases BACKWARDS from the observer at infinity.
// Integrating the second-order potential equations passes turning points
// continuously; no ad-hoc sign flips or square-root reflection is required.
#include <cmath>
#include <cstdint>
#include <algorithm>
#include <thread>
#include <vector>

struct State { double u, v, z, w; };
struct Constants { double A, K, a2q, invB2; };
static inline State derivative(const State& y,const Constants& c,double a2) {
    return {y.v,(c.A*y.u+3*c.K*y.u*y.u-2*c.a2q*y.u*y.u*y.u)*c.invB2,
            y.w,(c.A*y.z-2*a2*y.z*y.z*y.z)*c.invB2};
}
static inline State plus(const State& y,const State& k,double h) {
    return {y.u+h*k.u,y.v+h*k.v,y.z+h*k.z,y.w+h*k.w};
}
static inline State step(const State& y,double h,const Constants& c,double a2) {
    State k1=derivative(y,c,a2), k2=derivative(plus(y,k1,h*.5),c,a2),
          k3=derivative(plus(y,k2,h*.5),c,a2), k4=derivative(plus(y,k3,h),c,a2);
    return {y.u+h/6*(k1.u+2*k2.u+2*k3.u+k4.u),
            y.v+h/6*(k1.v+2*k2.v+2*k3.v+k4.v),
            y.z+h/6*(k1.z+2*k2.z+2*k3.z+k4.z),
            y.w+h/6*(k1.w+2*k2.w+2*k3.w+k4.w)};
}
static inline double hermite(double q0,double p0,double q1,double p1,double h,double t) {
    double t2=t*t,t3=t2*t;
    return (2*t3-3*t2+1)*q0+(t3-2*t2+t)*h*p0+(-2*t3+3*t2)*q1+(t3-t2)*h*p1;
}
static inline double hermite_derivative(double q0,double p0,double q1,double p1,double h,double t) {
    double t2=t*t;
    return ((6*t2-6*t)*q0+(3*t2-4*t+1)*h*p0+(-6*t2+6*t)*q1+(3*t2-2*t)*h*p1)/h;
}

extern "C" void trace_kerr(const double* alpha,const double* beta,long long n,
                           double a,double inc,double rin,double rout,double h,
                           int max_steps,int nthreads,float* radius,float* shift,
                           float* emission_mu,int16_t* layer,uint8_t* status,
                           double* residual) {
    const double si=std::sin(inc),ci=std::cos(inc),a2=a*a;
    const double uh=1/(1+std::sqrt(1-a2));
    auto worker=[&](long long begin,long long end) {
      for(long long pixel=begin;pixel<end;++pixel) {
        double al=alpha[pixel],be=beta[pixel],lam=-al*si;
        double eta=be*be+(al*al-a2)*ci*ci;
        double B=std::sqrt(al*al+be*be+1),invB=1/B;
        Constants c{a2-lam*lam-eta,eta+(lam-a)*(lam-a),a2*eta,invB*invB};
        State y{0,invB,ci,be*si*invB};
        int crossings=0;
        radius[pixel]=0;shift[pixel]=1;emission_mu[pixel]=0;
        layer[pixel]=-1;status[pixel]=3;residual[pixel]=0;
        for(int k=0;k<max_steps;++k) {
          State next=step(y,h,c,a2);
          if(y.z*next.z<0 && eta>=0) {
            double lo=0,hi=1;
            // Dense cubic Hermite interpolation retains fourth-order global
            // radius accuracy while locating the equator within each step.
            for(int j=0;j<32;++j) {
              double mid=(lo+hi)*.5;
              double zm=hermite(y.z,y.w,next.z,next.w,h,mid);
              if((y.z>0)==(zm>0)) lo=mid; else hi=mid;
            }
            double t=(lo+hi)*.5;
            double u=hermite(y.u,y.v,next.u,next.v,h,t);
            double r=1/u;
            if(u>0 && r>=rin && r<=rout) {
              double sr=std::sqrt(r),r32=r*sr;
              double omega=1/(r32+a);
              double ut=(r32+a)/std::sqrt(r*r*r-3*r*r+2*a*r32);
              double g=1/(ut*(1-omega*lam));
              radius[pixel]=r;shift[pixel]=g;
              emission_mu[pixel]=g*std::sqrt(eta)/r;
              layer[pixel]=crossings;status[pixel]=1;
              double v=hermite_derivative(y.u,y.v,next.u,next.v,h,t)*B;
              double w=hermite_derivative(y.z,y.w,next.z,next.w,h,t)*B;
              double U=1+c.A*u*u+2*c.K*u*u*u-c.a2q*u*u*u*u;
              residual[pixel]=std::max(std::abs(v*v-U)/(1+std::abs(U)),std::abs(w*w-eta)/(1+std::abs(eta)));
              break;
            }
            ++crossings;
          }
          if(next.u>=uh) {status[pixel]=2;break;}
          if(next.u<=0 && k>0) {status[pixel]=0;break;}
          if(!std::isfinite(next.u) || !std::isfinite(next.z)) {status[pixel]=4;break;}
          y=next;
        }
      }
    };
    nthreads=std::max(1,std::min(nthreads,int(n)));
    std::vector<std::thread> threads;
    for(int i=0;i<nthreads;++i) threads.emplace_back(worker,n*i/nthreads,n*(i+1)/nthreads);
    for(auto& t:threads) t.join();
}

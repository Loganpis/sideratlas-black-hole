// Original Kerr-Newman null-ray tracer with azimuth and relative coordinate delays.
// Equations: Wang, Lee & Lin (2022), PRD 106 084048, (9)-(15),(21).
// s=B*tau increases BACKWARD from a static observer at infinity.
// A pole-regular angular unit vector N removes the coordinate singularity
// of lambda/sin(theta)^2. The radial frame-dragging azimuth is integrated
// separately. No photon is reflected at a polar or radial turning point.
#include <cmath>
#include <cstdint>
#include <algorithm>
#include <thread>
#include <vector>

struct State { double u,v,x,y,z,px,py,pz,q,phi; };
struct Constants { double a,a2,charge2,lam,eta,A,K,C4,invB,invB2; };

static inline State derivative(const State& y,const Constants& c) {
    const double u2=y.u*y.u,D=1-2*y.u+(c.a2+c.charge2)*u2;
    const double P=1+(c.a2-c.a*c.lam)*u2,H=1+c.a2*u2;
    const double vt=y.v/c.invB;
    // q = elapsed backward coordinate time + r*(r) - r*(observer).
    // P^2-vt^2 = D*K*u^2 regularizes the incoming infinity cancellation.
    const double qr = vt>=0 ? H*c.K/(P+vt) : H*(P-vt)/(u2*D);
    const double qd=(qr+c.a*c.lam-c.a2*(1-y.z*y.z))*c.invB;
    const double az=-c.a*y.u*(2-(c.charge2+c.a*c.lam)*y.u)/D*c.invB;
    const double angular=-(c.eta+c.lam*c.lam+2*c.a2*y.z*y.z)*c.invB2;
    return {y.v,(c.A*y.u+3*c.K*u2-2*c.C4*u2*y.u)*c.invB2,
            y.px,y.py,y.pz,angular*y.x,angular*y.y,
            angular*y.z+c.a2*y.z*c.invB2,qd,az};
}
static inline State plus(const State& y,const State& k,double h) {
    return {y.u+h*k.u,y.v+h*k.v,y.x+h*k.x,y.y+h*k.y,y.z+h*k.z,
            y.px+h*k.px,y.py+h*k.py,y.pz+h*k.pz,y.q+h*k.q,y.phi+h*k.phi};
}
static inline State step(const State& y,double h,const Constants& c) {
    const State k1=derivative(y,c),k2=derivative(plus(y,k1,h*.5),c),
                k3=derivative(plus(y,k2,h*.5),c),k4=derivative(plus(y,k3,h),c);
    const double f=h/6;
    return {y.u+f*(k1.u+2*k2.u+2*k3.u+k4.u),y.v+f*(k1.v+2*k2.v+2*k3.v+k4.v),
            y.x+f*(k1.x+2*k2.x+2*k3.x+k4.x),y.y+f*(k1.y+2*k2.y+2*k3.y+k4.y),
            y.z+f*(k1.z+2*k2.z+2*k3.z+k4.z),
            y.px+f*(k1.px+2*k2.px+2*k3.px+k4.px),
            y.py+f*(k1.py+2*k2.py+2*k3.py+k4.py),
            y.pz+f*(k1.pz+2*k2.pz+2*k3.pz+k4.pz),
            y.q+f*(k1.q+2*k2.q+2*k3.q+k4.q),
            y.phi+f*(k1.phi+2*k2.phi+2*k3.phi+k4.phi)};
}
static inline double hermite(double q0,double p0,double q1,double p1,double h,double t) {
    const double t2=t*t,t3=t2*t;
    return (2*t3-3*t2+1)*q0+(t3-2*t2+t)*h*p0+(-2*t3+3*t2)*q1+(t3-t2)*h*p1;
}
static inline double tortoise(double r,double rp,double rm,double charge2) {
    return r+(2*rp-charge2)/(rp-rm)*std::log((r-rp)/2)-(2*rm-charge2)/(rp-rm)*std::log((r-rm)/2);
}

// The integrator and regularized delay above are exactly the established
// Kerr-Newman equations. This variant records every emitting thin plunging
// layer, stopping at the first opaque outer-disk intersection. Layer count is
// a storage/convergence limit only; there is no preferred physical image order.
extern "C" void trace_kn_accretion(const double* alpha,const double* beta,long long n,
                           double a,double charge_q,double si,double ci,double rin,double rout,
                           double rmin,double h,double max_parameter,
                           int max_steps,int nthreads,int max_hits,
                           float* radius,float* azimuth,float* delay,float* radial_momentum,
                           int16_t* layer,uint8_t* kind,uint8_t* status,uint8_t* overflow,
                           double* residual) {
    const double a2=a*a,charge2=charge_q*charge_q,rp=1+std::sqrt(1-a2-charge2),rm=1-std::sqrt(1-a2-charge2),uh=1/rp;
    const double umin=1/rmin;
    auto worker=[&](long long begin,long long end) {
      for(long long pixel=begin;pixel<end;++pixel) {
        const double al=alpha[pixel],be=beta[pixel],lam=-al*si;
        const double eta=be*be+(al*al-a2)*ci*ci;
        const double B=std::sqrt(al*al+be*be+1),invB=1/B;
        const double K=eta+(lam-a)*(lam-a);
        const Constants c{a,a2,charge2,lam,eta,a2-lam*lam-eta,K,
                          a2*eta+charge2*K,invB,invB*invB};
        State y{0,invB,0,-si,ci,al*invB,be*ci*invB,be*si*invB,0,0};
        const bool coplanar=(ci==0 && be==0);
        int crossings=0,hits=0;
        double parameter=0;
        status[pixel]=3;overflow[pixel]=0;residual[pixel]=0;
        for(int slot=0;slot<max_hits;++slot) {
            const long long ix=pixel*max_hits+slot;
            radius[ix]=azimuth[ix]=delay[ix]=radial_momentum[ix]=0;
            layer[ix]=-1;kind[ix]=0;
        }
        for(int k=0;k<max_steps && parameter<max_parameter;++k) {
          // BL t and phi diverge at the horizon; all emitted intersections
          // are outside rmin. Limit approach steps to stay outside the metric
          // singularity, then classify a ray as captured at the explicit
          // near-horizon emissivity cutoff. No horizon surface emits light.
          double hs=std::min(h,max_parameter-parameter);
          if(y.v>0 && y.u>0) hs=std::min(hs,.15*(uh-y.u)/y.v);
          const State next=step(y,hs,c);
          const bool crossing=(y.z*next.z<0 && eta>=0);
          const bool rim=coplanar && y.u<1/rout && next.u>=1/rout;
          if(crossing || rim) {
            double lo=0,hi=1;
            for(int j=0;j<30;++j) {
              const double mid=(lo+hi)*.5;
              const double f=rim ? hermite(y.u,y.v,next.u,next.v,hs,mid)-1/rout :
                                    hermite(y.z,y.pz,next.z,next.pz,hs,mid);
              if(rim ? f<0 : ((y.z>0)==(f>0))) lo=mid; else hi=mid;
            }
            const State hit=step(y,hs*(lo+hi)*.5,c);
            const double r=rim ? rout : 1/hit.u;
            if(hit.u>0 && r>=rmin && r<=rout) {
              if(hits>=max_hits) {overflow[pixel]=1;status[pixel]=5;break;}
              const long long ix=pixel*max_hits+hits;
              const double u=hit.u,v=hit.v*B,w=hit.pz*B;
              const double U=1+c.A*u*u+2*c.K*u*u*u-c.C4*u*u*u*u;
              const double V=eta+c.A*hit.z*hit.z-a2*std::pow(hit.z,4);
              const double D=1-2*u+(a2+charge2)*u*u;
              radius[ix]=r;
              azimuth[ix]=std::remainder(std::atan2(hit.y,hit.x)+hit.phi,2*M_PI);
              delay[ix]=hit.q-tortoise(r,rp,rm,charge2);
              // Future-directed covariant k_r/E = sign(backward du/ds)
              // sqrt(R)/Delta = sign(v)*sqrt(U)/D. Conserved potential
              // supplies the on-shell magnitude at the independently located
              // event; diagnostics retain the RK4 constraint residual.
              radial_momentum[ix]=std::copysign(std::sqrt(std::max(0.,U)),v)/D;
              layer[ix]=crossings;kind[ix]=(r>=rin ? 1 : 2);
              residual[pixel]=std::max(residual[pixel],std::max({std::abs(v*v-U)/(1+std::abs(U)),
                  std::abs(w*w-V)/(1+std::abs(V)),
                  std::abs(hit.x*hit.x+hit.y*hit.y+hit.z*hit.z-1)}));
              ++hits;
              if(r>=rin) {status[pixel]=1;break;}
            }
            ++crossings;
          }
          if(next.u>=umin) {status[pixel]=2;break;}
          if(next.u<=0 && k>0) {status[pixel]=0;break;}
          if(!std::isfinite(next.u) || !std::isfinite(next.z) || !std::isfinite(next.q)) {status[pixel]=4;break;}
          y=next;parameter+=hs;
        }
      }
    };
    nthreads=std::max(1,std::min(nthreads,int(n)));
    std::vector<std::thread> threads;
    for(int i=0;i<nthreads;++i) threads.emplace_back(worker,n*i/nthreads,n*(i+1)/nthreads);
    for(auto& t:threads) t.join();
}

"""KAM analysis of the QZS-ADV oscillator in its conservative limit (zeta1 = zeta2 = 0).

Hamiltonian (primary X, absorber X2, mass ratio mu = 0.1, beta = 0.35):
  H = p^2/2 + p2^2/(2 mu) + U(X) + mu beta^2 (X - X2)^2 / 2 - f X cos t
  U(X) = 1.55 X^2 - 3 sqrt((1.5 eta)^2 + X^2)
(these reproduce the equations used in all maps: X'' = f cos t - 3.1 X + 3X/r - 0.01225 Z,
 Z = X - X2, Z'' = X'' - 0.1225 Z, when zeta1 = zeta2 = 0).
Reduced 1.5-DOF model: absorber frozen (Z = W = 0), H1 = p^2/2 + U(X) - f X cos t.
Integrable skeleton (f = 0): frequency omega(E), twist domega/dE, resonances omega/Omega = p/q.

Stages (cached in outputs/qzs_kam/cache):
  omega      unperturbed omega(E) by two methods (ODE event, weighted quadrature)
  sections   stroboscopic sections (t = 2 pi n), 64 orbits x 1500 periods, 4th-order
             symplectic Yoshida (200 steps/period), with FLI for every orbit
  fli        FLI maps over initial conditions (x0, v0), 240 x 240, 400 periods
  fraction   regular fraction vs f (120 x 120, 400 periods, H0 - Umin <= 0.5)
  absorber   2.5-DOF (with absorber) FLI maps, shallow wells, Z = W = 0 at t = 0
  damping    fate of the tori with primary damping zeta1 (RK4, 180 steps/period)
  checks     energy conservation, dt/2 and DOP853 comparisons, analytic limits
"""
import json,sys,time,os,math
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
R=Path(__file__).resolve().parents[1];OUT=R/'outputs/qzs_kam';CACHE=OUT/'cache';CACHE.mkdir(parents=True,exist_ok=True)
T=2*np.pi;MU=.1;BETA2=.1225
CASES={'monostable':dict(eta=2/3.1,xw=1.4,vw=1.4,x0max=1.25),
       'shallow_wells':dict(eta=.6,xw=1.4,vw=1.4,x0max=1.30),
       'deep_wells':dict(eta=.3,xw=1.9,vw=1.5,x0max=1.75)}
F_SECTIONS=[.002,.01,.05,.15];F_FRACTION=[.001,.002,.005,.01,.02,.05,.1,.15,.2,.3]
ECUT=.5;W1=1/(2-2**(1/3));W0=-2**(1/3)/(2-2**(1/3))
CY=[W1/2,(W0+W1)/2,(W0+W1)/2,W1/2];DY=[W1,W0,W1]

def U(x,eta):return 1.55*x*x-3*np.sqrt((1.5*eta)**2+x*x)
def dU(x,eta):return 3.1*x-3*x/np.sqrt((1.5*eta)**2+x*x)
def d2U(x,eta):a2=(1.5*eta)**2;return 3.1-3*a2/(a2+x*x)**1.5
def xmin(eta):return 1.5*np.sqrt(max((2/3.1)**2-eta**2,0))

# ------------------------------------------------------------------ omega(E)
def omega_curve(name):
 """Unperturbed frequency by two independent methods.
 The outer turning point b fixes the energy exactly, E = U(b). With A = 1.5 eta,
 r(x) = sqrt(A^2 + x^2) and r* = 6/3.1 - r(b), the energy gap factors as
   U(b) - U(x) = 1.55 (b^2 - x^2) (r(x) - r*) / (r(b) + r(x)),
 so no cancellation occurs near the turning points. If r* > A the orbit stays in
 a well with inner turning point a = sqrt(r*^2 - A^2) and
   r(x) - r* = (x^2 - a^2) / (r(x) + r*);
 otherwise it passes over the barrier and is symmetric (turning points -b, b).
 Method 1: weighted Gauss quadrature of T = oint dx / sqrt(2 (E - U)).
 Method 2: DOP853 integration from (b, 0) to the next turning point (event).
 """
 from scipy.integrate import solve_ivp,quad
 from scipy.optimize import brentq
 eta=CASES[name]['eta'];A=1.5*eta;xm=xmin(eta);Um=U(xm,eta);Ub=U(0.,eta);depth=Ub-Um
 r=lambda x:np.sqrt(A*A+x*x);C6=6/3.1
 bsep=np.sqrt(max((C6-A)**2-A*A,0.))  # outer turning point of the separatrix
 def T_ode(b):
  ev=lambda t,y:y[1];ev.direction=1   # v crosses zero upwards at the next turning point
  s=solve_ivp(lambda t,y:[y[1],-dU(y[0],eta)],[0,20000],[b,0.],method='DOP853',rtol=1e-12,atol=1e-14,events=ev)
  te=s.t_events[0];te=te[te>1e-9];return 2*te[0] if len(te) else np.nan
 def T_quad(b):
  rb=r(b);rs=C6-rb
  if rs>A:  # inside a well: turning points a < x < b
   a=np.sqrt(rs*rs-A*A)
   g=lambda x:1/np.sqrt(2*1.55*(b+x)*(x+a)/((rb+r(x))*(r(x)+rs)))
   val,_=quad(g,a,b,weight='alg',wvar=(-.5,-.5),limit=800);return 2*val,a,0
  g=lambda x:1/np.sqrt(2*1.55*(b+x)*(r(x)-rs)/(rb+r(x)))
  val,_=quad(g,0.,b,weight='alg',wvar=(0.,-.5),limit=800);return 4*val,-b,1
 targets=[]
 if depth>0:targets+=[(xm,bsep,depth*q) for q in np.r_[1e-6,1e-4,1e-2,np.linspace(.03,.97,30),1-np.logspace(-2,-6,25)]]
 targets+=[(max(bsep,1e-12),6.,depth+q) for q in np.logspace(-6,np.log10(1.2),90)]
 rows=[]
 for lo,hi,dE in targets:
  b=brentq(lambda x:U(x,eta)-Um-dE,lo,hi,xtol=1e-15)
  Tq,a,over=T_quad(b);To=T_ode(b);E=U(b,eta)
  rows.append((E-Um,E,a,b,2*np.pi/To,2*np.pi/Tq,over))
 arr=np.array(rows)
 # analytic limits
 K=d2U(xm,eta);lim=dict(omega0_sqrtK=float(np.sqrt(max(K,0))))
 if depth==0:
  a=1.5*eta;c=3/(8*a**3);I=math.gamma(.25)**2/(4*math.sqrt(2*math.pi))
  lim['quartic_c']=c;lim['quartic_omega_over_E14']=float(2*np.pi/(4*I*c**-.25/np.sqrt(2)))
 np.savez_compressed(CACHE/f'omega_{name}.npz',dE=arr[:,0],E=arr[:,1],a=arr[:,2],b=arr[:,3],omega_ode=arr[:,4],
  omega_quad=arr[:,5],over=arr[:,6],Umin=Um,Ubar=Ub,depth=depth,xmin=xm,**{k:v for k,v in lim.items()})
 return name,'omega',float(np.nanmax(np.abs(arr[:,4]-arr[:,5])/arr[:,5]))

# ------------------------------------------------------------------ symplectic integrators
def yoshida_1dof(x,v,eta,f,periods,ns,tangent=True,record=True,tag=''):
 """Stroboscopic samples and FLI (log10 growth, renormalized each period)."""
 n=x.size;h=T/ns;t=0.;dx=np.ones(n)/np.sqrt(2);dv=dx.copy();logg=np.zeros(n);t0=time.time()
 pts=np.empty((periods,2,n),np.float32) if record else None
 for c in range(periods):
  for s in range(ns):
   for k in range(3):
    x=x+CY[k]*h*v;t+=CY[k]*h
    if tangent:dx=dx+CY[k]*h*dv
    v=v+DY[k]*h*(f*np.cos(t)-dU(x,eta))
    if tangent:dv=dv-DY[k]*h*d2U(x,eta)*dx
   x=x+CY[3]*h*v;t+=CY[3]*h
   if tangent:dx=dx+CY[3]*h*dv
  t=(c+1)*T  # remove round-off drift of the clock
  if tangent:
   nrm=np.sqrt(dx*dx+dv*dv);logg+=np.log10(nrm);dx/=nrm;dv/=nrm
  if record:pts[c,0]=x;pts[c,1]=v
  if tag and c%200==199:print(f'  {tag}: {c+1}/{periods} {time.time()-t0:.0f}s',flush=True)
 return x,v,logg,pts

def yoshida_2dof(x,v,x2,v2,eta,f,periods,ns,tag=''):
 """Conservative 2.5-DOF system (p2 = mu v2), FLI of the 4D tangent map."""
 n=x.size;h=T/ns;t=0.;mb=MU*BETA2
 d=[np.ones(n)/2 for _ in range(4)];logg=np.zeros(n);t0=time.time()
 for c in range(periods):
  for s in range(ns):
   for k in range(4):
    x=x+CY[k]*h*v;x2=x2+CY[k]*h*v2;t+=CY[k]*h
    d[0]=d[0]+CY[k]*h*d[1];d[2]=d[2]+CY[k]*h*d[3]
    if k==3:break
    coup=mb*(x-x2)
    v=v+DY[k]*h*(f*np.cos(t)-dU(x,eta)-coup);v2=v2+DY[k]*h*(BETA2*(x-x2))
    k2=d2U(x,eta)
    d[1]=d[1]+DY[k]*h*(-(k2+mb)*d[0]+mb*d[2]);d[3]=d[3]+DY[k]*h*(BETA2*(d[0]-d[2]))
  t=(c+1)*T
  nrm=np.sqrt(d[0]**2+d[1]**2+d[2]**2+d[3]**2);logg+=np.log10(nrm);d=[q/nrm for q in d]
  if tag and c%100==99:print(f'  {tag}: {c+1}/{periods} {time.time()-t0:.0f}s',flush=True)
 return logg

# ------------------------------------------------------------------ tasks
def ics_line(name):
 m=CASES[name];xs=np.linspace(.03,m['x0max'],32);return np.r_[-xs[::-1],xs],np.zeros(64)

def task_section(name,f):
 eta=CASES[name]['eta'];x0,v0=ics_line(name)
 _,_,fli,pts=yoshida_1dof(x0.copy(),v0.copy(),eta,f,1500,200,tag=f'sec {name} f={f}')
 E0=.5*v0**2+U(x0,eta)-U(xmin(eta),eta)
 np.savez_compressed(CACHE/f'section_{name}_{f}.npz',x0=x0,v0=v0,E0=E0,fli=fli,pts=pts,f=f,eta=eta)
 return name,f'section f={f}',float(np.mean(fli>10))

def grid(name,nx):
 m=CASES[name];xs=np.linspace(-m['xw'],m['xw'],nx);vs=np.linspace(-m['vw'],m['vw'],nx);X,V=np.meshgrid(xs,vs);return xs,vs,X.ravel(),V.ravel()

def task_fli(name,f,nx,periods,ns,prefix):
 eta=CASES[name]['eta'];xs,vs,X,V=grid(name,nx)
 _,_,fli,_=yoshida_1dof(X.copy(),V.copy(),eta,f,periods,ns,record=False,tag=f'{prefix} {name} f={f}')
 H0=.5*V**2+U(X,eta)-U(xmin(eta),eta)
 np.savez_compressed(CACHE/f'{prefix}_{name}_{f}.npz',xs=xs,vs=vs,fli=fli.reshape(nx,nx),H0=H0.reshape(nx,nx),f=f,periods=periods)
 inside=H0<=ECUT;return name,f'{prefix} f={f}',float(np.mean(fli[inside]<=8))

def task_absorber(name,f,nx,periods,ns):
 eta=CASES[name]['eta'];xs,vs,X,V=grid(name,nx)
 fli=yoshida_2dof(X.copy(),V.copy(),X.copy(),V.copy(),eta,f,periods,ns,tag=f'abs {name} f={f}')
 H0=.5*V**2+U(X,eta)-U(xmin(eta),eta)
 np.savez_compressed(CACHE/f'absorber_{name}_{f}.npz',xs=xs,vs=vs,fli=fli.reshape(nx,nx),H0=H0.reshape(nx,nx),f=f,periods=periods)
 return name,f'absorber f={f}',float(np.mean(fli[H0<=ECUT]<=8))

def task_damping(name,f,zeta1):
 eta=CASES[name]['eta'];x,v=ics_line(name);ns=180;h=T/ns;cs=np.cos(np.arange(2*ns+1)*np.pi/ns)
 periods=3000;pts=np.empty((periods,2,x.size),np.float32);t0=time.time()
 acc=lambda x,v,force:force-2*zeta1*v-dU(x,eta)
 for c in range(periods):
  for s in range(ns):
   k1x=v;k1v=acc(x,v,f*cs[2*s]);k2x=v+.5*h*k1v;k2v=acc(x+.5*h*k1x,k2x,f*cs[2*s+1])
   k3x=v+.5*h*k2v;k3v=acc(x+.5*h*k2x,k3x,f*cs[2*s+1]);k4x=v+h*k3v;k4v=acc(x+h*k3x,k4x,f*cs[2*s+2])
   x=x+h/6*(k1x+2*k2x+2*k3x+k4x);v=v+h/6*(k1v+2*k2v+2*k3v+k4v)
  pts[c,0]=x;pts[c,1]=v
 E0=U(ics_line(name)[0],eta)-U(xmin(eta),eta)
 np.savez_compressed(CACHE/f'damping_{name}_{f}_{zeta1}.npz',pts=pts,E0=E0,f=f,zeta1=zeta1)
 return name,f'damping f={f} zeta1={zeta1}',time.time()-t0

def task_checks(name):
 from scipy.integrate import solve_ivp
 eta=CASES[name]['eta'];m=CASES[name];out={}
 # 1) energy conservation, f = 0, 2000 periods
 x0,v0=ics_line(name);x=x0.copy();v=v0.copy();H0=.5*v*v+U(x,eta)
 x,v,_,_=yoshida_1dof(x,v,eta,0.,2000,200,tangent=False,record=False)
 out['energy_rel_error_f0_2000T']=float(np.max(np.abs(.5*v*v+U(x,eta)-H0)/np.maximum(np.abs(H0),1e-12)))
 # 2) regular orbits at f = 0.002: dt vs dt/2 vs DOP853 over 300 periods (sections)
 f=.002;sel=np.array([8,40,60])  # interior orbits on both sides
 xa,va,_,pa=yoshida_1dof(x0[sel].copy(),v0[sel].copy(),eta,f,300,200)
 xb,vb,_,pb=yoshida_1dof(x0[sel].copy(),v0[sel].copy(),eta,f,300,400)
 out['section_diff_dt_vs_dt2_300T']=float(np.max(np.abs(pa-pb)))
 dd=[]
 for i,j in enumerate(sel):
  s=solve_ivp(lambda t,y:[y[1],f*np.cos(t)-dU(y[0],eta)],[0,300*T],[x0[j],v0[j]],method='DOP853',rtol=1e-12,atol=1e-14,t_eval=np.arange(1,301)*T)
  dd.append(np.max(np.abs(s.y.T-pa[:,:,i].astype(float))))
 out['section_diff_yoshida_vs_dop853_300T']=float(max(dd))
 # 3) 2.5-DOF energy conservation, f = 0
 X,V=x0[::8].copy(),v0[::8].copy();X2,V2=X.copy()+.05,V.copy()
 H=lambda x,v,x2,v2:.5*v*v+.5*MU*v2*v2+U(x,eta)+.5*MU*BETA2*(x-x2)**2
 Hs=H(X,V,X2,V2);n=X.size;h=T/200;t=0.
 for c in range(500):
  for s in range(200):
   for k in range(4):
    X=X+CY[k]*h*V;X2=X2+CY[k]*h*V2
    if k==3:break
    V=V+DY[k]*h*(-dU(X,eta)-MU*BETA2*(X-X2));V2=V2+DY[k]*h*(BETA2*(X-X2))
 out['energy_rel_error_absorber_f0_500T']=float(np.max(np.abs(H(X,V,X2,V2)-Hs)/np.abs(Hs)))
 (CACHE/f'checks_{name}.json').write_text(json.dumps(out,indent=2))
 return name,'checks',out

def main():
 workers=int(os.environ.get('KAM_WORKERS','4'));t0=time.time();jobs=[];queued=set()
 def want(path,fn,*args):
  if not Path(path).exists() and str(path) not in queued:jobs.append((fn,args));queued.add(str(path))
 for nm in CASES:
  want(CACHE/f'omega_{nm}.npz',omega_curve,nm);want(CACHE/f'checks_{nm}.json',task_checks,nm)
  for f in F_SECTIONS:
   want(CACHE/f'section_{nm}_{f}.npz',task_section,nm,f)
   want(CACHE/f'fli_{nm}_{f}.npz',task_fli,nm,f,200,400,100,'fli')
  # the four section forcings reuse the 200 x 200 FLI maps; the others use 100 x 100
  for f in F_FRACTION:
   if f not in F_SECTIONS:want(CACHE/f'frac_{nm}_{f}.npz',task_fli,nm,f,100,400,100,'frac')
 for f in [.01,.05,.15]:
  want(CACHE/f'absorber_shallow_wells_{f}.npz',task_absorber,'shallow_wells',f,160,400,100)
 for nm in ['shallow_wells','deep_wells']:
  for z in [0.,1e-4,1e-3,1e-2]:want(CACHE/f'damping_{nm}_0.05_{z}.npz',task_damping,nm,.05,z)
 print(f'{len(jobs)} tasks to run, {workers} workers',flush=True)
 with ProcessPoolExecutor(max_workers=workers) as pool:
  fut={pool.submit(fn,*args):(fn.__name__,args) for fn,args in jobs}
  for k,fu in enumerate(as_completed(fut),1):
   r=fu.result();print(f'done {k}/{len(jobs)} {r[0]} {r[1]}: {r[2]}  [{time.time()-t0:.0f}s]',flush=True)
 checks={nm:json.loads((CACHE/f'checks_{nm}.json').read_text()) for nm in CASES}
 for nm in CASES:
  d=np.load(CACHE/f'omega_{nm}.npz');checks[nm]['omega_ode_vs_quad_max_rel']=float(np.nanmax(np.abs(d['omega_ode']-d['omega_quad'])/d['omega_quad']))
  if float(d['depth'])>0:
   w=d['over']==0;i=np.argmin(np.where(w,d['dE'],np.inf))
   checks[nm]['omega_small_amplitude_vs_sqrtK']=[float(d['omega_quad'][i]),float(d['omega0_sqrtK']),float(d['dE'][i])]
  else:
   i=np.argmin(d['dE']);checks[nm]['quartic_omega_over_E14_numeric_vs_analytic']=[float(d['omega_quad'][i]/d['dE'][i]**.25),float(d['quartic_omega_over_E14'])]
 (OUT/'checks.json').write_text(json.dumps(checks,indent=2));print(json.dumps(checks,indent=1))
 print(f'ALL DONE {time.time()-t0:.0f}s',flush=True)

if __name__=='__main__':main()

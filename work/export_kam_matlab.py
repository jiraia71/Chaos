"""Export the KAM analysis data (work/kam_compute.py) to MATLAB."""
import json,hashlib
from pathlib import Path
import numpy as np
from scipy.io import savemat,loadmat
import kam_compute as K
import kam_render as KR
R=Path(__file__).resolve().parents[1];OUT=R/'outputs/matlab_qzs_kam';OUT.mkdir(exist_ok=True)
CASES=[('monostable','Monostable'),('shallow_wells','Shallow wells'),('deep_wells','Deep wells')]
D={};report={}
for nm,label in CASES:
 d=np.load(K.CACHE/f'omega_{nm}.npz');res=[r for r in KR.resonances(d) if r['dE']<=1.0 and r['q']<=3]
 m=K.CASES[nm]
 c=dict(name=nm,label=label,eta=m['eta'],xw=m['xw'],vw=m['vw'],Umin=float(d['Umin']),Ubar=float(d['Ubar']),
  depth=float(d['depth']),xmin=float(d['xmin']),omega0_sqrtK=float(d['omega0_sqrtK']),
  dE=d['dE'].reshape(-1,1),omega_quad=d['omega_quad'].reshape(-1,1),omega_ode=d['omega_ode'].reshape(-1,1),over=d['over'].reshape(-1,1),
  res_p=np.array([r['p'] for r in res],float).reshape(-1,1),res_q=np.array([r['q'] for r in res],float).reshape(-1,1),
  res_rho=np.array([r['rho'] for r in res]).reshape(-1,1),res_dE=np.array([r['dE'] for r in res]).reshape(-1,1),
  res_over=np.array([r['over'] for r in res],float).reshape(-1,1))
 sec=[np.load(K.CACHE/f'section_{nm}_{f}.npz') for f in K.F_SECTIONS]
 c.update(sec_f=np.array(K.F_SECTIONS).reshape(1,-1),sec_pts=np.stack([s['pts'] for s in sec]).astype(np.float32),
  sec_fli=np.stack([s['fli'] for s in sec]),sec_E0=sec[0]['E0'].reshape(-1,1))
 fli=[np.load(K.CACHE/f'fli_{nm}_{f}.npz') for f in K.F_SECTIONS]
 c.update(fli_f=np.array(K.F_SECTIONS).reshape(1,-1),fli_maps=np.stack([g['fli'] for g in fli]).astype(np.float32),
  fli_xs=fli[0]['xs'].reshape(1,-1),fli_vs=fli[0]['vs'].reshape(-1,1),fli_periods=400.)
 fr=[];lo=[];hi=[]
 for f in K.F_FRACTION:
  g=np.load(K.CACHE/(f'fli_{nm}_{f}.npz' if f in K.F_SECTIONS else f'frac_{nm}_{f}.npz'));v=g['fli'][g['H0']<=K.ECUT]
  fr.append(np.mean(v<=8));lo.append(np.mean(v<=6));hi.append(np.mean(v<=10))
 c.update(frac_f=np.array(K.F_FRACTION).reshape(1,-1),frac_regular=np.array(fr).reshape(1,-1),
  frac_fli6=np.array(lo).reshape(1,-1),frac_fli10=np.array(hi).reshape(1,-1))
 if nm=='shallow_wells':
  ab=[np.load(K.CACHE/f'absorber_shallow_wells_{f}.npz') for f in [.01,.05,.15]]
  c.update(abs_f=np.array([.01,.05,.15]).reshape(1,-1),abs_maps=np.stack([g['fli'] for g in ab]).astype(np.float32),
   abs_xs=ab[0]['xs'].reshape(1,-1),abs_vs=ab[0]['vs'].reshape(-1,1),
   abs_regular=np.array([float(np.mean(g['fli'][g['H0']<=K.ECUT]<=8)) for g in ab]).reshape(1,-1))
 if nm in ('shallow_wells','deep_wells'):
  zs=[0.,1e-4,1e-3,1e-2];dm=[np.load(K.CACHE/f'damping_{nm}_0.05_{z}.npz') for z in zs]
  c.update(damp_zeta1=np.array(zs).reshape(1,-1),damp_f=.05,damp_E0=dm[0]['E0'].reshape(-1,1),
   damp_early=np.stack([g['pts'][:150] for g in dm]).astype(np.float32),
   damp_late=np.stack([g['pts'][-300:] for g in dm]).astype(np.float32))
 D[nm]=c
 report[nm]=dict(resonances=len(res),chaotic_fraction_sections={str(f):float(np.mean(s['fli']>10)) for f,s in zip(K.F_SECTIONS,sec)},
  regular_fraction_fli={str(f):float(np.mean(g['fli'][g['H0']<=K.ECUT]<=8)) for f,g in zip(K.F_SECTIONS,fli)})
checks=json.loads((K.OUT/'checks.json').read_text())
mat=OUT/'dados_kam.mat'
savemat(mat,dict(D=D,ECUT=K.ECUT,MU=K.MU,BETA2=K.BETA2,fli_chaos_threshold=8.,section_chaos_threshold=10.),do_compression=True,oned_as='column')
back=loadmat(mat,squeeze_me=False)['D'][0,0]
for nm,_ in CASES:
 b=back[nm][0,0]
 for key in ['dE','omega_quad','sec_pts','fli_maps','frac_regular','res_dE']:
  assert np.array_equal(b[key],D[nm][key]),(nm,key)
(OUT/'export_checks.json').write_text(json.dumps(dict(mat=mat.name,roundtrip_exact=True,cases=report,
 numerical_checks=checks,source='outputs/qzs_kam/cache',
 sha256_kam_compute=hashlib.sha256((R/'work/kam_compute.py').read_bytes()).hexdigest()),indent=2))
print(mat,round(mat.stat().st_size/1e6,2),'MB; roundtrip exact')
print(json.dumps(report,indent=1))

"""Mapas KAM (FLI 2.5 GL conservativo) com frequência de EXCITAÇÃO Ω e sintonia β
arbitrárias. Generaliza os mapas anteriores (que fixavam Ω=1) para excitar o sistema
em outra frequência — aqui Ω=0.35 com β=0.35 (β=Ω, antirressonância para esse drive).

O drive passa a ser f·cos(Ω·t); o período estroboscópico é 2π/Ω e o passo dt é mantido
constante (ns = round(NS_BASE/Ω)) para preservar a precisão do integrador simplético.

Exporta dados_kam_omega{Ω}_beta{β}.mat com a struct B (mesmo formato de fig_kam_beta):
por caso, eta/xw/vw/Umin/depth/xs/vs/maps(1×nx×nx)/reg, além de B.omega, B.beta, B.f.
"""
import os
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
from scipy.io import savemat

import kam_compute as K

NAMES = ["monostable", "shallow_wells", "deep_wells"]
LABELS = ["Monostable", "Shallow wells", "Deep wells"]
OMEGA = float(os.environ.get("KAM_OMEGA", "0.35"))
BETA = float(os.environ.get("KAM_BETA", "0.35"))
FFORCE = float(os.environ.get("KAM_F", "0.05"))
NX = int(os.environ.get("KAM_NX", "200"))
NS_BASE = 100
SNAPS = (200, 400)   # períodos estroboscópicos do drive


def snap_map(nm):
    eta = K.CASES[nm]["eta"]; m = K.CASES[nm]; b2 = BETA ** 2; mb = K.MU * b2
    Tdrive = 2 * np.pi / OMEGA
    ns = max(NS_BASE, int(round(NS_BASE / OMEGA)))   # mantém dt ~ constante
    h = Tdrive / ns
    xs = np.linspace(-m["xw"], m["xw"], NX); vs = np.linspace(-m["vw"], m["vw"], NX)
    X, V = np.meshgrid(xs, vs); x = X.ravel(); v = V.ravel(); x2 = x.copy(); v2 = v.copy()
    n = x.size; t = 0.
    d = [np.ones(n) / 2 for _ in range(4)]; logg = np.zeros(n); snaps = {}
    for c in range(max(SNAPS)):
        for s in range(ns):
            for k in range(4):
                x = x + K.CY[k] * h * v; x2 = x2 + K.CY[k] * h * v2; t += K.CY[k] * h
                d[0] = d[0] + K.CY[k] * h * d[1]; d[2] = d[2] + K.CY[k] * h * d[3]
                if k == 3:
                    break
                coup = mb * (x - x2)
                v = v + K.DY[k] * h * (FFORCE * np.cos(OMEGA * t) - K.dU(x, eta) - coup)
                v2 = v2 + K.DY[k] * h * (b2 * (x - x2))
                k2 = K.d2U(x, eta)
                d[1] = d[1] + K.DY[k] * h * (-(k2 + mb) * d[0] + mb * d[2])
                d[3] = d[3] + K.DY[k] * h * (b2 * (d[0] - d[2]))
        t = (c + 1) * Tdrive
        nrm = np.sqrt(d[0]**2 + d[1]**2 + d[2]**2 + d[3]**2); logg += np.log10(nrm); d = [q / nrm for q in d]
        if (c + 1) in SNAPS:
            snaps[c + 1] = logg.copy()
    H0 = .5 * V.ravel()**2 + K.U(x, eta) - K.U(K.xmin(eta), eta); inside = H0 <= K.ECUT
    reg = float(np.mean((snaps[SNAPS[1]] - snaps[SNAPS[0]])[inside] <= 1.0))
    return nm, xs, vs, snaps[SNAPS[1]].reshape(NX, NX).astype(np.float32), reg


def main():
    res = {}
    with ProcessPoolExecutor(max_workers=int(os.environ.get("KAM_WORKERS", "4"))) as pool:
        futs = {pool.submit(snap_map, nm): nm for nm in NAMES}
        for k, fu in enumerate(as_completed(futs), 1):
            nm, xs, vs, mp, reg = fu.result(); res[nm] = (xs, vs, mp, reg)
            print(f"  {k}/{len(NAMES)} {nm} Ω={OMEGA} β={BETA}: regular={100*reg:.0f}%", flush=True)
    B = dict(omega=OMEGA, beta=BETA, f=FFORCE,
             names=np.array(NAMES, dtype=object).reshape(1, -1),
             labels=np.array(LABELS, dtype=object).reshape(1, -1))
    for nm, lab in zip(NAMES, LABELS):
        eta = K.CASES[nm]["eta"]; m = K.CASES[nm]; Um = K.U(K.xmin(eta), eta); depth = K.U(0., eta) - Um
        xs, vs, mp, reg = res[nm]
        B[nm] = dict(eta=eta, xw=m["xw"], vw=m["vw"], Umin=float(Um), depth=float(depth),
                     xs=xs.reshape(1, -1), vs=vs.reshape(1, -1),
                     maps=mp[None, :, :].astype(np.float32), reg=np.array([[reg]]), label=lab)
    tag = f"omega{OMEGA:g}_beta{BETA:g}".replace(".", "p")
    out = f"/home/user/Chaos/dados_kam_{tag}.mat"
    savemat(out, dict(B=B), do_compression=True, oned_as="column")
    print(f"salvo {out} ({os.path.getsize(out)/1e6:.2f} MB), nx={NX}, Ω={OMEGA}, β={BETA}")


if __name__ == "__main__":
    main()

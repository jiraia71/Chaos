"""Computa os mapas KAM (FLI 2.5 GL conservativo) para as sintonias β e exporta um .mat
para o MATLAB (fig_kam_beta.m) gerar a figura em alta definição.

Saída: dados_kam_beta.mat com a struct B:
  B.betas, B.f, B.names, B.labels e, por caso, um sub-struct com
  eta, xw, vw, Umin, depth, xs, vs, maps (nβ×nx×nx), reg (1×nβ).
"""
import os
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
from scipy.io import savemat

import kam_compute as K
CHAOS_ROOT = os.environ.get("CHAOS_ROOT") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORK = os.path.join(CHAOS_ROOT, "work")

NAMES = ["monostable", "shallow_wells", "deep_wells"]
LABELS = ["Monostable", "Shallow wells", "Deep wells"]
BETAS = [0.35, 0.70, 1.00]
FFORCE = 0.05
NX = int(os.environ.get("KAM_NX", "200"))
SNAPS = (200, 400)
NS = 100


def snap_map(nm, beta):
    eta = K.CASES[nm]["eta"]; m = K.CASES[nm]; b2 = beta ** 2; mb = K.MU * b2
    xs = np.linspace(-m["xw"], m["xw"], NX); vs = np.linspace(-m["vw"], m["vw"], NX)
    X, V = np.meshgrid(xs, vs); x = X.ravel(); v = V.ravel(); x2 = x.copy(); v2 = v.copy()
    n = x.size; h = K.T / NS; t = 0.
    d = [np.ones(n) / 2 for _ in range(4)]; logg = np.zeros(n); snaps = {}
    for c in range(max(SNAPS)):
        for s in range(NS):
            for k in range(4):
                x = x + K.CY[k] * h * v; x2 = x2 + K.CY[k] * h * v2; t += K.CY[k] * h
                d[0] = d[0] + K.CY[k] * h * d[1]; d[2] = d[2] + K.CY[k] * h * d[3]
                if k == 3:
                    break
                coup = mb * (x - x2)
                v = v + K.DY[k] * h * (FFORCE * np.cos(t) - K.dU(x, eta) - coup); v2 = v2 + K.DY[k] * h * (b2 * (x - x2))
                k2 = K.d2U(x, eta)
                d[1] = d[1] + K.DY[k] * h * (-(k2 + mb) * d[0] + mb * d[2]); d[3] = d[3] + K.DY[k] * h * (b2 * (d[0] - d[2]))
        t = (c + 1) * K.T
        nrm = np.sqrt(d[0]**2 + d[1]**2 + d[2]**2 + d[3]**2); logg += np.log10(nrm); d = [q / nrm for q in d]
        if (c + 1) in SNAPS:
            snaps[c + 1] = logg.copy()
    H0 = .5 * V.ravel()**2 + K.U(X.ravel(), eta) - K.U(K.xmin(eta), eta); inside = H0 <= K.ECUT  # posição INICIAL
    reg = float(np.mean((snaps[SNAPS[1]] - snaps[SNAPS[0]])[inside] <= 1.0))
    return nm, beta, xs, vs, snaps[SNAPS[1]].reshape(NX, NX).astype(np.float32), reg


def main():
    jobs = [(nm, b) for nm in NAMES for b in BETAS]
    res = {}
    with ProcessPoolExecutor(max_workers=int(os.environ.get("KAM_WORKERS", "4"))) as pool:
        futs = [pool.submit(snap_map, nm, b) for nm, b in jobs]
        for k, fu in enumerate(as_completed(futs), 1):
            nm, beta, xs, vs, mp, reg = fu.result()
            res[(nm, beta)] = (xs, vs, mp, reg)
            print(f"  {k}/{len(jobs)} {nm} β={beta}: regular={100*reg:.0f}%", flush=True)

    B = dict(betas=np.array(BETAS).reshape(1, -1), f=FFORCE,
             names=np.array(NAMES, dtype=object).reshape(1, -1),
             labels=np.array(LABELS, dtype=object).reshape(1, -1))
    for nm, lab in zip(NAMES, LABELS):
        eta = K.CASES[nm]["eta"]; m = K.CASES[nm]; Um = K.U(K.xmin(eta), eta); depth = K.U(0., eta) - Um
        xs, vs, _, _ = res[(nm, BETAS[0])]
        maps = np.stack([res[(nm, b)][2] for b in BETAS]).astype(np.float32)
        reg = np.array([res[(nm, b)][3] for b in BETAS]).reshape(1, -1)
        B[nm] = dict(eta=eta, xw=m["xw"], vw=m["vw"], Umin=float(Um), depth=float(depth),
                     xs=xs.reshape(1, -1), vs=vs.reshape(1, -1), maps=maps, reg=reg, label=lab)
    out = f"{CHAOS_ROOT}/dados_kam_beta.mat"
    savemat(out, dict(B=B), do_compression=True, oned_as="column")
    print(f"salvo {out} ({os.path.getsize(out)/1e6:.2f} MB), nx={NX}")


if __name__ == "__main__":
    main()

"""Robustez da sintonia do absorvedor a uma incerteza ±Δ em β, com a excitação em Ω fixo.

Computa os mapas KAM (FLI 2.5 GL conservativo) para β = β0-Δ, β0, β0+Δ (padrão β0=Ω=0.35,
Δ=0.10), nos três casos. Mostra se a proteção da antirressonância (β=Ω) sobrevive a uma
desintonização de ±Δ. Exporta dados_kam_omega{Ω}_robust.mat no formato de fig_kam_beta.
"""
import os
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
from scipy.io import savemat

import kam_compute as K

NAMES = ["monostable", "shallow_wells", "deep_wells"]
LABELS = ["Monostable", "Shallow wells", "Deep wells"]
OMEGA = float(os.environ.get("KAM_OMEGA", "0.35"))
B0 = float(os.environ.get("KAM_BETA0", "0.35"))
DELTA = float(os.environ.get("KAM_DELTA", "0.10"))
BETAS = [round(B0 - DELTA, 3), round(B0, 3), round(B0 + DELTA, 3)]
FFORCE = float(os.environ.get("KAM_F", "0.05"))
NX = int(os.environ.get("KAM_NX", "200"))
NS_BASE = 100
SNAPS = (200, 400)


def snap_map(nm, beta):
    eta = K.CASES[nm]["eta"]; m = K.CASES[nm]; b2 = beta ** 2; mb = K.MU * b2
    Tdrive = 2 * np.pi / OMEGA; ns = max(NS_BASE, int(round(NS_BASE / OMEGA))); h = Tdrive / ns
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
    return nm, beta, xs, vs, snaps[SNAPS[1]].reshape(NX, NX).astype(np.float32), reg


def main():
    jobs = [(nm, b) for nm in NAMES for b in BETAS]
    res = {}
    with ProcessPoolExecutor(max_workers=int(os.environ.get("KAM_WORKERS", "4"))) as pool:
        futs = [pool.submit(snap_map, nm, b) for nm, b in jobs]
        for k, fu in enumerate(as_completed(futs), 1):
            nm, beta, xs, vs, mp, reg = fu.result(); res[(nm, beta)] = (xs, vs, mp, reg)
            print(f"  {k}/{len(jobs)} {nm} Ω={OMEGA} β={beta}: regular={100*reg:.0f}%", flush=True)
    B = dict(omega=OMEGA, beta0=B0, delta=DELTA, f=FFORCE, betas=np.array(BETAS).reshape(1, -1),
             names=np.array(NAMES, dtype=object).reshape(1, -1),
             labels=np.array(LABELS, dtype=object).reshape(1, -1))
    for nm, lab in zip(NAMES, LABELS):
        eta = K.CASES[nm]["eta"]; m = K.CASES[nm]; Um = K.U(K.xmin(eta), eta); depth = K.U(0., eta) - Um
        xs, vs, _, _ = res[(nm, BETAS[0])]
        maps = np.stack([res[(nm, b)][2] for b in BETAS]).astype(np.float32)
        reg = np.array([res[(nm, b)][3] for b in BETAS]).reshape(1, -1)
        B[nm] = dict(eta=eta, xw=m["xw"], vw=m["vw"], Umin=float(Um), depth=float(depth),
                     xs=xs.reshape(1, -1), vs=vs.reshape(1, -1), maps=maps, reg=reg, label=lab)
    tag = f"omega{OMEGA:g}_robust".replace(".", "p")
    out = f"/home/user/Chaos/dados_kam_{tag}.mat"
    savemat(out, dict(B=B), do_compression=True, oned_as="column")
    print(f"salvo {out} ({os.path.getsize(out)/1e6:.2f} MB), nx={NX}, Ω={OMEGA}, β={BETAS}")


if __name__ == "__main__":
    main()

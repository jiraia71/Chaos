"""Sweep da sintonia β do absorvedor: qual β preserva melhor os toros KAM (conservativo)?

O caos extra do QZS-ADV vem da ressonância interna primário↔absorvedor em ω(E)=ω_a=β.
Este script varre β, e para cada β mede a fração regular do sistema conservativo 2.5 GL
(classificador por LEI DE CRESCIMENTO do FLI, independente da dimensão), na região
regular H0≤ECUT. Assim acha o β de afinação que MAXIMIZA os toros KAM e mostra como o
KAM muda com a sintonia.

μ é mantido em 0.1; só β varia (β²=BETA2). Forçamento fixo (f=0.05, o mais sensível),
ζ=0 (conservativo/simplético).

Uso:
  python beta_tuning.py            # preview (grade/β pequenos)
  python beta_tuning.py --full     # produção
"""
import argparse
import json
import os
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

import kam_compute as K

NAMES = ["monostable", "shallow_wells", "deep_wells"]
FFORCE = 0.05
TAU = 1.0
OUT = os.environ.get("KAM_OUT", "/home/user/Chaos/work/beta_tuning_result.json")

PRESET = {"preview": dict(nx=40, snaps=(60, 120), ns=40, betas=np.linspace(0.15, 1.8, 8)),
          "full":    dict(nx=84, snaps=(200, 400), ns=100, betas=np.round(np.linspace(0.15, 1.9, 24), 3))}


def snap_2dof_beta(x, v, x2, v2, eta, f, beta2, snaps, ns):
    """FLI 4D do 2.5 GL conservativo com β² arbitrário; snapshots de logg em `snaps`."""
    n = x.size; h = K.T / ns; t = 0.; mb = K.MU * beta2
    d = [np.ones(n) / 2 for _ in range(4)]; logg = np.zeros(n)
    out = {}; smax = max(snaps)
    for c in range(smax):
        for s in range(ns):
            for k in range(4):
                x = x + K.CY[k] * h * v; x2 = x2 + K.CY[k] * h * v2; t += K.CY[k] * h
                d[0] = d[0] + K.CY[k] * h * d[1]; d[2] = d[2] + K.CY[k] * h * d[3]
                if k == 3:
                    break
                coup = mb * (x - x2)
                v = v + K.DY[k] * h * (f * np.cos(t) - K.dU(x, eta) - coup); v2 = v2 + K.DY[k] * h * (beta2 * (x - x2))
                k2 = K.d2U(x, eta)
                d[1] = d[1] + K.DY[k] * h * (-(k2 + mb) * d[0] + mb * d[2]); d[3] = d[3] + K.DY[k] * h * (beta2 * (d[0] - d[2]))
        t = (c + 1) * K.T
        nrm = np.sqrt(d[0]**2 + d[1]**2 + d[2]**2 + d[3]**2); logg += np.log10(nrm); d = [q / nrm for q in d]
        if (c + 1) in snaps:
            out[c + 1] = logg.copy()
    return out


def job(spec):
    nm, beta, cfg = spec
    eta = K.CASES[nm]["eta"]; m = K.CASES[nm]
    xs = np.linspace(-m["xw"], m["xw"], cfg["nx"]); vs = np.linspace(-m["vw"], m["vw"], cfg["nx"])
    X, V = np.meshgrid(xs, vs); X = X.ravel(); V = V.ravel()
    snp = snap_2dof_beta(X.copy(), V.copy(), X.copy(), V.copy(), eta, FFORCE, beta**2, cfg["snaps"], cfg["ns"])
    H0 = .5 * V**2 + K.U(X, eta) - K.U(K.xmin(eta), eta); inside = H0 <= K.ECUT
    t1, t2 = cfg["snaps"]; growth = snp[t2] - snp[t1]
    reg = float(np.mean(growth[inside] <= TAU))
    return (nm, float(beta), reg)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", action="store_true")
    ap.add_argument("--workers", type=int, default=int(os.environ.get("KAM_WORKERS", "4")))
    args = ap.parse_args()
    cfg = PRESET["full" if args.full else "preview"]
    mode = "full" if args.full else "preview"
    betas = list(cfg["betas"])
    print(f"[{mode}] f={FFORCE} nx={cfg['nx']} snaps={cfg['snaps']} | {len(betas)} betas x {len(NAMES)} casos", flush=True)
    specs = [(nm, b, cfg) for nm in NAMES for b in betas]
    res = {nm: {} for nm in NAMES}
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futs = [pool.submit(job, s) for s in specs]
        for k, fu in enumerate(as_completed(futs), 1):
            nm, b, reg = fu.result(); res[nm][b] = reg
            print(f"  {k}/{len(specs)} {nm} β={b:.3f}: regular={100*reg:.0f}%", flush=True)
    print("\n===== FRAÇÃO REGULAR (KAM) vs β  —  f=0.05, conservativo =====")
    out = {"f": FFORCE, "tau": TAU, "mode": mode, "betas": betas, "beta_baseline": 0.35, "cases": {}}
    for nm in NAMES:
        bs = sorted(res[nm]); vals = [res[nm][b] for b in bs]
        bopt = bs[int(np.argmax(vals))]; ropt = max(vals)
        # baseline mais próximo de 0.35
        b035 = min(bs, key=lambda x: abs(x - 0.35)); r035 = res[nm][b035]
        print(f"\n{nm}:")
        print("  β    :", " ".join(f"{b:5.2f}" for b in bs))
        print("  reg% :", " ".join(f"{100*res[nm][b]:5.0f}" for b in bs))
        print(f"  → ótimo β={bopt:.2f} (reg={100*ropt:.0f}%)  vs  baseline β≈0.35 (reg={100*r035:.0f}%)  ganho {100*(ropt-r035):+.0f} p.p.")
        out["cases"][nm] = dict(betas=bs, regular=vals, beta_opt=float(bopt), reg_opt=float(ropt),
                                beta_baseline=float(b035), reg_baseline=float(r035))
    json.dump(out, open(OUT, "w"), indent=1)
    print(f"\nsalvo {OUT}")


if __name__ == "__main__":
    main()
